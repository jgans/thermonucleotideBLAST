#ifdef USE_MPI

#include <mpi.h>
#include <stdlib.h>

#include "tntblast.h"
#include "options.h"
#include "mpi_util.h"
#include "bitmask.h"

using namespace std;


// Global variables
extern int mpi_numtasks;
extern int mpi_rank;

#ifdef PROFILE
extern unsigned int num_plus_tm_eval;
extern unsigned int num_minus_tm_eval;
#endif // PROFILE

int worker(int argc, char *argv[])
{

	try{
	
		#ifdef MEMORY_CHECK
		test_memory(MEMORY_CHECK);
		#endif // MEMORY_CHECK
		
		int continue_exec;
		MPI_Bcast(&continue_exec, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		if(continue_exec == false){
			return EXIT_SUCCESS;
		}

		vector<hybrid_sig> sig_list; // The queries
		bitmask sig_match; // Has a given query matched at least one target (for OUTPUT_INVERSE_QUERY)
		list<hybrid_sig> results_list; // The query-target matches
		
		Options opt;

		// Receive the command-line options from the master
		broadcast(opt, mpi_rank, 0);
		
		// Allow different strand concentrations for the forward
		// and reverse primers (i.e. asymmetric PCR). When
		// opt.asymmetric_strand_ratio != 1, the opt.primer_strand
		// is assumed to be the concentration of the reverse primer.
		const float forward_primer_strand = opt.asymmetric_strand_ratio*
			opt.primer_strand;
		const float reverse_primer_strand = opt.primer_strand;
		
		DNAHash dbase(opt.hash_word_size);
		
		// The database file from the master. This may overwrite the filename
		// provided on the command line, depending on the format of the file
		receive(opt.dbase_filename);
		
		int sequence_file_format;

		MPI_Bcast(&sequence_file_format, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		// Get the queries from the master
		receive_queries(sig_list);
		
		const unsigned long int num_sig = sig_list.size();
		
		// If the user has selected OUTPUT_INVERSE_QUERY for their output format,
		// we will only track whether each query has matched one or more targets
		if(opt.output_format & OUTPUT_INVERSE_QUERY){
			sig_match.resize(num_sig, false);
		}
		
		// Worker can either directly read sequence data or get it from the master
		sequence_data seq_file;		
		
		// Suppress output from the workers
		seq_file.verbose(false);
		
		if(opt.dbase_filename != ""){
			
			if(seq_file.wants_indicies(sequence_file_format) == true){
				
				// Load the indicies from the master
				deque<file_index> tmp;
				
				receive(tmp);
				
				seq_file.indicies(tmp);
			}
			
			seq_file.open( (char*)opt.dbase_filename.c_str(), opt.blast_include, opt.blast_exclude);
		}
		
		// Initialize the melting engine. There is a fair amount of overhead involved in
		// initialization (handled by the constructor) so it is best to do it just once
		// per program invocation.
		NucCruc melt(opt.melting_param, opt.target_t);
	
		melt.Salt(opt.salt);

		melt.dangle(opt.allow_dangle_5, opt.allow_dangle_3);
		melt.dinkelbach(opt.use_dinkelbach);
		
		unordered_map<BindCacheKey, BindCacheValue> plus_strand_melt_cache;
		unordered_map<BindCacheKey, BindCacheValue> minus_strand_melt_cache;

		pair<string, SEQPTR> bio_seq = make_pair("", SEQPTR(NULL) );
		long int last_target = -1;
		unsigned int last_target_start = 0;
		unsigned int target_len = 0;
		
		#ifdef PROFILE
		// How much time are we spending doing work, and how much time are we spending
		// on communication.
		double profile_time[NUM_PROFILE] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

		double profile = 0.0;
		#endif // PROFILE
		
		// For adaptive query segmentation, we need to know the ratio of the time spent
		// searching a single query (Q_time) to the time spend loading and hashing a single
		// target (T_time).		
		double T_time = -1.0;
		
		// Wait for the master to dispatch work for us to do
		while(true){
			
			#ifdef PROFILE
			profile = MPI_Wtime();
			#endif // PROFILE
			
			MPI_Status status;
			MPI_Probe(0 /* listen to the master */, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			
			#ifdef PROFILE
			profile_time[PROFILE_IDLE] += MPI_Wtime() - profile;
			profile = MPI_Wtime();
			#endif // PROFILE
			
			if(status.MPI_TAG == SEARCH_COMPLETE){
			
				// Break out of the message loop and send our results back to the master
				int msg;
				
				if(MPI_Recv(&msg, 1, MPI_INT, status.MPI_SOURCE, 
					SEARCH_COMPLETE, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS){

					throw __FILE__ ":worker: Error receiving SEARCH_COMPLETE";
				}
				
				break;
			}
			
			// If we get here, the message had better be "SEARCH_QUERY"
			if(status.MPI_TAG != SEARCH_QUERY){
			
				cerr << "[" << mpi_rank << "] got an unknown message tag: " << status.MPI_TAG << endl;
				throw __FILE__ ":worker: Unknown message tag";
			}
			
			unsigned int buffer[SEARCH_QUERY_BUFFER_SIZE];
			
			if(MPI_Recv(buffer, SEARCH_QUERY_BUFFER_SIZE, MPI_UNSIGNED, status.MPI_SOURCE, 
				SEARCH_QUERY, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS){

				throw __FILE__ ":worker: Error receiving SEARCH_QUERY";
			}
			
			unsigned int cur_query = buffer[0];
			const unsigned int num_query = buffer[1];
			const unsigned int cur_target = buffer[2];
			const unsigned int cur_target_start = buffer[3];
			const unsigned int cur_target_stop = buffer[4];
			const unsigned int cur_target_max_stop = buffer[5];
			
			#ifdef PROFILE
			profile_time[PROFILE_COMM] += MPI_Wtime() - profile;
			#endif // PROFILE
			
			// Read a DNA sequence (and defline) from either the database or the master
			// node.
			/////////////////////////////////////////////////////////////////////////////////////
			// Load the database sequence if not already loaded. There is no need to test both 
			// cur_target_start and cur_target_stop (i.e. cur_target_start will suffice).
			const bool same_target = (last_target == (long int)cur_target) && 
				(last_target_start == cur_target_start);
			
			if(!same_target){

				#ifdef PROFILE
				profile = MPI_Wtime();
				#endif // PROFILE
				
				T_time = MPI_Wtime();
				
				// Free the memory used to store the sequence
				if(bio_seq.second != NULL){

					delete [] bio_seq.second;
					bio_seq.second = NULL;
				}
				
				// Free any existing cached binding data (which is only valid for a single target sequence).
				// There is no guarentee that unordered_map will free memory when calling unordered_map::clear().
				// To be on the same size, use the swap trick
				unordered_map<BindCacheKey, BindCacheValue>().swap(plus_strand_melt_cache);
				unordered_map<BindCacheKey, BindCacheValue>().swap(minus_strand_melt_cache);
				
				// Read a DNA sequence (and defline) from either the database or the master
				// node.
				target_len = seq_file.read_bio_seq(bio_seq, cur_target, cur_target_start, cur_target_stop);
				
				if( target_len < dbase.min_sequence_size() ){

					// This sequence is too small to hash! Genbank contains some very small
					// DNA sequences (4-6 nucleotides). These sequences include ligands
					// bound to proteins in the PDB.
					if(bio_seq.second != NULL){

						// Free the memory used to store the sequence
						delete [] bio_seq.second;
						bio_seq.second = NULL;
					}
			
					// Tell the master that we're ready for more work
					if(MPI_Send( (void*)&num_query, 1, MPI_UNSIGNED, 0, 
					   STATUS_UPDATE, MPI_COMM_WORLD) != MPI_SUCCESS){

					   throw __FILE__ ":worker: Error sending STATUS_UPDATE to master";
					}
					
					// There is no stored sequence
					last_target = -1;
				
					// Don't bother accumulating profile information for the
					// small sequences that are not searched.
					
					continue;
				}
				
				#ifdef PROFILE
				profile_time[PROFILE_LOAD] += MPI_Wtime() - profile;
				profile = MPI_Wtime();
				#endif // PROFILE
				
				// Hash the sequence for fast searching. Note that hash_dbase will
				// re-use previously allocated memory, so there is no need to deallocate
				// dbase after every sequence
				dbase.hash( bio_seq.second, SEQ_SIZE(bio_seq.second), 0, SEQ_SIZE(bio_seq.second) );
			
				#ifdef PROFILE
				profile_time[PROFILE_HASH] += MPI_Wtime() - profile;
				#endif // PROFILE
				
				last_target = cur_target;
				last_target_start = cur_target_start;
				
				// The time to load and hash a single target sequence
				T_time = MPI_Wtime() - T_time;
 			}
			
			#ifdef PROFILE
			profile = MPI_Wtime();
			#endif // PROFILE
			
			double Q_time = MPI_Wtime();
			
			for(unsigned int i = 0;i < num_query;i++, cur_query++){
			
				if(cur_query >= num_sig){
					throw __FILE__ ": cur_query >= num_sig";
				}
				
				list<hybrid_sig> local_results;
				const hybrid_sig &sig_ref = sig_list[cur_query];

				if(sig_ref.has_primers() == true){

					switch(opt.assay_format){

						case ASSAY_PCR:

							// What amplicons do these primers/probe produce?
							local_results = amplicon(dbase, bio_seq,
								sig_ref, melt, plus_strand_melt_cache, minus_strand_melt_cache,
								forward_primer_strand, reverse_primer_strand, opt.probe_strand,
								opt.min_primer_tm, opt.max_primer_tm, 
								opt.min_primer_dg, opt.max_primer_dg, 
								opt.min_probe_tm, opt.max_probe_tm, 
								opt.min_probe_dg, opt.max_probe_dg, 
								opt.primer_clamp, opt.min_max_primer_clamp,
								opt.probe_clamp_5, opt.probe_clamp_3, 
								opt.max_gap, opt.max_mismatch,
								opt.max_len,
								opt.single_primer_pcr);

							mask_binding_sites(local_results, opt.mask_options,
								opt.min_primer_tm, opt.min_probe_tm, melt,
								forward_primer_strand, reverse_primer_strand, opt.probe_strand);

							break;
						case ASSAY_PADLOCK:

							local_results = padlock(dbase, bio_seq,
								sig_ref, melt, plus_strand_melt_cache, minus_strand_melt_cache,
								forward_primer_strand, reverse_primer_strand, 
								opt.min_probe_tm, opt.max_probe_tm,
								opt.min_probe_dg, opt.max_probe_dg,
								opt.probe_clamp_5, opt.probe_clamp_3, 
								opt.max_gap, opt.max_mismatch,
								opt.target_strand);

							break;
					};
				}
				else{

					if(sig_ref.has_probe() == true){
						
						local_results = hybrid(dbase, bio_seq,
							sig_ref, melt, opt.probe_strand, 
							opt.min_probe_tm, opt.max_probe_tm, 
							opt.min_probe_dg, opt.max_probe_dg, 
							opt.probe_clamp_5, opt.probe_clamp_3,
							opt.max_gap, opt.max_mismatch, 
							opt.target_strand);
					}
				}
				
				// If the user has selected OUTPUT_INVERSE_QUERY for their output format,
				// we will only track whether each query has matched one or more targets
				if(opt.output_format & OUTPUT_INVERSE_QUERY){
					
					// Have we found one or more matches to this query?
					if(local_results.empty() == false){
						sig_match[cur_query] = true;
					}
					
					continue;
				}
				
				// Before we add the local_results to the main results list, add the 
				// sequence id so we know which sequence produced this match
				list<hybrid_sig>::iterator local_iter;
				
				list< list<hybrid_sig>::iterator > reaper;
				
				for(local_iter = local_results.begin();local_iter != local_results.end();local_iter++){
						
					// Is this match truncated due to target sequence fragmentation?
					if( (cur_target_start != 0) && local_iter->start_overlap(0) ){
					
						reaper.push_back(local_iter);
						continue;
					}
					
					if( (cur_target_stop != cur_target_max_stop) && local_iter->stop_overlap(target_len - 1) ){
					
						reaper.push_back(local_iter);
						continue;
					}
					
					local_iter->seq_id(cur_target);
					
					// Since we may be processing a subsequence, add the appropriate offset
					// so that match coordinates are correct
					local_iter->offset_ranges(cur_target_start);
					
					// Compute the dimer and hairpin temperatures for this assay
					if(local_iter->has_primers() == true){
					
						melt.set_duplex(local_iter->forward_oligo);
						melt.Strand(forward_primer_strand, forward_primer_strand);

						local_iter->forward_hairpin_tm = melt.approximate_tm_hairpin();
						local_iter->forward_dimer_tm = melt.approximate_tm_homodimer();

						melt.set_duplex(local_iter->reverse_oligo);
						melt.Strand(reverse_primer_strand, reverse_primer_strand);

						local_iter->reverse_hairpin_tm = melt.approximate_tm_hairpin();
						local_iter->reverse_dimer_tm = melt.approximate_tm_homodimer();

						melt.set_query(local_iter->forward_oligo);
						melt.set_target(local_iter->reverse_oligo);

						melt.Strand(forward_primer_strand, reverse_primer_strand);

						local_iter->primer_dimer_tm = melt.approximate_tm_heterodimer();
					}
					
					if(local_iter->has_probe() == true){
					
						melt.set_duplex(local_iter->probe_oligo);

						melt.strand(opt.probe_strand, opt.probe_strand);

						local_iter->probe_hairpin_tm = melt.approximate_tm_hairpin();
						local_iter->probe_dimer_tm = melt.approximate_tm_homodimer();
					}
				}
				
				// Remove truncated matches before we save these local results
				while(reaper.empty() == false){
					
					local_results.erase( reaper.back() );
					reaper.pop_back();
				}
				
				results_list.splice(results_list.end(), local_results);

				// If enabled, only keep the best matches between a given query and a given target
				if(opt.best_match == true){
					select_best_match(results_list);
				}
			}
			
			#ifdef PROFILE
			profile_time[PROFILE_WORK] += MPI_Wtime() - profile;
			profile = MPI_Wtime();
			#endif // PROFILE
			
			// The average time required to search the target sequence with a *single* query
			Q_time = (num_query > 0) ? (MPI_Wtime() - Q_time)/num_query : -1.0;
			
			float QT = ( (Q_time < 0.0) || (T_time < 0.0) ) ? -1.0f: float(Q_time/T_time);
			
			// Tell the master that we're ready for more work
			if(MPI_Send((void*)&QT, sizeof(float), MPI_BYTE, 0, 
			   STATUS_UPDATE, MPI_COMM_WORLD) != MPI_SUCCESS){

			   throw __FILE__ ":worker: Error sending STATUS_UPDATE to master";
			}
			
			#ifdef PROFILE
			profile_time[PROFILE_COMM] += MPI_Wtime() - profile;
			#endif // PROFILE
		}
		
		// Clean up any previously allocated sequence
		if(bio_seq.second != NULL){

			delete [] bio_seq.second;
			bio_seq.second = NULL;
		}
		
		// Remove the no-longer-needed input signatures
		sig_list.clear();
		
		#ifdef PROFILE
		profile_time[PROFILE_NUM_PLUS_TM_EVAL] = num_plus_tm_eval; 	// Copy from global var
		profile_time[PROFILE_NUM_MINUS_TM_EVAL] = num_minus_tm_eval;	// Copy from global var
		
		// Send all of the profile information back to the master
		if(MPI_Send(profile_time, NUM_PROFILE, MPI_DOUBLE, 0, 
		   PROFILE_INFO, MPI_COMM_WORLD) != MPI_SUCCESS){

		   throw __FILE__ ":worker: Error sending PROFILE_INFO to master";
		}
		#endif // PROFILE
		
		// The number of query results to send back in a single chunk
		unsigned int query_chunck_size = 0;
		int num_chunck = 0;
		
		// If the user has selected OUTPUT_INVERSE_QUERY for their output format,
		// we will only send back the bitmask of queries that *did* match a target
		if(opt.output_format & OUTPUT_INVERSE_QUERY){

			unsigned int buffer_size = sig_match.mpi_size();
			
			unsigned char* buffer = new unsigned char [buffer_size];
		
			if(buffer == NULL){
				throw __FILE__ ":worker: Error allocating send buffer for bitmask";
			}
			
			sig_match.mpi_pack(buffer);
			
			if(MPI_Ssend(buffer, buffer_size, MPI_BYTE, 0, 
			   SIGNATURE_RESULTS, MPI_COMM_WORLD) != MPI_SUCCESS){

			   throw __FILE__ ":worker: Error sending SIGNATURE_RESULTS (bitmask) to master";
			}
			
			delete [] buffer;
		}
		else{
		
			MPI_Bcast(&query_chunck_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
			num_chunck = num_sig/query_chunck_size + (num_sig%query_chunck_size != 0);
		}
		
		for(int sig_index = 0;sig_index < num_chunck;sig_index++){
			
			const int max_sig = (sig_index + 1)*query_chunck_size;
			
			// Return all of the results for query "sig_index" back to the master for output to the user
			unsigned int num_results = 0;
			
			// Track the amount of buffer space required to send the results back to the master
			unsigned int buffer_size = sizeof(unsigned int);
			
			typedef list<hybrid_sig>::iterator I;
			
			I iter;
			list<I> reaper;
			
			for(iter = results_list.begin();iter != results_list.end();iter++){
				
				if(iter->my_id() < max_sig){
					
					++num_results;
					buffer_size += mpi_size(*iter);
					reaper.push_back(iter);
				}
			}
			
			unsigned char* buffer = new unsigned char [buffer_size];
		
			if(buffer == NULL){
				throw __FILE__ ":worker: Error allocating send buffer";
			}

			unsigned char *ptr = buffer;

			memcpy( ptr, &num_results, sizeof(unsigned int) );
			ptr += sizeof(unsigned int);

			for(iter = results_list.begin();iter != results_list.end();iter++){
			
				if(iter->my_id() < max_sig){
					ptr = mpi_pack(ptr, *iter);
				}
			}
			
			// Send the results back to the master
			if(MPI_Ssend(buffer, buffer_size, MPI_BYTE, 0, 
			   SIGNATURE_RESULTS, MPI_COMM_WORLD) != MPI_SUCCESS){

			   throw __FILE__ ":worker: Error sending SIGNATURE_RESULTS to master";
			}

			delete [] buffer;		

			// Free results as we go
			while(reaper.empty() == false){
				
				results_list.erase( reaper.back() );
				reaper.pop_back();
			}
		}
		
		if(results_list.empty() == false){
			
			throw __FILE__ ":worker: results_list still has entries!";
		}
		
	}
	catch(const char *error){
	
		cerr << "Caught the worker [" << mpi_rank << "] error: " << error << endl;
		
		const char *name = getenv("hostname");
		
		if(name != NULL){
			cerr << "Worker rank [" << mpi_rank << "] -> " << name << endl;
		}
		
		return EXIT_FAILURE;
	}
	catch(std::exception &e){
		
		cerr << "Caught the std exception: " << e.what() << endl;
		return EXIT_FAILURE;
	} 
	catch(...){
		
		cerr << "Caught an unhandled worker [" << mpi_rank << "] error" << endl;
		
		const char *name = getenv("hostname");
		
		if(name != NULL){
			cerr << "Worker rank [" << mpi_rank << "] -> " << name << endl;
		}
		
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

#endif // USE_MPI
