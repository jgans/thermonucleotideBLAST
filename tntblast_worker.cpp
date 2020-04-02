// ThermonucleotideBLAST
// 
// Copyright (c) 2007, Los Alamos National Security, LLC
// All rights reserved.
// 
// Copyright 2007. Los Alamos National Security, LLC. This software was produced under U.S. Government 
// contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos 
// National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, 
// reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, 
// LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  
// If software is modified to produce derivative works, such modified software should be clearly marked, 
// so as not to confuse it with the version available from LANL.
// 
// Additionally, redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:
// 
//      * Redistributions of source code must retain the above copyright notice, this list of conditions 
//        and the following disclaimer.
//      * Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
//        and the following disclaimer in the documentation and/or other materials provided with the distribution.
//      * Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, 
//        the U.S. Government, nor the names of its contributors may be used to endorse or promote products 
//        derived from this software without specific prior written permission.
// 
// 
// THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY 
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY 
// AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC 
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
// OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
		
		int max_len;
		int primer_clamp;
		int min_max_primer_clamp;
		int max_gap;
		int max_mismatch;
		int probe_clamp_5;
		int probe_clamp_3;
		int target_strand;
		float min_primer_tm;
		float max_primer_tm;
		float min_primer_dg;
		float max_primer_dg;
		float min_probe_tm;
		float max_probe_tm;
		float min_probe_dg;
		float max_probe_dg;
		float salt;
		float forward_primer_strand;
		float reverse_primer_strand;
		float probe_strand;
		float target_t;
		int melting_param;
		int mask_options;
		int assay_format;
		unsigned int output_format;
		int hash_word_size;
		int allow_fasta_mmap;
		int sequence_file_format;
		int single_primer_pcr;
		int allow_dangle_5;
		int allow_dangle_3;
		int use_dinkelbach;
		int best_match;

		vector<hybrid_sig> sig_list; // The queries
		bitmask sig_match; // Has a given query matched at least one target (for OUTPUT_INVERSE_QUERY)
		list<hybrid_sig> results_list; // The query-target matches
		
		string dbase_filename;
		
		MPI_Bcast(&salt, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&forward_primer_strand, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&reverse_primer_strand, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&probe_strand, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
		MPI_Bcast(&min_primer_tm, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&max_primer_tm, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
		MPI_Bcast(&min_primer_dg, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&max_primer_dg, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
		MPI_Bcast(&min_probe_tm, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&max_probe_tm, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
		MPI_Bcast(&min_probe_dg, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&max_probe_dg, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
		MPI_Bcast(&target_t, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
		MPI_Bcast(&primer_clamp, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&min_max_primer_clamp, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		MPI_Bcast(&probe_clamp_5, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&probe_clamp_3, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&max_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&target_strand, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&melting_param, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&mask_options, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&assay_format, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&output_format, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
		MPI_Bcast(&hash_word_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&allow_fasta_mmap, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&single_primer_pcr, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&allow_dangle_5, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&allow_dangle_3, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&use_dinkelbach, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&best_match, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&max_gap, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&max_mismatch, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		DNAHash dbase(hash_word_size);
		
		// The database file from the master
		receive(dbase_filename);
		
		MPI_Bcast(&sequence_file_format, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		// Get the queries from the master
		receive_queries(sig_list);
		
		const unsigned long int num_sig = sig_list.size();
		
		// If the user has selected OUTPUT_INVERSE_QUERY for their output format,
		// we will only track whether each query has matched one or more targets
		if(output_format & OUTPUT_INVERSE_QUERY){
			sig_match.resize(num_sig, false);
		}
		
		// Worker can either directly read sequence data or get it from the master
		sequence_data seq_file;		
		
		// Suppress output from the workers
		seq_file.verbose(false);
		
		if(dbase_filename != ""){
			
			if(seq_file.wants_indicies(sequence_file_format) == true){
				
				// Load the indicies from the master
				deque<file_index> tmp;
				
				receive(tmp);
				
				seq_file.indicies(tmp);
			}
			
			seq_file.open( (char*)dbase_filename.c_str(), allow_fasta_mmap, false);
		}
		
		// Initialize the melting engine. There is a fair amount of overhead involved in
		// initialization (handled by the constructor) so it is best to do it just once
		// per program invocation.
		NucCruc melt(melting_param, target_t);
		
		// Boolean values are transmitted as ints (so we need to convert back)
		melt.dangle( allow_dangle_5 == int(true), allow_dangle_3 == int(true) );
		
		melt.dinkelbach( use_dinkelbach == int(true) );
		
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

					switch(assay_format){

						case ASSAY_PCR:

							// What amplicons do these primers/probe produce?
							local_results = amplicon(dbase, bio_seq,
								sig_ref, melt,
								salt, forward_primer_strand, reverse_primer_strand, probe_strand,
								min_primer_tm, max_primer_tm, 
								min_primer_dg, max_primer_dg, 
								min_probe_tm, max_probe_tm, 
								min_probe_dg, max_probe_dg, 
								primer_clamp, min_max_primer_clamp,
								probe_clamp_5, probe_clamp_3, 
								max_gap, max_mismatch,
								max_len,
								(single_primer_pcr == true) );

							mask_binding_sites(local_results, mask_options,
								min_primer_tm, min_probe_tm, melt,
								salt, forward_primer_strand, reverse_primer_strand, probe_strand);

							break;
						case ASSAY_PADLOCK:

							local_results = padlock(dbase, bio_seq,
								sig_ref, melt,
								salt, forward_primer_strand, reverse_primer_strand, 
								min_probe_tm, max_probe_tm,
								min_probe_dg, max_probe_dg,
								probe_clamp_5, probe_clamp_3, 
								max_gap, max_mismatch,
								target_strand);

							break;
					};
				}
				else{

					if(sig_ref.has_probe() == true){
						
						local_results = hybrid(dbase, bio_seq,
							sig_ref, melt, salt, probe_strand, 
							min_probe_tm, max_probe_tm, 
							min_probe_dg, max_probe_dg, 
							probe_clamp_5, probe_clamp_3,
							max_gap, max_mismatch, 
							target_strand);
					}
				}
				
				// If the user has selected OUTPUT_INVERSE_QUERY for their output format,
				// we will only track whether each query has matched one or more targets
				if(output_format & OUTPUT_INVERSE_QUERY){
					
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
						melt.Strand(forward_primer_strand);

						local_iter->forward_hairpin_tm = melt.approximate_tm_hairpin();
						local_iter->forward_dimer_tm = melt.approximate_tm_homodimer();

						melt.set_duplex(local_iter->reverse_oligo);
						melt.Strand(reverse_primer_strand);

						local_iter->reverse_hairpin_tm = melt.approximate_tm_hairpin();
						local_iter->reverse_dimer_tm = melt.approximate_tm_homodimer();

						melt.set_query(local_iter->forward_oligo);
						melt.set_target(local_iter->reverse_oligo);

						// The strand is the *total* concentration. For forward+reverse 
						// heterodimers we must take the sum of strands.
						melt.Strand(forward_primer_strand + reverse_primer_strand);

						local_iter->primer_dimer_tm = melt.approximate_tm_heterodimer();
					}
					
					if(local_iter->has_probe() == true){
					
						melt.set_duplex(local_iter->probe_oligo);

						melt.strand(probe_strand);

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
				if(best_match == true){
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
		if(output_format & OUTPUT_INVERSE_QUERY){

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
					
					num_results ++;
					buffer_size += iter->mpi_size();
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
					ptr = iter->mpi_pack(ptr);
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
