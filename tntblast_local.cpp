#include "tntblast.h"
#include "options.h"
#include "hybrid_sig.h"
#include "degenerate_na.h"
#include "primer.h"
#include "bitmask.h"

#include <stdlib.h>
#include <time.h>

#ifdef _OPENMP
// Under windows, we need to include omp.h to load
// vcomp.dll (which is required for openMP on windows)
#include <omp.h>
#endif // _OPENMP

#include <iostream>
#include <iomanip>
#include <map>
#include <sstream>

using namespace std;

int local_main(int argc, char *argv[])
{
		
	try{
		ofstream fout;
		ostream *ptr_out = NULL;
		
		// For writing network (i.e. cytoscape) output files
		ofstream fout_atr;
		ofstream fout_sif;
		
		Options opt;
		
		try{
			opt.parse(argc, argv);
		}
		catch(const char *error){
			
			cerr << "Input error: " << error << endl;			
			return EXIT_FAILURE;
		}
		catch(...){
		
			cerr << "Unhandled input error, please report to " 
				<< EMAIL_ADDRESS << endl;			
			return EXIT_FAILURE;
		}
		
		// Is the user printing the command line arguments and then exiting?
		if(opt.print_usage){
			return EXIT_FAILURE;
		}
		
		if(opt.input_filename != ""){
		
			if(opt.verbose){
				cout << "Reading assays from " << opt.input_filename << endl;
			}
			
			read_input_file(opt.input_filename, opt.sig_list, 
				opt.ignore_probe, (opt.assay_format == ASSAY_PROBE) );
		}
		
		// Bind either stdout of fout to ptr_out
		if(opt.output_filename == ""){
		
			ptr_out = &cout;
		}
		else{
		
			if(opt.one_output_file_per_query == false){
				
				if( (opt.output_format & OUTPUT_STANDARD) || 
				    (opt.output_format & OUTPUT_FASTA) ){
				    
					fout.open( opt.output_filename.c_str() );

					if(!fout){
						throw "Unable to open output file";
					}
				}
				
				if(opt.output_format & OUTPUT_NETWORK){
				    
				    	const string filename_sif =  opt.output_filename + ".sif";
				    
					fout_sif.open( filename_sif.c_str() );

					if(!fout_sif){
					
						cerr << "Unable to open " << filename_sif << endl;
						throw ":local_main: I/O error";
					}
				}
			}
			
			// There is only one attribute file per run -- even if the user
			// has selected one_output_file_per_query == true
			if(opt.output_format & OUTPUT_NETWORK){

				const string filename_atr =  opt.output_filename + ".atr";

				fout_atr.open( filename_atr.c_str() );

				if(!fout_atr){

					cerr << "Unable to open " << filename_atr << endl;
					throw ":local_main: I/O error";
				}

				// Write the attribute header
				fout_atr << "FunctionalCatagory" << endl;
			}
			
			if( (opt.output_format & OUTPUT_INVERSE_TARGET) ||
			    (opt.output_format & OUTPUT_INVERSE_QUERY) ){
			
				fout.open( opt.output_filename.c_str() );

				if(!fout){
					throw "Unable to open output file";
				}
			}

			ptr_out = &fout;
		}
		
		// Consider all multiplex combinations of primers and probes
		if(opt.multiplex){
			opt.sig_list = multiplex_expansion(opt.sig_list, opt.assay_format);
		}
		
		// Expand the primers/probes (if needed);
		opt.sig_list = expand_degenerate_signatures(opt.sig_list, opt.degen_rescale_ct);
		
		if(opt.dump_query){
			
			// Write all queries to stdout
			opt.write_queries(cout);
		}
		
		// Make sure that the user has provided the search contraints (in Tm or Delta G)
		// that are appropriate for the given assays
		opt.validate_search_threshold();
		
		const size_t num_sig = opt.sig_list.size();
		
		if(num_sig == 0){
			throw __FILE__ ":local_main: No primers or probes found!";
		}
		
		// Count the number of probe only queries in the list of queries
		const unsigned int num_probes = probe_only_count(opt.sig_list);
		
		const unsigned int max_product_length = opt.max_product_length() + 2; // Room for dangling end bases
		
		// Read the sequence data base
		sequence_data seq_file;
			
		if(opt.dbase_filename != ""){
			
			if(opt.verbose){
				cout << "Reading sequence database: " << opt.dbase_filename << endl;
			}
			
			seq_file.open(opt.dbase_filename, opt.blast_include, opt.blast_exclude);
		}
		else{
			
			if(opt.verbose){
				cout << "Reading sequence database: " << opt.local_dbase_filename << endl;
			}
			
			seq_file.open(opt.local_dbase_filename, opt.blast_include, opt.blast_exclude);
		}
		
		// How many sequences are in the database? The user can specify either a global
		// database (dbase) accessible to all nodes, or a local database(local_dbase) that
		// is only visible to the master node.
		const size_t num_seq = seq_file.size();
		
		if(num_seq == 0){
			throw __FILE__ ":local_main: Empty database -- no sequences found!";
		}
		
		// If we take into account target sequence fragmentation, how many sequences
		// are there?
		unsigned long int effective_num_seq = seq_file.effective_size(opt.fragment_target_threshold);
		
		if(opt.verbose){
		
			// Let the user know what parameters are being used
			cout << "Found " << num_seq << " database sequences";
			
			if(num_seq == effective_num_seq){
				cout << endl;
			}
			else{
				cout << " (" << effective_num_seq << " after fragmentation)" << endl;
			}
			
			// Write a summary of program options to the screen
			cout << opt;
		}
		
		// Track the total time required to perform the search
		time_t profile = time(NULL);
				
		// Allow different strand concentrations for the forward
		// and reverse primers (i.e. asymmetric PCR). When
		// opt.asymmetric_strand_ratio != 1, the opt.primer_strand
		// is assumed to be the concentration of the reverse primer.
		const float forward_primer_strand = opt.asymmetric_strand_ratio*
			opt.primer_strand;
		const float reverse_primer_strand = opt.primer_strand;
				
		// Output statistics
		pair<float, float> forward_tm_range = make_pair(9999.0f, -1.0f);
		pair<float, float> reverse_tm_range = make_pair(9999.0f, -1.0f);;
		pair<float, float> probe_tm_range = make_pair(9999.0f, -1.0f);
		
		pair<float, float> forward_dg_range = make_pair(9999.0f, -9999.0f);
		pair<float, float> reverse_dg_range = make_pair(9999.0f, -9999.0f);;
		pair<float, float> probe_dg_range = make_pair(9999.0f, -9999.0f);
		
		float max_primer_hairpin = -1.0f;
		float max_primer_homodimer = -1.0f;
		float max_primer_heterodimer = -1.0f;
		
		float max_probe_homodimer = -1.0f;
		float max_probe_hairpin = -1.0f;
		
		pair<float, float> forward_gc_range = make_pair(9999.0f, -1.0f);
		pair<float, float> reverse_gc_range = make_pair(9999.0f, -1.0f);
		pair<float, float> probe_gc_range = make_pair(9999.0f, -1.0f);
		
		pair<unsigned int, unsigned int> amplicon_size_range = make_pair(9999,0);
		
		pair<unsigned int, unsigned int> forward_size_range = make_pair(9999,0);
		pair<unsigned int, unsigned int> reverse_size_range = make_pair(9999,0);
		pair<unsigned int, unsigned int> probe_size_range = make_pair(9999,0);
		
		unsigned int num_primer = 0;
		unsigned int num_probe = 0;
		
		vector< list<hybrid_sig> > search_results(num_sig);
		
		bitmask query_matches; // Has a given query matched at least one target (for OUTPUT_INVERSE_QUERY)
		
		// If the user has selected OUTPUT_INVERSE_QUERY for their output format,
		// we will only track whether each query has matched one or more targets
		if(opt.output_format & OUTPUT_INVERSE_QUERY){
			query_matches.resize(num_sig, false);
		}
		
		const unsigned int update_buffer_size = 15;
		const double inv_total_num_of_comparisons = 1.0/( (double)(num_seq)*(double)(num_sig) );
		float search_display_every = 0.01f;
		unsigned int search_display_precision = 3;
		
		unsigned int cur_target = 0;
		
		unsigned int cur_target_len = seq_file.approx_seq_len(cur_target);
		unsigned int cur_target_max_stop = cur_target_len - 1;
		unsigned int cur_target_delta = seq_len_increment(cur_target_len, 
			opt.fragment_target_threshold).first;
		
		// Start and stop are inclusive
		unsigned int cur_target_start = 0;
		unsigned int cur_target_stop = cur_target_delta;
		
		// This flag is set of we have fragmented *any* target sequences. If we haven't
		// fragmented, then we can save some time by skipping the uniquify_results()
		// function.
		bool fragment_target = false;
		
		unsigned int cur_query = 0;
		bool segment_queries = false;
		
		// Track the ratio of time spent searching an individual query against a target
		// sequence to the time spent loading and hashing a target sequence
		float total_QT = 0.0f;
		unsigned int QT_count = 0;
		
		// To estimate the cost of query seqmentation we need the ratio of query search 
		// time to target sequence load and hash time. Scale the defined value according to
		// the assay format:
		// 1XDEFAULT_QT	for PROBE	query targeting a single strand
		// 2XDEFAULT_QT	for PROBE queries targeting both strands
		// 4XDEFAULT_QT	for PCR queries
		// 4XDEFAULT_QT	for PADLOCK queries
		const float default_qt = DEFAULT_QT*( 
			num_probes*( (opt.target_strand == Seq_strand_both) ?  2.0f: 1.0f) +
			(num_sig - num_probes)*4.0f
			)/num_sig;
			
#pragma omp parallel \
	shared(search_results, cur_target, cur_target_len, cur_target_max_stop, cur_target_delta, cur_target_start, \
		cur_target_stop, fragment_target, cur_query, segment_queries, effective_num_seq, total_QT, QT_count)
{
		
		// Parallelization strategy: 
		// 1) Always segment the database targets
		// 2) Segment the assay queries if the number of
		//    remaining targets is *less* than the number of
		//    workers (times a constant greater than 1).
		// 3) Fragment the target sequence into sub-sequences if the target sequence
		// 	 length is larger than opt.fragment_target_threshold.
		
		// When should we start sending individual queries out to the workers?
		// Too soon and we incur a high overhead penalty. Too late and we incur
		// a load balancing penalty.
		
		#ifdef _OPENMP
		const unsigned int num_worker = omp_get_num_threads();
		#else
		const unsigned int num_worker = 1;
		#endif // _OPENMP
		
		#pragma omp master
		{
			segment_queries = query_sched(effective_num_seq, num_sig, num_worker,
				(QT_count == 0 ? default_qt : total_QT/QT_count),
				opt.query_segmentation);
				
			cur_query = segment_queries ? 0 : num_sig;
		}
		
		DNAHash dbase(opt.hash_word_size);
		
		float last_search_status = 0.0f;
		
		#pragma omp master
		if(opt.verbose){
			
			cout << "Searching database: ";
			cout.flush();

			for(unsigned int i = 0;i < update_buffer_size;i++){
				cout << ' ';
			}
		}
		
		// Initialize the melting engine. There is a fair amount of overhead involved in
		// initialization (handled by the contructor) so it is best to do it just once
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
		
		tnt_time T_time;
		
		while(true){

			unsigned int local_target;
			unsigned int local_query;
			
			unsigned int local_target_start;
			unsigned int local_target_stop;
			unsigned int local_target_max_stop;

			// Use a critical section to avoid race conditions while allowing
			// each thread to update the global sequence counter for 
			// load balancing (rather than use a parallel for-loop)
			#pragma omp critical (UpdateIndex)
			{
				local_target = cur_target;
				local_query = cur_query;
				
				local_target_start = cur_target_start;
				local_target_stop = cur_target_stop;
				local_target_max_stop = cur_target_max_stop;
				
				bool increment_target = false;
				
				if(segment_queries){ // We're segmenting queries. Once we start, we don't stop!
				
					// In order to minimize network load, the assay queries are the 
					// inner-loop, while the database targets are the outer-loop
					cur_query ++;

					if(cur_query == num_sig){
						
						increment_target = true;
						cur_query = 0;
					}
				}
				else{
					
					// We're *not* segmenting queries. Every database target is
					// screened against every query assay by a single worker
					increment_target = true;
					
					// Update the query segmentation status (that is, do we need to
					// start segmenting queries?)
					segment_queries = query_sched(effective_num_seq, num_sig, num_worker,
						(QT_count == 0 ? default_qt : total_QT/QT_count),
						opt.query_segmentation);
					
					if(segment_queries){
					
						// We're about to start segmenting queries for the first time,
						// so reset the query counter to 0.
						cur_query = 0;
					}
				}
				
				if(increment_target == true){
					
					// Don't accidentally force an unsigned int below zero!
					effective_num_seq -= (effective_num_seq == 0) ? 0 : 1;
					
					if( cur_target_stop == cur_target_max_stop ){

						cur_target ++;
						
						cur_target_len = seq_file.approx_seq_len(cur_target);
						cur_target_max_stop = cur_target_len - 1;
						cur_target_delta = seq_len_increment(cur_target_len, 
							opt.fragment_target_threshold).first;

						cur_target_start = 0;
						cur_target_stop = cur_target_delta;
					}
					else{

						cur_target_start = cur_target_stop + 1;
						cur_target_stop = min(cur_target_stop + cur_target_delta, cur_target_max_stop);

						// If we get here, then we have fragmented the target sequence into 
						// two or more peices.
						fragment_target = true;
					}
				}
			}
			
			if(local_target >= num_seq){
			
				// Free memory before we leave
				if(bio_seq.second != NULL){

					// Free the memory used to store the sequence
					delete [] bio_seq.second;
					bio_seq.second = NULL;
				}

				break;
			}
						
			/////////////////////////////////////////////////////////////////////////////////////
			// Load the database sequence if not already loaded. There is no need to test both 
			// cur_target_start and cur_target_stop (i.e. cur_target_start will suffice).
			const bool same_target = (last_target == (long int)local_target) && 
				(last_target_start == local_target_start);
				
			if(!same_target){
				
				T_time = get_time();
				
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
				target_len = seq_file.read_bio_seq(bio_seq, local_target, 
					local_target_start, local_target_stop + max_product_length);
				
				if( target_len < dbase.min_sequence_size() ){

					// This sequence is too small to hash! Genbank contains some very small
					// DNA sequences (4-6 nucleotides). These sequences include ligands
					// bound to proteins in the PDB.
					if(bio_seq.second != NULL){

						// Free the memory used to store the sequence
						delete [] bio_seq.second;
						bio_seq.second = NULL;
					}
					
					// There is no stored sequence
					last_target = -1;
					
					continue;
				}
				
				// Hash the sequence for fast searching. Note that hash_dbase will
				// re-use previously allocated memory, so there is no need to deallocate
				// dbase after every sequence
				dbase.hash( bio_seq.second, SEQ_SIZE(bio_seq.second), 0, SEQ_SIZE(bio_seq.second) );
			
				last_target = local_target;
				last_target_start = local_target_start;
				
				// The time to load and hash a single target sequence
				T_time = get_time() - T_time;				
 			}
			
			// Are we searching with *all* queries, or just a single, specified query?
			const bool single_query = (local_query < num_sig);
			
			if(single_query == false){
				
				// Since we're searching will all queries, reset the query counter to 0
				local_query = 0;
			}
			
			tnt_time Q_time = get_time();
			
			while(true){
			
				list<hybrid_sig> local_results;
				const hybrid_sig &sig_ref = opt.sig_list[local_query];
				
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

				// Merge the results of this local calculation into the global list
				if(local_results.empty() == false){

					list< list<hybrid_sig>::iterator > reaper;
					
					// Before we add the local_results to the main results list, add the 
					// sequence id so we know which sequence produced this match
					list<hybrid_sig>::iterator local_iter;

					for(local_iter = local_results.begin();local_iter != local_results.end();local_iter++){
					
						// Is this match truncated due to target sequence fragmentation?
						if( (local_target_start != 0) && local_iter->start_overlap(0) ){
							
							reaper.push_back(local_iter);
							continue;
						}

						if( (local_target_stop != local_target_max_stop) && local_iter->stop_overlap(target_len - 1) ){
							
							reaper.push_back(local_iter);
							continue;
						}

						local_iter->seq_id(local_target);
						
						// Since we may be processing a subsequence, add the appropriate offset
						// so that match coordinates are correct
						local_iter->offset_ranges(local_target_start);
						
						// Compute the dimer and hairpin temperatures for this assay
						if(local_iter->has_primers() == true){

							melt.set_duplex(local_iter->forward_oligo);
							melt.strand(forward_primer_strand, forward_primer_strand);

							local_iter->forward_hairpin_tm = melt.approximate_tm_hairpin();
							local_iter->forward_dimer_tm = melt.approximate_tm_homodimer();

							melt.set_duplex(local_iter->reverse_oligo);
							melt.strand(reverse_primer_strand, reverse_primer_strand);

							local_iter->reverse_hairpin_tm = melt.approximate_tm_hairpin();
							local_iter->reverse_dimer_tm = melt.approximate_tm_homodimer();

							melt.set_query(local_iter->forward_oligo);
							melt.set_target(local_iter->reverse_oligo);

							melt.strand(forward_primer_strand, reverse_primer_strand);
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
				}
				
				// Update the local_query here so we can test for exit conditions in
				// the critical section below (for updating Q_time).
				local_query ++;
				
				// Use a omp critical section to avoid a race condition when modifying
				// the search_results vector.
				#pragma omp critical (SpliceResults)
				{
					// If the user has selected OUTPUT_INVERSE_QUERY for their output format,
					// we will only track whether each query has matched one or more targets
					if(opt.output_format & OUTPUT_INVERSE_QUERY){

						// Have we found one or more matches to this query?
						if(local_results.empty() == false){
							query_matches[sig_ref.my_id()] = true;
						}
					}
					else{
					
						list<hybrid_sig> &results_ref = search_results[sig_ref.my_id()];

						results_ref.splice(results_ref.begin(), local_results);
					}
					
					if(local_query >= num_sig){
						
						// The time required to search a single query
						Q_time = get_time() - Q_time;
						
						// Always normalize by num_sig (if we're already in single_query
						// mode we don't care about Q_time any more!).
						const float Q_tmp = Q_time.seconds()/num_sig;
						const float T_tmp = T_time.seconds();
						
						if(T_tmp > 0.0f){
						
							total_QT += Q_tmp/T_tmp;
							QT_count ++;
						}
					}
				}

				
				if(single_query || (local_query >= num_sig) ){
					break;
				}
			}
			
			#pragma omp master
			if(opt.verbose){

				const float search_status = (local_target*num_sig + local_query)*inv_total_num_of_comparisons;

				// Only update the search status if we've progressed at least 1% 
				// towards completion.
				if(search_status - last_search_status > search_display_every){

					stringstream ssout;

					ssout << setprecision(search_display_precision) << 100*search_status << '%';

					if(segment_queries){
						ssout << " [qs]" ;
					}

					int i;

					for(i = 0;i < int(update_buffer_size);i++){
						cout << '\b';
					}

					cout << ssout.str();

					const int num_space = (int)update_buffer_size - (int)ssout.str().size();

					for(i = 0;i < num_space;i++){
						cout << ' ';
					}

					// We need an explicit flush to make sure that MPI
					// updates the terminal
					cout.flush();

					last_search_status = search_status;

					if(search_status > 0.9f){

						search_display_every = 0.001f;
						search_display_precision = 4;

						if(search_status > 0.99f){

							search_display_every = 0.0001f;
							search_display_precision = 5;
						}
					}
				}
			}
		}
} // pragma omp parallel

		if(opt.verbose){
			
			stringstream ssout;

			ssout << setprecision(3) << 100 << '%';

			int i;

			for(i = 0;i < int(update_buffer_size);i++){
				cout << '\b';
			}

			cout << ssout.str();

			const int num_space = (int)update_buffer_size - (int)ssout.str().size();

			for(i = 0;i < num_space;i++){
				cout << ' ';
			}

			cout << endl;
		}
		
		if(opt.output_format & OUTPUT_INVERSE_QUERY){
			
			// Since individual assays may have been expanded into multiple assays (due to degenerate bases),
			// only output the list of query *names*
			set<string> query_set;
			set<string> match_set;
			
			// Write the name of all queries that did *not* match any targets
			for(unsigned int i = 0;i < num_sig;i++){
				
				query_set.insert( opt.sig_list[i].assay_string() );
				
				if(query_matches[i] == true){
					match_set.insert( opt.sig_list[i].assay_string() );
				}
			}
			
			set<string> diff;
			
			set_difference( query_set.begin(), query_set.end(), 
				match_set.begin(), match_set.end(), 
				insert_iterator< set<string> >( diff, diff.begin() ) );
			
			for(set<string>::const_iterator i = diff.begin();i != diff.end();i++){
				(*ptr_out) << *i << endl;
			}
		}
		
		// Track the number of unique targets "detected" by each signature
		set<string> total_unique_targets;
		vector<unsigned int> match_count(num_sig);
		
		// Format all of the results for output.
		for(vector< list<hybrid_sig> >::iterator result_iter = search_results.begin();result_iter != search_results.end();result_iter++){
		
			list<hybrid_sig> &tmp = *result_iter;
			
			if(tmp.empty() == true){
				continue;
			}
			
			// If enabled, only keep the best matches between a given query and a given target
			if(opt.best_match == true){
				select_best_match(tmp);
			}

			if(fragment_target == true){
			
				// Make the output unique (since different workers can be sent overlapping
				// target sequences) and save only the highest scoring exactly overlapping matches.
				uniquify_results(tmp);
			}
			
			// Sort the output by the lowest melting temperature
			tmp.sort();
			
			if(opt.one_output_file_per_query == true){
			
				const string filename = opt.output_filename + "." + tmp.front().name;
				
				if( (opt.output_format & OUTPUT_STANDARD) || 
				    (opt.output_format & OUTPUT_FASTA) ){
				    
					fout.close();

					fout.open( filename.c_str() );

					if(!fout){
					
						cerr << "Unable to open " << filename << endl;
						throw ":master: I/O error";
					}
				}
				
				if(opt.output_format & OUTPUT_NETWORK){
				    
				    	fout_sif.close();
				    
				    	const string filename_sif = filename + ".sif";
					
					fout_sif.open( filename_sif.c_str() );

					if(!fout_sif){
					
						cerr << "Unable to open " << filename_sif << endl;
						throw ":master: I/O error";
					}
				}
			}
			
			if(opt.output_format & OUTPUT_STANDARD){
			
				(*ptr_out) << "#####################################################################################" 
					<< endl;
			}
			
			set<string> unique_targets;
			
			list<hybrid_sig>::const_iterator iter;
			
			for(iter = tmp.begin();iter != tmp.end();iter++){
			
				if(opt.output_format & OUTPUT_STANDARD){
					(*ptr_out) << "name = " << iter->name << endl;
				}
				
				// What should we call the primers (if present)?
				string fp;
				string rp;
				
				if(iter->has_primers() == true){
				
					// Print the actual primers used (since this will catch the cases when
					// the forward or the reverse primer alone produces an amplicon)
						
					num_primer++;
					
					// What should we call the primers?
					fp = (opt.assay_format == ASSAY_PCR) ? "forward primer" : "5' probe";
					rp = (opt.assay_format == ASSAY_PCR) ? "reverse primer" : "3' probe";
					
					if(opt.output_format & OUTPUT_STANDARD){
					
						(*ptr_out) << fp << " = 5' " << iter->forward_oligo << " 3'" << endl;
						(*ptr_out) << rp << " = 5' " << iter->reverse_oligo << " 3'" << endl;
					}
					
					const float forward_dg = iter->forward_dH - opt.target_t*iter->forward_dS;
					const float reverse_dg = iter->reverse_dH - opt.target_t*iter->reverse_dS;

					
					// Track Tm bounds for final summary to the user
					max_primer_hairpin = max(max_primer_hairpin, iter->forward_hairpin_tm);
					max_primer_hairpin = max(max_primer_hairpin, iter->reverse_hairpin_tm);

					max_primer_homodimer = max(max_primer_homodimer, iter->forward_dimer_tm);
					max_primer_homodimer = max(max_primer_homodimer, iter->reverse_dimer_tm);

					max_primer_heterodimer = max(max_primer_heterodimer, iter->primer_dimer_tm);
					
					if(opt.output_format & OUTPUT_STANDARD){
					
						(*ptr_out) << fp << " tm = " << iter->forward_tm << endl;
						(*ptr_out) << rp << " tm = " << iter->reverse_tm << endl;

						(*ptr_out) << fp << " hairpin tm = " << iter->forward_hairpin_tm
							<< endl;
						(*ptr_out) << rp << " hairpin tm = " << iter->reverse_hairpin_tm
							<< endl;

						(*ptr_out) << fp << " homodimer tm = " << iter->forward_dimer_tm
							<< endl;
						(*ptr_out) << rp << " homodimer tm = " << iter->reverse_dimer_tm
							<< endl;

						(*ptr_out) << "heterodimer tm = " << iter->primer_dimer_tm
							<< endl;
						
						(*ptr_out) << fp << " dG[" << forward_dg << "] = "
							<< "dH[" << iter->forward_dH << "] - T*dS[" 
							<< iter->forward_dS << "]" << endl;	
						
						(*ptr_out) << rp << " dG[" << reverse_dg << "] = "
							<< "dH[" << iter->reverse_dH << "] - T*dS[" 
							<< iter->reverse_dS << "]" << endl;	
						
						(*ptr_out) << fp << " mismatches = " << iter->forward_mm << endl;
						(*ptr_out) << rp << " mismatches = " << iter->reverse_mm << endl;
						
						(*ptr_out) << fp << " gaps = " << iter->forward_gap << endl;
						(*ptr_out) << rp << " gaps = " << iter->reverse_gap << endl;
						
						if(opt.assay_format == ASSAY_PCR){
						
							(*ptr_out) << "min 3' clamp = " << iter->min_primer_clamp() << endl;
							(*ptr_out) << "max 3' clamp = " << iter->max_primer_clamp() << endl;
						}
						
						if(opt.assay_format == ASSAY_PADLOCK){
								
							(*ptr_out) << "5' probe 3' ligation clamp = " << iter->forward_primer_clamp 
								<< endl;
							(*ptr_out) << "3' probe 5' ligation clamp = " << iter->reverse_primer_clamp 
								<< endl;
						}
					}
					
					size_t len = iter->forward_oligo.size();
					
					forward_size_range.first = 
						min(forward_size_range.first, (unsigned int)len);
					forward_size_range.second = 
						max(forward_size_range.second, (unsigned int)len);
					
					len = iter->reverse_oligo.size();
					 
					reverse_size_range.first = 
						min(reverse_size_range.first, (unsigned int)len);
					reverse_size_range.second = 
						max(reverse_size_range.second, (unsigned int)len);
						
					forward_tm_range.first = 
						min(forward_tm_range.first, iter->forward_tm);
					forward_tm_range.second = 
						max(forward_tm_range.second, iter->forward_tm);
					
					forward_dg_range.first = 
						min(forward_dg_range.first, forward_dg);
					forward_dg_range.second = 
						max(forward_dg_range.second, forward_dg);
						
					reverse_tm_range.first = 
						min(reverse_tm_range.first, iter->reverse_tm);
					reverse_tm_range.second = 
						max(reverse_tm_range.second, iter->reverse_tm);
					
					reverse_dg_range.first = 
						min(reverse_dg_range.first, reverse_dg);
					reverse_dg_range.second = 
						max(reverse_dg_range.second, reverse_dg);
						
					float gc = 100.0f*gc_content(iter->forward_oligo);
					
					if(opt.output_format & OUTPUT_STANDARD){
					
						(*ptr_out) << fp << " %GC = "
							<< gc << endl;
					}
					
					forward_gc_range.first = 
						min(forward_gc_range.first, gc);
					forward_gc_range.second = 
						max(forward_gc_range.second, gc);
						
					gc = 100.0f*gc_content(iter->reverse_oligo);
					
					if(opt.output_format & OUTPUT_STANDARD){
					
						(*ptr_out) << rp << " %GC = "
							<< gc << endl;
					}
					
					reverse_gc_range.first = 
						min(reverse_gc_range.first, gc);
					reverse_gc_range.second = 
						max(reverse_gc_range.second, gc);
					
					if(opt.output_format & OUTPUT_STANDARD){
					
						(*ptr_out) << fp << " heuristics = " 
							<< primer_heuristics(iter->forward_oligo) << endl;
						(*ptr_out) << rp << " heuristics = " 
							<< primer_heuristics(iter->reverse_oligo) << endl;

						if(opt.assay_format == ASSAY_PCR){
						
							(*ptr_out) << "amplicon range = " << iter->amplicon_range.first 
								<< " .. " << iter->amplicon_range.second  << endl;
							(*ptr_out) << "amplicon length = " << iter->amplicon.size() << endl;
						}
						else{
							
							if(opt.assay_format == ASSAY_PADLOCK){
								(*ptr_out) << "product range = " << iter->amplicon_range.first 
									<< " .. " << iter->amplicon_range.second  << endl;
								(*ptr_out) << "product length = " << iter->amplicon.size() << endl;
							}
						}
						
						if(iter->primer_strand == hybrid_sig::PLUS){
							(*ptr_out) << "Forward primer is contained in the target plus strand (+)" << endl;
						}
						else{
							(*ptr_out) << "Forward primer is contained in the target minus strand (-)" << endl;
						}
					}
					
					amplicon_size_range.first = (unsigned int)min( (size_t)amplicon_size_range.first, 
						iter->amplicon.size() );
						 
					amplicon_size_range.second = (unsigned int)max( (size_t)amplicon_size_range.second, 
						iter->amplicon.size() );
				}
				
				if(iter->has_probe() == true){
				
					num_probe ++;
					
					const float gc = 100.0f*gc_content(iter->probe_oligo);
					
					probe_gc_range.first = 
						min(probe_gc_range.first, gc);
					probe_gc_range.second = 
						max(probe_gc_range.second, gc);					
					
					const float probe_dg = iter->probe_dH - opt.target_t*iter->probe_dS;
					
					max_probe_hairpin = max(max_probe_hairpin, iter->probe_hairpin_tm);
					max_probe_homodimer = max(max_probe_homodimer, iter->probe_dimer_tm);

					if(opt.output_format & OUTPUT_STANDARD){
					
						(*ptr_out) << "probe = 5' " << iter->probe_oligo << " 3'" << endl;
						(*ptr_out) << "probe tm = " << iter->probe_tm << endl;
						
						(*ptr_out) << "probe hairpin tm = " << iter->probe_hairpin_tm << endl;
						(*ptr_out) << "probe homodimer tm = " << iter->probe_dimer_tm << endl;
						
						(*ptr_out) << "probe dG[" << probe_dg << "] = "
							<< "dH[" << iter->probe_dH << "] - T*dS[" 
							<< iter->probe_dS << "]" << endl;	
						
						(*ptr_out) << "probe mismatches = " << iter->probe_mm << endl;
						(*ptr_out) << "probe gaps = " << iter->probe_gap << endl;
						
						(*ptr_out) << "probe %GC = " 
							<< gc << endl;

						(*ptr_out) << "probe range = " << iter->probe_range.first 
							<< " .. " << iter->probe_range.second  << endl;
						
						if(iter->probe_strand != iter->primer_strand){
							(*ptr_out) << "probe contained in forward strand (+)" << endl;
						}
						else{
							(*ptr_out) << "probe contained in reverse strand (-)" << endl;
						}
					}
					
					probe_tm_range.first = 
						min(probe_tm_range.first, iter->probe_tm);
					probe_tm_range.second = 
						max(probe_tm_range.second, iter->probe_tm);
						
					probe_dg_range.first = 
						min(probe_dg_range.first, probe_dg);
					probe_dg_range.second = 
						max(probe_dg_range.second, probe_dg);
						
					size_t len = iter->probe_oligo.size();
					
					probe_size_range.first = 
						min(probe_size_range.first, (unsigned int)len);
					probe_size_range.second = 
						max(probe_size_range.second, (unsigned int)len);
						
				}
				
				if(opt.output_format & OUTPUT_STANDARD){
				    
				    if(opt.output_format & OUTPUT_ALIGNMENTS){
						
						write_alignment(*ptr_out, fp + " align ", iter->forward_align);
						write_alignment(*ptr_out, rp + " align ", iter->reverse_align);
						write_alignment(*ptr_out, "probe align ", iter->probe_align);
					}
					
					if( seq_file.is_annot_format() ){
						write_annotation(*ptr_out, *iter, seq_file);
					}
				}
				
				if( (opt.output_format & OUTPUT_STANDARD) || 
				    (opt.output_format & OUTPUT_FASTA) ){
					
					(*ptr_out) << '>' << iter->amplicon_def;

					if(opt.append_name_to_defline){

						// The user has request that the assay name
						// be appended to the defline
						(*ptr_out) << ' ' << iter->name;
					}

					(*ptr_out) << endl;
					
					if(opt.output_format & OUTPUT_SEQ_MATCH){
						(*ptr_out) << iter->amplicon << endl;
					}
				}
				
				if(opt.output_format & OUTPUT_STANDARD){
					(*ptr_out) << endl;
				}
				
				if(opt.output_format & OUTPUT_NETWORK){
					
					fout_sif << mask_white_space(iter->name) << " binds " 
						    << mask_white_space(iter->amplicon_def) << endl;
				}
				
				// Track the number of unique targets (for the current signature)
				unique_targets.insert(iter->amplicon_def);
				
				// Track the number of unique targets
				total_unique_targets.insert(iter->amplicon_def);
			}
			
			match_count[result_iter - search_results.begin()] = unique_targets.size();
		}
		
		if(opt.output_format & OUTPUT_NETWORK){
			
			// All of the assays listed in opt.sig_list are parents
			vector<hybrid_sig>::const_iterator p_iter;

			for(p_iter = opt.sig_list.begin();p_iter != opt.sig_list.end();p_iter++){
				
				fout_atr << mask_white_space(p_iter->name) << " = parent" << endl;
			}
			
			set<string>::const_iterator c_iter;

			for(c_iter = total_unique_targets.begin();c_iter != total_unique_targets.end();c_iter++){

				fout_atr << mask_white_space(*c_iter) << " = child" << endl;
			}
		}
		
		if(opt.output_format & OUTPUT_INVERSE_TARGET){		
			
			// Write out the deflines of the target sequences that *don't* match any query!
			const unsigned int count = 
				write_inverse_matches(*ptr_out, seq_file, total_unique_targets);
			
			if(opt.verbose){
				cout << "Wrote " << count 
					<< " inverse target matches (that did not match any query!)" << endl;
			}
		}
		
		if(opt.verbose && !(opt.output_format & OUTPUT_INVERSE_QUERY) ){
			cout << "Found " << total_unique_targets.size() << " (total) target sequence matches" << endl;
		}
		
		if(opt.verbose && (num_primer > 0) && !(opt.output_format & OUTPUT_INVERSE_QUERY) ){
		
			cout << "Amplicon:" << endl
				<< "\t" << amplicon_size_range.first << " <= Amplicon length <= " 
				<< amplicon_size_range.second << endl;
				
			cout << "Forward primer:" << endl
				<< "\t" << forward_tm_range.first << " <= Tm (C) <= " << forward_tm_range.second  << endl
				<< "\t" << forward_dg_range.first << " <= Delta G (Kcal/Mol) <= " << forward_dg_range.second  << endl
				<< "\t" << forward_gc_range.first << " <= %GC <= " << forward_gc_range.second << endl
				<< "\t" << forward_size_range.first << " <= length <= " << forward_size_range.second << endl;
				
			cout << "Reverse primer:" << endl
				<< "\t" << reverse_tm_range.first << " <= Tm (C) <= " << reverse_tm_range.second << endl
				<< "\t" << reverse_dg_range.first << " <= Delta G (Kcal/Mol) <= " << reverse_dg_range.second << endl
				<< "\t" << reverse_gc_range.first << " <= %GC <= " << reverse_gc_range.second << endl
				<< "\t" << reverse_size_range.first << " <= length <= " << reverse_size_range.second << endl;
				
			cout << "Max primer hairpin Tm = " << max_primer_hairpin << endl;
			cout << "Max primer heterodimer Tm = " << max_primer_heterodimer << endl;
			cout << "Max primer homodimer Tm = " << max_primer_homodimer << endl;
		}
		
		if(opt.verbose && (num_probe > 0) && !(opt.output_format & OUTPUT_INVERSE_QUERY) ){
		
			cout << "Probe:" << endl
				<< "\t" << probe_tm_range.first << " <= Tm (C) <= " << probe_tm_range.second << endl
				<< "\t" << probe_dg_range.first << " <= Delta G (Kcal/Mol) <= " << probe_dg_range.second << endl
				<< "\t" << probe_gc_range.first << " <= %GC <= " << probe_gc_range.second << endl
				<< "\t" << probe_size_range.first << " <= length <= " << probe_size_range.second << endl;
				
			cout << "Max probe hairpin Tm = " << max_probe_hairpin << endl;
			cout << "Max probe homodimer Tm = " << max_probe_homodimer << endl;
		}
		
		if(opt.assay_summary && !(opt.output_format & OUTPUT_INVERSE_QUERY) ){
			
			cout << "*** Assay Summary ***" << endl;
			
			// Print the number of matches found for each assay
			for(vector<hybrid_sig>::const_iterator i = opt.sig_list.begin();i != opt.sig_list.end();i++){
				
				cout  << i->name << " matched " << match_count[ i->my_id() ] << " sequences" << endl;
				
				if( (i->forward_oligo != "") && (i->reverse_oligo != "") ){
					cout << "\tF::R = " << i->forward_oligo << " :: " 
						<< i->reverse_oligo <<  endl;
				}
								
				if(i->probe_oligo != ""){
					cout << "\tP = " << i->probe_oligo << endl;
				}
			}
		}

		profile = time(NULL) - profile;
		
		if(opt.verbose){
			cout << "Search completed in " << profile << " sec" << endl;
		}
		
		// All done!
		
	}
	catch(const char *error){
	
		cerr << "Caught the error: " << error << endl;
		return EXIT_FAILURE;
	}
	catch(std::exception &e){
		
		cerr << "Caught the std exception: " << e.what() << endl;
		return EXIT_FAILURE;
	} 
	catch(...){
	
		cerr << "Caught an unhandled error. Please report to " << EMAIL_ADDRESS << endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
