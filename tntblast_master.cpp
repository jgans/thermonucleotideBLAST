#ifdef USE_MPI

#include <mpi.h>
#include "tntblast.h"
#include "options.h"
#include "hybrid_sig.h"
#include "degenerate_na.h"
#include "primer.h"
#include "mpi_util.h"
#include "bitmask.h"

#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>

using namespace std;

// Global variables
extern int mpi_numtasks;
extern int mpi_rank;

int master(int argc, char *argv[])
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
			
			int continue_exec = false;
			MPI_Bcast(&continue_exec, 1, MPI_INT, 0, MPI_COMM_WORLD);
			
			return EXIT_FAILURE;
		}
		catch(...){
		
			cerr << "Unhandled input error, please report to " 
				<< EMAIL_ADDRESS << endl;
			
			int continue_exec = false;
			MPI_Bcast(&continue_exec, 1, MPI_INT, 0, MPI_COMM_WORLD);
			
			return EXIT_FAILURE;
		}

		// Is the user printing the command line arguments and then exiting?
		if(opt.print_usage){
			
			int continue_exec = false;
			MPI_Bcast(&continue_exec, 1, MPI_INT, 0, MPI_COMM_WORLD);
			
			return EXIT_SUCCESS;
		}
		else{
		
			int continue_exec = true;
			MPI_Bcast(&continue_exec, 1, MPI_INT, 0, MPI_COMM_WORLD);
		}
		
		if(opt.input_filename != ""){
		
			if(opt.verbose){
				cout << "Reading assays from " << opt.input_filename << endl;
			}
			
			// Note that using the assay format ASSAY_PROBE forces all oligos to be 
			// treated as probes!
			read_input_file( opt.input_filename, opt.sig_list, 
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
						cerr << "Unable to open " << opt.output_filename << endl;
						throw "Unable to open output file";
					}
				}
				
				if(opt.output_format & OUTPUT_NETWORK){
				    
				    	const string filename_sif =  opt.output_filename + ".sif";
				    
					fout_sif.open( filename_sif.c_str() );

					if(!fout_sif){
					
						cerr << "Unable to open " << filename_sif << endl;
						throw ":master: I/O error";
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
					throw ":master: I/O error";
				}

				// Write the attribute header
				fout_atr << "FunctionalCatagory" << endl;
			}
			
			if( (opt.output_format & OUTPUT_INVERSE_TARGET) ||
			    (opt.output_format & OUTPUT_INVERSE_QUERY) ){
			
				fout.open( opt.output_filename.c_str() );

				if(!fout){
					cerr << "Unable to open " << opt.output_filename << endl;
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
		
		const unsigned long int num_sig = opt.sig_list.size();
		
		if(num_sig == 0){
			throw __FILE__ ":master: No assay oligos found!";
		}
		
		// Count the number of probe only queries in the list of queries
		const unsigned int num_probes = probe_only_count(opt.sig_list);
		
		const unsigned int max_product_length = opt.max_product_length() + 2; // Room for dangling end bases
				
		const unsigned int num_worker = (unsigned int)(mpi_numtasks - 1);
		
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
		const unsigned long int num_seq = seq_file.size();
		
		if(num_seq == 0){
			throw __FILE__ ":master: Empty database -- no sequences found!";
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
		
		#ifdef PROFILE
		// Track the total time required to perform the search
		double profile = MPI_Wtime();
		#endif // PROFILE
		
		// Broadcast the command line options to the worker nodes
		broadcast(opt, mpi_rank, 0);
		
		// Send the database file to the worker nodes
		if(seq_file.is_annot_format() == true){
		
			// Force the worker nodes to read sequence data from the master
			// (not directly from the disk). This will avoid the overhead of
			// workers parsing and storing annotation information.
			send("");
		}
		else{
		
			send(opt.dbase_filename);
		}
		
		int sequence_file_format = seq_file.file_format();
		
		MPI_Bcast(&sequence_file_format, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		// Send the format of the data base file to the workers. If the format
		// is sequence_data::FASTA_SLOW, then the master will send the fasta
		// record indicies to the workers (so we don't have to parse the fasta
		// file multiple times)
				
		// Distribute the signature queries
		distribute_queries(opt.sig_list);
		
		// Do we need to send indicies to the workers?
		if( (opt.dbase_filename != "") && seq_file.wants_indicies(sequence_file_format) ){
			
			cout << "Broadcasting fasta indicies to the workers" << endl;
			send( seq_file.indicies() );
		}
		
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
				
		list<int> idle;
		
		// All nodes start off as idle
		for(int i = 1;i < mpi_numtasks;i++){
			idle.push_back(i);
		}
				
		const double inv_total_num_of_comparisons = 1.0/( (double)(num_seq)*(double)(num_sig) );
		unsigned long int num_of_comparisons = 0;
		
		float last_search_status = 0.0f;
		const unsigned int update_buffer_size = 15;
		float search_display_every = 0.01f;
		unsigned int search_display_precision = 3;
		
		if(opt.verbose){
			
			cout << "Searching database: ";
			
			for(unsigned int i = 0;i < update_buffer_size;i++){
				cout << ' ';
			}
			
			// We need an explicit flush to make sure that MPI
			// updates the terminal
			cout.flush();
		}
		
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
		
		// Parallelization strategy: 
		// 1) Always segment the database targets
		// 2) Segment the assay queries if the number of
		//    remaining targets is *less* than the number of
		//    workers (times a constant greater than 1).
		// 3) Fragment the target sequence into sub-sequences if the target sequence
		// 	 length is larger than opt.fragment_target_threshold.
		
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
		
		// When should we segmenting queries?
		// Too soon and we incur a high overhead penalty. Too late and we incur
		// a load balancing penalty.
		bool segment_queries = query_sched(effective_num_seq, num_sig, num_worker, 
			(QT_count == 0 ? default_qt : total_QT/QT_count), opt.query_segmentation);
			
		const long unsigned int delta_query = max(1UL, num_sig/num_worker);
		
		unsigned int cur_query = 0;
		
		// Scratch variables to store results from the workers
		MPI_Status status;
		
		// ret_buffer is a flexible data structure that can store either
		// 1) Requests for a sequence from a worker (an unsigned int)
		// 2) Timing information for adaptive query segmentation (a float)
		union ref_buffer{
			int int_value;
			unsigned int uint_value;
			float float_value;
		} ret;
			
		#ifdef PROFILE
		// The time spent searching
		double search_time = MPI_Wtime();
		
		// The time spent collecting data from the workers and outputting the results
		double output_time = 0.0;
		#endif // PROFILE
		
		// Keep looping until we have searched everything!
		while(true){
		
			if( (cur_target >= num_seq) && (idle.size() == num_worker) ){
				
				// We're done!
				break;
			}
				
			// Dispatch work to any idle nodes
			while( (idle.empty() == false) && (cur_target < num_seq) ){
			
				
				const int node = idle.front();
				idle.pop_front();
				
				const unsigned int num_query = segment_queries ? min( delta_query, (num_sig - cur_query) ) :
					num_sig;
								
				unsigned int buffer[SEARCH_QUERY_BUFFER_SIZE];
				
				// The query assay to search
				buffer[0] = cur_query;
				
				// The number of queries to search
				buffer[1] = num_query;
				
				// The database sequence to search.
				buffer[2] = cur_target;
				
				buffer[3] = cur_target_start;
				
				// Add an overlap of max_product length so that fragmenting targets does
				// not remove potential matches.
				buffer[4] = min(cur_target_stop + max_product_length, cur_target_max_stop);
				
				buffer[5] = cur_target_max_stop;
								
				if(MPI_Send(buffer, SEARCH_QUERY_BUFFER_SIZE, MPI_UNSIGNED, node, 
				   SEARCH_QUERY, MPI_COMM_WORLD) != MPI_SUCCESS){

				   throw __FILE__ ":master: Error sending SEARCH_QUERY";
				}
				
				// DEBUG
				//cerr << "Sending target " << buffer[2] << " (index " << cur_target << ") and query " << buffer[0] 
				//	<< " to [" << node << "] (segment query = " 
				//	<< (segment_queries ? "true" : "false") << ")" << endl;
					
				// In order to minimize network load, the assay queries are the 
				// inner-loop, while the database targets are the outer-loop
				cur_query += num_query;
				
				if(cur_query == num_sig){

					cur_query = 0;
					
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
						
						num_of_comparisons += num_sig;
					}
					else{
						
						cur_target_start = cur_target_stop + 1;
						cur_target_stop = min(cur_target_stop + cur_target_delta, cur_target_max_stop);
						
						// If we get here, then we have fragmented the target sequence into 
						// two or more peices.
						fragment_target = true;
					}
				}

				if(segment_queries == false){
					
					// Update the query segmentation status (that is, do we need to
					// start segmenting queries?)
					segment_queries = query_sched(effective_num_seq, num_sig, num_worker,
						(QT_count == 0 ? default_qt : total_QT/QT_count), 
						opt.query_segmentation);
				}
			}

			// Every message below this point has an unsigned int payload
			if(MPI_Recv(&ret, sizeof(ret), MPI_BYTE, MPI_ANY_SOURCE, 
				MPI_ANY_TAG, MPI_COMM_WORLD, &status) != MPI_SUCCESS){

				throw __FILE__ ":master: Error receiving worker info";
			}
						
			if(status.MPI_TAG == STATUS_UPDATE){
				
				total_QT += ret.float_value;
				QT_count ++;
				
				if(opt.verbose){
				
					const float search_status = num_of_comparisons*inv_total_num_of_comparisons;

					// Only update the search status if we've progressed at least 1% 
					// towards completion.
					if(search_status - last_search_status > search_display_every){

						stringstream ssout;

						ssout << setprecision(search_display_precision) << 100*search_status << '%';
						
						if(segment_queries){
							ssout << " [qs]" ;
						}
						
						for(unsigned int i = 0;i < update_buffer_size;i++){
							cout << '\b';
						}

						cout << ssout.str();

						const int num_space = (int)update_buffer_size - (int)ssout.str().size();

						for(int i = 0;i < num_space;i++){
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
				
				// This worker is ready for more work
				idle.push_back(status.MPI_SOURCE);
				
				continue;
			}
			
			if(status.MPI_TAG != SEQ_REQUEST){
				throw __FILE__ ":master: Unexpected message";
			}
						
			// Send the sequence "ret" to worker "status.MPI_SOURCE"
			serve_sequence(status.MPI_SOURCE, ret.uint_value, seq_file);
		}

		if(opt.verbose){
						
			// Start a new line
			for(unsigned int i = 0;i < update_buffer_size;i++){
				cout << '\b';
			}
		
			cout << "search complete" << endl;
		}
		
		// Tell the workers to start sending back their results
		int msg = SEARCH_COMPLETE;
		
		for(int j = 1;j < mpi_numtasks;++j){

			// I have observed LAM MPI hanging with the error:
			// MPI_Ssend: internal MPI error: Bad address (rank 0, MPI_COMM_WORLD)
			// when MPI_Ssend is used. This may indicate:
			// a) memory issue (nothing found using valgrind)
			// b) conflict with NCBI tool kit
			// The problem only occurs for certain queries and while using 
			// -ssi rpi lamd
			if(MPI_Send(&msg, 1, MPI_INT, j, 
			   SEARCH_COMPLETE, MPI_COMM_WORLD) != MPI_SUCCESS){

			   throw __FILE__ ":master: Error sending SEARCH_COMPLETE";
			}
		}
		
		#ifdef PROFILE
		search_time = MPI_Wtime() - search_time;

		output_time = MPI_Wtime();
		#endif // PROFILE
		
		if(opt.verbose){
		
			cout << "Collecting results: ";
			
			for(unsigned int i = 0;i < update_buffer_size;i++){
				cout << ' ';
			}
			
			cout.flush();
		}
		
		#ifdef PROFILE
		// How much time are we spending doing work, and how much time are we spending
		// on communication.
		double work_time = 0.0;			// Time spent searching
		double seq_load_time = 0.0;		// Time spent loading sequence data
		double seq_hash_time = 0.0;		// Time spent hashing sequence data
		double communication_time = 0.0; 	// Time spent on communication
		double idle_time = 0.0; 		// Time spent idle
		double num_tm_eval_plus = 0.0; 		// Number of plus strand tm evaluations
		double num_tm_eval_minus = 0.0; 	// Number of minus strand tm evaluations
		
		// Collect all of the profile information from the workers
		for(unsigned int worker = 1;worker <= num_worker;worker++){
			
			double buffer[NUM_PROFILE];
			
			if(MPI_Recv(buffer, NUM_PROFILE, MPI_DOUBLE, worker, 
				PROFILE_INFO, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS){

				throw __FILE__ ":master: Error receiving buffer";
			}
			
			work_time += buffer[PROFILE_WORK];
			seq_load_time += buffer[PROFILE_LOAD];
			seq_hash_time += buffer[PROFILE_HASH];
			communication_time += buffer[PROFILE_COMM];
			idle_time += buffer[PROFILE_IDLE];
			num_tm_eval_plus += buffer[PROFILE_NUM_PLUS_TM_EVAL];
			num_tm_eval_minus += buffer[PROFILE_NUM_MINUS_TM_EVAL];
		}
		
		#endif // PROFILE
		
		unsigned int query_chunck_size = 0;
		int num_chunck = 0;
		
		// If our output format is OUTPUT_INVERSE_QUERY, we do not need to transmit the 
		// results in chunks.
		if(opt.output_format & OUTPUT_INVERSE_QUERY){
		
			bitmask query_matches(num_sig, false);
			bitmask tmp_matches(num_sig, false);
			
			const unsigned int buffer_size = query_matches.mpi_size();
			unsigned char* buffer = new unsigned char [buffer_size];
			
			if(buffer == NULL){
				throw __FILE__ ":master: Unable to allocate bitmask buffer";
			}

			// Download the results from each worker. 
			for(int worker = 1;worker <= int(num_worker);worker++){
				
				if(MPI_Recv(buffer, buffer_size, MPI_BYTE, worker, 
					SIGNATURE_RESULTS, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS){

					throw __FILE__ ":master: Error receiving bitmask buffer";
				}
				
				tmp_matches.mpi_unpack(buffer);
				
				// Compute the set union
				query_matches += tmp_matches;
			}
			
			delete [] buffer;
			
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
		else{
			
			// Retrieve data from the workers in blocks of query_chunck_size queries.
			// This is to reduce the memory constraints on the master (which needs to
			// sort the results for output).
			// Tell the workers how many query results to send back in a single chunk
			query_chunck_size = min(1000UL, max(1UL, num_sig/10) );
			num_chunck = num_sig/query_chunck_size + (num_sig%query_chunck_size != 0);
			
			MPI_Bcast(&query_chunck_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
		}
		
		// Track the number of unique targets "detected" by each signature
		set<string> total_unique_targets;
		vector<unsigned int> match_count(num_sig);
		
		last_search_status = 0.0f;
		
		// Collect the results from the workers for each query
		for(int sig_index = 0;sig_index < num_chunck;sig_index++){
		
			if(opt.verbose){
			
				const float search_status = sig_index/float(num_chunck);

				if(search_status - last_search_status > 0.01f){

					stringstream ssout;

					ssout << setprecision(3) << 100*search_status << '%';

					for(unsigned int i = 0;i < update_buffer_size;i++){
						cout << '\b';
					}

					cout << ssout.str();

					const int num_space = (int)update_buffer_size - (int)ssout.str().size();

					for(int i = 0;i < num_space;i++){
						cout << ' ';
					}

					// We need an explicit flush to make sure that MPI
					// updates the terminal
					cout.flush();

					last_search_status = search_status;
				}
			}
			
			list<hybrid_sig> search_results;
			
			for(int worker = 1;worker <= int(num_worker);worker++){
				
				MPI_Status status;

				MPI_Probe(worker, SIGNATURE_RESULTS, MPI_COMM_WORLD, &status);

				int buffer_size;

				MPI_Get_count(&status, MPI_BYTE, &buffer_size);

				if(buffer_size == 0){
					throw __FILE__ ":master: Error buffer size is 0";
				}

				unsigned char* buffer = new unsigned char [buffer_size];

				if(buffer == NULL){
					throw __FILE__ ":master: Error allocating receive buffer";
				}

				if(MPI_Recv(buffer, buffer_size, MPI_BYTE, status.MPI_SOURCE, 
					SIGNATURE_RESULTS, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS){

					throw __FILE__ ":master: Error receiving buffer";
				}

				unsigned char *ptr = buffer;

				// Unpack the number of results
				unsigned int num_results;

				memcpy(&num_results, ptr, sizeof(unsigned int) );
				ptr += sizeof(unsigned int);
				
				for(unsigned int i = 0;i < num_results;i++){

					search_results.push_back( hybrid_sig() );
					//ptr = search_results.back().mpi_unpack(ptr);
					ptr = mpi_unpack( ptr, search_results.back() );
				}

				delete [] buffer;
			}

			// Format all of the results for output.
			if(search_results.empty() == true){
				continue;
			}
			
			// If enabled, only keep the best matches between a given query and a given target
			if(opt.best_match == true){
				select_best_match(search_results);
			}

			// We only need to make the results unique if we have fragmented a target sequence
			if(fragment_target == true){
			
				// Make the output unique (since different workers can be sent overlapping
				// target sequences) and save only the highest scoring exactly overlapping matches.
				uniquify_results(search_results);
			}
			
			// Sort the output by the lowest melting temperature
			search_results.sort();
			
			set<string> unique_targets;

			list<hybrid_sig>::const_iterator iter;

			int last_id = -1;
			
			for(iter = search_results.begin();iter != search_results.end();iter++){

				if( last_id != iter->my_id() ){
				
					if(opt.one_output_file_per_query == true){

						const string filename = opt.output_filename + "." + search_results.front().name;

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
				}
				
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

					// Track Tm and delta G bounds for final summary to the user
					max_primer_hairpin = max(max_primer_hairpin, iter->forward_hairpin_tm);
					max_primer_hairpin = max(max_primer_hairpin, iter->reverse_hairpin_tm);

					max_primer_homodimer = max(max_primer_homodimer, iter->forward_dimer_tm);
					max_primer_homodimer = max(max_primer_homodimer, iter->reverse_dimer_tm);

					max_primer_heterodimer = max(max_primer_heterodimer, iter->primer_dimer_tm);

					const float forward_dg = iter->forward_dH - opt.target_t*iter->forward_dS;
					const float reverse_dg = iter->reverse_dH - opt.target_t*iter->reverse_dS;

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

					unsigned int len = iter->forward_oligo.size();

					forward_size_range.first = 
						min(forward_size_range.first, len);
					forward_size_range.second = 
						max(forward_size_range.second, len);

					len = iter->reverse_oligo.size();

					reverse_size_range.first = 
						min(reverse_size_range.first, len);
					reverse_size_range.second = 
						max(reverse_size_range.second, len);

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

					amplicon_size_range.first = min( (size_t)amplicon_size_range.first, 
						iter->amplicon.size() );

					amplicon_size_range.second = max( (size_t)amplicon_size_range.second, 
						iter->amplicon.size() );
				}

				if(iter->has_probe() == true){

					num_probe ++;

					const float gc = 100.0f*gc_content(iter->probe_oligo);

					probe_gc_range.first = 
						min(probe_gc_range.first, gc);
					probe_gc_range.second = 
						max(probe_gc_range.second, gc);

					max_probe_hairpin = max(max_probe_hairpin, iter->probe_hairpin_tm);
					max_probe_homodimer = max(max_probe_homodimer, iter->probe_dimer_tm);

					const float probe_dg = iter->probe_dH - opt.target_t*iter->probe_dS;

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

					unsigned int len = iter->probe_oligo.size();

					probe_size_range.first = 
						min(probe_size_range.first, len);
					probe_size_range.second = 
						max(probe_size_range.second, len);

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
				
				if( last_id != iter->my_id() ){
				
					if(last_id >= 0){
			
						// Save the number of matches for this query     
						match_count[last_id] = unique_targets.size();
						
						unique_targets.clear();
					}
				}
				
				// Track the number of unique targets (for the current signature)
				unique_targets.insert(iter->amplicon_def);
				
				// Track the number of unique targets
				total_unique_targets.insert(iter->amplicon_def);
				
				last_id = iter->my_id();
			}
			
			if(last_id >= 0){
			
				// Save the number of matches for this query     
				match_count[last_id] = unique_targets.size();
			}
		}
		
		if(opt.verbose){
						
			// Start a new line
			for(unsigned int i = 0;i < update_buffer_size;i++){
				cout << '\b';
			}
		
			cout << "collection complete" << endl;
		}
			
		if(opt.output_format & OUTPUT_NETWORK){

			// All of the assays listed in opt.sig_list are parents
			for(vector<hybrid_sig>::const_iterator i = opt.sig_list.begin();i != opt.sig_list.end();i++){

				fout_atr << mask_white_space(i->name) << " = parent" << endl;
			}

			for(set<string>::const_iterator i = total_unique_targets.begin();
				i != total_unique_targets.end();i++){

				fout_atr << mask_white_space(*i) << " = child" << endl;
			}
		}

		if(opt.output_format & OUTPUT_INVERSE_TARGET){

			// Write out the deflines of the sequences that *don't* match any query!
			const unsigned int count = 
				write_inverse_matches(*ptr_out, seq_file, total_unique_targets);

			if(opt.verbose){
				cout << "Wrote " << count 
					<< " inverse target matches (that did not match any query!)" << endl;
			}
		}

		if(opt.verbose && !(opt.output_format & OUTPUT_INVERSE_QUERY) ){
			cout << "Found " << total_unique_targets.size() << " (total) target sequence  matches" << endl;
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
				
				cout  << i->name << " matched " << match_count[ i->my_id() ] 
					<< " sequences" << endl;
				
				if( (i->forward_oligo != "") && (i->reverse_oligo != "") ){
					cout << "\tF::R = " << i->forward_oligo << " :: " 
						<< i->reverse_oligo <<  endl;
				}
								
				if(i->probe_oligo != ""){
					cout << "\tP = " << i->probe_oligo << endl;
				}
			}
		}
		
		#ifdef PROFILE
		profile = MPI_Wtime() - profile;
		output_time = MPI_Wtime() - output_time;
		#endif // PROFILE
		
		if(opt.verbose){
			
			#ifdef PROFILE
			cout << "Completed in " << profile << " sec" << endl;
			
			cout << "Time spent searching was " << search_time << " sec" << endl;
			
			double time_norm = 100.0/(work_time + seq_load_time + 
				seq_hash_time + communication_time + idle_time);
			
			cout << "\t" << time_norm*work_time 
				<< "% of time spent comparing" << endl;
			cout << "\t" << time_norm*seq_load_time 
				<< "% of time spent loading database sequences" << endl;
				cout << "\t" << time_norm*seq_hash_time 
				<< "% of time spent hashing database sequences" << endl;
			cout << "\t" << time_norm*communication_time 
				<< "% of time spent in interprocess communication" << endl;
			cout << "\t" << time_norm*idle_time 
				<< "% of time spent idle" << endl;
			
			cout << "Time required for output was " << output_time << " sec" << endl;
			
			cout << "Performed " << int(num_tm_eval_plus) << " (+) strand and " 
				<< int(num_tm_eval_minus) << " (-) strand Tm eval; " 
				<< int(num_tm_eval_plus + num_tm_eval_minus) << " total" << endl;
				
			// Uncomment this code to output the average ratio of query search time to 
			// target sequence load and hash time. This value (computed for relevant test
			// sets with query segmentation always on) is used to set DEFAULT_QT when 
			// algorithms and hardware change significantly.
			if(QT_count != 0){
			
				// To estimate the cost of query seqmentation we need the ratio of query search 
				// time to target sequence load and hash time. Scale the defined value according to
				// the assay format:
				// 1XDEFAULT_QT	for PROBE	query targeting a single strand
				// 2XDEFAULT_QT	for PROBE queries targeting both strands
				// 4XDEFAULT_QT	for PCR queries
				// 4XDEFAULT_QT	for PADLOCK queries
				const float benchmark_qt = num_sig*(total_QT/QT_count)/
					(num_probes*( (opt.target_strand == Seq_strand_both) ?  2.0f: 1.0f) +
					(num_sig - num_probes)*4.0f);

				
				//cout << "scaled <query search time to target sequence load and hash time> = " 
				//	<< benchmark_qt << endl;
			}
			#endif // PROFILE
		}
		
		// All done!
		
	}
	catch(const char *error){
	
		cerr << "Caught the master error: " << error << endl;
		return EXIT_FAILURE;
	}
	catch(...){
	
		cerr << "Caught an unhandled master error. Please report to " << EMAIL_ADDRESS << endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

#endif // USE_MPI

