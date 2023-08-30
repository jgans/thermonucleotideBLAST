#include "options.h"

#ifdef WIN32
extern "C" {
#endif // WIN32

	#include "getopt.h"

#ifdef WIN32
}
#endif // WIN32
	
#include <sstream>

using namespace std;

void Options::parse_command_line(int argc, char *argv[])
{
	// Command line options:
	// -i <input file of query oligos>
	// -o <output file>
	// -d <database of target sequences to search against>
	// -D <local database of target sequences to search against>
	// -l <maximum amplicon length>
	// -e <min primer Tm>
	// -E <min probe Tm>
	// -z <min primer delta G>
	// -Z <min probe delta G>
	// -x <max primer Tm>
	// -X <max probe Tm>
	// -g <max primer delta G>
	// -G <max probe delta G>
	// -s <salt concentration (in MOL)>
	// -t <primer strand concentration (in MOL)>
	// -T <Probe strand concentration (in MOL)>
	// -y <ratio of forward/reverse strand concentrations>
	// -A <assay format>
	// -W <hash word length>
	// -m <output format>
	// -a <T|F> (show alignments)
	// -M <T|F> (show matching sequence)
	// -k <T|F> (Mask primer binding sites)
	// -K <T|F> (Mask probe binding sites)
	// -r <T|F> (Replace primer binding sites w/ primer sequence)
	// -v <T|F> (Disable verbose output)
	// -p <T|F> (Ignore probe oligos in inputfile)
	// -n <T|F> (One output file per query)
	// -L <T|F> (Append assay name to output defline)
	// -S <T|F> (Print assay summary after search)
	// -? help
	// -h help
	// --help
	// --primer-clamp <number of exact 3' primer matches requried>
	// --probe-clamp5 <number of exact 5' probe matches requried>
	// --probe-clamp3 <number of exact 3' probe matches requried>
	// --dangle5 <T|F> (Allow dangling bases on the 5' query side of an alignment)
	// --dangle3 <T|F> (Allow dangling bases on the 3' query side of an alignment)
	// --target-strand <sense|antisense|both>
	// --plex <T|F> (All input assays in a single multiple reaction)
	// --single-primer-pcr <T|F> (Allow amplicons produced by a single PCR primer
	// binding in both forward and reverse orientation?)
	// --temperature <temperature for computing Delta G (in Kelvin)>
	// --hash-size <max size> (max number of hash elements)
	// --max-target-len <max len> (max sequence length, in bases, before targets are split)
	// --query-seg <always | never | adaptive> (query segmentation algorithm)
	// --dump-query <T|F> (write queries to stdout)
	// --dinkelbach <T|F> (Use the Dinkelbach fractional programming algorithm)
	// --max-gap <number of gaps> (maximum number of gaps allowed in a DNA duplex)
	// --max-mismatch <number of mismatches> (maximum number of mismatches allowed in a DNA duplex)
	// --rescale-ct <T|F> (Use of degenerate bases results in rescaling of oligo concentration)
	// --best-match (Only report the best match, in Tm, between an assay and a target)
	// --blast-include <Accession or NCBI TaxId to include from BLAST database>
	// --blast-exclude <Accession or NCBI TaxId to exclude from BLAST database>

	const char* options = "-i:o:d:D:l:e:E:z:Z:x:X:g:G:s:t:T:y:A:W:m:a:M:k:K:r:v:p:n:L:S:?h";
	int config_opt = 0;
	int long_index = 0;

	struct option long_opts[] = {
		{"help", false, &config_opt, 1},
		{"primer-clamp", true, &config_opt, 3},
		{"probe-clamp5", true, &config_opt, 4},
		{"probe-clamp3", true, &config_opt, 5},
		{"plex", true, &config_opt, 6},
		{"single-primer-pcr", true, &config_opt, 7},
		{"hash-size", true, &config_opt, 8},
		{"target-strand", true, &config_opt, 9},
		{"temperature", true, &config_opt, 10},
		{"max-target-len", true, &config_opt, 11},
		{"query-seg", true, &config_opt, 12},
		{"dump-query", true, &config_opt, 13},
		{"dangle5", true, &config_opt, 14},
		{"dangle3", true, &config_opt, 15},
		{"min-max-primer-clamp", true, &config_opt, 16},
		{"dinkelbach", true, &config_opt, 17},
		{"max-gap", true, &config_opt, 18},
		{"max-mismatch", true, &config_opt, 19},
		{"rescale-ct", true, &config_opt, 20},
		{"best-match", false, &config_opt, 21},
		{"blast-include", true, &config_opt, 22},
		{"blast-exclude", true, &config_opt, 23},
		{0,0,0,0} // Terminate options list
	};

	int opt_code;
	opterr = 0;
	
	threshold_format = THRESHOLD_NONE;
	
	// If no command line arguments are specified, then print the command line usage
	print_usage = (argc == 1);

	// Read the command line options
	while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){
	
		switch( opt_code ){
			case 0:
				// --help
				if(config_opt == 1){
					print_usage = true;
					break;
				}
				
				// --primer-clamp
				if(config_opt == 3){
					primer_clamp = atoi(optarg);
					break;
				}
				
				// --probe-clamp5
				if(config_opt == 4){
					probe_clamp_5 = atoi(optarg);
					break;
				}
				
				// --probe-clamp3
				if(config_opt == 5){
					probe_clamp_3 = atoi(optarg);
					break;
				}
				
				// --plex
				if(config_opt == 6){
					multiplex = parse_bool(optarg);
					break;
				}
				
				// --single-primer-pcr
				if(config_opt == 7){
					single_primer_pcr = parse_bool(optarg);
					break;
				}
				
				// --target-strand
				if(config_opt == 9){
				
					target_strand = parse_strand(optarg);
					break;
				}
				
				// --temperature
				if(config_opt == 10){
				
					target_t = atof(optarg);
					
					if(target_t < 0.0f){
						cerr << "Warning: --temperature is less than zero!" << endl;
					}
					
					break;
				}
				
				// --max-target-len
				if(config_opt == 11){
				
					fragment_target_threshold = atoi(optarg);
					
					if(fragment_target_threshold <= 1){
						throw "Error: --max-target-len is <= 1";
					}
					
					break;
				}
				
				// --query-seg
				if(config_opt == 12){
				
					query_segmentation = parse_query_seg(optarg);
					break;
				}
				
				// --dump-query
				if(config_opt == 13){
				
					dump_query = parse_bool(optarg);
					break;
				}
				
				// --dangle5
				if(config_opt == 14){
				
					allow_dangle_5 = parse_bool(optarg);
					break;
				}
				
				// --dangle3
				if(config_opt == 15){
				
					allow_dangle_3 = parse_bool(optarg);
					break;
				}
				
				// --min-max-primer-clamp
				if(config_opt == 16){
					min_max_primer_clamp = atoi(optarg);
					break;
				}
				
				// --dinkelbach
				if(config_opt == 17){
					use_dinkelbach = parse_bool(optarg);
					break;
				}
				
				// --max-gap
				if(config_opt == 18){
				
					max_gap = atoi(optarg);
					break;
				}
				
				// --max-mismatch
				if(config_opt == 19){
				
					max_mismatch = atoi(optarg);
					break;
				}
				
				// --rescale-ct
				if(config_opt == 20){
				
					degen_rescale_ct = parse_bool(optarg);
					break;
				}

				// --best-match
				if(config_opt == 21){
				
					best_match = true;
					break;
				}

				// --blast-include
				if(config_opt == 22){
				
					blast_include.push_back(optarg);
					break;
				}

				// --blast-exclude
				if(config_opt == 23){
				
					blast_exclude.push_back(optarg);
					break;
				}

				cerr << "Unknown flag!" << endl;
				break;
			case 'i':
				input_filename = optarg;
				break;
			case 'o':
				output_filename = optarg;
				break;
			case 'd':
				dbase_filename = optarg;
				break;
			case 'D':
				local_dbase_filename = optarg;
				break;
			case 'l':
				max_len = atoi(optarg);
				break;
			case 'e':
				min_primer_tm = float( atof(optarg) );
				threshold_format |= THRESHOLD_PRIMER_TM;
				break;
			case 'E':
				min_probe_tm = float( atof(optarg) );
				threshold_format |= THRESHOLD_PROBE_TM;
				break;
			case 'z':
				min_primer_dg = float( atof(optarg) );
				threshold_format |= THRESHOLD_PRIMER_DELTA_G;
				break;
			case 'Z':
				min_probe_dg = float( atof(optarg) );
				threshold_format |= THRESHOLD_PROBE_DELTA_G;
				break;
			case 'x':
				max_primer_tm = float( atof(optarg) );
				threshold_format |= THRESHOLD_PRIMER_TM;
				break;
			case 'X':
				max_probe_tm = float( atof(optarg) );
				threshold_format |= THRESHOLD_PROBE_TM;
				break;
			case 'g':
				max_primer_dg = float( atof(optarg) );
				threshold_format |= THRESHOLD_PRIMER_DELTA_G;
				break;
			case 'G':
				max_probe_dg = float( atof(optarg) );
				threshold_format |= THRESHOLD_PROBE_DELTA_G;
				break;
			case 's':
				salt = float( atof(optarg) );
				break;
			case 't':
				primer_strand = float( atof(optarg) );
				break;
			case 'T':
				probe_strand = float( atof(optarg) );
				break;
			case 'y':
				asymmetric_strand_ratio = float( atof(optarg) );
				break;
			case 'A':
				assay_format = parse_assay_format(optarg);
				break;
			case 'W':
				hash_word_size = atoi(optarg);
				break;
			case 'm':
				parse_output_file(optarg);
				break;
			case 'a':
				if(parse_bool(optarg) == true){
				
					// Set the bit
					output_format |= OUTPUT_ALIGNMENTS;
				}
				else{
					// Unset the bit
					output_format &= ~OUTPUT_ALIGNMENTS;
				}
				
				break;
			case 'M':
				if(parse_bool(optarg) == true){
				
					// Set the bit
					output_format |= OUTPUT_SEQ_MATCH;
				}
				else{
					// Unset the bit
					output_format &= ~OUTPUT_SEQ_MATCH;
				}
				
				break;
			case 'k':
				if(parse_bool(optarg) == true){
					mask_options |= MASK_PRIMERS;
				}
				else{
					mask_options &= ~MASK_PRIMERS;
				}
				break;
			case 'K':
				if(parse_bool(optarg) == true){
					mask_options |= MASK_PROBE;
				}
				else{
					mask_options &= ~MASK_PROBE;
				}
				break;
			case 'r':
				if(parse_bool(optarg) == true){
					mask_options |= REPLACE_PRIMERS;
				}
				else{
					mask_options &= ~REPLACE_PRIMERS;
				}
				break;
			case 'v':
				verbose = parse_bool(optarg);
				break;
			case 'p':
				ignore_probe = parse_bool(optarg);
				break;
			case 'n':
				one_output_file_per_query = parse_bool(optarg);
				break;
			case 'L':
				append_name_to_defline = parse_bool(optarg);
				break;
			case 'S':
				assay_summary = parse_bool(optarg);
				break;
			case 'h':
			case '?':
				print_usage = true;
				break;
			default:
				cerr << '\"' << (char)opt_code << "\" is not a valid option!" << endl;
				break;
		};
	}
	
	if(print_usage){
	
		cerr << "thermonucleotideBLAST v." << TNTBLAST_VERSION << endl;
		cerr << "Options:" << endl;
		cerr << "\t-i <input file of query oligos>" << endl;
		cerr << "\t-o <output file> (default is stdout)" << endl;
		cerr << "\t-d <database of target sequences to search against>" << endl;
		cerr << "\t[-D <local database of target sequences to search against>]" << endl;
		cerr << "\t[-l <maximum amplicon length> (default is " << DEFAULT_MAX_LEN << " bases)" << endl;
		cerr << "\t-e <minimum primer Tm>" << endl;
		cerr << "\t-E <minimum probe Tm>" << endl;
		cerr << "\t[-z <minimum primer delta G (in Kcal/Mol)>] (default is no limit)" << endl;
		cerr << "\t[-Z <minimum probe delta G (in Kcal/Mol)>] (default is no limit)" << endl;
		cerr << "\t[-x <maximum primer Tm>] (default is no limit)" << endl;
		cerr << "\t[-X <maximum probe Tm>] (default is no limit)" << endl;
		cerr << "\t[-g <maximum primer delta G (in Kcal/Mol)>] (default is no limit)" << endl;
		cerr << "\t[-G <maximum probe delta G(in Kcal/Mol)>] (default is no limit)" << endl;
		cerr << "\t[-s <salt concentration (in MOL)>] (default is " << DEFAULT_SALT << " M)" << endl;
		cerr << "\t[-t <primer strand concentration (in MOL)>] (default is " << DEFAULT_PRIMER_STRAND << " M)" << endl;
		cerr << "\t[-T <Probe strand concentration (in MOL)>] (default is " << DEFAULT_PROBE_STRAND << " M)" << endl;
		cerr << "\t[-y <ratio of forward/reverse strand concentrations>] (default is 1, i.e. symmetric PCR)" << endl;
		cerr << "\t[-A <PCR | PROBE | PADLOCK | AFFY>] (assay format, default is PCR)" << endl;
		cerr << "\t[-W <2-8>] (hash word length, default is " << DEFAULT_HASH_WORD_SIZE << ")" << endl;
		cerr << "\t[-m <output format>] " << endl;
		cerr << "\t\t0 = verbose output file (default)" << endl;
		cerr << "\t\t1 = fasta output file" << endl;
		cerr << "\t\t2 = network output files (*.atr and *.sif)" << endl;
		cerr << "\t\t3 = \"inverse target\" (targets that *don't* match any query)" << endl;
		cerr << "\t\t4 = \"inverse query\" (queries that *don't* match any target)" << endl;
		cerr << "\t[-a <T|F>] (show alignments, default is T)" << endl;
		cerr << "\t[-M <T|F>] (show matching sequence, default is T)" << endl;
		cerr << "\t[-k <T|F>] (Mask primer binding sites, default is F)" << endl;
		cerr << "\t[-K <T|F>] (Mask probe binding sites, default is F)" << endl;
		cerr << "\t[-r <T|F>] (Replace primer binding sites w/ primer sequence, default is F)" << endl;
		cerr << "\t[-v <T|F>] (Disable verbose terminal output, default is T)" << endl;
		cerr << "\t[-p <T|F>] (Ignore all probe oligos in inputfile, default is F)" << endl;
		cerr << "\t[-n <T|F>] (One output file per query, default is F)" << endl;
		cerr << "\t[-L <T|F>] (Append assay name to output defline, default is F)" << endl;
		cerr << "\t[-S <T|F>] (Ouput assay summary after searching, default is F)" << endl;
		cerr << "\t[-h|-?] (Command-line usage)" << endl;
		cerr << "\t[--primer-clamp <number of exact 3' primer matches requried>] (default is " 
			<< DEFAULT_PRIMER_CLAMP << " bases)" << endl;
		cerr << "\t[--min-max-primer-clamp <the minimum max number of exact 3' primer matches requried>] (default is no limit)" << endl;
		cerr << "\t[--probe-clamp5 <number of exact 5' probe matches requried>] (default is " 
			<< DEFAULT_PROBE_CLAMP_5 << " bases)" << endl;
		cerr << "\t[--probe-clamp3 <number of exact 3' probe matches requried>] (default is " 
			<< DEFAULT_PROBE_CLAMP_3 << " bases)" << endl;
		cerr << "\t[--dangle5 <T|F>] (Allow dangling bases on the 5' query side of an alignment, default is "
			<< (DEFAULT_DANGLE_5 ? "T" : "F") << ")" << endl;
		cerr << "\t[--dangle3 <T|F>] (Allow dangling bases on the 3' query side of an alignment, default is "
			 << (DEFAULT_DANGLE_3 ? "T" : "F") << ")" << endl;
		cerr << "\t[--plex <T|F>] (All input assays in a single multiple reaction, default is F)" << endl;
		cerr << "\t[--temperature <temperature for computing Delta G (in Kelvin)>] (default is " 
			<< DEFAULT_TARGET_T << " K)" << endl;
		cerr << "\t[--single-primer-pcr <T|F>] (Allow amplicons produced by a single PCR primer binding in both forward and reverse orientation, default is T)" << endl;
		cerr << "\t[--target-strand <plus|minus|both>] (which strand to target with probes, default is \"both\")" << endl;
		cerr << "\t[--max-target-len <max len>] (max sequence length before targets are split, default is " 
			<< DEFAULT_FRAGMENT_TARGET_LENGTH << " bases)" << endl;
		cerr << "\t[--query-seg <always | never | adaptive>] (query segmentation algorithm, default is \"adaptive\")" << endl;
		cerr << "\t[--dump-query <T|F>] (write queries to stdout, default is F)" << endl;
		cerr << "\t[--dinkelbach <T|F>] (Use the Dinkelbach fractional programming algorithm, default is F)" << endl;
		cerr << "\t[--max-gap <number of gaps>] (Max number of allowed gaps in a DNA duplex, default is " 
			<< DEFAULT_MAX_GAP << ")" << endl;
		cerr << "\t[--max-mismatch <number of mismatches>] (Max number of allowed mismatches in a DNA duplex, default is "
			<< DEFAULT_MAX_MISMATCH << ")" << endl;
		cerr << "\t[--rescale-ct <T|F>] (Use of degenerate bases results in rescaling of oligo concentration, default is "
			<< (DEFAULT_RESCALE_CT ? "T" : "F") << ")" << endl;
		cerr << "\t[--best-match] (Only save the best match, in Tm, between a query and target)" << endl;

		#ifdef USE_BLAST_DB
		cerr << "\t[--blast-include <Limit search to include accessions or NCBI TaxIds from a BLAST database>] (may be repeated)" << endl;
		cerr << "\t[--blast-exclude <Limit search to exclude accessions or NCBI TaxId from a BLAST database>] (may be repeated)" << endl;
		#endif // USE_BLAST_DB
	}
}

unsigned int Options::parse_assay_format(string m_opt)
{

	for(string::iterator i = m_opt.begin();i != m_opt.end();i++){

		*i = toupper(*i);
	}

	if(m_opt == "PCR"){
		return ASSAY_PCR;
	}
	
	if(m_opt == "PROBE"){
		return ASSAY_PROBE;
	}
	
	if(m_opt == "PADLOCK"){
		return ASSAY_PADLOCK;
	}
	
	if( (m_opt == "AFFYMETRIX") || (m_opt == "AFFY") ){
		return ASSAY_AFFYMETRIX;
	}
	
	return ASSAY_NONE;
}

void Options::validate_parameters()
{
	if( (dbase_filename == "") && (local_dbase_filename == "") ){
		throw "Unable to read either dbase or local_dbase";
	}
	
	if( (dbase_filename != "") && (local_dbase_filename != "") ){
		throw "Please specify either dbase or local_dbase (but not both)";
	}
	
	if(ignore_probe){
	
		if(assay_format == ASSAY_PROBE){
			throw "Error: Ignore probes (i.e. -p T) can not be used with the probe assay format";
		}
		
		if(verbose){
			cout << "** Ignoring all probe sequences **" << endl;
		}
	}
	
	if(salt <= 0.0f){
		throw "[Na+] (i.e. \"salt\") is less than zero";
	}

	if(salt >= 1.0f){
		throw "[Na+] (i.e. \"salt\") is greater than 1M";
	}
	
	if(primer_strand <= 0.0f){
		throw "[Ct] (i.e. \"primer_strand\") is less than zero";
	}

	if(primer_strand > 10.0f){
		throw "[Ct] (i.e. \"primer_strand\") is greater than 10M";
	}
	
	if(probe_strand < 0.0f){
	
		if(verbose){
			cout << "Setting probe strand concentration equal to primer strand concentration" << endl;
		}
		
		probe_strand = primer_strand;
	}
	
	if(probe_strand <= 0.0f){
		throw "[Ct] (i.e. \"probe_strand\") is less than zero";
	}

	if(probe_strand > 10.0f){
		throw "[Ct] (i.e. \"probe_strand\") is greater than 10M";
	}
	
	if(asymmetric_strand_ratio <= 0.0){
		throw "The ratio of forward to reverse primer [Ct] is <= 0";
	}
	 
	if(min_primer_tm < 0.0f){
		throw "min_primer_tm is less than zero";
	}

	if(min_primer_tm > 200.0f){
		throw "min_primer_tm is greater than 200 C -- that's too hot!";
	}
	
	if(max_primer_tm < 0.0f){
		throw "max_primer_tm is less than zero";
	}
	
	if(min_primer_tm > max_primer_tm){
		throw "min_primer_tm > max_primer_tm. Please use consistent values!";
	}
	
	if(min_probe_tm < 0.0f){
		throw "min_probe_tm is less than zero";
	}

	if(min_probe_tm > 200.0f){
		throw "min_probe_tm is greater than 200 C -- that's too hot!";
	}
	
	if(max_probe_tm < 0.0f){
		throw "max_probe_tm is less than zero";
	}

	if(min_probe_tm > max_probe_tm){
		throw "min_probe_tm > max_probe_tm. Please use consistent values!";
	}
	
	if(max_len <= 0){
		throw "max_len is less than 1 base -- too small!";
	}
	
	if(primer_clamp < 0){
		throw "primer_clamp is less than 0 -- too small!";
	}
	
	if(probe_clamp_5 < 0){
		throw "probe_clamp_5 is less than 0 -- too small!";
	}
	
	if(probe_clamp_3 < 0){
		throw "probe_clamp_3 is less than 0 -- too small!";
	}
	
	if(assay_format == ASSAY_NONE){
		throw "Please specify a valid assay format";
	}
	
	if( (hash_word_size < 3) || (hash_word_size > 8) ){
		throw "Please specify a valid hash word size";
	}
	
	if(output_format & OUTPUT_NETWORK){
		
		if(output_filename == ""){
			throw "Please specify an output filename when writing "
				"network files";
		}
	}
	
	if(max_gap < 0){
		throw "Error: --max-gap < 0";
	}
	
	if(max_mismatch < 0){
		throw "Error: --max-mismatch < 0";
	}

	if(verbose){
		
		switch(query_segmentation){
			case QUERY_SEGMENTATION_ON:
				cout << "Query segmentation: always on" << endl;
				break;
			case QUERY_SEGMENTATION_OFF:
				cout << "Query segmentation: disabled" << endl;
				break;
			case QUERY_SEGMENTATION_ADAPTIVE:
				cout << "Query segmentation: adaptive" << endl;
				break;
			default:
				throw "Unknown option for query segmentation";
				break;
		};
	}
}

void Options::parse_output_file(const string &m_format)
{
	const int opt = atoi( m_format.c_str() );

	// Unset any previous output file format
	output_format &= ~OUTPUT_STANDARD;
	output_format &= ~OUTPUT_FASTA;
	output_format &= ~OUTPUT_NETWORK;
	output_format &= ~OUTPUT_INVERSE_TARGET;
	output_format &= ~OUTPUT_INVERSE_QUERY;
	
	switch(opt){
		case 0: // 0 = standard verbose output file
			output_format |= OUTPUT_STANDARD;
			break;
		case 1: // 1 = fasta output file
			output_format |= OUTPUT_FASTA;
			break;
		case 2: // 2 = network files (*.atr and *.sif)
			output_format |= OUTPUT_NETWORK;
			break;
		case 3: // 3 = targets (deflines) that *don't* match a query
			output_format |= OUTPUT_INVERSE_TARGET;
			break;
		case 4: // 4 = queries that *don't* match a target
			output_format |= OUTPUT_INVERSE_QUERY;
			break;
		default:
			throw "Unknown output format. Please specify a number between 0-3";
			break;
	};
}

bool Options::parse_bool(string m_opt)
{
	// Make the input string upper case
	for(string::iterator i = m_opt.begin();i != m_opt.end();i++){
		*i = toupper(*i);
	}

	if( (m_opt == "T") || (m_opt == "TRUE") ){
		return true;
	}
	
	if( (m_opt == "F") || (m_opt == "FALSE") ){
		return false;
	}
	
	throw "Unknown boolean options -- please use \"T\" or \"F\"";
	
	return false;
}

int Options::parse_strand(string m_opt)
{
	// Make the input string upper case
	for(string::iterator i = m_opt.begin();i != m_opt.end();i++){
		*i = toupper(*i);
	}
	
	if( (m_opt == "PLUS") || (m_opt == "+") || (m_opt == "SENSE") ){
		return Seq_strand_plus;
	}
	
	if( (m_opt == "MINUS") || (m_opt == "-") || (m_opt == "ANTISENSE") ){
		return Seq_strand_minus;
	}
	
	if(m_opt == "BOTH"){
		return Seq_strand_both;
	}
	
	cerr << "Unknown target-strand option." << endl;
	cerr << "Use \"+\", \"plus\" or \"sense\" for the sense strand" << endl;
	cerr << "Use \"-\", \"minus\" or \"antisense\" for the antisense strand" << endl;
	cerr << "Use \"both\" for both sense and antisense strands" << endl;
	
	// If we get here, we could not identify the users input
	throw "Unknown target-strand option";
}

int Options::parse_query_seg(string m_opt)
{
	// Make the input string upper case
	for(string::iterator i = m_opt.begin();i != m_opt.end();i++){
		*i = toupper(*i);
	}
	
	if(m_opt == "ALWAYS"){
		
		return QUERY_SEGMENTATION_ON;
	}
	
	if(m_opt == "NEVER"){
		
		return QUERY_SEGMENTATION_OFF;
	}
	
	if(m_opt == "ADAPTIVE"){
		
		return QUERY_SEGMENTATION_ADAPTIVE;
	}
	
	cerr << "Unknown query segmentation option." << endl;
	cerr << "Use \"always\" to force segmentation" << endl;
	cerr << "Use \"never\" to disable segmentation" << endl;
	cerr << "Use \"adaptive\" for an adaptive algorithm" << endl;
	
	// If we get here, we could not identify the users input
	throw "Unknown query segmentation option";

}

size_t Options::max_product_length() const
{
	size_t ret = 0;
	
	if(assay_format == ASSAY_PCR){
	
		vector<hybrid_sig>::const_iterator iter;
	
		for(iter = sig_list.begin();iter != sig_list.end();iter++){
		
			if(iter->has_primers() == true){
			
				ret = max_len;
				break;
			}
			
			ret = max( ret, iter->probe_oligo.size() );
		}
		
		return ret;
	}
	
	if(assay_format == ASSAY_PADLOCK){
		
		vector<hybrid_sig>::const_iterator iter;
	
		for(iter = sig_list.begin();iter != sig_list.end();iter++){
			
			ret = max( ret, iter->forward_oligo.size() + iter->reverse_oligo.size() );
		}
		
		return ret;
	}
	
	vector<hybrid_sig>::const_iterator iter;
	
	for(iter = sig_list.begin();iter != sig_list.end();iter++){

		ret = max( ret, iter->probe_oligo.size() );
	}

	return ret;
}

void Options::validate_search_threshold()
{	
	// Make sure that the user has specific a threshold that is appropriate to 
	// the assay format
	switch(assay_format){
		case ASSAY_PCR:
			
			// Since the ASSAY_PCR format can be either PCR primers, or PCR primers and a probe,
			// or just a probe, we need to check every assay that the user has provided
			for(vector<hybrid_sig>::const_iterator iter = sig_list.begin();iter != sig_list.end();iter++){
				
				if(iter->has_primers() == true){
					
					if( !(threshold_format & THRESHOLD_PRIMER_DELTA_G) &&
					    !(threshold_format & THRESHOLD_PRIMER_TM) ){

						throw "Please specify primer search bounds in Tm and/or Delta G";
					}
				}
				
				if(iter->has_probe() == true){
					
					if( !(threshold_format & THRESHOLD_PROBE_DELTA_G) &&
					    !(threshold_format & THRESHOLD_PROBE_TM) ){

						throw "Please specify probe search bounds in Tm and/or Delta G";
					}
				}
			}
			
			break;
		case ASSAY_AFFYMETRIX:
		case ASSAY_PROBE:
			
			if( !(threshold_format & THRESHOLD_PROBE_DELTA_G) &&
			    !(threshold_format & THRESHOLD_PROBE_TM) ){
			    	
				throw "Please specify probe search bounds in Tm and/or Delta G";
			}
			
			break;
		case ASSAY_PADLOCK:
			
			if( !(threshold_format & THRESHOLD_PROBE_DELTA_G) &&
			    !(threshold_format & THRESHOLD_PROBE_TM) ){
			    	
				throw "Please specify probe search bounds in Tm and/or Delta G";
			}
			
			break;
		case ASSAY_NONE:
			
			throw "No assay format has been specified!";
			break;
		default:
			throw __FILE__ ":validate_search_threshold: Unknown assay format!";
			break;
	};
}

void Options::write_queries(std::ostream &s)
{

	vector<hybrid_sig>::const_iterator iter;
	
	for(iter = sig_list.begin();iter != sig_list.end();iter ++){
		
		s << iter->name;
		
		if(iter->has_primers() == true){
			
			s << '\t' << iter->forward_oligo << '\t' << iter->reverse_oligo;
		}
			
		if(iter->has_probe() == true){
			s << '\t' << iter->probe_oligo;
		}
		
		s << endl;
	}
}

ostream& operator << (ostream &s, const Options &m_opt)
{
	
	s << "Found " << m_opt.sig_list.size() << " query assays" << endl;
	s << "Search parameters:" << endl;
	s << "\tOutput = " << m_opt.output_filename << endl;
	s << "\t[Na+] = " << m_opt.salt << " M" << endl;
	s << "\tmax gap = " << m_opt.max_gap << endl;
	s << "\tmax mismatch = " << m_opt.max_mismatch << endl;
	
	if( m_opt.has_primers() ){

		if(m_opt.asymmetric_strand_ratio != 1.0f){

			s << "\t[reverse primer Ct] = " << m_opt.primer_strand 
				<< " M" << endl;
			s << "\t[forward primer Ct]/[reverse primer Ct] = " 
				<< m_opt.asymmetric_strand_ratio << endl;
		}
		else{
			s << "\t[primer Ct] = " << m_opt.primer_strand << " M" << endl;
		}
	}

	if( m_opt.has_probe() ){
		s << "\t[probe Ct] = " << m_opt.probe_strand << " M" << endl;
	}

	if( m_opt.has_primers() ){

		if(m_opt.assay_format == ASSAY_PCR){ // These are PCR primers
		
			s << '\t' << m_opt.min_primer_tm << " <= Primer Tm (C) <= " << m_opt.max_primer_tm << endl;
			s << '\t' << m_opt.min_primer_dg << " <= Primer Delta G (Kcal/Mol) <= " 
				<< m_opt.max_primer_dg << endl;
		}
		else{ // These are Padlock probes (not actually primers)
		
			s << '\t' << m_opt.min_primer_tm << " <= Padlock Tm (C) <= " << m_opt.max_primer_tm << endl;
			s << '\t' << m_opt.min_primer_dg << " <= Padlock Delta G (Kcal/Mol) <= " 
				<< m_opt.max_primer_dg << endl;
		}
	}

	if( m_opt.has_probe() ){

		s << '\t' << m_opt.min_probe_tm << " <= Probe Tm (C) <= " << m_opt.max_probe_tm << endl;
		s << '\t' << m_opt.min_probe_dg << " <= Probe Delta G (Kcal/Mol) <= " 
			<< m_opt.max_probe_dg << endl;
	}

	if(m_opt.assay_format == ASSAY_PADLOCK){

		s << "\t5' Ligation clamp = " << m_opt.probe_clamp_5 << endl;
		s << "\t3' Ligation clamp = " << m_opt.probe_clamp_3 << endl;
		s << "Assay format is PADLOCK/MOL-PCR" << endl;
	}
	else{

		if( m_opt.has_primers() ){

			s << "\t3' Primer clamp = " << m_opt.primer_clamp << endl;
			
			if(m_opt.min_max_primer_clamp >= 0){
				s << "\tThe minimum, maximum 3' Primer clamp = " << m_opt.min_max_primer_clamp << endl;
			}
		}

		if( m_opt.has_probe() ){

			s << "\t5' Probe clamp = " << m_opt.probe_clamp_5 << endl;
			s << "\t3' Probe clamp = " << m_opt.probe_clamp_3 << endl;
		}

		if( m_opt.has_primers() ){

			s << "\tMax amplicon len = " << m_opt.max_len << endl;
			s << "Assay format is PCR and/or PROBE" << endl;

			if(m_opt.single_primer_pcr == false){				
				s << "Single primers will *not* be tested for amplicon generation" << endl;
			}
		}

		if(m_opt.assay_format == ASSAY_AFFYMETRIX){
			s << "Assay format is Affymetrix PROBE" << endl;
		}
	}
	
	return s;
}

