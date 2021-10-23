#ifndef __TNTBLAST_OPTIONS
#define __TNTBLAST_OPTIONS

#include "tntblast.h"

// Enumerate all output formats
#define	OUTPUT_STANDARD				(1 << 0)
#define	OUTPUT_FASTA				(1 << 1)
#define	OUTPUT_NETWORK				(1 << 2)
#define	OUTPUT_INVERSE_TARGET		(1 << 3)
#define	OUTPUT_INVERSE_QUERY		(1 << 4)
#define	OUTPUT_ALIGNMENTS			(1 << 5)
#define	OUTPUT_SEQ_MATCH			(1 << 6)

#define	THRESHOLD_NONE 				0
#define	THRESHOLD_PRIMER_DELTA_G 	(1 << 0)
#define	THRESHOLD_PRIMER_TM 		(1 << 1)
#define	THRESHOLD_PROBE_DELTA_G 	(1 << 2)
#define	THRESHOLD_PROBE_TM 			(1 << 3)

struct Options
{
	// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
	// ensure that structure variable are correctly serialized.
	#define OPTIONS_MEMBERS \
		VARIABLE(std::string, dbase_filename) \
		VARIABLE(std::string, local_dbase_filename) \
		VARIABLE(std::string, output_filename) \
		VARIABLE(std::string, input_filename) \
		VARIABLE(std::vector<hybrid_sig>, sig_list) \
		VARIABLE(std::vector<std::string>, blast_include) \
		VARIABLE(std::vector<std::string>, blast_exclude) \
		VARIABLE(int, max_len) \
		VARIABLE(int, primer_clamp) \
		VARIABLE(int, min_max_primer_clamp) \
		VARIABLE(int, probe_clamp_5) \
		VARIABLE(int, probe_clamp_3) \
		VARIABLE(int, max_gap) \
		VARIABLE(int, max_mismatch) \
		VARIABLE(int, target_strand) \
		VARIABLE(int, melting_param) \
		VARIABLE(int, mask_options) \
		VARIABLE(int, assay_format) \
		VARIABLE(int, query_segmentation) \
		VARIABLE(int, hash_word_size) \
		VARIABLE(int, fragment_target_threshold) \
		VARIABLE(float, min_primer_tm) \
		VARIABLE(float, max_primer_tm) \
		VARIABLE(float, min_primer_dg) \
		VARIABLE(float, max_primer_dg) \
		VARIABLE(float, min_probe_tm) \
		VARIABLE(float, max_probe_tm) \
		VARIABLE(float, min_probe_dg) \
		VARIABLE(float, max_probe_dg) \
		VARIABLE(float, salt) \
		VARIABLE(float, primer_strand) \
		VARIABLE(float, probe_strand) \
		VARIABLE(float, asymmetric_strand_ratio) \
		VARIABLE(float, target_t) \
		VARIABLE(unsigned int, output_format) \
		VARIABLE(unsigned int, threshold_format) \
		VARIABLE(bool, ignore_probe) \
		VARIABLE(bool, one_output_file_per_query) \
		VARIABLE(bool, append_name_to_defline) \
		VARIABLE(bool, verbose) \
		VARIABLE(bool, print_usage) \
		VARIABLE(bool, assay_summary) \
		VARIABLE(bool, multiplex) \
		VARIABLE(bool, single_primer_pcr) \
		VARIABLE(bool, dump_query) \
		VARIABLE(bool, use_dinkelbach) \
		VARIABLE(bool, allow_dangle_5) \
		VARIABLE(bool, allow_dangle_3) \
		VARIABLE(bool, degen_rescale_ct) \
		VARIABLE(bool, best_match)
	
	public:

		#define VARIABLE(A, B) A B;
			OPTIONS_MEMBERS
		#undef VARIABLE
	
		// Set sensible default values for program parameters
		void default_values()
		{
			dbase_filename = "";
			local_dbase_filename = "";
			output_filename = "";
			input_filename = "";

			// Maximum amplicon length
			max_len = DEFAULT_MAX_LEN;
			
			// 3' primer clamp (i.e. number of exact matches required)
			primer_clamp = DEFAULT_PRIMER_CLAMP;
			min_max_primer_clamp = DEFAULT_MIN_MAX_PRIMER_CLAMP;
			probe_clamp_5 = DEFAULT_PROBE_CLAMP_5;
			probe_clamp_3 = DEFAULT_PROBE_CLAMP_3;
			
			max_gap = DEFAULT_MAX_GAP;
			max_mismatch = DEFAULT_MAX_MISMATCH;
		
			target_strand = Seq_strand_both;
			
			min_primer_tm = DEFAULT_MIN_PRIMER_TM;
			max_primer_tm = DEFAULT_MAX_PRIMER_TM;
			
			min_primer_dg = DEFAULT_MIN_PRIMER_DG;
			max_primer_dg = DEFAULT_MAX_PRIMER_DG;
			
			min_probe_tm = DEFAULT_MIN_PROBE_TM;
			max_probe_tm = DEFAULT_MAX_PROBE_TM;
			
			min_probe_dg = DEFAULT_MIN_PROBE_DG;
			max_probe_dg = DEFAULT_MAX_PROBE_DG;

			salt = DEFAULT_SALT;
			primer_strand = DEFAULT_PRIMER_STRAND;
			probe_strand = DEFAULT_PROBE_STRAND;
			
			target_t = DEFAULT_TARGET_T;
			
			// By default the forward and reverse primers are present in equal abundance
			asymmetric_strand_ratio = 1.0f;
			
			print_usage = false;
			
			output_format = OUTPUT_STANDARD | OUTPUT_ALIGNMENTS | OUTPUT_SEQ_MATCH;
			
			mask_options = NO_MASK;
			verbose = true;
			ignore_probe = false;
			one_output_file_per_query = false;
			append_name_to_defline = false;
			assay_summary = false;
			multiplex = false;
			dump_query = false;
			
			use_dinkelbach = false;
			
			allow_dangle_5 = DEFAULT_DANGLE_5;
			allow_dangle_3 = DEFAULT_DANGLE_3;
		
			degen_rescale_ct = DEFAULT_RESCALE_CT;

			best_match = false;

			// Allow amplicons produced by a single PCR primer
			// binding in both forward and reverse orientation?
			single_primer_pcr = true;
			
			// Use the adaptive algorithm (in query_sched() -- see tntblast_util.cpp)
			// to decide when to start segmenting queries
			query_segmentation = QUERY_SEGMENTATION_ADAPTIVE;
			
			// Use the SantaLucia hybridization parameters
			melting_param = DEFAULT_MELT_PARAM;
			assay_format = DEFAULT_ASSAY_FORMAT;

			// How should we identify sequences for hybridization?
			hash_word_size = DEFAULT_HASH_WORD_SIZE;
			
			fragment_target_threshold = DEFAULT_FRAGMENT_TARGET_LENGTH;
			
			threshold_format = THRESHOLD_NONE;
		};
		
		Options()
		{
			// Initialize the parameters
			default_values();
		};
		
		Options(int argc, char *argv[])
		{
			// Initialize the parameters
			default_values();
			
			parse(argc, argv);
		}
		
		void parse(int argc, char *argv[])
		{
			// The -c flag controls reading of a configuration file
			// (and parse_config_file is called by parse_command_line).
			parse_command_line(argc, argv);
			
			// Don't validate parameters, if the user is only
			// printing the command line options
			if(print_usage == false){
				validate_parameters();
			}
		};
		
		//unsigned int parse_config_file(const std::string &m_filename);
		void parse_command_line(int argc, char *argv[]);
		unsigned int parse_assay_format(std::string m_opt);
		void parse_output_file(const std::string &m_format);
		bool parse_bool(std::string m_opt);
		int parse_strand(std::string m_opt);
		int parse_query_seg(std::string m_opt);
		
		void validate_parameters();
		
		void validate_search_threshold();
		
		inline bool has_probe() const
		{
			return (assay_format == ASSAY_PROBE) || 
				(assay_format == ASSAY_PCR) ||
				(assay_format == ASSAY_AFFYMETRIX);
		};
		
		inline bool has_primers() const
		{
			return (assay_format == ASSAY_PCR) ||
				(assay_format == ASSAY_PADLOCK);
		};
		
		size_t max_product_length() const;
		
		void write_queries(std::ostream &s);
		
		friend std::ostream& operator << (std::ostream &s, const Options &m_opt);

		template<class T> friend size_t mpi_size(const T &m_obj);
		template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, 
			const T &m_obj);
		template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr,
			const T &m_obj);
};

template<> size_t mpi_size(const Options &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const Options &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, Options &m_obj);

#endif // __TNTBLAST_OPTIONS
