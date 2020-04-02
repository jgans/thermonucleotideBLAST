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

struct program_options
{
	private:
	
	// The threshold format allows us to insure that the
	// user has specified an appropriate search threshold
	// (i.e. primer and/or probe free energy and/or Tm).
	unsigned int threshold_format;
	
	public:
	
	std::string dbase_filename;
	std::string local_dbase_filename;
	std::string output_filename;
	std::string input_filename;

	int max_len;
	int primer_clamp;
	int min_max_primer_clamp;
	int probe_clamp_5;
	int probe_clamp_3;
	
	// Maximum number of gaps allowed in a DNA duplex
	int max_gap;
	
	// Maximum number of mismatches allowed in a DNA duplex
	int max_mismatch;
	
	// For the target sequences, which strand(s) are present?
	// Choices are: Seq_strand_plus, Seq_strand_minus or Seq_strand_both
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
	float primer_strand;
	float probe_strand;
	float asymmetric_strand_ratio;
	
	// The target temperature for computing delta G
	float target_t;
	
	unsigned int output_format;
	
	bool ignore_probe;
	bool one_output_file_per_query;
	bool append_name_to_defline;
	bool verbose;
	bool print_usage;
	bool assay_summary;
	bool allow_fasta_mmap;
	bool multiplex;
	bool single_primer_pcr;
	bool dump_query;
	
	bool use_dinkelbach;
	
	// Allow dangling bases on the 5' and/or 3' query side of the alignment
	bool allow_dangle_5;
	bool allow_dangle_3;
	
	// Should the use of degenerate bases result in a rescaling of oligo concentration?
	bool degen_rescale_ct;

	bool best_match;

	int melting_param;
	int mask_options;
	int assay_format;

	int query_segmentation;

	int hash_word_size;
	
	int fragment_target_threshold;
	
	std::vector<hybrid_sig> sig_list;
	
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
		probe_strand = -1.0f;
		
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
		allow_fasta_mmap = true;
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
	
	program_options()
	{
		// Initialize the parameters
		default_values();
	};
	
	program_options(int argc, char *argv[])
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
	
	unsigned int parse_config_file(const std::string &m_filename);
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
	
	friend std::ostream& operator << (std::ostream &s, const program_options &m_opt);
};

#endif // __TNTBLAST_OPTIONS
