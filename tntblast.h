#ifndef __AMPLICON
#define __AMPLICON

#include "hybrid_sig.h"
#include "seq.h"
#include "sequence_data.h"
#include "nuc_cruc.h"

#include "tntblast_timer.h"

#include <list>
#include <vector>
#include <set>
#include <unordered_map>

/////////////////////////////////////////////////////////////////////////
// Default parameter values (applied in options.h)

#define	DEFAULT_MAX_LEN							2000
		
// 3' primer clamp (i.e. number of exact matches required)
#define	DEFAULT_PRIMER_CLAMP					0

// The smallest allowed maximum primer clamp (like the lowest, highest point)
// (a value of -1 deactivates this test)
#define	DEFAULT_MIN_MAX_PRIMER_CLAMP			-1

#define	DEFAULT_PROBE_CLAMP_5					0
#define	DEFAULT_PROBE_CLAMP_3					0

#define	DEFAULT_MIN_PRIMER_TM					0.0f
#define	DEFAULT_MAX_PRIMER_TM					9999.0f

#define	DEFAULT_MIN_PROBE_TM					0.0f
#define	DEFAULT_MAX_PROBE_TM					9999.0f

#define	DEFAULT_MIN_PRIMER_DG					-9999.0f
#define	DEFAULT_MAX_PRIMER_DG					0.0f

#define	DEFAULT_MIN_PROBE_DG					-9999.0f
#define	DEFAULT_MAX_PROBE_DG					0.0f

#define	DEFAULT_MELT_PARAM 						NucCruc::SANTA_LUCIA
#define	DEFAULT_SALT							50.0e-3f	// Moles
#define	DEFAULT_PRIMER_STRAND					9.0e-7f // Moles
#define	DEFAULT_PROBE_STRAND					2.5e-7f // Moles
#define	DEFAULT_TARGET_T						310.15f	// Kelvin

#define	DEFAULT_ASSAY_FORMAT					ASSAY_PCR

#define	DEFAULT_HASH_WORD_SIZE					7

// Should dangling bases be allowed on the 5' and/or 3' query side of an alignment?
// (Note that the hairpin calculation always uses a dangling end)
#define	DEFAULT_DANGLE_5						false
#define	DEFAULT_DANGLE_3						false

// The largest sequence size (in bases) that will be searched
// without fragmentation (i.e. spliting up into overlapping sub-sequences
// that are all smaller than or equal to DEFAULT_FRAGMENT_TARGET_LENGTH).
#define	DEFAULT_FRAGMENT_TARGET_LENGTH			500000

// The maximum number of gaps allowed in DNA duplex
#define	DEFAULT_MAX_GAP							999

// The maximum number of mismatched allowed in DNA duplex
#define	DEFAULT_MAX_MISMATCH					999

// By default, rescale the strand concentration to reflect the use of degenerate bases
#define	DEFAULT_RESCALE_CT						true

// The number of additional flanking target bases to include on the 5' and 3'
// of the target when computing the target-query alignment
#define	NUM_FLANK_BASE							4

/////////////////////////////////////////////////////////////////////////

// The default ratio of query search time to target sequence load and hash time
// (for deciding when to segment queries between workers). Note that this value
// depends on a number of factors, including choice of alignment algorithm. 
// As the computational cost of the alignment algorithm increases, so should 
// the vallue of DEFAULT_QT. Note that this value is scaled depending
// on the assay type of the query:
// 1X	for PROBE	query targeting a single strand
// 2X	for PROBE queries targeting both strands
// 4X	for PCR queries
// 4X	for PADLOCK queries
#define	DEFAULT_QT						0.25f

/////////////////////////////////////////////////////////////////////////

// Version control
#define		TNTBLAST_VERSION		"2.61 (September 1, 2023)"

// Email address for send complaints, questions, kudos, rants, etc.
#define		EMAIL_ADDRESS			"jgans@lanl.gov"

// MPI message tags
#define		SIGNATURE			102
#define		SIGNATURE_COMPLETE	103
#define		SIGNATURE_RESULTS	104
#define		STATUS_UPDATE		105
#define		SEQ_REQUEST			106
#define		SEARCH_QUERY		107
#define		SEARCH_COMPLETE		108
#define		PROFILE_INFO		109

// The size of the SEARCH_QUERY message buffer
#define		SEARCH_QUERY_BUFFER_SIZE	6

// Sequence masking options
#define		NO_MASK			0
#define		MASK_PRIMERS	1	// lower-case primer binding sites
#define		MASK_PROBE		2	// lower-case probe binding sites
#define		REPLACE_PRIMERS	4	// primer sequence replaces primer binding sites

#define		DEFLINE_BUFFER_SIZE		512

// How should we perform query segmentation?
#define	QUERY_SEGMENTATION_ON			0
#define	QUERY_SEGMENTATION_OFF			1
#define	QUERY_SEGMENTATION_ADAPTIVE		2

// Best match rules
enum {BEST_MATCH_NOT_SET, BEST_MATCH_TM, BEST_MATCH_DG};

#ifdef PROFILE

enum {
	PROFILE_WORK, 			// Time spent searching
	PROFILE_LOAD,			// Time spent loading sequence data
	PROFILE_HASH,			// Time spent hashing sequence data
	PROFILE_COMM,			// Time spent on communication
	PROFILE_IDLE,			// Time spent idle
	PROFILE_NUM_PLUS_TM_EVAL,	// Number of + strand Tm evaluations
	PROFILE_NUM_MINUS_TM_EVAL,	// Number of - strand Tm evaluations
	NUM_PROFILE
};


#endif // PROFILE

struct oligo_info{

	enum {
		F = (1 << 0),
		R = (1 << 1),
		P = (1 << 2),
		PLUS_STRAND = (1 << 3),
		MINUS_STRAND = (1 << 4),
		VALID = (1 << 5)
	};
	
	oligo_info(const int &m_loc_5, const int &m_loc_3, const float &m_tm, 
		const float &m_dH, const float &m_dS, 
		const unsigned int &m_anchor_5, const unsigned int &m_anchor_3,
		const unsigned int &m_num_mm, const unsigned int &m_num_gap,
		const std::string &m_align):
		loc_5(m_loc_5), loc_3(m_loc_3), tm(m_tm), dH(m_dH), dS(m_dS), 
		anchor_5(m_anchor_5), anchor_3(m_anchor_3), 
		num_mm(m_num_mm), num_gap (m_num_gap),
		alignment(m_align)
	{
		query_loc = 0;
		target_loc = 0;
		mask = 0;
	};
	
	oligo_info(const unsigned int &m_query_loc,
		const unsigned int &m_target_loc,
		const unsigned char &m_mask) : query_loc(m_query_loc), target_loc(m_target_loc), mask(m_mask)
	{
		loc_5 = 0;
		loc_3 = 0;
		tm = -1.0f;
		dH = -1.0f;
		dS = -1.0f;
		anchor_5 = 0;
		anchor_3 = 0;
		num_mm = 0;
		num_gap = 0;
	};
	
	~oligo_info()
	{
		// Do nothing
	};
	
	int loc_5;
	int loc_3;
	float tm;
	float dH;
	float dS;
	unsigned int anchor_5;
	unsigned int anchor_3;
	unsigned int num_mm;
	unsigned int num_gap;
	std::string alignment;
	
	unsigned int query_loc;
	unsigned int target_loc;
	unsigned char mask;
	
	inline bool operator==(const oligo_info &m_oligo)
	{
		return (loc_5 == m_oligo.loc_5) && (loc_3 == m_oligo.loc_3);
	};
	
	inline bool operator!=(const oligo_info &m_oligo)
	{
		return !(*this == m_oligo);
	};

	inline bool operator>(const oligo_info &m_rhs) const
	{
		if(loc_5 != m_rhs.loc_5){
			return (loc_5 > m_rhs.loc_5);
		}
		
		if(loc_3 != m_rhs.loc_3){
			return (loc_3 > m_rhs.loc_3);
		}
		
		// Reverse the tm comparison to sort by descending tm
		return (tm < m_rhs.tm);
	}

	inline bool operator<(const oligo_info &m_rhs) const
	{
		if(loc_5 != m_rhs.loc_5){
			return (loc_5 < m_rhs.loc_5);
		}
		
		if(loc_3 != m_rhs.loc_3){
			return (loc_3 < m_rhs.loc_3);
		}
		
		// Reverse the tm comparison to sort by descending tm
		return (tm > m_rhs.tm);
	}
};

// Cache the melting temperature calculation results for a given oligo and target sequence range
// sequence location
struct BindCacheKey
{
	std::string oligo;
	unsigned int target_start;
	unsigned int target_stop;

	BindCacheKey()
	{
		// Do nothing
	};

	BindCacheKey(const std::string &m_oligo, const unsigned int &m_start,
		const unsigned int& m_stop) : oligo(m_oligo), target_start(m_start), target_stop(m_stop)
	{
	};

	inline bool operator==(const BindCacheKey &m_rhs) const
	{
		return (oligo == m_rhs.oligo) && 
			(target_start == m_rhs.target_start) && 
			(target_stop == m_rhs.target_stop);
	};
};

// Provide a hash function so we can store the cache results in an unordered_map
namespace std{

	template <>
	struct hash<BindCacheKey>
	{
		size_t operator()(const BindCacheKey& m_key) const
		{
			// Pack the target_start and target_stop values into a single 64 bit value
			return hash<string>()(m_key.oligo) ^ ( (size_t(m_key.target_start) << 32) | size_t(m_key.target_stop) );
		};
	};
}

struct BindCacheValue
{
	float tm;
	float dg;
	float dH;
	float dS;
	unsigned int anchor_5;
	unsigned int anchor_3;
	unsigned int target_5;
	unsigned int target_3;
	unsigned int num_mismatch;
	unsigned int num_gap;
	float fraction_real_base_pairs;
	std::string seq_align;

	BindCacheValue()
	{
		// Do nothing
	};

	BindCacheValue(const float &m_tm,
		const float &m_dg,
		const float &m_dH,
		const float &m_dS,
		const unsigned int &m_anchor_5,
		const unsigned int &m_anchor_3,
		const unsigned int &m_target_5,
		const unsigned int &m_target_3,
		const unsigned int &m_num_mismatch,
		const unsigned int &m_num_gap,
		const float &m_fraction_real_base_pairs,
		const std::string &m_seq_align) :
			tm(m_tm), dg(m_dg), dH(m_dH), dS(m_dS), anchor_5(m_anchor_5),
			anchor_3(m_anchor_3), target_5(m_target_5), target_3(m_target_3),
			num_mismatch(m_num_mismatch), num_gap(m_num_gap), fraction_real_base_pairs(m_fraction_real_base_pairs),
			seq_align(m_seq_align)
	{

	};
};

#ifdef USE_MPI
//////////////////////////////////////////////////////////////////////////////////////////////
// In tntblast_master.cpp
int worker(int argc, char *argv[]);

//////////////////////////////////////////////////////////////////////////////////////////////
// In tntblast_worker.cpp
int master(int argc, char *argv[]);

void distribute_queries(const std::vector<hybrid_sig> &m_sig);
void receive_queries(std::vector<hybrid_sig> &m_sig);
void serve_sequence(const int &m_dest, 
	const unsigned int &m_index, const sequence_data &m_data);
#endif // USE_MPI

//////////////////////////////////////////////////////////////////////////////////////////////
// Im tntblast_local.cpp
int local_main(int argc, char *argv[]);

//////////////////////////////////////////////////////////////////////////////////////////////
// In tntblast_util.cpp

// Show all primer binding sites in lower case (if requested by the user).
void mask_binding_sites(std::list<hybrid_sig> &m_sig, const int &m_mask, 
	const float &m_min_primer_tm, const float &m_min_probe_tm, NucCruc &m_melt,
	const float &m_forward_primer_strand, 
	const float &m_reverse_primer_strand, const float &m_probe_strand);
	
void mask_primer_5(std::string &m_amp, const std::string &m_oligo, NucCruc &m_melt,
	const bool &m_mask, const bool &m_replace);
void mask_primer_3(std::string &m_amp, const std::string &m_oligo, NucCruc &m_melt,
	const bool &m_mask, const bool &m_replace);
void mask_probe(std::string &m_amp, const std::string &m_primer, NucCruc &m_melt,
	const float &m_min_tm);
	
std::vector<hybrid_sig> expand_degenerate_signatures(const std::vector<hybrid_sig> &m_sig,
													 const bool &m_degen_rescale_ct);

std::vector<hybrid_sig> multiplex_expansion(const std::vector<hybrid_sig> &m_sig, 
	const unsigned int &m_format);

std::string primer_heuristics(const std::string &m_primer);

char base_complement(const char &m_base);

float gc_content(const std::string &m_seq);

unsigned int write_inverse_matches(std::ostream &m_fout, const sequence_data &m_data,
	std::set<std::string> &m_targets);

std::string mask_white_space(const std::string &m_str);

void write_alignment(std::ostream &fout, const std::string &m_prefix, 
	const std::string &m_alignment);

void write_annotation(std::ostream &fout, const hybrid_sig &m_sig, 
	const sequence_data &m_annot);

void test_memory(const int &m_num_mb);

void uniquify_results(std::list<hybrid_sig> &m_results);

void select_best_match(std::list<hybrid_sig> &m_results);

unsigned int probe_only_count(const std::vector<hybrid_sig> &m_queries);

bool query_sched(const unsigned int &m_num_target, 
	const unsigned int &m_num_query, const unsigned int &m_num_worker,
	const float &m_S_div_H, const int &m_mode);
	
//////////////////////////////////////////////////////////////////////////////////////////////
// In amplicon_search.cpp
std::list<hybrid_sig> amplicon(DNAHash &m_hash, const std::pair<std::string, SEQPTR> &m_seq, 
	const hybrid_sig &m_sig, NucCruc &m_melt,
	std::unordered_map<BindCacheKey, BindCacheValue> &m_plus_strand_melt_cache, 
	std::unordered_map<BindCacheKey, BindCacheValue> &m_minus_strand_melt_cache,
	const float &m_forward_primer_strand, 
	const float &m_reverse_primer_strand, const float &m_probe_strand, 
	const float &m_min_primer_tm, const float &m_max_primer_tm,
	const float &m_min_primer_dg, const float &m_max_primer_dg,
	const float &m_min_probe_tm, const float &m_max_probe_tm, 
	const float &m_min_probe_dg, const float &m_max_probe_dg, 
	const unsigned int &m_primer_clamp, 
	const int &m_min_primer_clamp, 
	const unsigned int &m_probe_clamp_5,
	const unsigned int &m_probe_clamp_3,
	const unsigned int &m_max_gap,
	const unsigned int &m_max_mismatch,
	const unsigned int &m_max_amplicon_len,
	const bool &m_single_primer_pcr);

//////////////////////////////////////////////////////////////////////////////////////////////
// In padlock_search.cpp
std::list<hybrid_sig> padlock(DNAHash &m_hash, const std::pair<std::string, SEQPTR> &m_seq, 
	const hybrid_sig &m_sig, NucCruc &m_melt,
	std::unordered_map<BindCacheKey, BindCacheValue> &m_plus_strand_melt_cache, 
	std::unordered_map<BindCacheKey, BindCacheValue> &m_minus_strand_melt_cache,
	const float &m_forward_primer_strand, 
	const float &m_reverse_primer_strand, 
	const float &m_min_primer_tm, const float &m_max_primer_tm,
	const float &m_min_primer_dg, const float &m_max_primer_dg,
	const unsigned int &m_probe_clamp_5, const unsigned int &m_probe_clamp_3, 
	const unsigned int &m_max_gap, const unsigned int &m_max_mismatch,
	const int &m_target_strand);

//////////////////////////////////////////////////////////////////////////////////////////////
// In probe_search.cpp
bool hybrid(
	const std::list<oligo_info> &m_bind_minus,
	const std::list<oligo_info> &m_bind_plus,
	hybrid_sig &m_ret, 
	const hybrid_sig &m_sig, 
	const int &m_amp_start, const int &m_amp_stop);

std::list<hybrid_sig> hybrid(DNAHash &m_hash, const std::pair<std::string, SEQPTR> &m_seq, 
	const hybrid_sig &m_sig, NucCruc &m_melt,
	const float &m_probe_strand, 
	const float &m_min_probe_tm, const float &m_max_probe_tm, 
	const float &m_min_probe_dg, const float &m_max_probe_dg, 
	const unsigned int &m_probe_clamp_5,
	const unsigned int &m_probe_clamp_3,
	const unsigned int &m_max_gap,
	const unsigned int &m_max_mismatch,
	const int &m_target_strand);
	

//////////////////////////////////////////////////////////////////////////////////////////////
// In bind_oligo.cpp

// The dangling ends of oligos are included in the match region
void bind_oligo_to_minus_strand(std::list<oligo_info> &info_list, 
		DNAHash &m_hash, SEQPTR m_seq, 
		const std::string &m_oligo,
		NucCruc &m_melt, std::unordered_map<BindCacheKey, BindCacheValue> &m_melt_cache,
		const float &m_min_oligo_tm, const float &m_max_oligo_tm,
		const float &m_min_oligo_dg, const float &m_max_oligo_dg,
		const unsigned int &m_clamp_5,
		const unsigned int &m_clamp_3,
		const unsigned int &m_max_gap,
		const unsigned int &m_max_mismatch);

void bind_oligo_to_plus_strand(std::list<oligo_info> &info_list, 
		DNAHash &m_hash, SEQPTR m_seq, 
		const std::string &m_oligo,
		NucCruc &m_melt, std::unordered_map<BindCacheKey, BindCacheValue> &m_melt_cache,
		const float &m_min_oligo_tm, const float &m_max_oligo_tm,
		const float &m_min_oligo_dg, const float &m_max_oligo_dg,
		const unsigned int &m_clamp_5,
		const unsigned int &m_clamp_3,
		const unsigned int &m_max_gap,
		const unsigned int &m_max_mismatch);

void bind_oligo_to_minus_strand(std::list<oligo_info> &info_list, 
		const unsigned char &m_oligo_mask, SEQPTR m_seq, 
		const std::string &m_oligo,
		NucCruc &m_melt, std::unordered_map<BindCacheKey, BindCacheValue> &m_melt_cache,
		const float &m_min_oligo_tm, const float &m_max_oligo_tm,
		const float &m_min_oligo_dg, const float &m_max_oligo_dg,
		const unsigned int &m_clamp_5,
		const unsigned int &m_clamp_3,
		const unsigned int &m_max_gap,
		const unsigned int &m_max_mismatch);

void bind_oligo_to_plus_strand(std::list<oligo_info> &info_list, 
		const unsigned char &m_oligo_mask, SEQPTR m_seq, 
		const std::string &m_oligo,
		NucCruc &m_melt, std::unordered_map<BindCacheKey, BindCacheValue> &m_melt_cache,
		const float &m_min_oligo_tm, const float &m_max_oligo_tm,
		const float &m_min_oligo_dg, const float &m_max_oligo_dg,
		const unsigned int &m_clamp_5,
		const unsigned int &m_clamp_3,
		const unsigned int &m_max_gap,
		const unsigned int &m_max_mismatch);
		
void match_oligo_to_minus_strand(std::list<oligo_info> &info_list, 
		DNAHash &m_hash, const std::string &m_oligo, const unsigned char &m_mask);
		
void match_oligo_to_plus_strand(std::list<oligo_info> &info_list, 
		DNAHash &m_hash, const std::string &m_oligo, const unsigned char &m_mask);
#endif // __AMPLICON
