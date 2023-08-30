// Nucleic Crucible
// J. D. Gans
// Los Alamos National Laboratory
// 6/9/2006

// version 2.1
// - integer DP matrix score for fast alignments
// - Added is_watson_and_crick(); does the computed alignment
//   contain one or more non-watson and crick base alignments?
// version 2.2
// - Fixed anchor computation and added target anchor functions
// version 2.3
// - Added debug code to dump parameter tables (#define MELT_DEBUG to turn on)
// - Added functions to output the internal query and target buffers
// - The following perl script can be used to generate paramter entries for the
//	nearest neighbor parameter sets:
// version 2.4
// - Changed the way dangling end virtual bases are handled:
//		Old version: pushed and popped on to and off of the query and target stacks
//		New version: Virtual bases are stored separately and manipulated with dedicated
//				   member functions.
// version 2.5
// - Added gapless alignment
// version 2.6
// - Moved NC_Elem into the NucCruc class (private).
// - MAX_SEQUENCE_LENGTH can now be specified via the makefile
// - Added near-neighbor dynamic programming (woot!). This is only approximate and
//   does not take into account:
//		1) Base pair dependent initiation terms
//		2) Length dependance of bulges and loops
//		3) Asymmetric loop term
// - Todo: Allow the user to specify a variable folding temperature (rather than the default T = 37.0C).
//
// version 2.7
// - Added a set_query_reverse_complement (similar to the set_target_reverse_complement)
// - Added boolean flags to indicate whether dangling base parameters have been incorporated into 
//   the alignment.
// - Added a #undef for _I, which is a predefined macro for OSX
//
// version 2.8 
// - Modified is_watson_and_crick() to require that 
//		1) all bases in the target and query are contained in
//		   the alignment
//		2) Gaps in the alignment cause the function to return false
// - The anchor functions now count from the ends of the query and
//   target oligos, not from the first aligned base of the query and
//   target oligos.
// - Added functions to test for exact matches at the 5' and 3 terminal bases
//   of the query and target oligos.
//
// version 2.9
// - Increased the sizes of MAX_LOOP_LENGTH, MAX_BULGE_LENGTH and MAX_HAIRPIN_LENGTH.
//
// version 3.0
// - Added functions to provide \delta G, \delta H and \delta S from the most recent Tm calculation.
// - Added a default constructor that defaults to SANTA_LUCIA melting parameters
//
// version 5.0
// - Almost complete re-write implementing the near neighbor dynamic programming algorithm of Leber et. al.
//
// version 5.1
// - Fixed the specification and definition of strand concentration!
//
// version 5.2
// - Added a function for computing fast (no gap, diagonal only) alignments.
// 
// version 5.3
// - Fixed the hairpin closure term. We now use the unpublished param_hairpin_terminal values from Santa Lucia 
//  (extracted from the unafold program) when computing hairpin structures, melting temperatures and other
//  parameters.
// - Added the parameters from the unafold tstackh.DHD file.
// - Fixed the estimation of the loop entropies. Now using interpolation between known values and 
//   Jacobson-Stockmayer entropy extrapolation for values greater than 30. Fixed the sign error in the
//   Jacobson-Stockmayer extrapolation implementation.
// - Added additional steps to the hairpin evaluation code to test a larger number of potential structures. This adds
//   computational complexity, but improves the agreement with unafold.
//
// version 5.4 (April 13, 2020)
//	- Added the ability to match target sequences that have degenerate nucleotides (finally)! We assume the most
//	  optimistic base pairing (i.e. that, when possible, the degenerate nucleotide is perfect match to the assay
//	  nucleotide).
//		-- Please note that we still require short, perfect, non-degenerate matches to "seed" longer alignments.
//		   If degenerate nucleotides are spaced close-enough together (i.e a long run of poly-N), then we will
//		   *not* initiate an alignment due to the absence of a perfect match seed k-mer to anchor an alignment.
//	- Removed some legacy code.
//	- Cleaned-up and stream-lined several parts of the code, including:
//		-- Handling conversions from different nucleotide encoding schemes.
//		-- Display of pairwise sequence alignments (alignments are no longer trimmed to remove dangling mismatches).
//	- Added an explicit version string in the NUC_CRUC_VERSION macro
//
// version 5.5 (August 4, 2023)
//	- Fixed test that was intended to remove spurious matches to target sequences with a large fraction of degerate bases.
//    This poorly implemented test had the side effect of removing matches to primers with more than 50% unalignable bases.
//    This test has been fixed and the `num_read_bases()` member function in `nuc_cruc.h` has been replaced by the 
//    `fraction_aligned_real_base_pairs()` member function.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __NUC_CRUC
#define __NUC_CRUC

#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <vector>
#include <deque>

#include "circle_buffer.h"

#define	NUC_CRUC_VERSION	"5.5"

// Define ENUMERATE_PATH to enable enumeration of multiple, equally high scoring
//  paths through the dynamic programming matrix
#define ENUMERATE_PATH

// Define UNAFOLD_COMPATIBILITY to better match the parameter values
// used by UNAFOLD.
#define UNAFOLD_COMPATIBILITY

// Physical constants
#define	NC_ZERO_C	(273.15f) 			// K
#define	NC_DEFAULT_T	(NC_ZERO_C + 37.0f)		// 37 C
#define	NC_R		(1.9872e-3f) 			// Kcal/(Mol . K) -- McQuarrie Stat. Mech.


// Compute the entropy from the free energy, enthalpy @ 37.0 C
#define	ENTROPY(DG, DH)		( ( (DH) - (DG) )/(310.15f) )

#define	BASE_PAIR(X, Y)		( (X)*NUM_BASE + (Y) )
#define	NON_VIRTUAL_BASE_PAIR(X)	( ( (X)%NUM_BASE < BASE::E ) && ( (X)/NUM_BASE < BASE::E ) )
#define	VIRTUAL_BASE_PAIR(X)	( ( (X)%NUM_BASE >= BASE::E ) || ( (X)/NUM_BASE >= BASE::E ) )
#define	HAS_GAP(X)			( ( (X)%NUM_BASE == BASE::GAP ) || ( (X)/NUM_BASE >= BASE::GAP ) )
#define	QUERY_BASE(X)			( (X)/NUM_BASE )
#define	TARGET_BASE(X)			( (X)%NUM_BASE )
#define	IS_REAL_BASE(X)			( (X) <= BASE::I )
#define	IS_DEGENERATE_BASE(X)	( (X) >= NUM_BASE )
#define	IS_VIRTUAL_BASE(X)	( ( (X) == BASE::E ) || ( (X) == BASE::GAP) )

// Trace back locations:
// 
// im1 = i - 1
// jm1 = j - 1
//
// im1_jm1	im1_j
// i_jm1		x

#define	im1_jm1		(1 << 0)	// match or mismatch
#define	im1_j		(1 << 1)	// gap the query (i.e. target matches gap)
#define	i_jm1		(1 << 2)	// gap the target (i.e. query matches gap)
#define	invalid_trace	(1 << 3)

// Match states
#define	query_target	im1_jm1
#define	query_gap		im1_j
#define	gap_target	i_jm1
#define	invalid_match	invalid_trace

#define PATH_SPLIT(X)	( ( ( (X) & im1_jm1 ) + ( ( (X) & im1_j ) >> 1) + ( ( (X) & i_jm1 ) >> 2) ) > 1 )

// Allow the user to override MAX_SEQUENCE_LENGTH in a makefile
#ifndef MAX_SEQUENCE_LENGTH

// For hardware efficiency, MAX_SEQUENCE_LENGTH should be
// an integer power of 2
#define	MAX_SEQUENCE_LENGTH		1024

#endif // MAX_SEQUENCE_LENGTH

// old version NUM_BASE_PAIR = |{A, C, G, T, E, GAP}|*|{A, C, G, T, E, GAP}| = 36
// NUM_BASE_PAIR = |{A, C, G, T, I, E, GAP}|*|{A, C, G, T, I, E, GAP}| = 49
#define	NUM_BASE				7 // The number of real (5) and virtual (2) bases
#define	NUM_BASE_PAIR			(NUM_BASE*NUM_BASE)

// The orginial maximum loop, buldge and hairpin lengths were 30.
// The new lengths are set from the MAX_SEQUENCE_LENGTH but must
// always be >= 30.
#define MAX_LOOP_LENGTH		(MAX_SEQUENCE_LENGTH/2)
#define	MAX_BULGE_LENGTH		(MAX_SEQUENCE_LENGTH/2)
#define	MAX_HAIRPIN_LENGTH		(MAX_SEQUENCE_LENGTH/2)

#if (MAX_LOOP_LENGTH < 30)
#define 	MAX_LOOP_LENGTH	30
#endif

#if (MAX_BULGE_LENGTH < 30)
#define 	MAX_BULGE_LENGTH	30
#endif

#if (MAX_HAIRPIN_LENGTH < 30)
#define 	MAX_HAIRPIN_LENGTH	30
#endif

// Include code below if we need to protect our constants
// OS X has macros named _A, _C, _I and _G
#ifdef _A
#undef _A
#endif

#ifdef _C
#undef _C
#endif

#ifdef _G
#undef _G
#endif

#ifdef _I
#undef _I
#endif

// Windows has a macro named _T
#ifdef _T
#undef _T
#endif

// Older versions of the Visual Studio C++ compiler do
// not povide std::min or std::max. _MSC_VER < 1300 should
// capture all of these older version (like 6.0).
#if defined(WIN32) && (_MSC_VER < 1300)

#ifdef min
#undef min
#endif 

template <class A>
inline A min(const A &m_a, const A &m_b)
{
	return (m_a < m_b) ? m_a : m_b;
}

#ifdef max
#undef max
#endif

template <class A>
inline A max(const A &m_a, const A &m_b)
{
	return (m_a > m_b) ? m_a : m_b;
}

#endif // WIN32

/////////////////////////////////////////////////////////
// By default, use an integer based alignment
#ifdef FLOATING_POINT_ALIGNMENT
	typedef float NC_Score;

	// We're using float's -- no scaling is required!
	#define	NC_SCORE_SCALE(X)		(X)
	#define	NC_INV_SCORE_SCALE(X)	(X)
#else // Integer alignment
	typedef int NC_Score;

	// Rescale the floating point delta G values to store
	// as integers
	#define	NC_SCORE_SCALE(X)		( NC_Score( (X)*10000.0f ) )
	#define	NC_INV_SCORE_SCALE(X)	( float(X)/10000.0f )
#endif // FLOATING_POINT_ALIGNMENT

// Allowed bases, no 'N' for now
// These values have been chosen to
// match nuc_cruc.h and hash_db.h (up to 'T').
//	A = adenine
//	C = cytosine
//	G = guanine
//	T = thymine
//	I = inosine
//	E = dangling end virtual base
//	GAP = no base
namespace BASE {

	typedef enum {
		// Real bases
		A, C, G, T, I, 
		// "Virtual bases"
		E, GAP,
		// IUPAC degenerate bases
		M, R, S, V, W, Y, H, K, D, B, N
		} nucleic_acid;

	inline nucleic_acid char_to_nucleic_acid(char m_base)
	{
		switch(m_base){
			case 'A': case 'a':
				return BASE::A;
			case 'T': case 't':
				return BASE::T;
			case 'G': case 'g':
				return BASE::G;
			case 'C': case 'c':
				return BASE::C;
			case 'I': case 'i':
				return BASE::I;
			// IUPAC degenerate bases
			case 'M': case 'm':// A or C
				return BASE::M;
			case 'R': case 'r': // G or A
				return BASE::R;
			case 'S': case 's': // G or C
				return BASE::S;
			case 'V': case 'v': // G or C or A
				return BASE::V;
			case 'W': case 'w': // A or T
				return BASE::W;
			case 'Y': case 'y': // T or C
				return BASE::Y;
			case 'H': case 'h': // A or C or T
				return BASE::H;
			case 'K': case 'k': // G or T
				return BASE::K;
			case 'D': case 'd': // G or A or T
				return BASE::D;
			case 'B': case 'b': // G or T or C
				return BASE::B;
			case 'N': case 'n': // A or T or G or C
				return BASE::N;
			default:
				throw __FILE__ ":char_to_nucleic_acid: Illegal base";
		};

		return BASE::N; // We should never get here!
	}

	inline nucleic_acid char_to_complement_nucleic_acid(char m_base)
	{
		switch(m_base){
			case 'A': case 'a':
				return BASE::T;
			case 'T': case 't':
				return BASE::A;
			case 'G': case 'g':
				return BASE::C;
			case 'C': case 'c':
				return BASE::G;
			case 'I': case 'i': // Inosine is defined to be it's own complement
				return BASE::I;
			// IUPAC degenerate bases
			case 'M': case 'm':// A or C
				return BASE::K;
			case 'R': case 'r': // G or A
				return BASE::Y;
			case 'S': case 's':// G or C
				return BASE::S;
			case 'V': case 'v': // G or C or A
				return BASE::B;
			case 'W': case 'w': // A or T
				return BASE::W;
			case 'Y': case 'y': // T or C
				return BASE::R;
			case 'H': case 'h': // A or C or T
				return BASE::D;
			case 'K': case 'k': // G or T
				return BASE::M;
			case 'D': case 'd': // G or A or T
				return BASE::H;
			case 'B': case 'b': // G or T or C
				return BASE::V;
			case 'N': case 'n': // A or T or G or C
				return BASE::N;
			default:
				throw __FILE__ ":char_to_complement_nucleic_acid: Illegal base";
		};

		return BASE::N; // We should never get here!
	}

	bool is_complemetary_base(const nucleic_acid &m_query, const nucleic_acid &m_target);
}

class trace_branch{

	private:
		unsigned char *mask_ptr;
		unsigned char curr_trace;
	
	public:
	
		trace_branch(unsigned char &m_trace)
		{
			mask_ptr = &m_trace;
			
			if(mask_ptr != NULL){
			
				if(*mask_ptr & im1_jm1){
					curr_trace = im1_jm1;
				}
				else{
					
					if(*mask_ptr & im1_j){
						curr_trace = im1_j;
					}
					else{
						curr_trace = i_jm1;
					}
				}
			}
		};

		unsigned char trace() const
		{
			return curr_trace;
		};
		
		inline bool next_trace()
		{
			if(mask_ptr == NULL){
				throw "mask_ptr == NULL";
			}
			
			while( (curr_trace = curr_trace << 1) < invalid_trace){
				
				// Is this trace in the mask?
				if(curr_trace & *mask_ptr){
				
					return true;
				}
			}
			
			// We have exhausted the mask
			return false;
		};
		
		inline bool operator==(const trace_branch &m_tb) const
		{
			return (mask_ptr == m_tb.mask_ptr);
		};
		
		inline bool operator==(const unsigned char &m_ptr) const
		{
			return (mask_ptr == &m_ptr);
		};
};

struct alignment {
	
	alignment()
	{
		valid = false;
		
		dH = 0.0f;
		dS = 0.0f;
		tm = 0.0f;
		dp_dg = 0.0f;
	};
	
	~alignment()
	{
		// Do nothing
	};
	
	void clear()
	{
		valid = false;
		
		dH = 0.0f;
		dS = 0.0f;
		tm = 0.0f;
		dp_dg = 0.0f;
		
		query_align.clear();
		target_align.clear();		
	};
	
	bool valid;
	
	float dH;
	float dS;
	float tm;
	
	// The value of the high scoring cell in the dynamic programming 
	// matrix (times -1 and scaled by a constant to convert to a float)
	float dp_dg;
	
	std::deque<BASE::nucleic_acid> query_align;
	std::deque<BASE::nucleic_acid> target_align;
	
	std::pair<int, int> first_match;
	std::pair<int, int> last_match;
	
	inline unsigned int num_gap() const
	{
		std::deque<BASE::nucleic_acid>::const_iterator iter;
		
		unsigned int gap_count = 0;
		
		// Counting the number of gaps in the aligned query and target sequences
		// assumes that a pair of gaps never occurs.
		for(iter = query_align.begin();iter != query_align.end();iter++){
		
			gap_count += (*iter == BASE::GAP) ? 1 : 0;
		}
		
		for(iter = target_align.begin();iter != target_align.end();iter++){
		
			gap_count += (*iter == BASE::GAP) ? 1 : 0;
		}
		
		return gap_count;
	};
	
	inline unsigned int num_mismatch_by_query(const unsigned int &m_query_len) const
	{
		std::deque<BASE::nucleic_acid>::const_iterator q = query_align.begin();
		std::deque<BASE::nucleic_acid>::const_iterator t = target_align.begin();
		
		unsigned int mismatch_count = 0;
		unsigned int num_aligned_query_bases = 0;
		
		while( q != query_align.end() ){

			if( !IS_VIRTUAL_BASE(*q) ){

				if( !IS_VIRTUAL_BASE(*t) && !is_complemetary_base(*q, *t) ){
					++mismatch_count;
				}
				
				++num_aligned_query_bases;
			}

			++q;
			++t;
		}
		
		if(m_query_len < num_aligned_query_bases){
			throw __FILE__ ":num_mismatch_by_query: m_query_len < num_aligned_query_bases";
		}
		
		mismatch_count += m_query_len - num_aligned_query_bases;
		
		return mismatch_count;
	};

	// Count the fraction of real base pairs in the alignment
	inline unsigned int fraction_aligned_real_base_pairs() const
	{
		std::deque<BASE::nucleic_acid>::const_iterator q = query_align.begin();
		std::deque<BASE::nucleic_acid>::const_iterator t = target_align.begin();
		
		unsigned int num_real = 0;
		unsigned int num_aligned = 0;

		while( q != query_align.end() ){

			if( IS_REAL_BASE(*q) && IS_REAL_BASE(*t) ){
				++num_real;
			}

			++num_aligned;
			++q;
			++t;
		}
		
		return (num_aligned == 0) ? 0.0f : float(num_real)/num_aligned;
	};
};

class NucCruc{

	public:

		// Allowed parameter sets		
		typedef enum {SANTA_LUCIA = 0} parameter_set;

		// Allowed comparisons
		typedef enum {HOMO_DIMER = 0, HETERO_DIMER, HAIRPIN, INVALID} mode;
		
		enum {
			// A pair of mismatches i.e.:
			// 5' AT 3'
			// 3' CG 5'
			LOOP_H = 0, LOOP_S, 
			
			// A pair of gaps i.e.:
			// 5' AT 3'
			// 3' -- 5'
			BULGE_H, BULGE_S,
			
			// A matched pair (AT or GC) with a gap i.e.:
			// 5' AT 3' or 5' -C 3'
			// 3' T- 5'    3' TG 5'
			TERMINAL_MATCH_AT_H, TERMINAL_MATCH_AT_S,
			TERMINAL_MATCH_GC_H, TERMINAL_MATCH_GC_S,
			TERMINAL_MATCH_I_H, TERMINAL_MATCH_I_S,
			
			// A mimatched pair with a gap i.e.:
			// 5' GT 3' or 5' -C 3'
			// 3' G- 5'    3' TC 5'
			TERMINAL_MISMATCH_H, TERMINAL_MISMATCH_S,
			
			NUM_SUPP_PARAM};
		
		enum {
			LOOP_SALT,
			BULGE_SALT,
			TERMINAL_MATCH_SALT,
			TERMINAL_MISMATCH_SALT,			
			NUM_SALT_PARAM};
	private:
		
		struct NC_Elem{

			NC_Elem()
			{
				M = I_query = I_target = -1;
				
				M_trace = I_query_trace = I_target_trace = invalid_trace;
			}

			~NC_Elem()
			{
				// Do nothing
			};

			// Use the notation of "Biological sequence analysis" by Durbin, Eddy, Krogh and Mitchison
			NC_Score M;
			NC_Score I_query;  // insertion in query
			NC_Score I_target; // insertion in target

			// The trace-back pointers
			unsigned char M_trace;
			unsigned char I_query_trace;
			unsigned char I_target_trace;
		};
		
		// Enumerate all possible nearest-neighboor base stacking arrangements
		enum {
			AA = 0, AC, AG, AT, AI, AE, A_, 
			CA, CC, CG, CT, CI, CE, C_, 
			GA, GC, GG, GT, GI, GE, G_, 
			TA, TC, TG, TT, TI, TE, T_, 
			IA, IC, IG, IT, II, IE, I_,
			EA, EC, EG, ET, EI, EE, E_,
			_A, _C, _G, _T, _I, _E, __
		};

		// Enumerate all 3 and 4 base "special" hairpin loops
		// Supplemental Material: Annu.Rev.Biophs.Biomol.Struct.33:415-40
		// doi: 10.1146/annurev.biophys.32.110601.141800
		// The Termodynamicso f DNA Structural Motifs
		// SantaLucia and Hicks, 2004

		enum {AAAAAT = 0, AAAACT, AAACAT, ACTTGT, AGAAAT, AGAAT, AGAGAT, AGATAT,
			AGCAAT, AGCAT, AGCGAT, AGCTTT, AGGAAT, AGGAT, AGGGAT, AGGGGT, AGTAAT,
			AGTAT, AGTGAT, AGTTCT, ATTCGT, ATTTGT, ATTTTT, CAAAAG, CAAACG,
			CAACAG, CAACCG, CCTTGG, CGAAAG, CGAAG, CGAGAG, CGATAG, CGCAAG,
			CGCAG, CGCGAG, CGCTTG, CGGAAG, CGGAG, CGGGAG, CGGGGG, CGTAAG, CGTAG,
			CGTGAG, CGTTCG, CTTCGG, CTTTGG, CTTTTG, GAAAAC, GAAAAT, GAAACC,
			GAAACT, GAACAC, GAACAT, GCTTGC, GCTTGT, GGAAAC, GGAAAT, GGAAC, 
			GGAGAC, GGAGAT, GGATAC, GGATAT, GGCAAC, GGCAAT, GGCAC, GGCGAC,
			GGCGAT, GGCTTC, GGCTTT, GGGAAC, GGGAAT, GGGAC, GGGGAC, GGGGAT,
			GGGGGC, GGGGGT, GGTAAC, GGTAAT, GGTAC, GGTGAC, GGTGAT, GGTTCC,
			GTATAT, GTTCGC, GTTCGT, GTTTGC, GTTTGT, GTTTTC, GTTTTT, TAAAAA,
			TAAAAG, TAAACA, TAAACG, TAACAA, TAACAG, TCTTGA, TCTTGG, TGAAA,
			TGAAAA, TGAAAG, TGAGAA, TGAGAG, TGATAA, TGATAG, TGCAA, TGCAAA,
			TGCAAG, TGCGAA, TGCGAG, TGCTTA, TGCTTG, TGGAA, TGGAAA, TGGAAG,
			TGGGAA, TGGGAG, TGGGGA, TGGGGG, TGTAA, TGTAAA, TGTAAG, TGTGAA,
			TGTGAG, TGTTCA, TTTCGA, TTTCGG, TTTTAG, TTTTGA, TTTTGG, TTTTTA,
			TTTTTG, NUM_SPECIAL_HAIRPIN_LOOP
		};
			
		// Save the last type of alignment performed (for printing)
		mode tm_mode;
		parameter_set tm_param;
		
		// The effective temperature used to align sequences
		float target_T;
		
		// Use the Dinkelbach algorithm to compute Tm?
		bool use_dinkelbach;
		
		// If false, compute standard alignments. If false, only compute alignments along the
		// diagonal of the dynamic programming matrix (no gaps allowed).
		bool diagonal_alignment;
		
		// Pre-compute the delta G values corresponding to the current temperature
		// and parameter set
		NC_Score delta_g[NUM_BASE_PAIR][NUM_BASE_PAIR];
		
		float param_H[NUM_BASE_PAIR][NUM_BASE_PAIR];
		float param_S[NUM_BASE_PAIR][NUM_BASE_PAIR];
		
		float param_init_H;
		float param_init_S;
		
		float param_loop_S[MAX_LOOP_LENGTH + 1];
		float param_asymmetric_loop_dS;
		
		float param_loop_terminal_H[NUM_BASE_PAIR][NUM_BASE_PAIR];
		float param_loop_terminal_S[NUM_BASE_PAIR][NUM_BASE_PAIR];
		
		float param_hairpin_terminal_H[NUM_BASE_PAIR][NUM_BASE_PAIR];
		float param_hairpin_terminal_S[NUM_BASE_PAIR][NUM_BASE_PAIR];
		
		float param_bulge_S[MAX_BULGE_LENGTH + 1];
		float param_bulge_AT_closing_S;
		
		float param_hairpin_S[MAX_HAIRPIN_LENGTH + 1];
		
		float param_hairpin_special_H[NUM_SPECIAL_HAIRPIN_LOOP];
		float param_hairpin_special_S[NUM_SPECIAL_HAIRPIN_LOOP];
		
		float param_AT_closing_H;
		float param_AT_closing_S;
		
		float param_symmetry_S;
		
		float param_SALT;
		
		// Supplementary parameters are needed to completely specify
		// the dynamic programming solution
		float param_supp[NUM_SUPP_PARAM];
		
		// The salt correction for supplementary parameters
		float param_supp_salt[NUM_SALT_PARAM];
		
		bool watson_and_crick[NUM_BASE_PAIR];
		
		float na_concentration;
		float strand_concentration;
		
		NC_Elem *dp_matrix;
		
		// The input sequences in 5'-3' orientation
		CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH> query;
		CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH> target;
		
		// The current alignment
		alignment curr_align;
		
		// The maximum number of equally high scoring paths to explore
		// through the dynamic programming matrix. If max_dp_path_enum == 0,
		// then *all* high scoring paths will be explored.
		unsigned int max_dp_path_enum;
		
		// Should we allow dangling bases?
		// enable_dangle.first ==  5' query side of alignment
		// enable_dangle.second ==  3' query side of alignment
		std::pair<bool, bool> enable_dangle;
						
		// A scratch buffer for fast target sequence access
		BASE::nucleic_acid target_buffer[MAX_SEQUENCE_LENGTH];
		
		// Pointers to the highest scoring element(s) in the DP matrix
		// (according to gprof, it is faster to use a vector than a deque).
		std::vector<NC_Elem*> max_ptr;
		
		int find_loop_index(const CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q,
			const unsigned int &m_start, const unsigned int &m_len);
		
		// Align the query sequence to iteself to check for potential hairpin
		// secondary structures.
		NC_Score align_hairpin(const CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q);
		
		NC_Score align_dimer(const CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q,
			const CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH> &m_t);
		
		NC_Score align_dimer_diagonal(const CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q,
			const CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH> &m_t);
			
		void enumerate_dimer_alignments(NC_Elem *m_dp_matrix, 
			NC_Elem *m_max_ptr,
			alignment &m_best_align,
			const CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q,
			const CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH> &m_t,
			const mode &m_mode);
		void enumerate_hairpin_alignments(NC_Elem *m_dp_matrix, 
			NC_Elem *m_max_ptr, 
			alignment &m_best_align, 
			const CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q);
			
		void trace_back(NC_Elem *m_dp_matrix, NC_Elem *m_cell_ptr, 
			std::deque<trace_branch> &m_trace_stack, int &m_zero_count, 
			alignment &m_local_align,
			const CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q,
			const CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH> &m_t);
		
		bool evaluate_alignment(alignment &local_align, const mode &m_mode);
		bool evaluate_hairpin_alignment(alignment &local_align);
		
		// Given a query and target sequence, compute an alignment
		inline NC_Score align_heterodimer()
		{
			tm_mode = HETERO_DIMER;
			
			return ( diagonal_alignment ? align_dimer_diagonal(query, target) : align_dimer(query, target) );
		};
		
		// Given a query sequence, compute an alignment of the query to itself
		inline NC_Score align_homodimer()
		{
			tm_mode = HOMO_DIMER;
			
			return ( diagonal_alignment ? align_dimer_diagonal(query, query) : align_dimer(query, query) );
		};
		
		// Given a query and target sequence, compute the melting temperature
		float tm_dimer(const CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q,
			const CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH> &m_t,
			const mode &m_mode);
		
		// Different parameter sets:
		// The hybridization parameters from the SantaLucia lab
		void init_param_Santa_Lucia();
		
		// Recompute the dynamic programming parameters
		void update_dp_param();
			
		bool is_internal_to_loop( std::deque<BASE::nucleic_acid>::const_iterator m_q, 
			const std::deque<BASE::nucleic_acid>::const_iterator &m_query_end, 
			std::deque<BASE::nucleic_acid>::const_iterator m_t, 
			const std::deque<BASE::nucleic_acid>::const_iterator &m_target_end );
		
		bool has_AT_initiation( std::deque<BASE::nucleic_acid>::const_iterator m_q, 
			const std::deque<BASE::nucleic_acid>::const_iterator &m_query_begin, 
			std::deque<BASE::nucleic_acid>::const_iterator m_t, 
			const std::deque<BASE::nucleic_acid>::const_iterator &m_target_begin );
	public:
		
		NucCruc(const unsigned int &m_param_set = SANTA_LUCIA,
			const float &m_t = NC_DEFAULT_T);
			
		~NucCruc()
		{
			if(dp_matrix != NULL){
				delete [] dp_matrix;
			}			
		};
		
		inline bool dinkelbach() const
		{
			return use_dinkelbach;
		};
		
		inline void dinkelbach(const bool &m_dink)
		{
			use_dinkelbach = m_dink;
		};
		
		// Given a query and target sequence, compute the melting temperature
		float approximate_tm_heterodimer();
		
		// Given a query, compute the melting temperature of the query to itself
		// (as a homodimer)
		float approximate_tm_homodimer();
		
		// Given *aligned* query and target sequences, compute the melting temperature
		inline float tm_from_align(const std::deque<BASE::nucleic_acid> &m_query_align,
			const std::deque<BASE::nucleic_acid> &m_target_align)
		{
			
			if( m_query_align.size() != m_target_align.size() ){
				throw "query target size mismatch";
			}
			
			curr_align.clear();
			
			curr_align.query_align = m_query_align;
			curr_align.target_align = m_target_align;
	
			evaluate_alignment(curr_align, HETERO_DIMER);
			
			return curr_align.tm;
		};
		
		// Given a query sequence, compute the perfect match duplex melting temperature
		inline float tm_pm_duplex(const std::string &m_query)
		{
			
			curr_align.clear();
			
			for(std::string::const_iterator i = m_query.begin();i != m_query.end();i++){
			
				switch(*i){
					case 'A':
					case 'a':
						curr_align.query_align.push_back(BASE::A);
						curr_align.target_align.push_back(BASE::T);
						break;
					case 'T':
					case 't':
						curr_align.query_align.push_back(BASE::T);
						curr_align.target_align.push_back(BASE::A);
						break;
					case 'G':
					case 'g':
						curr_align.query_align.push_back(BASE::G);
						curr_align.target_align.push_back(BASE::C);
						break;
					case 'C':
					case 'c':
						curr_align.query_align.push_back(BASE::C);
						curr_align.target_align.push_back(BASE::G);
						break;
					default:
						throw "Unknown base in tm_pm_duplex";
				};
			}
			
			evaluate_alignment(curr_align, HETERO_DIMER);
			
			return curr_align.tm;
		};
		
		// Given a query, compute the melting temperature of the query to itself
		// (as a hairpin).
		float approximate_tm_hairpin();
		
		inline void fast_alignment(const bool &m_fast_align)
		{
			diagonal_alignment = m_fast_align;
		};
		
		inline bool fast_alignment()
		{
			return diagonal_alignment;
		};
		
		inline float salt() const
		{
			return na_concentration;
		};

		inline void salt(const float &m_na_concentration)
		{
			if(m_na_concentration < 1.0e-6f){
				throw ":salt: [Na+] < 1.0e-6f";
			}
			
			if(m_na_concentration > 1.0f){
				throw ":salt: [Na+] > 1.0f";
			}
			
			na_concentration = m_na_concentration;
			
			// Recompute the dynamic programming parameters (which depend
			// on the salt concentration).
			update_dp_param();
		};
		
		inline void Salt(const float &m_na_concentration)
		{
			salt(m_na_concentration);
		};
		
		inline float strand() const
		{
			return strand_concentration;
		};

		inline void strand(const float &m_strand_concentration)
		{
			if(m_strand_concentration < 0.0f){
				throw ":strand: strand_concentration < 0.0f";
			}
			
			strand_concentration = m_strand_concentration;
		};
		
		// Set the total strand concentration from the strand
		// concentrations of each strand
		inline void strand(const float &m_c_a, const float &m_c_b)
		{
			if(m_c_a < 0.0f){
				throw ":strand: m_c_a < 0.0f";
			}
			
			if(m_c_b < 0.0f){
				throw ":strand: m_c_b < 0.0f";
			}
			
			// For A + B -> D,
			// Ct = C_excess - C_limit/2
			// If C_excess == C_limit, Ct = C_excess/2
			// If C_excess >> C_limit, Ct = C_excess
			if(m_c_a > m_c_b){
				strand_concentration = m_c_a - 0.5f*m_c_b;
			}
			else{
				strand_concentration = m_c_b - 0.5f*m_c_a;
			}
		};
		
		inline void Strand(const float &m_strand_concentration)
		{
			strand(m_strand_concentration);
		};
		
		inline void Strand(const float &m_c_a, const float &m_c_b)
		{
			strand(m_c_a, m_c_b);
		};
		
		inline void clear()
		{
			query.clear();
			target.clear();
						
			// Invalidate any existing alignment
			tm_mode = INVALID;
		};
		
		inline void clear_query()
		{
			query.clear();
						
			// Invalidate any existing alignment
			tm_mode = INVALID;
		};
		
		inline void clear_target()
		{
			target.clear();
			
			// Invalidate any existing alignment
			tm_mode = INVALID;
		};
		
		inline void set_query(const std::string &m_query)
		{
			const size_t query_len = m_query.size();
			
			if(query_len > MAX_SEQUENCE_LENGTH){
				throw __FILE__ ":set_query: Query size out of bounds";
			}
			
			query.clear();
			
			for(unsigned int i = 0;i < query_len;i++){
				query.push_back( BASE::char_to_nucleic_acid( m_query[i]) );
			}
		};
		
		inline void set_query_reverse_complement(const std::string &m_query)
		{
			const size_t query_len = m_query.size();
			
			if(query_len > MAX_SEQUENCE_LENGTH){
				throw __FILE__ ":set_query_reverse_complement: Query size out of bounds";
			}
			
			query.clear();
			
			for(unsigned int i = 0;i < query_len;i++){
				query.push_front( BASE::char_to_complement_nucleic_acid(m_query[i]) );
			}
		};
		
		inline void set_target(const std::string &m_target)
		{
			const size_t target_len = m_target.size();
			
			if(target_len > MAX_SEQUENCE_LENGTH){
				throw __FILE__ ":set_target: Target size out of bounds";
			}
			
			target.clear();
			
			for(unsigned int i = 0;i < target_len;i++){
				target.push_back( BASE::char_to_nucleic_acid(m_target[i]) );
			}
		};
		
		inline void set_target_reverse_complement(const std::string &m_target)
		{
			const size_t target_len = m_target.size();
			
			if(target_len > MAX_SEQUENCE_LENGTH){
				throw __FILE__ ":set_target_reverse_complement: Target size out of bounds";
			}
			
			target.clear();
			
			for(unsigned int i = 0;i < target_len;i++){
				target.push_front( BASE::char_to_complement_nucleic_acid(m_target[i]) );
			}
		};
		
		// Set query to be m_seq and target to be the reverse complement
		// of m_seq
		inline void set_duplex(const std::string &m_seq)
		{
			const size_t seq_len = m_seq.size();
			
			if(seq_len > MAX_SEQUENCE_LENGTH){
				throw __FILE__ ":set_duplex: m_seq size out of bounds";
			}
			
			query.clear();
			target.clear();
						
			for(unsigned int i = 0;i < seq_len;i++){
			
				query.push_back( BASE::char_to_nucleic_acid(m_seq[i]) );
				target.push_front( BASE::char_to_complement_nucleic_acid(m_seq[i]) );
			}
		};
		
		inline void push_back_query(const BASE::nucleic_acid &m_base)
		{
			// Since we changed the method from adding dangling bases
			// test for any illegal attempts to add a dangling base here
			// (this check can be removed after all dependent code has been corrected)
			if(m_base == BASE::E){
				throw __FILE__ ":push_back_query: Illegal push of '$'";
			}
			
			query.push_back(m_base);
		};
		
		inline void push_back_query(const char &m_base)
		{
			query.push_back( BASE::char_to_nucleic_acid(m_base) );
		};
		
		inline void pop_back_query()
		{
			query.pop_back();
		};
		
		inline void push_front_query(const BASE::nucleic_acid &m_base)
		{
			// Since we changed the method from adding dangling bases
			// test for any illegal attempts to add a dangling base here
			// (this check can be removed after all dependent code has been corrected)
			if(m_base == BASE::E){
				throw __FILE__ ":push_front_query: Illegal push of '$'";
			}
			
			query.push_front(m_base);
		};
		
		inline void push_front_query(const char &m_base)
		{
			query.push_front( BASE::char_to_nucleic_acid(m_base) );
		};
		
		inline void pop_front_query()
		{
			query.pop_front();
		};
		
		inline void push_back_target(const BASE::nucleic_acid &m_base)
		{
			// Since we changed the method from adding dangling bases
			// test for any illegal attempts to add a dangling base here
			// (this check can be removed after all dependent code has been corrected)
			if(m_base == BASE::E){
				throw __FILE__ ":push_back_target: Illegal push of '$'";
			}
			
			target.push_back(m_base);
		};
		
		inline void push_back_target(const char &m_base)
		{
			target.push_back( BASE::char_to_nucleic_acid(m_base) );
		};
		
		inline void pop_back_target()
		{
			target.pop_back();
		};
		
		inline void push_front_target(const BASE::nucleic_acid &m_base)
		{
			// Since we changed the method from adding dangling bases
			// test for any illegal attempts to add a dangling base here
			// (this check can be removed after all dependent code has been corrected)
			if(m_base == BASE::E){
				throw __FILE__ ":push_front_target: Illegal push of '$'";
			}
			
			target.push_front(m_base);
		};
		
		inline void push_front_target(const char &m_base)
		{
			target.push_front( BASE::char_to_nucleic_acid(m_base) );
		};
		
		inline void pop_front_target()
		{
			target.pop_front();
		};
		
		inline unsigned int size_query() const
		{
			return query.size();
		};
		
		inline unsigned int size_target() const
		{
			return target.size();
		};
			
		// The number of 5' query bases that exactly complement the
		// target. Note that the alignment must already be computed!
		unsigned int anchor5_query() const;
		unsigned int anchor5_target() const;
		
		// The number of 3' query bases that exactly complement the
		// target. Note that the alignment must already be computed!
		unsigned int anchor3_query() const;
		unsigned int anchor3_target() const;
		
		// The coordinates of the first and last aligned base in the query
		std::pair<unsigned int, unsigned int> alignment_range_query() const;
		
		// The coordinates of the first and last aligned base in the query
		std::pair<unsigned int, unsigned int> alignment_range_target() const;
		
		void alignment_range(std::pair<unsigned int, unsigned int> &m_query_range,
			std::pair<unsigned int, unsigned int> &m_target_range) const;

		// Does the 5'-most base of the query have an exact match to the
		// target?
		bool match_terminal5_query() const;
		
		// Does the 3'-most base of the query have an exact match to the
		// target?
		bool match_terminal3_query() const;
		
		// Does the 5'-most base of the target have an exact match to the
		// query?
		bool match_terminal5_target() const;
			
		// Does the 3'-most base of the target have an exact match to the
		// query?
		bool match_terminal3_target() const;
			
		// Does the alignment contain any non-Watson and Crick base pairs?
		// If so, this function returns false. Gaps do *not* effect the return value!
		bool is_watson_and_crick() const;
		
		inline unsigned int num_gap() const
		{
			return curr_align.num_gap();
		};
		
		inline unsigned int num_mismatch() const
		{
			return curr_align.num_mismatch_by_query( query.size() );
		};
		
		inline unsigned int fraction_aligned_real_base_pairs() const
		{
			return curr_align.fraction_aligned_real_base_pairs();
		};

		inline void dangle(const bool &m_dangle_5, const bool &m_dangle_3)
		{
			enable_dangle = std::make_pair(m_dangle_5, m_dangle_3);
		};
		
		inline size_t alignment_size() const
		{
			return curr_align.query_align.size();
		};
		
		// Print an alignment (for debugging purposes)
		friend std::ostream& operator << (std::ostream &s, const NucCruc &m_melt);
		
		// Return the query and target sequences
		std::string query_seq() const;
		std::string target_seq() const;
		
		// Return the change in energy on binding for the most recent melting temperature
		// calculation
		inline float delta_H() const
		{
			return curr_align.dH;
		};
		
		// Return the change in entropy on binding for the most recent melting temperature
		// calculation
		inline float delta_S() const
		{
			return curr_align.dS;
		};
		
		// Return the change in free energy on binding for the most recent melting temperature
		// calculation
		inline float delta_G() const
		{
			return curr_align.dH - target_T*curr_align.dS;
		};
		
		// Return the change in free energy on binding for the most recent melting temperature
		// calculation as computing by the best scoring cell found in the dynamic programming
		// matrix
		float delta_G_dp() const;

		inline float temperature() const
		{
			return target_T;
		};

		inline void temperature(const float &m_tm)
		{
			if(m_tm < 0.0f){
				throw "temperature: tm < 0";
			}

			target_T = m_tm;
			
			// Recompute the dynamic programming parameters
			update_dp_param();
		};
		
		// These functions are for optimizing the dynamic programming
		// parameters
		void set_supp_param(double m_param[NUM_SUPP_PARAM]);
		void set_supp_salt_param(double m_param[NUM_SALT_PARAM]);
		
		#ifdef __DEBUG
		void dump_tables();
		#endif // __DEBUG
};

// For debugging
std::string print_alignment(std::deque<BASE::nucleic_acid> &m_query_align, 
	std::deque<BASE::nucleic_acid> &m_target_align);
	
#endif // __NUC_CRUC

