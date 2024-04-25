#ifndef __SEQ_ALIGN
#define __SEQ_ALIGN

#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <vector>
#include <deque>
#include <algorithm>

#define	SEQ_ALIGN_VERSION	"1.1 (sse/avx)"

// Change log
//
// 1.0 (June 12, 2019)
//	* Added arguments to the pack_query and pack_target functions to
//	  enable reverse complementing the input sequence.
//	* Added "Overlap" and "Global" sequence alignments.
//	* Added a reverse_complement function for sequence alignments.
//	* Fixed the calculation of max_loc in align_smith_waterman to
//	  only accept new max_loc values if the alignment score is
// 	  greater than the old score (as opposed to greater or equal to).
//	  The old ('>') version was allowing masked 'N' characters to be
// 	  appended to the sequence alignment (when 'N' is treated as a masking
// 	  character, it is assigned a score of zero).
// 1.1 (June 19, 2019)
//	* Added support for AVX2-based SIMD alignment.

/////////////////////////////////////////////////////////////////////////////////////////////
// SSE & AVX instruction defines, macros and datatypes. Note that using SSE4.1 requires
// g++ version >= 4.4 and the compiler flag -msse4.1. Using AVX2 requires the compiler 
// flag -mavx2.
/////////////////////////////////////////////////////////////////////////////////////////////
#include <immintrin.h>
	
#ifdef ATA32x8

	// Using __m256i requires memory aligned on 32 byte boundaries
	#define	MEMORY_ALIGNMENT		32

	// Use a 32 bit integer data type for dynamic programming. Since AVX provides
	// a 256-bit integer data type, we can process 8 elements in parallel.
	#define		SA_LEN				8u
	#define		SIMD_DATA			__m256i
	
	#define		SIMD_SET1			_mm256_set1_epi32
	#define		SIMD_ADD			_mm256_add_epi32
	#define		SIMD_SUB			_mm256_sub_epi32
	#define		SIMD_CMPEQ			_mm256_cmpeq_epi32
	#define		SIMD_CMPGT			_mm256_cmpgt_epi32
	#define		SIMD_MAX			_mm256_max_epi32
	
	#define		SIMD_BITWISE_AND		_mm256_and_si256
	#define		SIMD_BITWISE_OR			_mm256_or_si256
	#define		SIMD_BITWISE_AND_NOT		_mm256_andnot_si256
	
	#define		SIMD_MALLOC			_mm_malloc
	#define		SIMD_FREE			_mm_free
	
	typedef 	int				SA_Score;
	
#elif ATA32x4
	
	// 32-bit OSes require allocated memory be aligned to 16 byte boundaries. This requires use of
	// the _mm_malloc/_mm_free functions (instead of new and delete). This dependance does not
	// appear to be necessary for 64-bit OSes (but is important for backwards compatibility).
	#define	MEMORY_ALIGNMENT		16

	// Use a 32 bit integer data type for dynamic programming. Since SSE provides
	// a 128-bit integer data type, we can process 4 elements in parallel.
	#define		SA_LEN				4u
	#define		SIMD_DATA			__m128i
	
	#define		SIMD_SET1			_mm_set1_epi32
	#define		SIMD_ADD			_mm_add_epi32
	#define		SIMD_SUB			_mm_sub_epi32
	#define		SIMD_CMPEQ			_mm_cmpeq_epi32
	#define		SIMD_CMPGT			_mm_cmpgt_epi32
	#define		SIMD_MAX			_mm_max_epi32
	
	#define		SIMD_BITWISE_AND		_mm_and_si128
	#define		SIMD_BITWISE_OR			_mm_or_si128
	#define		SIMD_BITWISE_AND_NOT		_mm_andnot_si128
	
	#define		SIMD_MALLOC			_mm_malloc
	#define		SIMD_FREE			_mm_free
	
	typedef 	int				SA_Score;
#elif ATA32x2
	
	// 32-bit OSes require allocated memory be aligned to 16 byte boundaries. This requires use of
	// the _mm_malloc/_mm_free functions (instead of new and delete). This dependance does not
	// appear to be necessary for 64-bit OSes (but is important for backwards compatibility).
	#define	MEMORY_ALIGNMENT		16

	// Use a 32 bit integer data type for dynamic programming. Since MMX provides
	// a 64-bit integer data type, we can process 2 elements in parallel.
	#define		SA_LEN				2u
	#define		SIMD_DATA			__m64
	
	#define		SIMD_SET1			_mm_set1_pi32
	#define		SIMD_ADD			_mm_add_pi32
	#define		SIMD_SUB			_mm_sub_pi32
	#define		SIMD_CMPEQ			_mm_cmpeq_pi32
	#define		SIMD_CMPGT			_mm_cmpgt_pi32
	
	// Note that the MMX instruction set does not provide a built-in max function 
	// We need to provide this functionality ourselves (see _mm_max_pi32 defined below).
	#define		SIMD_MAX			_mm_max_pi32
	
	#define		SIMD_BITWISE_AND		_mm_and_si64
	#define		SIMD_BITWISE_OR			_mm_or_si64
	#define		SIMD_BITWISE_AND_NOT		_mm_andnot_si64
	
	#define		SIMD_MALLOC			_mm_malloc
	#define		SIMD_FREE			_mm_free
	
	typedef 	int				SA_Score;
	
	inline __m64 _mm_max_pi32(__m64 A, __m64 B)
	{
		__m64 mask = _mm_cmpgt_pi32(A, B);
		
		return _mm_or_si64( _mm_and_si64(mask, A), _mm_andnot_si64(mask, B) );
		
	};
	
#elif ATA32x1 // No SIMD!
	
	#define		MEMORY_ALIGNMENT				0

	#define		SA_LEN							1u
	#define		SIMD_DATA						int
	
	#define		SIMD_SET1(_X)					(_X)
	#define		SIMD_ADD(_X, _Y)				( (_X) + (_Y) )
	#define		SIMD_SUB(_X, _Y)				( (_X) - (_Y) )
	#define		SIMD_CMPEQ(_X, _Y)				( ( (_X) == (_Y) ) ? 0xFFFFFFFF : 0x00000000 )
	#define		SIMD_CMPGT(_X, _Y)				( ( (_X) > (_Y) ) ? 0xFFFFFFFF : 0x00000000 )
	
	#define		SIMD_MAX(_X, _Y)				( ((_X) > (_Y)) ? (_X) : (_Y) )
	
	#define		SIMD_BITWISE_AND(_X, _Y)		( (_X) & (_Y) )
	#define		SIMD_BITWISE_OR(_X, _Y)			( (_X) | (_Y) )
	#define		SIMD_BITWISE_AND_NOT(_X, _Y)	( ~(_X) & (_Y) )
	
	#define		SIMD_MALLOC(_X,_Y)				malloc(_X)
	#define		SIMD_FREE						free
	
	typedef 	int								SA_Score;
	
#elif ATA16x8
	
	// 32-bit OSes require allocated memory be aligned to 16 byte boundaries. This requires use of
	// the _mm_malloc/_mm_free functions (instead of new and delete). This dependance does not
	// appear to be necessary for 64-bit OSes (but is important for backwards compatibility).
	#define	MEMORY_ALIGNMENT		16
	
	// Use a 16 bit integer data type for dynamic programming. Since SSE provides
	// a 128bit integer data type, we can process 8 elements in parallel.
	#define		SA_LEN				8u
	#define		SIMD_DATA			__m128i
	
	#define		SIMD_SET1			_mm_set1_epi16
	#define		SIMD_ADD			_mm_add_epi16
	#define		SIMD_SUB			_mm_sub_epi16
	#define		SIMD_CMPEQ			_mm_cmpeq_epi16
	#define		SIMD_CMPGT			_mm_cmpgt_epi16
	#define		SIMD_MAX			_mm_max_epi16
	
	#define		SIMD_BITWISE_AND		_mm_and_si128
	#define		SIMD_BITWISE_OR			_mm_or_si128
	#define		SIMD_BITWISE_AND_NOT		_mm_andnot_si128
	
	#define		SIMD_MALLOC			_mm_malloc
	#define		SIMD_FREE			_mm_free
	
	typedef 	short int			SA_Score;

#elif ATA16x16

	// Using __m256i requires memory aligned on 32 byte boundaries
	#define	MEMORY_ALIGNMENT		32

	// Use a 32 bit integer data type for dynamic programming. Since AVX provides
	// a 256-bit integer data type, we can process 8 elements in parallel.
	#define		SA_LEN				16u
	#define		SIMD_DATA			__m256i
	
	#define		SIMD_SET1			_mm256_set1_epi16
	#define		SIMD_ADD			_mm256_add_epi16
	#define		SIMD_SUB			_mm256_sub_epi16
	#define		SIMD_CMPEQ			_mm256_cmpeq_epi16
	#define		SIMD_CMPGT			_mm256_cmpgt_epi16
	#define		SIMD_MAX			_mm256_max_epi16
	
	#define		SIMD_BITWISE_AND		_mm256_and_si256
	#define		SIMD_BITWISE_OR			_mm256_or_si256
	#define		SIMD_BITWISE_AND_NOT	_mm256_andnot_si256
	
	#define		SIMD_MALLOC			_mm_malloc
	#define		SIMD_FREE			_mm_free
	
	typedef 	short int				SA_Score;
#else
	#error "Please specificy a SIMD data size!"
#endif

// Default nucleic acid alignment scores
#define	DEFAULT_NA_ALIGN_MATCH		2
#define	DEFAULT_NA_ALIGN_MISMATCH	-3
#define	DEFAULT_NA_ALIGN_MASK		0
#define	DEFAULT_NA_ALIGN_GAP_OPEN	-5
#define	DEFAULT_NA_ALIGN_GAP_EXTEND	-2
//#define	DEFAULT_NA_ALIGN_GAP_EXTEND	0

// Default amino acid alignment scores
#define	DEFAULT_AA_ALIGN_GAP_OPEN	-11
#define	DEFAULT_AA_ALIGN_GAP_EXTEND	-1

union simd_elem {
	
	simd_elem()
	{
		// Do nothing
	};
	
	simd_elem(const SA_Score &m_scalar)
	{
		simd = SIMD_SET1(m_scalar);
	};
	
	#ifndef ATA32x1
	simd_elem(const SIMD_DATA &m_vector)
	{
		simd = m_vector;
	};
	#endif // ATA32x1
	
	~simd_elem()
	{
		// Do nothing
	};
	
	SA_Score max() const
	{
		SA_Score ret = v[0];
		
		for(unsigned int i = 1;i < SA_LEN;++i){
			ret = (ret < v[i]) ? v[i] : ret;
		}
		
		return ret;
	};
	
	SIMD_DATA simd;
	SA_Score v[SA_LEN];
};

namespace SA {

// Trace back locations:
// 
// im1 = i - 1
// jm1 = j - 1
//
//--------------------
// | im1_jm1 | im1_j |
//--------------------
// | i_jm1   |   x   |
//--------------------

#define	im1_jm1			(1 << 0)	// match or mismatch
#define	im1_j			(1 << 1)	// gap & query (i.e. query matches gap)
#define	i_jm1			(1 << 2)	// gap & target (i.e.  target matches gap)
#define	invalid_trace		(1 << 3)

// Match states
#define	query_target		im1_jm1
#define	query_gap			im1_j
#define	gap_target			i_jm1
#define	invalid_match		invalid_trace

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

// The definition of the amino acids comes *before* the definition of the nucleic acids
// to allow the nucleic acid gap value to be the same as the amino acid gap value.
namespace AA {
	
	// Allowed amino acid bases
	// The actual values must match the score matrix arrangement
	typedef enum {
		A = 0, // Alanine, Ala
		R, // Arginine, Arg
		N, // Asparagine, Asn
		D, // Aspartic acid, Asp
		C, // Cysteine, Cys
		Q, // Glutamine, Gln
		E, // Glutamic acid, Glu
		G, // Glycine, Gly
		H, // Histidine, His
		I, // Isoleucine, Ile
		L, // Leucine, Leu
		K, // Lysine, Lys
		M, // Methionine, Met
		F, // Phenylalanine, Phe
		P, // Proline, Pro
		S, // Serine, Ser
		T, // Threonine, Thr
		W, // Tryptophan, Trp
		Y, // Tyrosine, Tyr
		V, // Valine, Val
		B, // Aspartic acid or Asparagine, Asx
		Z, // Glutamine or Glutamic acid, Glx
		X, // Any amino acid, Xaa
		GAP = (1 << 5) // Make sure that no bits overlap with any nucleic acid

	} amino_acid;
	
	const size_t NUM_AA = (X + 1); // (20 aa + B, Z and X)
	const size_t MATRIX_SIZE = NUM_AA*NUM_AA; // (20 aa + B, Z and X) by (20 aa + B, Z and X)
}

namespace NA {
	
	// Allowed nucleic acid bases
	//	A = adenine
	//	C = cytosine
	//	G = guanine
	//	T = thymine
	//	M = A or C
	//	R = G or A
	//	S = G or C
	//	V = G or C or A
	//	W = A or T
	//	Y = T or C
	//	H = A or C or T
	//	K = G or T
	//	D = G or A or T
	//	B = G or T or C
	//	N = A or T or C or G
	typedef enum {
		A = (1 << 0), 
		C = (1 << 1), 
		G = (1 << 2), 
		T = (1 << 3), 
		M = (A | C),
		R = (G | A),
		S = (G | C),
		V = (G | C | A),
		W = (A | T),
		Y = (T | C),
		H = (A | C | T),
		K = (G | T),
		D = (G | A | T),
		B = (G | T | C),
		N = (A | T | C | G),
		// It make life easier to share the same value for a GAP between
		// nucleic acids and amino acids
		GAP = AA::GAP
	} nucleic_acid;
}

// Both base types (nucleic acid and amino acid) must fit
// into a variable of size base_type
typedef unsigned char base_type;

// Since the arguments to these functions are just a single byte, just make
// a copy of the byte (rather than passing as a const reference). Ideally, these
// should be inlined functions.
char bits_to_na(base_type m_bits);
base_type na_to_bits(char m_base);
char bits_to_aa(base_type m_bits);
base_type aa_to_bits(char m_base);
base_type complement(base_type m_base);

std::deque<base_type> reverse_complement(const std::deque<base_type> &m_seq);

class Alignment
{	
	private:
		bool is_na; // If is_na == false, then we are using amino acids
		bool mask_N_na;
		bool valid;

		SA_Score score;
		
		// The relative score is defined as:
		// 	  (observed score) - (worst possible score)
		//	----------------------------------------------
		//	(best possible score) - (worst possible score)
		float relative_score;
		
		unsigned int num_gap;
		unsigned int num_match;
		unsigned int num_mismatch;

		std::deque<base_type> query_align;
		std::deque<base_type> target_align;

		unsigned int query_start;
		unsigned int query_stop;
		
		unsigned int target_start;
		unsigned int target_stop;
		
	public:	
		Alignment()
		{
			is_na = true;
			mask_N_na = false;
			valid = false;

			score = 0;
			relative_score = 0.0f;
			
			num_gap = 0;
			num_match = 0;
			num_mismatch = 0;
			
			query_start = 0;
			query_stop = 0;
			target_start = 0;
			target_stop = 0;
		};

		~Alignment()
		{
			// Do nothing
		};

		void clear()
		{
			valid = false;

			num_gap = 0;
			num_match = 0;
			num_mismatch = 0;

			query_align.clear();
			target_align.clear();		
		};

		void push_front(const base_type &m_query_base, const base_type &m_target_base)
		{
			query_align.push_front(m_query_base);
			target_align.push_front(m_target_base);

			if( (m_query_base == NA::GAP) || (m_target_base == NA::GAP) ){
				++num_gap;
			}
			else{
				if(is_na){

					// Don't count matches/mismatches if we are treating 'N' as a masking character
					// and one or both of the query and target is an 'N'
					if(!mask_N_na || ( (m_query_base != NA::N) && (m_target_base != NA::N) ) ){
						if(m_query_base & m_target_base){
							++num_match;
						}
						else{
							++num_mismatch;
						}
					}
				}
				else{ // AA
					if(m_query_base == m_target_base){
						++num_match;
					}
					else{
						++num_mismatch;
					}
				}
			}
		};

		inline void enable_mask_N_na()
		{
			mask_N_na = true;
		};

		inline void disable_mask_N_na()
		{
			mask_N_na = false;
		};

		inline void set_na(const bool &m_is_na)
		{
			is_na = m_is_na;
		};
		
		inline bool get_na() const
		{
			return is_na;
		};
		
		inline void set_valid()
		{
			valid = true;
		};
		
		inline bool is_valid() const
		{
			return valid;
		};
		
		inline SA_Score get_score() const
		{
			return score;
		};
		
		inline float get_relative_score() const
		{
			return relative_score;
		};
		
		inline void set_score(const SA_Score &m_score)
		{
			score = m_score;
		};
		
		inline void set_relative_score(const float &m_relative_score)
		{
			relative_score = m_relative_score;
		};
		
		inline unsigned int get_num_match() const
		{
			return num_match;
		};
		
		inline unsigned int get_num_mismatch() const
		{
			return num_mismatch;
		};
		
		inline unsigned int get_num_gap() const
		{
			return num_gap;
		};
		
		inline const std::deque<base_type>& get_query_align() const
		{
			return query_align;
		};
		
		inline const std::deque<base_type>& get_target_align() const
		{
			return target_align;
		};
		
		inline unsigned int get_query_start() const
		{
			return query_start;
		};
		
		inline void set_query_start(const unsigned int &m_start)
		{
			query_start = m_start;
		};
		
		inline unsigned int get_query_stop() const
		{
			return query_stop;
		};
		
		inline void set_query_stop(const unsigned int &m_stop)
		{
			query_stop = m_stop;
		};
		
		inline unsigned int get_target_start() const
		{
			return target_start;
		};
		
		inline void set_target_start(const unsigned int &m_start)
		{
			target_start = m_start;
		};
		
		inline unsigned int get_target_stop() const
		{
			return target_stop;
		};
		
		inline void set_target_stop(const unsigned int &m_stop)
		{
			target_stop = m_stop;
		};
		
		std::string print() const;
		friend std::ostream& operator << (std::ostream &m_s, const Alignment &m_align);
};

class SeqAlign{

	public:

	typedef enum {SmithWaterman, Overlap, Global} AlignmentMode;

	private:
		
		struct SA_Elem{

			SA_Elem()
			{
				M.simd = I_query.simd = I_target.simd = SIMD_SET1(-1);
				
				M_trace.simd = I_query_trace.simd = I_target_trace.simd = SIMD_SET1(invalid_trace);
			}

			~SA_Elem()
			{
				// Do nothing
			};

			// Use the notation of "Biological sequence analysis" by Durbin, Eddy, Krogh and Mitchison
			simd_elem M;
			simd_elem I_query;  // insertion in query
			simd_elem I_target; // insertion in target

			// The trace-back pointers
			simd_elem M_trace;
			simd_elem I_query_trace;
			simd_elem I_target_trace;			
		};
		
		AlignmentMode mode;
		
		// The dynamic programming matrix is (max query len + 1) by (max target len + 1)
		SA_Elem *dp_matrix;
		simd_elem max_loc;
		simd_elem query_len;
		simd_elem target_len;

		simd_elem match;
		simd_elem mask;
		simd_elem mismatch;
		simd_elem gap_existance;
		simd_elem gap_extension;

		// The input sequences in 5'-3' orientation
		simd_elem *query;
		simd_elem *target;
				
		bool is_na; // If is_na == false, then we are using amino acids
		
		// If mask_N_na == true, then treat 'N' as a masking character that gets a score of mask
		bool mask_N_na;
		
		std::vector<SA_Score> aa_score_matrix;
		
		// The current alignments
		Alignment curr_align[SA_LEN];
		
		void init_BLOSUM62(std::vector<SA_Score> &m_matrix);
		
		void translate_na_seq(const std::string &m_input, std::vector<base_type> &m_output)
		{
			#ifdef _DEBUG
			if(!is_na){
				throw __FILE__ ":translate_na_seq: Only valid for NA sequences";
			}
			
			if( m_input.size() != m_output.size() ){
				throw __FILE__ ":translate_na_seq: input/output size mismatch";
			}
			#endif // _DEBUG
			
			std::string::const_iterator i = m_input.begin();
			std::vector<base_type>::iterator o = m_output.begin();

			for(;i != m_input.end();i++, o++){
			
				switch(*i){
					case 'A': case 'a':
						*o = NA::A;
						break;
					case 'T': case 't':
						*o = NA::T;
						break;
					case 'G': case 'g':
						*o = NA::G;
						break;
					case 'C': case 'c':
						*o = NA::C;
						break;
					case 'M': case 'm':
						*o = NA::M;
						break;
					case 'R': case 'r':
						*o = NA::R;
						break;
					case 'S': case 's':
						*o = NA::S;
						break;
					case 'V': case 'v':
						*o = NA::V;
						break;
					case 'W': case 'w':
						*o = NA::W;
						break;
					case 'Y': case 'y':
						*o = NA::Y;
						break;
					case 'H': case 'h':
						*o = NA::H;
						break;
					case 'K': case 'k':
						*o = NA::K;
						break;
					case 'D': case 'd':
						*o = NA::D;
						break;
					case 'B': case 'b':
						*o = NA::B;
						break;
					case 'N': case 'n':
					case 'I': case 'i': // For now, treat inosine as an 'N'
						*o = NA::N;
						break;
					default:
						throw __FILE__ ":translate_na_seq: Illegal base";
						break;
				};
			}
		};

		void translate_aa_seq(const std::string &m_input, std::vector<base_type> &m_output)
		{
			#ifdef _DEBUG
			if(is_na){
				throw __FILE__ ":translate_aa_seq: Only valid for AA sequences";
			}
			
			if( m_input.size() != m_output.size() ){
				throw __FILE__ ":translate_aa_seq: input/output size mismatch";
			}
			#endif // _DEBUG
			
			std::string::const_iterator i = m_input.begin();
			std::vector<base_type>::iterator o = m_output.begin();

			for(;i != m_input.end();i++, o++){
			
				switch(*i){
					case 'A': case 'a':
						*o = AA::A;
						break;
					case 'R': case 'r':
						*o = AA::R;
						break;
					case 'N': case 'n':
						*o = AA::N;
						break;
					case 'D': case 'd':
						*o = AA::D;
						break;
					case 'C': case 'c':
						*o = AA::C;
						break;
					case 'Q': case 'q':
						*o = AA::Q;
						break;
					case 'E': case 'e':
						*o = AA::E;
						break;
					case 'G': case 'g':
						*o = AA::G;
						break;
					case 'H': case 'h':
						*o = AA::H;
						break;
					// J = I or L. Since there is no 'J' in the
					// BLOSUM62 matrix, approximate J = I
					case 'J': case 'j':
					case 'I': case 'i':
						*o = AA::I;
						break;
					case 'L': case 'l':
						*o = AA::L;
						break;
					case 'K': case 'k':
						*o = AA::K;
						break;
					case 'M': case 'm':
						*o = AA::M;
						break;
					case 'F': case 'f':
						*o = AA::F;
						break;
					case 'P': case 'p':
						*o = AA::P;
						break;
					case 'S': case 's':
						*o = AA::S;
						break;
					case 'T': case 't':
						*o = AA::T;
						break;
					case 'W': case 'w':
						*o = AA::W;
						break;
					case 'Y': case 'y':
						*o = AA::Y;
						break;
					case 'V': case 'v':
						*o = AA::V;
						break;
					case 'B': case 'b':
						*o = AA::B;
						break;
					case 'Z': case 'z':
						*o = AA::Z;
						break;
					case 'X': case 'x':
					case 'U': case 'u': // Treat selenocysteine as 'X' (same as BLAST)
						*o = AA::X;
						break;
					default:
						throw __FILE__ ":translate_aa_seq: Illegal base";
						break;
				};
			}
		};

		void align_overlap();
		void align_global();
		void align_smith_waterman();

		void compute_score(const unsigned int &m_index);
		void trace_back(const unsigned int &m_index);
		
		simd_elem *resize_sequence(const size_t &new_len, simd_elem *ptr_old_seq, const size_t &old_len)
		{
			if(new_len == 0){
				
				if(ptr_old_seq != NULL){
					SIMD_FREE(ptr_old_seq);
				}
				
				return NULL;
			}
			
			simd_elem *ptr_new_seq = (simd_elem*)SIMD_MALLOC(new_len*sizeof(simd_elem), MEMORY_ALIGNMENT);

			if(ptr_new_seq == NULL){
				throw __FILE__ ":SeqAlign::resize_sequence: Unable to allocate sequence memory";
			}				

			// Copy any existing query data
			for(size_t i = 0;i < std::min(old_len, new_len);i++){
				ptr_new_seq[i].simd = ptr_old_seq[i].simd;
			}

			// If needed (i.e. new_len > old_len) initialize the 
			// remaining new query memory to be all gaps
			for(size_t i = old_len;i < new_len;i++){
				ptr_new_seq[i].simd = SIMD_SET1( SA_Score(NA::GAP) );
			}

			// Free the old sequence to avoid memory leaks
			if(ptr_old_seq != NULL){
				SIMD_FREE(ptr_old_seq);
			}

			return ptr_new_seq;
		};

	public:
		
		SeqAlign(const AlignmentMode &m_mode, const bool &m_is_na);
			
		~SeqAlign()
		{
			// Free the dynamic programming matrix
			if(dp_matrix != NULL){
			
				SIMD_FREE(dp_matrix);
				dp_matrix = NULL;
			}
			
			// Free the query sequence
			if(query != NULL){
			
				SIMD_FREE(query);
				query = NULL;
			}
			
			// Free the target sequence
			if(target != NULL){
			
				SIMD_FREE(target);
				target = NULL;
			}
		};

		void align()
		{
			switch(mode){
				case Overlap:
					align_overlap();
					break;
				case Global:
					align_global();
					break;
				case SmithWaterman:
					align_smith_waterman();
					break;
				default:
					throw __FILE__ ":align: Unknown mode";
			};
			
			for(unsigned int i = 0;i < SA_LEN;i++){
				compute_score(i);
			}
		}

		inline void set_mode(const AlignmentMode &m_mode)
		{
			mode = m_mode;
		};

		inline void enable_mask_N_na()
		{
			mask_N_na = true;
			
			for(unsigned int i = 0;i < SA_LEN;i++){
				curr_align[i].enable_mask_N_na();
			}
		};
		
		inline void disable_mask_N_na()
		{
			mask_N_na = false;
			
			for(unsigned int i = 0;i < SA_LEN;i++){
				curr_align[i].disable_mask_N_na();
			}
		};
		
		// Pack the supplied query into slot m_index
		inline void pack_query(const unsigned int &m_index, const std::string &m_query,
			const bool &m_reverse_complement = false)
		{
			#ifdef _DEBUG
			if(m_index > SA_LEN){
				throw __FILE__ ":pack_query: m_index > SA_LEN";
			}
			#endif // _DEBUG
			
			const unsigned int len = m_query.size();

			const size_t old_max_query_len = query_len.max();
			
			query_len.v[m_index] = len;
			
			const size_t new_max_query_len = query_len.max();
			
			if(old_max_query_len != new_max_query_len){
				
				// Grow or shrink the query sequence array as needed
				query = resize_sequence(new_max_query_len, query, old_max_query_len);
			}
						
			if(is_na){
				if(m_reverse_complement){
					for(unsigned int i = 0;i < len;i++){
						query[i].v[m_index] = complement( na_to_bits(m_query[len - (i + 1)]) );
					}
				}
				else{
					for(unsigned int i = 0;i < len;i++){
						query[i].v[m_index] = na_to_bits(m_query[i]);
					}
				}
			}
			else{ // aa
				for(unsigned int i = 0;i < len;i++){
					query[i].v[m_index] = aa_to_bits(m_query[i]);
				}
			}

			// If needed, pad the current sequence with a GAP symbol 
			// that will not match any other base
			for(unsigned int i = len;i < new_max_query_len;i++){
				query[i].v[m_index] = SA_Score(NA::GAP);
			}
		};
		
		inline void pack_query(const unsigned int &m_index, const char* m_query,
			const bool &m_reverse_complement = false)
		{
			pack_query(m_index, std::string(m_query), m_reverse_complement);
		};
		
		// Pack the supplied query into *all* slots
		inline void pack_query(const std::string &m_query, const bool &m_reverse_complement = false)
		{
			const unsigned int len = m_query.size();

			if(query != NULL){
				
				SIMD_FREE(query);
				query = NULL;
			}
			
			query = (simd_elem*)SIMD_MALLOC(len*sizeof(simd_elem), MEMORY_ALIGNMENT);
			
			if(query == NULL){
				throw __FILE__ ":SeqAlign::pack_query: Unable to allocate sequence memory";
			}
			
			// Set the query lengths
			query_len = SIMD_SET1(len);
			
			if(is_na){
				if(m_reverse_complement){
					for(unsigned int i = 0;i < len;i++){
						query[i].simd = SIMD_SET1( complement(na_to_bits(m_query[len - (i + 1)]) ) );
					}
				}
				else{
					for(unsigned int i = 0;i < len;i++){
						query[i].simd = SIMD_SET1( na_to_bits(m_query[i]) );
					}
				}
			}
			else{ // aa
				for(unsigned int i = 0;i < len;i++){
					query[i].simd = SIMD_SET1( aa_to_bits(m_query[i]) );
				}
			}
		};
		
		// Clear the query in slot m_index bby packing an empty string
		inline void clear_query(const unsigned int &m_index)
		{
			pack_query( m_index, std::string() );
		};
		
		// Clear all queries in all slots
		inline void clear_query()
		{
			query_len = SIMD_SET1(0);
			
			if(query != NULL){
			
				SIMD_FREE(query);
				query = NULL;
			}
		};
		
		// Pack the supplied target into slot m_index
		inline void pack_target(const unsigned int &m_index, const std::string &m_target,
			const bool &m_reverse_complement = false)
		{
			#ifdef _DEBUG
			if(m_index > SA_LEN){
				throw __FILE__ ":pack_target: m_index > SA_LEN";
			}
			#endif // _DEBUG
			
			const unsigned int len = m_target.size();

			const size_t old_max_target_len = target_len.max();
			
			target_len.v[m_index] = len;
			
			const size_t new_max_target_len = target_len.max();
			
			if(old_max_target_len != new_max_target_len){
				
				// Grow or shrink the target sequence array as needed
				target = resize_sequence(new_max_target_len, target, old_max_target_len);
			}
			
			if(is_na){
				
				if(m_reverse_complement){
					for(unsigned int i = 0;i < len;i++){
						target[i].v[m_index] = complement( na_to_bits(m_target[len - (i + 1)]) );
					}
				}
				else{
					for(unsigned int i = 0;i < len;i++){
						target[i].v[m_index] = na_to_bits(m_target[i]);
					}
				}
			}
			else{ // aa
				for(unsigned int i = 0;i < len;i++){
					target[i].v[m_index] = aa_to_bits(m_target[i]);
				}
			}

			// If needed, pad the sequence with a GAP symbol
			// that will not match any other base
			for(unsigned int i = len;i < new_max_target_len;i++){
				target[i].v[m_index] = SA_Score(NA::GAP);
			}
		};
		
		inline void pack_target(const unsigned int &m_index, const char* m_target,
			const bool &m_reverse_complement = false)
		{
			pack_target(m_index, std::string(m_target), m_reverse_complement);
		};
		
		// Pack the supplied target into *all* slots
		inline void pack_target(const std::string &m_target,
			const bool &m_reverse_complement = false)
		{
			const unsigned int len = m_target.size();

			if(target != NULL){
				
				SIMD_FREE(target);
				target = NULL;
			}
			
			target = (simd_elem*)SIMD_MALLOC(len*sizeof(simd_elem), MEMORY_ALIGNMENT);
			
			if(target == NULL){
				throw __FILE__ ":SeqAlign::pack_target: Unable to allocate sequence memory";
			}
			
			// Set the target lengths
			target_len = SIMD_SET1(len);
			
			if(is_na){
				if(m_reverse_complement){
					for(unsigned int i = 0;i < len;i++){
						target[i].simd = SIMD_SET1( complement( na_to_bits(m_target[len - (i + 1)]) ) );
					}
				}
				else{
					for(unsigned int i = 0;i < len;i++){
						target[i].simd = SIMD_SET1( na_to_bits(m_target[i]) );
					}
				}
			}
			else{ // aa
				for(unsigned int i = 0;i < len;i++){
					target[i].simd = SIMD_SET1( aa_to_bits(m_target[i]) );
				}
			}
		};
		
		// Clear the target in slot m_index by packing an empty target
		inline void clear_target(const unsigned int &m_index)
		{
			pack_target( m_index, std::string() );
		};
		
		// Clear all targets in all slots
		inline void clear_target()
		{
			target_len = SIMD_SET1(0);
			
			if(target != NULL){
			
				SIMD_FREE(target);
				target = NULL;
			}
		};
		
		// Return the query and target sequences as std::strings
		std::string query_seq(const unsigned int &m_index) const
		{
			#ifdef _DEBUG
			if(m_index > SA_LEN){
				throw __FILE__ ":query_seq: Index out of bounds";
			}
			#endif // _DEBUG
			
			std::string ret;
			
			if(is_na){
				for(SA_Score i = 0;i < query_len.v[m_index];i++){
					
					if(query[i].v[m_index] == NA::GAP){
						break;
					}
					
					ret.push_back( bits_to_na(query[i].v[m_index]) );
				}
			}
			else{
				for(SA_Score i = 0;i < query_len.v[m_index];i++){
				
					if(query[i].v[m_index] == NA::GAP){
						break;
					}
					
					ret.push_back( bits_to_aa(query[i].v[m_index]) );
				}
			}			
			
			return ret;
		};
		
		std::string target_seq(const unsigned int &m_index) const
		{
			#ifdef _DEBUG
			if(m_index > SA_LEN){
				throw __FILE__ ":target_seq: Index out of bounds";
			}
			#endif // _DEBUG
			
			std::string ret;
			
			if(is_na){
				for(SA_Score i = 0;i < target_len.v[m_index];i++){
				
					if(target[i].v[m_index] == NA::GAP){
						break;
					}
					
					ret.push_back( bits_to_na(target[i].v[m_index]) );
				}
			}
			else{
				for(SA_Score i = 0;i < target_len.v[m_index];i++){
				
					if(target[i].v[m_index] == NA::GAP){
						break;
					}
					
					ret.push_back( bits_to_aa(target[i].v[m_index]) );
				}
			}		
			
			return ret;
		};
		
		const Alignment& get_alignment(const unsigned int &m_index)
		{
			#ifdef _DEBUG
			if(m_index > SA_LEN){
				throw __FILE__ ":get_alignment: Index out of bounds";
			}
			#endif // _DEBUG
			
			if(curr_align[m_index].is_valid() == false){
				trace_back(m_index);
			}
			
			return curr_align[m_index];
		};
		
		inline void clear_alignments()
		{
			for(unsigned int i = 0;i < SA_LEN;i++){
				curr_align[i].clear();
			}
		};
			
		inline SA_Score score(const unsigned int &m_index) const
		{
			#ifdef _DEBUG
			if(m_index > SA_LEN){
				throw __FILE__ ":score: Index out of bounds";
			}
			#endif // _DEBUG
			
			return curr_align[m_index].get_score();
		};
		
		inline SA_Score get_match_score() const
		{
			return match.v[0];
		};
		
		inline void set_match_score(const SA_Score &m_score)
		{
			match.simd = SIMD_SET1(m_score);
		};

		inline SA_Score get_mismatch_score() const
		{
			return mismatch.v[0];
		};
		
		inline void set_mismatch_score(const SA_Score &m_score)
		{
			mismatch.simd = SIMD_SET1(m_score);
		};

		inline SA_Score get_gap_open_score() const
		{
			return gap_existance.v[0];
		};
		
		inline void set_gap_open_score(const SA_Score &m_score)
		{
			gap_existance.simd = SIMD_SET1(m_score);
		};

		inline SA_Score get_gap_extension_score() const
		{
			return gap_extension.v[0];
		};

		inline void set_gap_extension_score(const SA_Score &m_score)
		{
			gap_extension.simd = SIMD_SET1(m_score);
		};
};

std::ostream& operator << (std::ostream &m_s, const SeqAlign &m_align);

} // namespace::SA

#endif // __SEQ_ALIGN

