#ifndef __SIGNATURE_CLASS
#define __SIGNATURE_CLASS

#include <string.h> // for memcpy

#include <string>
#include <algorithm> // for std::pair
#include <vector>
#include <sstream>

#include "mpi_util.h"

// The symbol defined by COMMENT_SYMBOL is used to comment out
// lines in the assya input file
#define	COMMENT_SYMBOL						'#'

// What assay are we simulating?
enum {ASSAY_PCR, ASSAY_PROBE, ASSAY_PADLOCK, ASSAY_AFFYMETRIX, ASSAY_NONE};

// We need to protect any commas that appear in template variables
// if we are going to combine them with X Macros
#define SINGLE_ARG(...) __VA_ARGS__

// A class for storing hybridization signatures
class hybrid_sig {

	private:
		
		inline int def_to_gi(const std::string &m_defline) const
		{
			
			size_t start = m_defline.find("gi|");
			
			if(start == std::string::npos){
				return -1;
			}
		
			start += 3; // strlen("gi|")
			
			const size_t len = m_defline.size();
			
			size_t stop = start;
			
			while( (stop < len) && isdigit(m_defline[stop]) ){
				stop ++;
			}
			
			return atoi( m_defline.substr(start, stop - start).c_str() );
		};

	public:
		// Enumerate the DNA strand orientation
		// A value of PLUS for probe_strand means the probe binds to the 
		// forward strand, while a value of MINUS means the probe binds to
		// the reverse strand
		enum {PLUS = 0, MINUS};
		
		// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
		// ensure that structure variable are correctly serialized.
		#define HYBRID_SIG_MEMBERS \
			VARIABLE(std::string, name) \
			VARIABLE(std::string, forward_oligo) \
			VARIABLE(std::string, reverse_oligo) \
			VARIABLE(std::string, probe_oligo) \
			VARIABLE(std::string, amplicon_def) \
			VARIABLE(std::string, amplicon) \
			VARIABLE(std::string, forward_align) \
			VARIABLE(std::string, reverse_align) \
			VARIABLE(std::string, probe_align) \
			VARIABLE(SINGLE_ARG(std::pair<int, int>), amplicon_range) \
			VARIABLE(SINGLE_ARG(std::pair<int, int>), probe_range) \
			VARIABLE(int, id) \
			VARIABLE(int, degen_id) \
			VARIABLE(int, seq_index) \
			VARIABLE(float, forward_tm) \
			VARIABLE(float, reverse_tm) \
			VARIABLE(float, probe_tm) \
			VARIABLE(float, forward_hairpin_tm) \
			VARIABLE(float, reverse_hairpin_tm) \
			VARIABLE(float, forward_dimer_tm) \
			VARIABLE(float, reverse_dimer_tm) \
			VARIABLE(float, primer_dimer_tm) \
			VARIABLE(float, probe_hairpin_tm) \
			VARIABLE(float, probe_dimer_tm) \
			VARIABLE(float, forward_dH) \
			VARIABLE(float, forward_dS) \
			VARIABLE(float, reverse_dH) \
			VARIABLE(float, reverse_dS) \
			VARIABLE(float, probe_dH) \
			VARIABLE(float, probe_dS) \
			VARIABLE(int, primer_strand) \
			VARIABLE(int, probe_strand) \
			VARIABLE(int, forward_primer_clamp) \
			VARIABLE(int, reverse_primer_clamp) \
			VARIABLE(int, forward_degen) \
			VARIABLE(int, reverse_degen) \
			VARIABLE(int, probe_degen) \
			VARIABLE(int, forward_mm) \
			VARIABLE(int, forward_gap) \
			VARIABLE(int, reverse_mm) \
			VARIABLE(int, reverse_gap) \
			VARIABLE(int, probe_mm) \
			VARIABLE(int, probe_gap) \
			VARIABLE(float, ct)

		#define VARIABLE(A, B) A B;
			HYBRID_SIG_MEMBERS
		#undef VARIABLE

		hybrid_sig() // default constructor
		{
			id = degen_id = -1;
			seq_index = -1;
			
			amplicon_range = std::make_pair(0, 0);
			probe_range = std::make_pair(0, 0);
			
			forward_tm = -1.0f;
			reverse_tm = -1.0f;
			probe_tm = -1.0f;
			
			forward_hairpin_tm = -1.0f;
			reverse_hairpin_tm = -1.0f;
			forward_dimer_tm = -1.0f;
			reverse_dimer_tm = -1.0f;
			primer_dimer_tm = -1.0f;
		
			probe_hairpin_tm = -1.0f;
			probe_dimer_tm = -1.0f;
		
			forward_dH = 100.0f;
			forward_dS = 0.0f;
		
			reverse_dH = 100.0f;
			reverse_dS = 0.0f;
		
			probe_dH = 100.0f;
			probe_dS = 0.0f;
		
			primer_strand = PLUS;
			probe_strand = PLUS;

			forward_degen = 1;
			reverse_degen = 1;
			probe_degen = 1;
			
			forward_mm = -1;
			forward_gap = -1;
			reverse_mm = -1;
			reverse_gap = -1;
			probe_mm = -1;
			probe_gap = -1;
			
			ct = -1.0f;
		};
		
		// forward + reverse + probe
		hybrid_sig(const std::string &m_name, const std::string &m_forward, 
			const std::string &m_reverse, const std::string &m_probe, 
			const unsigned int &m_id)
		{
			id = degen_id = m_id;
			seq_index = -1;
			
			name = m_name;
			forward_oligo = m_forward;
			reverse_oligo = m_reverse;
			probe_oligo = m_probe;
			
			amplicon_range = std::make_pair(0, 0);
			probe_range = std::make_pair(0, 0);
			
			forward_tm = -1.0f;
			reverse_tm = -1.0f;
			probe_tm = -1.0f;
			
			forward_hairpin_tm = -1.0f;
			reverse_hairpin_tm = -1.0f;
			forward_dimer_tm = -1.0f;
			reverse_dimer_tm = -1.0f;
			primer_dimer_tm = -1.0f;
		
			probe_hairpin_tm = -1.0f;
			probe_dimer_tm = -1.0f;
			
			forward_dH = 100.0f;
			forward_dS = 0.0f;
		
			reverse_dH = 100.0f;
			reverse_dS = 0.0f;
		
			probe_dH = 100.0f;
			probe_dS = 0.0f;
			
			primer_strand = PLUS;
			probe_strand = PLUS;
			
			forward_primer_clamp = -1;
			reverse_primer_clamp = -1;

			forward_degen = 1;
			reverse_degen = 1;
			probe_degen = 1;
			
			forward_mm = -1;
			forward_gap = -1;
			reverse_mm = -1;
			reverse_gap = -1;
			probe_mm = -1;
			probe_gap = -1;
			
			ct = -1.0f;
		}
		
		// forward + reverse
		hybrid_sig(const std::string &m_name, const std::string &m_forward, 
			const std::string &m_reverse, const unsigned int &m_id)
		{
			id = degen_id = m_id;
			seq_index = -1;
			
			name = m_name;
			forward_oligo = m_forward;
			reverse_oligo = m_reverse;
			
			amplicon_range = std::make_pair(0, 0);
			probe_range = std::make_pair(0, 0);
			
			forward_tm = -1.0f;
			reverse_tm = -1.0f;
			probe_tm = -1.0f;
			
			forward_hairpin_tm = -1.0f;
			reverse_hairpin_tm = -1.0f;
			forward_dimer_tm = -1.0f;
			reverse_dimer_tm = -1.0f;
			primer_dimer_tm = -1.0f;
		
			probe_hairpin_tm = -1.0f;
			probe_dimer_tm = -1.0f;
			
			forward_dH = 100.0f;
			forward_dS = 0.0f;
		
			reverse_dH = 100.0f;
			reverse_dS = 0.0f;
		
			probe_dH = 100.0f;
			probe_dS = 0.0f;
			
			primer_strand = PLUS;
			probe_strand = PLUS;
			
			forward_primer_clamp = -1;
			reverse_primer_clamp = -1;

			forward_degen = 1;
			reverse_degen = 1;
			probe_degen = 1;
			
			forward_mm = -1;
			forward_gap = -1;
			reverse_mm = -1;
			reverse_gap = -1;
			probe_mm = -1;
			probe_gap = -1;
			
			ct = -1.0f;
		}
		
		// probe
		hybrid_sig(const std::string &m_name, const std::string &m_probe, 
			const unsigned int &m_id)
		{
			id = degen_id = m_id;
			seq_index = -1;
			
			name = m_name;
			probe_oligo = m_probe;
			
			amplicon_range = std::make_pair(0, 0);
			probe_range = std::make_pair(0, 0);
			
			forward_tm = -1.0f;
			reverse_tm = -1.0f;
			probe_tm = -1.0f;
			
			forward_hairpin_tm = -1.0f;
			reverse_hairpin_tm = -1.0f;
			forward_dimer_tm = -1.0f;
			reverse_dimer_tm = -1.0f;
			primer_dimer_tm = -1.0f;
		
			probe_hairpin_tm = -1.0f;
			probe_dimer_tm = -1.0f;
			
			forward_dH = 100.0f;
			forward_dS = 0.0f;
		
			reverse_dH = 100.0f;
			reverse_dS = 0.0f;
		
			probe_dH = 100.0f;
			probe_dS = 0.0f;
			
			primer_strand = PLUS;
			probe_strand = PLUS;
			
			forward_primer_clamp = -1;
			reverse_primer_clamp = -1;

			forward_degen = 1;
			reverse_degen = 1;
			probe_degen = 1;
			
			forward_mm = -1;
			forward_gap = -1;
			reverse_mm = -1;
			reverse_gap = -1;
			probe_mm = -1;
			probe_gap = -1;
			
			ct = -1.0f;
		}
		
		~hybrid_sig()
		{
			// Do nothing!
		};
		
		inline int my_id() const
		{
			return id;
		};
		
		inline int my_degen_id() const
		{
			return degen_id;
		};

		inline void my_id(const int &m_id)
		{
			// Id's must be >= 0
			if(m_id < 0){
				throw __FILE__ ":my_id: m_id < 0";
			}

			id = m_id;
		};

		inline void my_degen_id(const int &m_id)
		{
			// Id's must be >= 0
			if(m_id < 0){
				throw __FILE__ ":my_degen_id: m_id < 0";
			}
			
			degen_id = m_id;
		};
		
		inline void seq_id(const int &m_id)
		{
			seq_index = m_id;
		};
		
		inline int seq_id() const
		{
			return seq_index;
		};
		
		inline bool has_primers() const
		{
			return ( !forward_oligo.empty() && !reverse_oligo.empty() );
		};
		
		inline bool has_probe() const
		{
			return ( !probe_oligo.empty() );
		};
		
		inline bool probe_overlap(const hybrid_sig &m_sig) const
		{			
			// Overlapping probes must be on the same strand
			if(probe_strand != m_sig.probe_strand){
				return false;
			}
			
			// Is the first edge of this probe within the bounds of m_sig?
			if( (probe_range.first >= m_sig.probe_range.first) && 
			    (probe_range.first <= m_sig.probe_range.second) ){
			    
			    return true;
			}
			
			// Is the second edge of this probe within the bounds of m_sig?
			if( (probe_range.second >= m_sig.probe_range.first) && 
			    (probe_range.second <= m_sig.probe_range.second) ){
			    
			    return true;
			}
			
			// Does this probe contain m_sig?
			if( (probe_range.first <= m_sig.probe_range.first) && 
			    (probe_range.second >= m_sig.probe_range.second) ){
			    
			    return true;
			}
			
			// If we get here, the probes don't overlap
			return false;
		};

		inline bool operator>(const hybrid_sig &m_rhs) const
		{
			if(id == m_rhs.id){
			
				// Sort (in descending order) by:
				// 1) minimum primer melting temperature
				// 2) probe melting temperature
				// 3) maximum primer melting temperature
				// 4) target sequence index
				if( min_primer_tm() == m_rhs.min_primer_tm() ){
					
					if(probe_tm == m_rhs.probe_tm){
						
						if( max_primer_tm() == m_rhs.max_primer_tm() ){
						
							// Sort by target sequence
							return seq_index > m_rhs.seq_index;
						}
						
						return ( max_primer_tm() < m_rhs.max_primer_tm() );
					}
					
					return (probe_tm < m_rhs.probe_tm);
				}
				
				return ( min_primer_tm() < m_rhs.min_primer_tm() );
			}
			
			return (id > m_rhs.id);
		};

		inline bool operator<(const hybrid_sig &m_rhs) const
		{
			if(id == m_rhs.id){
			
				// Sort (in descending order) by:
				// 1) minimum primer melting temperature
				// 2) probe melting temperature
				// 3) maximum primer melting temperature
				// 4) target sequence index
				if( min_primer_tm() == m_rhs.min_primer_tm() ){
				
					if(probe_tm == m_rhs.probe_tm){
						
						if( max_primer_tm() == m_rhs.max_primer_tm() ){
						
							// Sort by target sequence
							return seq_index < m_rhs.seq_index;
						}
						
						return ( max_primer_tm() > m_rhs.max_primer_tm() );
					}
					
					return (probe_tm > m_rhs.probe_tm);
				}
				
				return ( min_primer_tm() > m_rhs.min_primer_tm() );
			}
			
			return (id < m_rhs.id);
		};
		
		inline int min_primer_clamp() const
		{
			return std::min(forward_primer_clamp, reverse_primer_clamp);
		};
		
		inline int max_primer_clamp() const
		{
			return std::max(forward_primer_clamp, reverse_primer_clamp);
		};
		
		inline float min_primer_tm() const
		{
			// Clamp the min Tm at zero
			return std::max(0.0f, std::min(forward_tm, reverse_tm) );
		};
		
		inline float max_primer_tm() const
		{
			return std::max(forward_tm, reverse_tm);
		};
		
		inline void offset_ranges(const int &m_off)
		{
		
			if(has_primers() == true){
			
				amplicon_range.first += m_off;
				amplicon_range.second += m_off;
			}
			
			if(has_probe() == true){
			
				probe_range.first += m_off;
				probe_range.second += m_off;
			}
		};
		
		inline bool start_overlap(const int &m_start) const
		{
			if(has_primers() == true){
			
				return (amplicon_range.first <= m_start);
			}
			else{ // has_probe == true
			
				return (probe_range.first <= m_start);
			}
		};
		
		inline bool stop_overlap(const int &m_stop) const
		{
			if(has_primers() == true){
			
				return (amplicon_range.second >= m_stop);
			}
			else{ // has_probe == true
			
				return (probe_range.second >= m_stop);
			}
		};
		
		inline std::string assay_string() const
		{
			return name;
			
			//std::stringstream ssout;
			//
			//ssout << name << '\t';
			//
			//if(forward_oligo.empty() == false){
			//	ssout << forward_oligo << '\t';
			//}
			//
			//if(reverse_oligo.empty() == false){
			//	ssout << reverse_oligo << '\t';
			//}
			//
			//if(probe_oligo.empty() == false){
			//	ssout << probe_oligo << '\t';
			//}
			//
			//return ssout.str();
		};
};

template<> size_t mpi_size(const hybrid_sig &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const hybrid_sig &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, hybrid_sig &m_obj);

class sort_by_match{

	public:
		sort_by_match() {};
		~sort_by_match() {};

		inline bool operator()(const hybrid_sig &m_a, const hybrid_sig &m_b)
		{
			if( m_a.my_id() == m_b.my_id() ){
				return ( m_a.seq_id() < m_b.seq_id() );
			}
			
			return ( m_a.my_id() < m_b.my_id() );
		};
};

class sort_by_loc{

	public:
	
		sort_by_loc() {};
		~sort_by_loc() {};
		
		inline bool operator()(const hybrid_sig &m_a, const hybrid_sig &m_b)
		{
			// First sort by id
			if( m_a.my_id() == m_b.my_id() ){
			
				// Then by target sequence
				if( m_a.seq_id() == m_b.seq_id() ){
					
					// Then by range
					if( m_a.has_primers() == true){
						
						return m_a.amplicon_range < m_b.amplicon_range;
					}
					else{ // m_a.has_probe() == true
						
						return m_a.probe_range < m_b.probe_range;
					}
				}
				
				return ( m_a.seq_id() < m_b.seq_id() );
			}
			
			return ( m_a.my_id() < m_b.my_id() );
		};
};

class unique_by_loc{

	public:
	
		unique_by_loc() {};
		~unique_by_loc() {};
		
		inline bool operator()(const hybrid_sig &m_a, const hybrid_sig &m_b)
		{
			// First sort by id
			if( m_a.my_id() == m_b.my_id() ){
			
				// Then by target sequence
				if( m_a.seq_id() == m_b.seq_id() ){
					
					// Then by range
					if( m_a.has_primers() == true){
						
						return m_a.amplicon_range == m_b.amplicon_range;
					}
					else{ // m_a.has_probe() == true
						
						return m_a.probe_range == m_b.probe_range;
					}
				}
				
				return false;
			}
			
			return false;
		};
};

size_t read_input_file(const std::string &m_file, std::vector<hybrid_sig> &m_sig_list, 
	const bool &m_ignore_probe, const bool &m_force_probe);

#endif // __SIGNATURE_CLASS
