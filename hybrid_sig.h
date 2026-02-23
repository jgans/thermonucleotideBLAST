#ifndef __SIGNATURE_CLASS
#define __SIGNATURE_CLASS

#include <string.h> // for memcpy

#include <string>
#include <algorithm> // for std::pair
#include <vector>
#include <sstream>

#include "mpi_util.h"
#include "throw.h"

// The symbol defined by COMMENT_SYMBOL is used to comment out
// lines in the assya input file
#define	COMMENT_SYMBOL						'#'

// What assay are we simulating?
enum {ASSAY_PCR, ASSAY_PROBE, ASSAY_PADLOCK, ASSAY_MIPS, ASSAY_AFFYMETRIX, ASSAY_NONE};

// We need to protect any commas that appear in template variables
// if we are going to combine them with X Macros
#define SINGLE_ARG(...) __VA_ARGS__

#define		INVALID_INDEX		0xFFFFFFFFFFFFFFFF

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

		inline void init()
		{
			id = degen_id = -1;
			seq_index = -1;
			
			name_str_index = INVALID_INDEX;
			forward_oligo_str_index = INVALID_INDEX;
			reverse_oligo_str_index = INVALID_INDEX;
			probe_oligo_str_index = INVALID_INDEX;
			amplicon_def_str_index = INVALID_INDEX;
			amplicon_str_index = INVALID_INDEX;
			forward_align_str_index = INVALID_INDEX;
			reverse_align_str_index = INVALID_INDEX;
			probe_align_str_index = INVALID_INDEX;

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
			VARIABLE(size_t, name_str_index) \
			VARIABLE(size_t, forward_oligo_str_index) \
			VARIABLE(size_t, reverse_oligo_str_index) \
			VARIABLE(size_t, probe_oligo_str_index) \
			VARIABLE(size_t, amplicon_def_str_index) \
			VARIABLE(size_t, amplicon_str_index) \
			VARIABLE(size_t, forward_align_str_index) \
			VARIABLE(size_t, reverse_align_str_index) \
			VARIABLE(size_t, probe_align_str_index) \
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
			VARIABLE(int8_t, primer_strand) \
			VARIABLE(int8_t, probe_strand) \
			VARIABLE(int8_t, forward_primer_clamp) \
			VARIABLE(int8_t, reverse_primer_clamp) \
			VARIABLE(int, forward_degen) \
			VARIABLE(int, reverse_degen) \
			VARIABLE(int, probe_degen) \
			VARIABLE(int8_t, forward_mm) \
			VARIABLE(int8_t, forward_gap) \
			VARIABLE(int8_t, reverse_mm) \
			VARIABLE(int8_t, reverse_gap) \
			VARIABLE(int8_t, probe_mm) \
			VARIABLE(int8_t, probe_gap)

		#define VARIABLE(A, B) A B;
			HYBRID_SIG_MEMBERS
		#undef VARIABLE

		hybrid_sig() // default constructor
		{
			init();
		};
		
		// forward + reverse + probe
		hybrid_sig(const size_t &m_name_str_index, const size_t &m_forward_str_index, 
			const size_t &m_reverse_str_index, const size_t &m_probe_str_index, 
			const unsigned int &m_id)
		{
			init();

			id = degen_id = m_id;
			
			name_str_index = m_name_str_index;
			forward_oligo_str_index = m_forward_str_index;
			reverse_oligo_str_index = m_reverse_str_index;
			probe_oligo_str_index = m_probe_str_index;
		}
		
		// forward + reverse
		hybrid_sig(const size_t &m_name_str_index, const size_t &m_forward_str_index, 
			const size_t &m_reverse_str_index, const unsigned int &m_id)
		{
			init();

			id = degen_id = m_id;
			
			name_str_index = m_name_str_index;
			forward_oligo_str_index = m_forward_str_index;
			reverse_oligo_str_index = m_reverse_str_index;
		}
		
		// probe
		hybrid_sig(const size_t &m_name_str_index, const size_t &m_probe_str_index, 
			const unsigned int &m_id)
		{
			init();

			id = degen_id = m_id;
			
			name_str_index = m_name_str_index;
			probe_oligo_str_index = m_probe_str_index;
		}
				
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
				THROW(__FILE__ ":my_id: m_id < 0");
			}

			id = m_id;
		};

		inline void my_degen_id(const int &m_id)
		{
			// Id's must be >= 0
			if(m_id < 0){
				THROW(__FILE__ ":my_degen_id: m_id < 0");
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
			return ( (forward_oligo_str_index != INVALID_INDEX) && (reverse_oligo_str_index != INVALID_INDEX) );
		};
		
		inline bool has_probe() const
		{
			return (probe_oligo_str_index != INVALID_INDEX);
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
		
		inline size_t assay_string_index() const
		{
			return name_str_index;
		};

		void reindex_str(const std::unordered_map<size_t, size_t> &m_old_to_new)
		{
			std::unordered_map<size_t, size_t>::const_iterator iter;

			#define REINDEX(VAR) \
				iter = m_old_to_new.find(VAR); \
				if(iter == m_old_to_new.end()){ \
					THROW(__FILE__ ":hybrid_sig::reindex_str: Unable to reindex " # VAR); \
				} \
				VAR = iter->second;

			REINDEX(name_str_index);
			REINDEX(forward_oligo_str_index);
			REINDEX(reverse_oligo_str_index);
			REINDEX(probe_oligo_str_index);
			REINDEX(amplicon_def_str_index);
			REINDEX(amplicon_str_index);
			REINDEX(forward_align_str_index);
			REINDEX(reverse_align_str_index);
			REINDEX(probe_align_str_index);
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

#include <iostream>

inline std::string index_to_str(const size_t &m_index, const std::vector<std::string> &m_str_table)
{
	if(m_index == INVALID_INDEX){
		// For backwards compatibility with the original codebase, return an empty string()
		return std::string();
	}

	if(m_index >= m_str_table.size() ){

		// DEBUG
		std::cerr << "m_index = " << m_index << std::endl;
		std::cerr << "|table| = " << m_str_table.size() << std::endl;

		THROW(__FILE__ ":index_to_str: Index out of bounds");
	}

	return m_str_table[m_index];
}

inline size_t str_to_index(const std::string &m_str, std::unordered_map<std::string, size_t> &m_str_table)
{
	std::unordered_map<std::string, size_t>::const_iterator iter = m_str_table.find(m_str);

	if( iter == m_str_table.end() ){

		// Add the string to the database
		const size_t ret = m_str_table.size();
		
		m_str_table[m_str] = ret;

		return ret;
	}

	return iter->second;
}

inline std::vector<std::string> ordered_keys(const std::unordered_map<std::string, size_t> &m_str_table)
{
	std::vector<std::string> ret( m_str_table.size() );

	for(std::unordered_map<std::string, size_t>::const_iterator i = m_str_table.begin();i != m_str_table.end();++i){
		
		if( i->second >= ret.size() ){
			THROW(__FILE__ ":ordered_keys: Key index is out of bounds");
		}

		ret[i->second] = i->first;
	}

	return ret;
}

size_t read_input_file(const std::string &m_file, std::vector<hybrid_sig> &m_sig_list, 
	const bool &m_ignore_probe, const bool &m_force_probe, std::unordered_map<std::string, size_t> &m_str_table);

#endif // __SIGNATURE_CLASS
