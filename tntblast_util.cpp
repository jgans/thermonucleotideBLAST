#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include "tntblast.h"

#include "degenerate_na.h"
#include "primer.h"
#include "throw.h"

using namespace std;

// Global variables
extern int mpi_numtasks;
extern int mpi_rank;

inline string toupper(const string &m_str)
{
	string ret(m_str);

	for(string::iterator i = ret.begin();i != ret.end();++i){
		*i = toupper(*i);
	}

	return ret;
}

// Allow sorting objects of type hybrid_sig by assay oligos to remove redundant assays that
// may have been created by multiplex assay expansion. Keep in mind that the forward and reverse
// primer oligos can be reversed and still yeild the *same* assay!
struct sort_by_seq // < using oligo sequences
{
	inline bool operator()(const hybrid_sig &m_a, const hybrid_sig &m_b) const
	{

		// Impose a consistent order on all assay oligo
		size_t F_a = m_a.forward_oligo_str_index;
		size_t R_a = m_a.reverse_oligo_str_index;

		if(F_a < R_a){
			swap(F_a, R_a);
		}

		size_t F_b = m_b.forward_oligo_str_index;
		size_t R_b = m_b.reverse_oligo_str_index;

		if(F_b < R_b){
			swap(F_b, R_b);
		}

		if(F_a == F_b){
			
			if(R_a == R_b){
				return m_a.probe_oligo_str_index < m_b.probe_oligo_str_index;
			}

			return R_a < R_b;
		}
		
		return F_a < F_b;
	};	
};

struct compare_by_seq // == using oligo sequences
{
	inline bool operator()(const hybrid_sig &m_a, const hybrid_sig &m_b) const
	{
		// Impose a consistent order on all assay oligo
		size_t F_a = m_a.forward_oligo_str_index;
		size_t R_a = m_a.reverse_oligo_str_index;

		if(F_a < R_a){
			swap(F_a, R_a);
		}

		size_t F_b = m_b.forward_oligo_str_index;
		size_t R_b = m_b.reverse_oligo_str_index;

		if(F_b < R_b){
			swap(F_b, R_b);
		}

		if( F_a != F_b){
			return false;
		}

		if( R_a != R_b){
			return false;
		}

		return m_a.probe_oligo_str_index == m_b.probe_oligo_str_index;
	};	
};

string top_strand(const string &m_align);

void mask_binding_sites(string &m_amplicon, const hybrid_sig &m_sig, const int &m_mask, 
	const float &m_min_primer_tm, const float &m_min_probe_tm, NucCruc &m_melt,
	const float &m_forward_primer_strand, 
	const float &m_reverse_primer_strand, const float &m_probe_strand,
	const vector<string> &m_oligo_table)
{
	if(m_mask == NO_MASK){
		return;
	}
	
	if( m_sig.has_primers() ){
	
		m_melt.strand(m_forward_primer_strand, 0.0f);
		
		mask_primer_5( 
			m_amplicon, 
			index_to_str(m_sig.forward_oligo_str_index, m_oligo_table), 
			m_melt, 
			(m_mask & MASK_PRIMERS) == 0 ? false : true, 
			(m_mask & REPLACE_PRIMERS) == 0 ? false : true);
		
		m_melt.strand(m_reverse_primer_strand, 0.0f);
		
		mask_primer_3(
			m_amplicon, 
			index_to_str(m_sig.reverse_oligo_str_index, m_oligo_table), 
			m_melt, 
			(m_mask & MASK_PRIMERS) == 0 ? false : true, 
			(m_mask & REPLACE_PRIMERS) == 0 ? false : true);
	}
	
	if( (m_mask & MASK_PROBE) && m_sig.has_probe() ){
	
		m_melt.strand(m_probe_strand, 0.0f);
		
		mask_probe(
			m_amplicon, 
			index_to_str(m_sig.probe_oligo_str_index, m_oligo_table), 
			m_melt, m_min_probe_tm);
	}
}

void mask_primer_5(string &m_amp, const string &m_oligo, NucCruc &m_melt,
	const bool &m_mask, const bool &m_replace)
{
	// If we don't need to mask or replace, return right away!
	if(!m_mask && !m_replace){
		return;
	}
	
	const size_t len = m_amp.size();
	const size_t oligo_len = m_oligo.size();
	
	// If the oligo sequence binds to the amplicon, mask the binding site with lower
	// case letters
	m_melt.set_query(m_oligo);
	
	m_melt.clear_target();
	
	bool valid_base = true;
	unsigned int gap_offset = 0;

	// Load the reverse complement of the target (while avoiding any illegal bases)
	for(unsigned int i = 0;valid_base && (i < oligo_len);i++){
		
		switch(toupper(m_amp[i])){
			case 'A':
				m_melt.push_front_target(BASE::T);
				break;
			case 'T':
				m_melt.push_front_target(BASE::A);
				break;
			case 'G':
				m_melt.push_front_target(BASE::C);
				break;
			case 'C':
				m_melt.push_front_target(BASE::G);
				break;
			case 'I':
				m_melt.push_front_target(BASE::I);
				break;
			// IUPAC degenerate bases
			case 'M': // A or C
				m_melt.push_front_target(BASE::K);
				break;
			case 'R': // G or A
				m_melt.push_front_target(BASE::Y);
				break;
			case 'S': // G or C
				m_melt.push_front_target(BASE::S);
				break;
			case 'V': // G or C or A
				m_melt.push_front_target(BASE::B);
				break;
			case 'W': // A or T
				m_melt.push_front_target(BASE::W);
				break;
			case 'Y': // T or C
				m_melt.push_front_target(BASE::R);
				break;
			case 'H': // A or C or T
				m_melt.push_front_target(BASE::D);
				break;
			case 'K': // G or T
				m_melt.push_front_target(BASE::M);
				break;
			case 'D': // G or A or T
				m_melt.push_front_target(BASE::H);
				break;
			case 'B': // G or T or C
				m_melt.push_front_target(BASE::V);
				break;
			case 'N': // A or T or G or C
				m_melt.push_front_target(BASE::N);
				break;
			case '-':
				// Count any terminal gap symbols added to pad the amplicon
				++gap_offset;
				break;
			default:
				// Stop adding bases if we encounter a bad base
				valid_base = false;
				break;
				
		};
	}
	
	const unsigned int target_len = m_melt.size_target();
	
	m_melt.approximate_tm_heterodimer();
	
	// What is the range of target bases that have been bound? The 
	// maximum range is [0, oligo_len].
	pair<unsigned int, unsigned int> range = m_melt.alignment_range_target();

	range.first = gap_offset + target_len - range.first - 1;
	range.second = gap_offset + target_len - range.second - 1;
	
	if(m_replace){
		
		// Append the primer sequence to the 5' tail of the amplicon
		m_amp = m_oligo + m_amp.substr(range.first + 1, len - range.first - 1);
		
		if(m_mask){
			
			for(unsigned int j = 0;j < oligo_len;j++){
				m_amp[j] = tolower(m_amp[j]);
			}
		}
	}
	else{
		if(m_mask){

			// Mask this region of the amplicon
			for(int j = (int)range.second;j <= (int)range.first;j++){
				m_amp[j] = tolower(m_amp[j]);
			}
		}
	}
}

void mask_primer_3(string &m_amp, const string &m_oligo, NucCruc &m_melt,
	const bool &m_mask, const bool &m_replace)
{
	// If we don't need to mask or replace, return right away!
	if(!m_mask && !m_replace){
		return;
	}
	
	const size_t len = m_amp.size();
	const size_t oligo_len = m_oligo.size();

	// If the oligo sequence binds to the amplicon, mask the binding site with lower
	// case letters
	m_melt.set_query(m_oligo);
	
	m_melt.clear_target();

	size_t i;

	unsigned int gap_offset = 0;

	// Load the target (while avoiding any illegal bases)
	for(i = len - oligo_len;i < len;i++){
		
		switch( toupper(m_amp[i]) ){
			case 'A':
				m_melt.push_back_target(BASE::A);
				break;
			case 'T':
				m_melt.push_back_target(BASE::T);
				break;
			case 'G':
				m_melt.push_back_target(BASE::G);
				break;
			case 'C':
				m_melt.push_back_target(BASE::C);
				break;
			case 'I':
				m_melt.push_back_target(BASE::I);
				break;
			// IUPAC degenerate bases
			case 'M': // A or C
				m_melt.push_back_target(BASE::K);
				break;
			case 'R': // G or A
				m_melt.push_back_target(BASE::Y);
				break;
			case 'S': // G or C
				m_melt.push_back_target(BASE::S);
				break;
			case 'V': // G or C or A
				m_melt.push_back_target(BASE::B);
				break;
			case 'W': // A or T
				m_melt.push_back_target(BASE::W);
				break;
			case 'Y': // T or C
				m_melt.push_back_target(BASE::R);
				break;
			case 'H': // A or C or T
				m_melt.push_back_target(BASE::D);
				break;
			case 'K': // G or T
				m_melt.push_back_target(BASE::M);
				break;
			case 'D': // G or A or T
				m_melt.push_back_target(BASE::H);
				break;
			case 'B': // G or T or C
				m_melt.push_back_target(BASE::V);
				break;
			case 'N': // A or T or G or C
				m_melt.push_back_target(BASE::N);
				break;
			case '-':
				// Count any terminal gap symbols added to pad the amplicon
				++gap_offset;
				break;
			default:
				m_melt.clear_target();
				break;
		};
	}
	
	const unsigned int target_len = m_melt.size_target();
	
	m_melt.approximate_tm_heterodimer();
	
	// What is the range of target bases that have been bound? The 
	// maximum range is [0, oligo_len].
	pair<int, int> range = m_melt.alignment_range_target();
	
	range.first -= gap_offset;
	range.second -= gap_offset;

	if(m_replace){
		
		// Compute the reverse complement of m_oligo
		string oligo_complement(oligo_len, 'A');
		
		i = 0;
		for(string::const_reverse_iterator iter = m_oligo.rbegin();iter != m_oligo.rend();iter++, i++){
		
			oligo_complement[i] = base_complement(*iter);
		}

		// Append the primer sequence to the 3' tail of the amplicon
		m_amp = m_amp.substr(0, len - target_len + range.first) + oligo_complement;
		
		if(m_mask){
			
			// Compute the new amplicon length
			const size_t new_len = m_amp.size();
		
			for(size_t j = new_len - oligo_len;j < new_len;j++){
				m_amp[j] = tolower(m_amp[j]);
			}
		}
	}
	else{
		if(m_mask){
			
			// Mask this region of the amplicon
			const size_t stop = (len + range.second + 1) - target_len;

			for(size_t j = len - target_len + range.first;j < stop;j++){
				m_amp[j] = tolower(m_amp[j]);
			}
		}
	}
}


void mask_probe(string &m_amp, const string &m_oligo, NucCruc &m_melt, const float &m_min_tm)
{
	const size_t len = m_amp.size();
	const size_t oligo_len = m_oligo.size();
	const size_t padded_oligo_len = oligo_len + 2; // Oligo length + dangling bases
		
	// If the oligo sequence binds to the amplicon, mask the binding site with lower
	// case letters
	m_melt.set_query(m_oligo);
	
	m_melt.clear_target();
	
	size_t i;

	// Test the oligo against the plus strand of the amplicon
	for(i = 0;i < len;i++){
	
		switch(toupper(m_amp[i])){
			case 'A':
				m_melt.push_back_target(BASE::A);
				break;
			case 'T':
				m_melt.push_back_target(BASE::T);
				break;
			case 'G':
				m_melt.push_back_target(BASE::G);
				break;
			case 'C':
				m_melt.push_back_target(BASE::C);
				break;
			case 'I':
				m_melt.push_back_target(BASE::I);
				break;
			// IUPAC degenerate bases
			case 'M': // A or C
				m_melt.push_back_target(BASE::K);
				break;
			case 'R': // G or A
				m_melt.push_back_target(BASE::Y);
				break;
			case 'S': // G or C
				m_melt.push_back_target(BASE::S);
				break;
			case 'V': // G or C or A
				m_melt.push_back_target(BASE::B);
				break;
			case 'W': // A or T
				m_melt.push_back_target(BASE::W);
				break;
			case 'Y': // T or C
				m_melt.push_back_target(BASE::R);
				break;
			case 'H': // A or C or T
				m_melt.push_back_target(BASE::D);
				break;
			case 'K': // G or T
				m_melt.push_back_target(BASE::M);
				break;
			case 'D': // G or A or T
				m_melt.push_back_target(BASE::H);
				break;
			case 'B': // G or T or C
				m_melt.push_back_target(BASE::V);
				break;
			case 'N': // A or T or G or C
				m_melt.push_back_target(BASE::N);
				break;
			default:
				// We have encountered a gap or unknown base
				m_melt.clear_target();
				break;
		};		
		
		// Compute the melting temperatures
		// Tm for the plus strand
		const float tm = m_melt.approximate_tm_heterodimer();
		
		if(tm >= m_min_tm){
		
			// What is the range of target bases that have been bound? The 
			// maximum range is [0, m_melt.size_target()].
			pair<unsigned int, unsigned int> range = m_melt.alignment_range_target();
			
			const unsigned int target_len = m_melt.size_target();
			
			// Shift the range into the target coordinate system
			range.first += (unsigned int)( i - (target_len - 1) );
			range.second += (unsigned int)( i - (target_len - 1) );
			
			// Mask this region of the amplicon
			for(int j = (int)range.second;j >= (int)range.first;j--){
				m_amp[j] = tolower(m_amp[j]);
			}
		}
		
		if(m_melt.size_target() == padded_oligo_len){
			m_melt.pop_front_target();
		}
	}
	
	m_melt.clear_target();
	
	// Test the oligo against the minus strand of the amplicon
	for(i = 0;i < len;i++){

		switch( toupper(m_amp[i]) ){
			case 'A':
				m_melt.push_front_target(BASE::T);
				break;
			case 'T':
				m_melt.push_front_target(BASE::A);
				break;
			case 'G':
				m_melt.push_front_target(BASE::C);
				break;
			case 'C':
				m_melt.push_front_target(BASE::G);
				break;
			case 'I':
				m_melt.push_front_target(BASE::I);
				break;
			// IUPAC degenerate bases
			case 'M': // A or C
				m_melt.push_front_target(BASE::K);
				break;
			case 'R': // G or A
				m_melt.push_front_target(BASE::Y);
				break;
			case 'S': // G or C
				m_melt.push_front_target(BASE::S);
				break;
			case 'V': // G or C or A
				m_melt.push_front_target(BASE::B);
				break;
			case 'W': // A or T
				m_melt.push_front_target(BASE::W);
				break;
			case 'Y': // T or C
				m_melt.push_front_target(BASE::R);
				break;
			case 'H': // A or C or T
				m_melt.push_front_target(BASE::D);
				break;
			case 'K': // G or T
				m_melt.push_front_target(BASE::M);
				break;
			case 'D': // G or A or T
				m_melt.push_front_target(BASE::H);
				break;
			case 'B': // G or T or C
				m_melt.push_front_target(BASE::V);
				break;
			case 'N': // A or T or G or C
				m_melt.push_front_target(BASE::N);
				break;
			default:
				// We have encountered a gap or unknown base
				m_melt.clear_target();
				break;
		};		
		
		// Compute the melting temperatures
		// Tm for the plus strand
		const float tm = m_melt.approximate_tm_heterodimer();
					
		if(tm >= m_min_tm){

			// What is the range of target bases that have been bound? The 
			// maximum range is [0, m_melt.size_target()].
			pair<unsigned int, unsigned int> range = m_melt.alignment_range_target();
			
			// Shift the range into the target coordinate system
			range.first = (unsigned int)(i - range.first);
			range.second = (unsigned int)(i - range.second);
			
			// Mask this region of the amplicon
			for(int j = (int)range.second;j <= (int)range.first;j++){
				m_amp[j] = tolower(m_amp[j]);
			}
		}
		
		if(m_melt.size_target() == padded_oligo_len){
			m_melt.pop_back_target();
		}
	}
}

// Expand degenerate NA bases as needed
vector<hybrid_sig> expand_degenerate_signatures(const vector<hybrid_sig> &m_sig, 
	const bool &m_degen_rescale_ct, 
	const vector<string> &m_oligo_table, 
	unordered_map<string, size_t> &m_str_table)
{
	vector<hybrid_sig> ret;
	
	// Renumber all of the assay id's to ensure that 
	// a) they start at 0
	// b) are contiguous
	int id = 0;
		
	for(vector<hybrid_sig>::const_iterator iter = m_sig.begin();iter != m_sig.end();iter++){
	
		list< pair<string, string> > primers;
		
		if( (iter->forward_oligo_str_index != INVALID_INDEX) && (iter->reverse_oligo_str_index != INVALID_INDEX) ){
			primers.push_back( make_pair( 
				index_to_str(iter->forward_oligo_str_index, m_oligo_table), 
				index_to_str(iter->reverse_oligo_str_index, m_oligo_table) ) );
		}

		try{
			primers = expand_nucleic_acid(primers);
		}
		catch(const char* error){
		
			cerr << "Error expanding primers in assay: " << index_to_str(iter->name_str_index, m_oligo_table) << endl;
			THROW(error);
		}
		
		list<string> probes;
		
		try{
			if(iter->probe_oligo_str_index != INVALID_INDEX){
				probes = expand_nucleic_acid( index_to_str(iter->probe_oligo_str_index, m_oligo_table) );
			}
		}
		catch(const char* error){
		
			cerr << "Error expanding probe in assay: " << index_to_str(iter->name_str_index, m_oligo_table) << endl;
			THROW(error);
		}
		
		const size_t num_expanded_assays = primers.size() * probes.size();
		
		if(num_expanded_assays > 1){
			cout << "Expanded degenerate bases in " 
				<< index_to_str(iter->name_str_index, m_oligo_table) << " to make " 
				<< num_expanded_assays << " non-degenerate assays" << endl;
		}
		
		const int degen_forward = m_degen_rescale_ct ? 
			degeneracy( index_to_str(iter->forward_oligo_str_index, m_oligo_table) ) : 1;
		const int degen_reverse = m_degen_rescale_ct ? 
			degeneracy( index_to_str(iter->reverse_oligo_str_index, m_oligo_table) ) : 1;
		const int degen_probe = m_degen_rescale_ct ? 
			degeneracy( index_to_str(iter->probe_oligo_str_index, m_oligo_table) ) : 1;

		list< pair<string, string> >::const_iterator primer_iter;
		list<string>::const_iterator probe_iter;
		
		if( primers.empty() ){

			for(probe_iter = probes.begin();probe_iter != probes.end();probe_iter++){
				
				hybrid_sig tmp( iter->name_str_index, 
					str_to_index(*probe_iter, m_str_table), 
					iter->my_id() );

				tmp.probe_degen = degen_probe;

				// "Real" assays derived from degenerate, "virtual" assays, get an additional identifier.
				// This additional identifier is needed to make the search results unique when using
				// target fragmentation.
				tmp.my_degen_id(id++);

				ret.push_back(tmp);
			}
		}
		else{

			for(primer_iter = primers.begin();primer_iter != primers.end();primer_iter++){

				if( probes.empty() ){

					hybrid_sig tmp( iter->name_str_index, 
						str_to_index(primer_iter->first, m_str_table), 
						str_to_index(primer_iter->second, m_str_table),
						iter->my_id() );

					tmp.forward_degen = degen_forward;
					tmp.reverse_degen = degen_reverse;

					// "Real" assays derived from degenerate, "virtual" assays, get an additional identifier.
					// This additional identifier is needed to make the search results unique when using
					// target fragmentation.
					tmp.my_degen_id(id++);

					ret.push_back(tmp);
				}
				else{
					for(probe_iter = probes.begin();probe_iter != probes.end();probe_iter++){
					
						// This was the first attempt (every "real" assay created from a degenerate assay
						// got a new unique assay id).
						//hybrid_sig tmp(iter->name, primer_iter->first, primer_iter->second,
						//	*probe_iter, id++);
						
						// This is the second attempt (every "real" assay created from a degenerate assay
						// uses the assay id of its degenerate parent).
						hybrid_sig tmp( iter->name_str_index, 
							str_to_index(primer_iter->first, m_str_table), 
							str_to_index(primer_iter->second, m_str_table),
							str_to_index(*probe_iter, m_str_table), 
							iter->my_id() );

						tmp.forward_degen = degen_forward;
						tmp.reverse_degen = degen_reverse;
						tmp.probe_degen = degen_probe;

						// "Real" assays derived from degenerate, "virtual" assays, get an additional identifier.
						// This additional identifier is needed to make the search results unique when using
						// target fragmentation.
						tmp.my_degen_id(id++);

						ret.push_back(tmp);
					}
				}
			}
		}
	}
	
	return ret;
}

vector<hybrid_sig> multiplex_expansion(const vector<hybrid_sig> &m_sig, 
	const unsigned int &m_format, vector<string> &m_index_table, 
	unordered_map<string, size_t> &m_str_table)
{
	vector<hybrid_sig> ret;

	// Renumber all of the assay id's to insure that 
	// a) they start at 0
	// b) are contiguous
	int id = 0;
	
	if( (m_format == ASSAY_PADLOCK) || (m_format == ASSAY_MIPS) ){
		
		// Assume that all Padlock/MOLpcr probes are stored as forward-reverse pairs in the input file.
		//
		// MIP assay oligos are tethered togther by a linker, so no multiplexing is needed.
		// It is possible to ligate *different* MIPS to form a single DNA molecule:
		//	- When a linear DNA molecule is created by ligating two MIPS together, this molecule will not be sequenced
		//	  and could mask an intended MIPS product (false negative)
		//	- When a pair of *nested* MIPS are ligated to form a circular DNA molecule, the expected target sequence
		//	  of the outer MIP can be masked (partial negative).
		for(vector<hybrid_sig>::const_iterator i = m_sig.begin();i != m_sig.end();i++){
		
			for(vector<hybrid_sig>::const_iterator j = m_sig.begin();j != m_sig.end();j++){
		
				const string assay_name = (i == j) ?  index_to_str(i->name_str_index, m_index_table)
					: index_to_str(i->name_str_index, m_index_table) + "(5')/" + 
					  index_to_str(j->name_str_index, m_index_table) + "(3')";

				hybrid_sig tmp( 
					str_to_index(assay_name, m_str_table), 
					i->forward_oligo_str_index, j->reverse_oligo_str_index, id++);
				
				ret.push_back(tmp);
			}
		}
	}
		
	if(m_format == ASSAY_PCR){
		
		// Is there at least one  probe in the collection of input assays?
		bool has_probes = false;
		
		for(vector<hybrid_sig>::const_iterator i = m_sig.begin();i != m_sig.end();i++){
		
			// Is this a probe-only assay?
			if(i->forward_oligo_str_index == INVALID_INDEX){
			
				hybrid_sig tmp(i->name_str_index, i->probe_oligo_str_index, id++);
				continue;
			}
			
			// Count the number of PCR primer pairs that have an associated probe
			if(i->probe_oligo_str_index != INVALID_INDEX){
				has_probes = true;
			}
			
			for(vector<hybrid_sig>::const_iterator j = m_sig.begin();j != m_sig.end();j++){
		
				// Since tntblast automatically tests for PCR with the *same* primers
				// (i.e. two forward or two reverse), there is no need to enumerate
				// an assay with the same primer in both positions
				if(i->forward_oligo_str_index == j->reverse_oligo_str_index){
					continue;
				}
				
				const string assay_name = index_to_str(i->name_str_index, m_index_table) + "(F)/" + 
						index_to_str(j->name_str_index, m_index_table) + "(R)";

				hybrid_sig tmp( 
					str_to_index(assay_name, m_str_table), 
					i->forward_oligo_str_index, j->reverse_oligo_str_index, id++);
				
				ret.push_back(tmp);
			}
		}
		
		for(vector<hybrid_sig>::const_iterator i = m_sig.begin();i != m_sig.end();i++){
		
			// Is this a probe-only assay?
			if(i->forward_oligo_str_index == INVALID_INDEX){
				continue;
			}
			
			for(vector<hybrid_sig>::const_iterator j = m_sig.begin();j != m_sig.end();j++){
		
				// Since tntblast automatically tests for PCR with the *same* primers
				// (i.e. two forward or two reverse), there is no need to enumerate
				// an assay with the same primer in both positions
				if(i->forward_oligo_str_index == j->forward_oligo_str_index){
					continue;
				}
				
				const string assay_name = index_to_str(i->name_str_index, m_index_table) + "(F)/" + 
						index_to_str(j->name_str_index, m_index_table) + "(F)";

				hybrid_sig tmp(
					str_to_index(assay_name, m_str_table), 
					i->forward_oligo_str_index, j->forward_oligo_str_index, id++);
				
				ret.push_back(tmp);
			}
		}
		
		for(vector<hybrid_sig>::const_iterator i = m_sig.begin();i != m_sig.end();i++){
		
			// Is this a probe-only assay?
			if(i->forward_oligo_str_index == INVALID_INDEX){
				continue;
			}
			
			for(vector<hybrid_sig>::const_iterator j = m_sig.begin();j != m_sig.end();j++){
		
				// Since tntblast automatically tests for PCR with the *same* primers
				// (i.e. two forward or two reverse), there is no need to enumerate
				// an assay with the same primer in both positions
				if(i->reverse_oligo_str_index == j->reverse_oligo_str_index){
					continue;
				}
				
				const string assay_name = index_to_str(i->name_str_index, m_index_table) + "(R)/" + 
						index_to_str(j->name_str_index, m_index_table) + "(R)";

				hybrid_sig tmp(
					str_to_index(assay_name, m_str_table),
					i->reverse_oligo_str_index, j->reverse_oligo_str_index, id++);
				
				ret.push_back(tmp);
			}
		}
		
		if(has_probes){
			
			// Add probes to all of the assays!
			vector<hybrid_sig> ret_with_probe;
			
			// Reset the id
			id = 0;
			
			// Update the oligo table to account for the new names referenced in ret
			m_index_table = ordered_keys(m_str_table);

			for(vector<hybrid_sig>::const_iterator i = ret.begin();i != ret.end();i++){
				
				for(vector<hybrid_sig>::const_iterator j = m_sig.begin();j != m_sig.end();j++){
					
					// Allow a mixture of assay with and without probes.
					// If one or more probes is present, then *all* of the
					// resulting assays will have probes.
					if(j->probe_oligo_str_index == INVALID_INDEX){
						continue;
					}
					
					const string assay_name = index_to_str(i->name_str_index, m_index_table) + "+" + 
							index_to_str(j->name_str_index, m_index_table) + "(P)";

					hybrid_sig tmp(
						str_to_index(assay_name, m_str_table),
						i->forward_oligo_str_index, i->reverse_oligo_str_index, 
						j->probe_oligo_str_index, id++);
					
					ret_with_probe.push_back(tmp);
				}
			}
			
			// Update the oligo table (again!) to account for the assay names that now include probes
			m_index_table = ordered_keys(m_str_table);

			ret = ret_with_probe;
		}
	}

	if(m_format == ASSAY_AFFYMETRIX){

		// This is a probe-only assay, no multiplexing is needed
		return m_sig;
	}
	
	// If our input assays shared oligos in common, than we have created additional
	// redundant assay
	sort(ret.begin(), ret.end(), sort_by_seq() );

	// Remove the redundant assays
	ret.erase( unique( ret.begin(), ret.end(), compare_by_seq() ), ret.end() );

	// Reset the id
	id = 0;

	for(vector<hybrid_sig>::iterator i = ret.begin();i != ret.end();++i){

		i->my_id(id);
		i->my_degen_id(id);

		++id;
	}

	cerr << "Multiplexing has created " << ret.size() << " assays from " << m_sig.size() 
		<< " input assays" << endl;
	
	return ret;
}
	
string primer_heuristics(const string &m_primer)
{
	ASSY_HEURISTIC::PCRPrimer pcr_test(POLY_3_GC | MULTI_5_GC | NO_POLY_RUNS | NO_3_T);
	
	// Set the run length -- we should make this a user parameter!
	pcr_test.run(5);
		
	return pcr_test.error( pcr_test(m_primer) );
}

char base_complement(const char &m_base)
{
	switch( toupper(m_base) ){
		case 'A':
			return 'T';
		case 'T':
			return 'A';
		case 'G':
			return 'C';
		case 'C':
			return 'G';
	};

	return 'N';
}

#ifdef USE_MPI

void distribute_queries(const vector<hybrid_sig> &m_sig)
{
	vector<hybrid_sig>::const_iterator iter;
	
	unsigned int num_query = m_sig.size();
	
	unsigned int buffer_size = sizeof(unsigned int);
		
	// Send all queries to all nodes
	for(iter = m_sig.begin();iter != m_sig.end();iter++){
		buffer_size += mpi_size(*iter);
	}
	
	// Tell the workers the amount of memory to allocate
	MPI_Bcast(&buffer_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	
	unsigned char *buffer = new unsigned char [buffer_size];

	if(buffer == NULL){
		THROW(__FILE__ "distribute_queries: Unable to allocate send buffer");
	}
	
	unsigned char *ptr = buffer;
	
	memcpy( ptr, &num_query, sizeof(unsigned int) );
	ptr += sizeof(unsigned int);
	
	for(iter = m_sig.begin();iter != m_sig.end();iter++){
		ptr = mpi_pack(ptr, *iter);
	}
	
	MPI_Bcast(buffer, buffer_size, MPI_BYTE, 0, MPI_COMM_WORLD);
	
	delete [] buffer;
}

void receive_queries(vector<hybrid_sig> &m_sig)
{
	unsigned int buffer_size = 0;
	unsigned int num_query = 0;
	
	MPI_Bcast(&buffer_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	
	unsigned char *buffer = new unsigned char [buffer_size];
	
	if(buffer == NULL){
		THROW(__FILE__ "receive_queries: Unable to allocate send buffer");
	}
	
	MPI_Bcast(buffer, buffer_size, MPI_BYTE, 0, MPI_COMM_WORLD);
	
	unsigned char *ptr = buffer;
	
	// The number of assays to read from the buffer that the master just
	// sent us.
	memcpy( &num_query, ptr, sizeof(unsigned int) );
	ptr += sizeof(unsigned int);
	
	m_sig = vector<hybrid_sig>(num_query);
	
	for(unsigned int i = 0;i < num_query;i++){
		ptr = mpi_unpack(ptr, m_sig[i]);
	}
	
	delete  [] buffer;
}

void serve_sequence(const int &m_dest, const unsigned int &m_index, const sequence_data &m_data)
{
	int range[2];
	
	// Read the [start, stop] range
	if(MPI_Recv(range, 2, MPI_INT, m_dest, 
	   SEQ_REQUEST, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS){

	   THROW(__FILE__ ":serve_sequence: Error receiving SEQ_REQUEST");
	}
	
	// Read this sequence from the database file
	pair<string, SEQPTR> bio_seq;
	
	// The size of the sequence (including terminator and two separators
	const unsigned int seq_len = m_data.read_bio_seq(bio_seq, m_index, range[0], range[1]) + 
		SEQ_HEADER_SIZE;
	
	const unsigned int defline_size = bio_seq.first.size() + 1;
	
	// Pack the defline and sequnece in a single buffer (with room for the '\0' after the
	// defline).
	unsigned int buffer_size = defline_size + seq_len;
	
	unsigned char *buffer = new unsigned char [buffer_size];
	
	if(!buffer){
		THROW(__FILE__ ":serve_sequence: Unable to allocate buffer");
	}
	
	unsigned char *ptr = buffer;
	
	// Pack the defline
	memcpy( ptr, bio_seq.first.c_str(), defline_size*sizeof(char) );
	ptr += defline_size*sizeof(char);
	
	// Pack the sequence
	memcpy(ptr, bio_seq.second, seq_len*sizeof(SEQBASE) );
	ptr += seq_len*sizeof(SEQBASE);
	
	// Send the defline and sequence data to the worker
	if(MPI_Send(buffer, buffer_size, MPI_BYTE, m_dest, 
	   SEQ_REQUEST, MPI_COMM_WORLD) != MPI_SUCCESS){

	   THROW(__FILE__ ":serve_sequence: Error sending sequence data");
	}
	
	delete [] buffer;
	
	delete [] bio_seq.second;
}

void distribute_string_table(const unordered_map<string, size_t> &m_table)
{
	const unsigned int num_str = m_table.size();

	unsigned int buffer_size = mpi_size(num_str);

	for(unordered_map<string, size_t>::const_iterator i = m_table.begin();i != m_table.end();++i){
		buffer_size += mpi_size(i->first) + mpi_size(i->second);
	}

	MPI_Bcast(&buffer_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	unsigned char* buffer = new unsigned char[buffer_size];

	if(buffer == NULL){
		THROW(__FILE__ ":distribute_string_table: Unable to allocate buffer");
	}

	unsigned char* ptr = buffer;

	ptr = mpi_pack(ptr, num_str);

	for(unordered_map<string, size_t>::const_iterator i = m_table.begin();i != m_table.end();++i){

		ptr = mpi_pack(ptr, i->first);
		ptr = mpi_pack(ptr, i->second);
	}

	MPI_Bcast(buffer, buffer_size, MPI_BYTE, 0, MPI_COMM_WORLD);

	delete [] buffer;
}

void receive_string_table(unordered_map<string, size_t> &m_table)
{
	unsigned int num_str;
	unsigned int buffer_size;

	MPI_Bcast(&buffer_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	unsigned char* buffer = new unsigned char[buffer_size];

	if(buffer == NULL){
		THROW(__FILE__ ":receive_string_table: Unable to allocate buffer");
	}

	MPI_Bcast(buffer, buffer_size, MPI_BYTE, 0, MPI_COMM_WORLD);

	unsigned char* ptr = buffer;

	ptr = mpi_unpack(ptr, num_str);

	m_table.clear();

	for(unsigned int i = 0;i < num_str;++i){

		string key;
		size_t value;

		ptr = mpi_unpack(ptr, key);
		ptr = mpi_unpack(ptr, value);

		m_table[key] = value;
	}

	delete [] buffer;
}

#endif // USE_MPI

vector<string> synchronize_keys(const unordered_map<string, size_t> &m_table)
{
	set<string> unique_keys;

	#ifdef USE_MPI
	// Optimization for a single rank
	if(mpi_numtasks == 1){
	#endif // USE_MPI

		for(unordered_map<string, size_t>::const_iterator i = m_table.begin();i != m_table.end();++i){
			unique_keys.insert(i->first);
		}

		return vector<string>( unique_keys.begin(), unique_keys.end() );

	#ifdef USE_MPI
	}
	#endif // USE_MPI

	#ifdef USE_MPI

	const size_t max_num_str = 10000;

	for(int i = 0;i < mpi_numtasks;++i){

		if(i == mpi_rank){

			unsigned int num_chunk = m_table.size()/max_num_str + (m_table.size()%max_num_str != 0);

			MPI_Bcast(&num_chunk, 1, MPI_UNSIGNED, i, MPI_COMM_WORLD);

			unordered_map<string, size_t>::const_iterator start_iter = m_table.begin();

			for(unsigned int chunk = 0;chunk < num_chunk;++chunk){

				unsigned int num_str = 0;
				unsigned int buffer_size = mpi_size(num_str);

				for(unordered_map<string, size_t>::const_iterator j = start_iter;(j != m_table.end()) && (num_str < max_num_str) ;++j, ++num_str){
					buffer_size += mpi_size(j->first);
				}

				MPI_Bcast(&buffer_size, 1, MPI_UNSIGNED, i, MPI_COMM_WORLD);

				unsigned char *buffer = new unsigned char [buffer_size];

				if(buffer == NULL){
					THROW(__FILE__ ":synchronize_keys: Unable to allocate send buffer");
				}

				unsigned char* ptr = buffer;

				ptr = mpi_pack(ptr, num_str);

				num_str = 0;

				for(;(start_iter != m_table.end()) && (num_str < max_num_str) ;++start_iter, ++num_str){
					
					ptr = mpi_pack(ptr, start_iter->first);
					unique_keys.insert(start_iter->first);
				}

				MPI_Bcast(buffer, buffer_size, MPI_BYTE, i, MPI_COMM_WORLD);

				delete [] buffer;
			}
		}
		else{

			unsigned int num_chunk;

			MPI_Bcast(&num_chunk, 1, MPI_UNSIGNED, i, MPI_COMM_WORLD);

			for(unsigned int chunk = 0;chunk < num_chunk;++chunk){

				unsigned int num_str;
				unsigned int buffer_size;

				MPI_Bcast(&buffer_size, 1, MPI_UNSIGNED, i, MPI_COMM_WORLD);

				unsigned char *buffer = new unsigned char [buffer_size];

				if(buffer == NULL){
					THROW(__FILE__ ":synchronize_keys: Unable to allocate recv buffer");
				}

				MPI_Bcast(buffer, buffer_size, MPI_BYTE, i, MPI_COMM_WORLD);

				unsigned char* ptr = buffer;

				ptr = mpi_unpack(ptr, num_str);

				for(unsigned int j = 0;j < num_str;++j){

					string key;

					ptr = mpi_unpack(ptr, key);
					unique_keys.insert(key);
				}

				delete [] buffer;
			}
		}
	}

	// Every rank gets the same set of strings in the same order
	return vector<string>(unique_keys.begin(), unique_keys.end() );
	#endif // USE_MPI
}

unsigned int write_inverse_matches(ostream &m_fout, const sequence_data &m_data,
	set<string> &m_targets)
{
	// How many sequences to we need to read?
	size_t num_seq = m_data.size();
	
	// How many deflines did we write?
	unsigned int count = 0;
	
	for(unsigned int i = 0;i < num_seq;i++){
	
		// Read this sequence from the database file
		pair<string, SEQPTR> bio_seq;
			
		m_data.read_bio_seq(bio_seq, i);
		
		// Write this defline if it is *not* found in m_targets
		if( m_targets.find(bio_seq.first) == m_targets.end() ){
		
			m_fout << bio_seq.first << endl;
			count ++;
		}
		
		// Free the memory used to store the sequence
		delete [] bio_seq.second;
	}
	
	return count;
}

float gc_content(const string &m_seq)
{
	if(m_seq.empty() == true){
		THROW(__FILE__ ":gc_content: empty sequence!");
	}
	
	size_t len = m_seq.size();
	float gc = 0.0f;
	
	for(unsigned int i = 0;i < len;i++){
		switch(m_seq[i]){
			case 'G':
			case 'g':
			case 'C':
			case 'c':
				gc ++;
				break;
		};
	}

	return (gc/len);
};

// Convert all white space to '_'
string mask_white_space(const string &m_str)
{
	string output(m_str);
	const char replacement = '_';
	
	for(string::iterator i = output.begin();i != output.end();i++){
		
		if(isspace(*i)){
			*i = replacement;
		}
	}
	
	return output;
}

void write_alignment(ostream &fout, const string &m_prefix, 
	const string &m_alignment)
{
	// If there is not alignment to write, return immediately
	if(m_alignment == ""){
		return;
	}

	fout << m_prefix;
	
	string::const_iterator i;
	
	for(i = m_alignment.begin();i != m_alignment.end();i++){
		
		fout << *i;
		
		if(*i == '\n'){
			fout << m_prefix;
		}
	}
	
	fout << endl;
}
	
void write_annotation(std::ostream &fout, const hybrid_sig &m_sig, 
	const sequence_data &m_annot)
{	
	if(m_sig.seq_id() < 0){
		return;
	}
	
	// Look up the corresponding molecule
	const DNAMol &mol = m_annot.annot( m_sig.seq_id() );
		
	// Check all of the annotations in mol, and print all of
	// annotations that overlap the start and stop of m_sig
	int start;
	int stop;
	
	if( (m_sig.amplicon_range.first == 0) && (m_sig.amplicon_range.second == 0) ){
	
		start = m_sig.probe_range.first;
		stop = m_sig.probe_range.second;
	}
	else{
	
		start = m_sig.amplicon_range.first;
		stop = m_sig.amplicon_range.second;
	}
	
	list<GeneAnnotation>::const_iterator iter;
		
	for(iter = mol.begin();iter != mol.end();iter++){
		
		if(int( iter->start() ) > stop){
			continue;
		}
		
		if( int( iter->stop() ) < start){
			continue;
		}
			
		// This annotation overlaps the signature
		fout << "annotation";
		
		char strand = (iter->strand() == Seq_strand_plus) ? '+' : '-';
		
		switch( iter->type() ){
			case GeneAnnotation::CDS:
				fout << "(CDS)[" << iter->start() << ".." << iter->stop() << "]" << strand << ":";
				break;
			case GeneAnnotation::GENE:
				fout << "(gene)[" << iter->start() << ".." << iter->stop() << "]" << strand << ":";
				break;
			case GeneAnnotation::PSEUDO_GENE:
				fout << "(pseudo-gene)[" << iter->start() << ".." << iter->stop() << "]" << strand << ":";
				break;
			case GeneAnnotation::RNA:
				fout << "(RNA)[" << iter->start() << ".." << iter->stop() << "]" << strand << ":";
				break;
			case GeneAnnotation::tRNA:
				fout << "(tRNA)[" << iter->start() << ".." << iter->stop() << "]" << strand << ":";
				break;
			case GeneAnnotation::TFBS:
				fout << "(TFBS)[" << iter->start() << ".." << iter->stop() << "]" << strand << ":";
				break;
			case GeneAnnotation::IMP:
				fout << "(misc)[" << iter->start() << ".." << iter->stop() << "]" << strand << ":";
				break;
			case GeneAnnotation::PRIMER:
				fout << "(primer)[" << iter->start() << ".." << iter->stop() << "]" << strand << ":";
				break;
			case GeneAnnotation::USER:
				fout << "(user)[" << iter->start() << ".." << iter->stop() << "]" << strand << ":";
				break;
			case GeneAnnotation::NONE:
				fout << "(intergenic space)[" << iter->start() << ".." << iter->stop() << "]";
				break;
			default:
				fout << "(unknown)[" << iter->start() << ".." << iter->stop() << "]" << strand << ":";
		}
		
		string tmp = iter->seq_id_str();
		
		if(tmp.size() != 0){
			fout << " " << tmp;
		}
		
		tmp = iter->name();
		
		if(tmp.size() != 0){
			fout << " " << tmp;
		}
		
		tmp = iter->info(GeneAnnotation::PRODUCT);
		
		if(tmp.size() != 0){
			fout << " " << tmp;
		}
		
		tmp = iter->info(GeneAnnotation::NOTE);

		if(tmp.size() != 0){

			fout << ' ';

			// Put the note on a single line (i.e. skip \n and \r)
			for(string::const_iterator i = tmp.begin();i != tmp.end();i++){
				if( (*i != '\n') && (*i != '\r') ){
					fout << *i;
				}
			}
		}
		
		fout << endl;
	}	
}

void test_memory(const int &m_num_mb)
{
	if(m_num_mb <= 0){
		return;
	}

	// The size of a megabyte (in bytes)
	const unsigned int mb = 1024*1024;
	const unsigned int len = m_num_mb*mb;
	const unsigned char c = 255;
	
	unsigned char* ptr = new unsigned char [len];
	
	if(ptr == NULL){
		THROW("Unable to allocate memory for memory check");
	}
	
	// Make sure that we can write to every byte of allocated memory
	memset(ptr, c, len);
	
	// Make sure that we can read the values that we just wrote
	for(unsigned int i = 0;i < len;i++){
		
		if(c != ptr[i]){
			THROW("Corrupted byte found in memory check");
		}
	}
	
	delete [] ptr;
}

void select_best_match(list<hybrid_sig> &m_results)
{
	if(m_results.empty() == true){
		return;
	}

	m_results.sort( sort_by_match() );

	list<hybrid_sig>::iterator best_iter = m_results.begin();
	list<hybrid_sig>::iterator curr_iter = best_iter;

	curr_iter ++;

	while( curr_iter != m_results.end() ){

		if( ( curr_iter->my_id() == best_iter->my_id() ) && ( curr_iter->seq_id() == best_iter->seq_id() ) ){

			bool delete_curr = false;

			if(curr_iter->has_primers() == true){

				if( curr_iter->min_primer_tm() == best_iter->min_primer_tm() ){
					
					// Test the probe
					if( curr_iter->probe_tm < best_iter->probe_tm ){
						delete_curr = true;
					}
					
					// Test the other (higher tm) primer
					if( curr_iter->max_primer_tm() < best_iter->max_primer_tm() ){
						delete_curr = true;
					}
				}
				else{

					if( curr_iter->min_primer_tm() < best_iter->min_primer_tm() ){
						delete_curr = true;
					}
				}
			}
			else{ // Probe only
				if( curr_iter->probe_tm < best_iter->probe_tm ){
					delete_curr = true;
				}
			}

			if(delete_curr){

				list<hybrid_sig>::iterator tmp_iter = curr_iter;
				curr_iter++;
				
				m_results.erase(tmp_iter);
			}
			else{ // delete best_iter

				m_results.erase(best_iter);
				best_iter = curr_iter;
				curr_iter++;
			}
		}
		else{
			best_iter = curr_iter;
			curr_iter ++;
		}
	}
}

// Make the output unique (since different workers can be sent overlapping
// target sequences) and save only the highest scoring exactly overlapping matches.
// Let me be the first to admit that this function is rather kludgy (especially the
// use of the top_strand() function -- which is way too easily broken by a change in
// output format). A better way to remove redundant matches induced by target sequence
// fragmentation would be welcome!
void uniquify_results(list<hybrid_sig> &m_results, const vector<string> &m_index_table)
{
	// If we fewer than two elements, then the results list is already unique
	if(m_results.size() < 2){
		return;
	}
	
	m_results.sort( sort_by_loc() );
	
	typedef list<hybrid_sig>::iterator I;
	
	I start_iter = m_results.begin();
	I stop_iter = start_iter;
	
	list<I> reaper;
	
	// There are there possible outcomes for a comparison of two match elements, A & B:
	typedef enum {NO_MATCH, A_CONTAINS_B, B_CONTAINS_A} MatchState;
	
	while(true){
		
		if( ( stop_iter != m_results.end() ) && 
		    ( start_iter->my_id() == stop_iter->my_id() ) && 
			( start_iter->my_degen_id() == stop_iter->my_degen_id() ) &&
		    ( start_iter->seq_id() == stop_iter->seq_id() ) ){
		    
		    stop_iter ++;
		    
		    continue;
		}

		// start_iter != stop_iter
		
		list<I> valid;
		
		for(I iter = start_iter;iter != stop_iter;iter++){
			
			if(valid.empty() == true){
			
				valid.push_back(iter);
				continue;
			}
			
			// The size of the primers divided by 2 (which will be zero if only probes are present)
			const int forward_primer_len = index_to_str(iter->forward_oligo_str_index, m_index_table).size()/2;
			const int reverse_primer_len = index_to_str(iter->reverse_oligo_str_index, m_index_table).size()/2;
			
			// Inflate the alignment strings once
			const string forward_alignment = iter->has_primers() ? inflate_dna_seq( index_to_str(iter->forward_align_str_index, m_index_table) ) : string();
			const string reverse_alignment = iter->has_primers() ? inflate_dna_seq( index_to_str(iter->reverse_align_str_index, m_index_table) ) : string();
			const string probe_alignment = iter->has_primers() ? string() : inflate_dna_seq( index_to_str(iter->probe_align_str_index, m_index_table) );

			list<I> valid_update;
			
			MatchState match_status = NO_MATCH;
			
			// Check iter against all elements of valid
			for(list<I>::iterator valid_iter = valid.begin();valid_iter != valid.end();valid_iter++){
				
				MatchState same_match = NO_MATCH;
				
				// Compare by ...
				if(iter->has_primers() == true){ // Primers
					
					// The definition of "overlap" for an amplicon is a little more complicated than for a 
					// probe. To avoid the removal of valid *nested* amplicons; i.e.
					// ->..................<-
					// ->......................<-
					// ->...........................<-
					// which have been observed in some viral genomes, the primer binding sites must
					// also overlap (in addition to the amplicon sequences themselves).
					const bool primers_overlap = (abs( (int)(iter->amplicon_range.first) - 
						(int)( (*valid_iter)->amplicon_range.first) ) < forward_primer_len) &&
						(abs( (int)(iter->amplicon_range.second) - 
						(int)( (*valid_iter)->amplicon_range.second) ) < reverse_primer_len);
					
					if(primers_overlap){
					
						const string valid_forward_alignment = inflate_dna_seq( 
								index_to_str( (*valid_iter)->forward_align_str_index, m_index_table) );

						const string valid_reverse_alignment = inflate_dna_seq( 
								index_to_str( (*valid_iter)->reverse_align_str_index, m_index_table) );

						// Does iter contain (*valid_iter)?
						if( (iter->amplicon_range.first <= (*valid_iter)->amplicon_range.first) && 
						    (iter->amplicon_range.second >= (*valid_iter)->amplicon_range.second) &&
						    ( top_strand(forward_alignment).find( top_strand(valid_forward_alignment) ) != string::npos) &&
						    ( top_strand(reverse_alignment).find( top_strand(valid_reverse_alignment) ) != string::npos) ){

						    same_match = A_CONTAINS_B;
						}
						else{

							// Does (*valid_iter) contain iter?
							if( ( (*valid_iter)->amplicon_range.first <= iter->amplicon_range.first) &&
							    ( (*valid_iter)->amplicon_range.second >= iter->amplicon_range.second) &&
							    ( top_strand(valid_forward_alignment).find( top_strand(forward_alignment) ) != string::npos) &&
							    ( top_strand(valid_reverse_alignment).find( top_strand(reverse_alignment) ) != string::npos) ){

							    same_match = B_CONTAINS_A;
							}
						}
						
						// If there is also a probe associated with these assays, make sure that the probes bind to the same location
					    	// in the target sequence (otherwise, these matches are distinct).
					    	if( iter->has_probe() && (*valid_iter)->has_probe() && (iter->probe_range != (*valid_iter)->probe_range) ){
								same_match = NO_MATCH;
					    	}
					}
				}
				else{ // Compare by probe
					
					const string valid_probe_alignment = 
						inflate_dna_seq( index_to_str( (*valid_iter)->probe_align_str_index, m_index_table) );

					// Does iter contain (*valid_iter)?
					if( (iter->probe_range.first <= (*valid_iter)->probe_range.first) &&
					    (iter->probe_range.second >= (*valid_iter)->probe_range.second) &&
					    ( top_strand(probe_alignment).find( top_strand(valid_probe_alignment) ) != string::npos) ){
					    
						same_match = A_CONTAINS_B;
					}
					else{
					
						// Does (*valid_iter) contain iter?
						if( ( (*valid_iter)->probe_range.first <= iter->probe_range.first) &&
						    ( (*valid_iter)->probe_range.second >= iter->probe_range.second) &&
						    ( top_strand(valid_probe_alignment).find( top_strand(probe_alignment) ) != string::npos) ){

						    same_match = B_CONTAINS_A;
						}
					}
				}
				
				// Do we have a match?
				if(same_match == NO_MATCH){
					
					// Keep searching
					continue;
				}
				
				if(same_match == A_CONTAINS_B){
					
					// Overwrite B with A and keep searching
					*valid_iter = iter;
					match_status = A_CONTAINS_B;
				}
				else{
					
					if(same_match == B_CONTAINS_A){
						
						if(match_status != NO_MATCH){
							
							// If we get here, we have most likely found a match in a repeat rich sequence.
							// Don't throw an error, but mark this match as B_CONTAINS_A (since it has already be placed
							// in the valid list, it will be returned to the user).
							//throw "illegal match status";
						}
						
						// A is *not* valid, so don't push it back
						match_status = B_CONTAINS_A;
						
						// Stop searching, since iter is not valid
						break;
					}
				}
			}
			
			if(match_status == NO_MATCH){
				
				// iter is valid
				valid_update.push_back(iter);
			}
			
			// splice valid_update onto valid
			valid.splice(valid.begin(), valid_update);
		}
		
		// Mark all non-valid elements for deletion
		for(I iter = start_iter;iter != stop_iter;iter++){
			
			if( find(valid.begin(), valid.end(), iter) == valid.end() ){
				
				reaper.push_back(iter);
			}
		}
		
		start_iter = stop_iter;
		
		if( stop_iter == m_results.end() ){
			break;
		}
	}
	
	while(reaper.empty() == false){
	
		m_results.erase( reaper.back() );
		reaper.pop_back();
	}
}

// Extract the top strand from a double stranded sequence alignment
string top_strand(const string &m_align)
{
	string::size_type start = m_align.find("5' ");
	
	if(start == string::npos){
		THROW(__FILE__ ":top_strand: Unable to parse alignment");
	}
	
	start += 3; // strlen("5' ") == 3
	
	const string::size_type stop = m_align.find(" 3'");
	
	if(stop == string::npos){
		THROW(__FILE__ ":top_strand: Unable to parse alignment");
	}
	
	return m_align.substr(start, stop - start);
}

unsigned int probe_only_count(const vector<hybrid_sig> &m_queries)
{
	vector<hybrid_sig>::const_iterator iter;
	unsigned int count = 0;
	
	for(iter = m_queries.begin();iter != m_queries.end();iter++){
		
		if( (iter->probe_oligo_str_index != INVALID_INDEX) && 
			(iter->forward_oligo_str_index == INVALID_INDEX) && (iter->reverse_oligo_str_index == INVALID_INDEX) ){
			count ++;
		}
	}
	
	return count;
}

bool query_sched(const size_t &m_num_target, 
	const size_t &m_num_query, const size_t &m_num_worker,
	const float &m_S_div_H, const int &m_mode)
{
	if(m_mode == QUERY_SEGMENTATION_ON){
		return true;
	}
	
	if(m_mode == QUERY_SEGMENTATION_OFF){
		return false;
	}
		
	//////////////////////////////////////////////////////////////////
	// Adaptive algorithm (m_mode == QUERY_SEGMENTATION_ADAPTIVE)
	//////////////////////////////////////////////////////////////////
	if(m_num_target == 0){
	
		//throw ":query_sched: m_num_target == 0";
		return false;
	}
	
	if(m_num_query == 0){
	
		//throw ":query_sched: m_num_query == 0";
		return false;
	}
	
	if(m_num_worker == 0){
		THROW(":query_sched: m_num_worker == 0");
	}
	
	// Don't segment queries if we only have a single worker
	// (just increases the amount of communications overhead).
	if(m_num_worker == 1){
		return false;
	}

	// The following calculations initially used unsigned int to store m_num_target, m_num_worker and m_num_query.
	// GenBank NT database is now large enough that m_num_target*(a smallish number of queries or workers) will exceed
	// the maximum value for an unsigned int. Wow! Like a Y2K bug for genomics!

	// The cost of segmenting queries
	const float cost_seq = float( m_num_target*min(m_num_query, m_num_worker)*
		(1.0f + m_S_div_H*max(size_t(1), m_num_query/m_num_worker) ) )/
		min(m_num_worker, m_num_target*m_num_query);
	
	// The cost of *not* segmenting queries
	const float cost_no_seq = float( m_num_target*(1.0 + m_S_div_H*m_num_query) )/
		min(m_num_worker, m_num_target);
	
	return (cost_seq < cost_no_seq);
}
