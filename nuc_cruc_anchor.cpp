#include "nuc_cruc.h"

using namespace std;
using namespace BASE;

// Which nucleic acids are "perfect" complement matches?
bool BASE::is_complemetary_base(const nucleic_acid &m_query, const nucleic_acid &m_target)
{
	const unsigned int MASK_A = 1;
	const unsigned int MASK_T = 1 << 1;
	const unsigned int MASK_G = 1 << 2;
	const unsigned int MASK_C = 1 << 3;

	unsigned int query_mask = 0;
	unsigned int target_mask = 0;

	switch(m_query){
		case A:
			query_mask = MASK_A;
			break;
		case C:
			query_mask = MASK_C;
			break;
		case G:
			query_mask = MASK_G;
			break;
		case T:
			query_mask = MASK_T;
			break;
		case I: // Treat Inosine like 'N'
			query_mask = MASK_A | MASK_T | MASK_G | MASK_C;
			break;
		// "Virtual bases"  -- they match nothing
		case E:
		case GAP:
			query_mask = 0;
			break;
		// IUPAC degenerate bases
		case M: // A or C
			query_mask = MASK_A | MASK_C;
			break;
		case R: // G or A
			query_mask = MASK_G | MASK_A;
			break;
		case S: // G or C
			query_mask = MASK_G | MASK_C;
			break;
		case V: // G or C or A
			query_mask = MASK_G | MASK_C | MASK_A;
			break;
		case W: // A or T
			query_mask = MASK_A | MASK_T;
			break;
		case Y: // T or C
			query_mask = MASK_T | MASK_C;
			break;
		case H: // A or C or T
			query_mask = MASK_A | MASK_C | MASK_T;
			break;
		case K: // G or T
			query_mask = MASK_G | MASK_T;
			break;
		case D: // G or A or T
			query_mask = MASK_G | MASK_A | MASK_T;
			break;
		case B: // G or T or C
			query_mask = MASK_G | MASK_T | MASK_C;
			break;
		case N: // A or T or G or C
			query_mask = MASK_A | MASK_T | MASK_G | MASK_C;
			break;
		default:
			throw __FILE__ ":is_complemetary_base: Unknown query base";
	};

	// Map the target base to its *complement*
	switch(m_target){
		case A:
			target_mask = MASK_T;
			break;
		case C:
			target_mask = MASK_G;
			break;
		case G:
			target_mask = MASK_C;
			break;
		case T:
			target_mask = MASK_A;
			break;
		case I: // Treat Inosine like 'N'
			target_mask = MASK_A | MASK_T | MASK_G | MASK_C;
			break;
		// "Virtual bases" -- they match nothing
		case E:
		case GAP:
			target_mask = 0;
			break;
		// IUPAC degenerate bases
		case M: // A or C
			target_mask = MASK_T | MASK_G;
			break;
		case R: // G or A
			target_mask = MASK_C | MASK_T;
			break;
		case S: // G or C
			target_mask = MASK_C | MASK_G;
			break;
		case V: // G or C or A
			target_mask = MASK_C | MASK_G | MASK_T;
			break;
		case W: // A or T
			target_mask = MASK_T | MASK_A;
			break;
		case Y: // T or C
			target_mask = MASK_A | MASK_G;
			break;
		case H: // A or C or T
			target_mask = MASK_T | MASK_G | MASK_A;
			break;
		case K: // G or T
			target_mask = MASK_C | MASK_A;
			break;
		case D: // G or A or T
			target_mask = MASK_C | MASK_T | MASK_A;
			break;
		case B: // G or T or C
			target_mask = MASK_C | MASK_A | MASK_G;
			break;
		case N: // A or T or G or C
			target_mask = MASK_A | MASK_T | MASK_G | MASK_C;
			break;
		default:
			throw __FILE__ ":is_complemetary_base: Unknown target base";
	};

	// If any bits overlap then the two bases are complementary
	return (query_mask & target_mask);
}

// The number of 5' query bases that exactly complement the
// target. Note that the alignment must already be computed!
unsigned int NucCruc::anchor5_query() const
{
	const int target_len = (int)( target.size() );
	const int query_len = (int)( query.size() );
	
	unsigned int anchor = 0;
		
	// The start (i.e. the 5' end) of the query sequence
	int query_index = 0;
	
	// The target base the corresponds to the 5' start of the query sequence (assumming a no
	// gap extension of the existing alignment)
	int target_index = curr_align.first_match.first + curr_align.first_match.second;
	
	// If there are dangling bases on this alignment, we will need to adjust the
	// target index
	if( !curr_align.target_align.empty() && curr_align.target_align.front() == E){
		return anchor;
	}
	
	if( !curr_align.query_align.empty() && curr_align.query_align.front() == E){
		target_index --;
	}
	
	if( target_index >= target_len ){
		return anchor;
	}
	
	while(true){
		
		// Have we reached the end of either the target or query sequence?
		if( (query_index >= query_len) || (target_index < 0) ){
			return anchor;
		}
		
		if( !is_complemetary_base(query[query_index], target[target_index]) ){
			return anchor;
		}

		anchor ++;
		
		query_index ++;
		target_index --;
	}
	
	// We should never get here
	throw "NucCruc::anchor5_query() logic error";
	
	return anchor;
}

// The number of 3' target bases that exactly complement the
// query. Note that the alignment must already be computed!
unsigned int NucCruc::anchor3_target() const
{
	const int target_len = (int)( target.size() );
	const int query_len = (int)( query.size() );
	
	unsigned int anchor = 0;
		
	// The end (i.e. the 3' end) of the target sequence
	int target_index = target_len - 1;
	
	// The query base the corresponds to the 3' end of the target sequence (assumming a no
	// gap extension of the existing alignment)
	int query_index = (int)(curr_align.first_match.first + curr_align.first_match.second + 1) - target_len;
		
	// If there are dangling bases on this alignment, we will need to adjust the
	// query index
	if( !curr_align.target_align.empty() && curr_align.target_align.front() == E){
		query_index ++;
	}
	
	if( !curr_align.query_align.empty() && curr_align.query_align.front() == E){
		return anchor;
	}
		
	if(query_index < 0){
		return anchor;
	}
	
	while(true){
		
		// Have we reached the end of either the target or query sequence?
		if( (target_index < 0) || (query_index >= query_len) ){
			return anchor;
		}
		
		if( !is_complemetary_base(query[query_index], target[target_index]) ){
			return anchor;
		}

		anchor ++;
		
		query_index ++;
		target_index --;
	}
	
	// We should never get here
	throw "NucCruc::anchor3_target() logic error";
	
	return anchor;
}

// The number of 5' target bases that exactly complement the
// query. Note that the alignment must already be computed!
unsigned int NucCruc::anchor3_query() const
{
	const int target_len = (int)( target.size() );
	const int query_len = (int)( query.size() );
	
	unsigned int anchor = 0;
		
	// The end (i.e. the 3' end) of the query sequence
	int query_index = query_len - 1;
	
	// The target base the corresponds to the 3' end of the query sequence (assumming a no
	// gap extension of the existing alignment)
	int target_index = (int)(curr_align.last_match.first + curr_align.last_match.second + 1) - query_len;
	
	// If there are dangling bases on this alignment, we will need to adjust the
	// target index
	if( !curr_align.target_align.empty() && curr_align.target_align.back() == E){
		return anchor;
	}
	
	if( !curr_align.query_align.empty() && curr_align.query_align.back() == E){
		target_index ++;
	}
	
	if( (target_index >= target_len) || (target_index < 0) ){
		return anchor;
	}
	
	while(true){
		
		// Have we reached the end of either the target or query sequence?
		if( (query_index < 0) || (target_index >= target_len) ){
			return anchor;
		}
		
		if( !is_complemetary_base(query[query_index], target[target_index]) ){
			return anchor;
		}

		anchor ++;
		
		query_index --;
		target_index ++;
	}
	
	// We should never get here
	throw "NucCruc::anchor3_query() logic error";
	
	return anchor;
}

// The number of 5' target bases that exactly complement the
// query. Note that the alignment must already be computed!
unsigned int NucCruc::anchor5_target() const
{
	const int target_len = (int)( target.size() );
	const int query_len = (int)( query.size() );
	
	unsigned int anchor = 0;
		
	// The start (i.e. the 5' end) of the target sequence
	int target_index = 0;
	
	// The query base the corresponds to the 5' end of the target sequence (assumming a no
	// gap extension of the existing alignment)
	int query_index = curr_align.last_match.first + curr_align.last_match.second;
	
	// If there are dangling bases on this alignment, we will need to adjust the
	// target index
	if( !curr_align.target_align.empty() && curr_align.target_align.back() == E){
		query_index --;
	}
	
	if( !curr_align.query_align.empty() && curr_align.query_align.back() == E){
		return anchor;
	}
	
	if( query_index >= query_len ){
		return anchor;
	}
	
	while(true){
		
		// Have we reached the end of either the target or query sequence?
		if( (query_index < 0) || (target_index >= target_len) ){
			return anchor;
		}
		
		if( !is_complemetary_base(query[query_index], target[target_index]) ){
			return anchor;
		}

		anchor ++;
		
		query_index --;
		target_index ++;
	}
	
	// We should never get here
	throw "NucCruc::anchor5_target() logic error";
	
	return anchor;
}

// Does the alignment contain any non-Watson and Crick base pairs?
// If so, this function returns false. Gaps force a return value!
// Note that the alignment must already be computed!
bool NucCruc::is_watson_and_crick() const
{
	deque<nucleic_acid>::const_iterator q_iter = curr_align.query_align.begin();
	deque<nucleic_acid>::const_iterator t_iter = curr_align.target_align.begin();
	
	while( q_iter != curr_align.query_align.end() ){
		
		if( (*q_iter != E) && (*t_iter != E) ){
			
			if( !is_complemetary_base(*q_iter, *t_iter) ){
				return false;
			}
		}
		
		q_iter ++;
		t_iter ++;
	}

	return true;
}

// Return the coordinates of the first and last aligned query base
// The range is zero based (i.e. [0, N_query - 1]) and runs from 5' -> 3'.
pair<unsigned int, unsigned int> NucCruc::alignment_range_query() const
{
	return make_pair(curr_align.first_match.first, curr_align.last_match.first);
}

// Return the coordinates of the first and last aligned target base, 
// The range is zero based (i.e. [0, N_target - 1]) and runs from 5' -> 3'.
pair<unsigned int, unsigned int> NucCruc::alignment_range_target() const
{
	return make_pair(curr_align.last_match.second, curr_align.first_match.second);
}

void NucCruc::alignment_range(pair<unsigned int, unsigned int> &m_query_range,
			pair<unsigned int, unsigned int> &m_target_range) const
{
	m_query_range = alignment_range_query();
	m_target_range = alignment_range_target();
}

// Does the 5'-most base of the query have an exact match to the
// target?
bool NucCruc::match_terminal5_query() const
{
	// Compute the alignment range of both the query and target sequences
	pair<unsigned int, unsigned int> query_range;
	pair<unsigned int, unsigned int> target_range;
	
	alignment_range(query_range, target_range);
	
	// Is the target base opposite the 5'-terminal query base an exact match?
	// If there is *no* base opposite the 5'-terminal query base, then return false.	
	
	const unsigned int terminal_target_3 = target_range.second + query_range.first;
	
	// Is there a target base alignable to the 5'-terminal query base (assuming a
	// gappless extrapolation)?
	return ( terminal_target_3 >= target.size() ) ? false :
		is_complemetary_base(query.front(), target[terminal_target_3]);
}

// Does the 3'-most base of the query have an exact match to the
// target?
bool NucCruc::match_terminal3_query() const
{
	// Compute the alignment range of both the query and target sequences
	pair<unsigned int, unsigned int> query_range;
	pair<unsigned int, unsigned int> target_range;
	
	alignment_range(query_range, target_range);
	
	// Is the target base opposite the 3'-terminal query base an exact match?
	// If there is *no* base opposite the 3'-terminal query base, then return false.	
	const int terminal_target_5 = int(target_range.first) - int(query.size() 
		- query_range.second) + 1;
		
	// Is there a target base alignable to the 5'-terminal query base (assuming a
	// gappless extrapolation)?
	return (terminal_target_5 < 0) ? false :
		is_complemetary_base(query.back(), target[terminal_target_5]);
}

// Does the 5'-most base of the target have an exact match to the
// query?
bool NucCruc::match_terminal5_target() const
{
	// Compute the alignment range of both the query and target sequences
	pair<unsigned int, unsigned int> query_range;
	pair<unsigned int, unsigned int> target_range;
	
	alignment_range(query_range, target_range);
	
	// Is the query base opposite the 5'-terminal target base an exact match?
	// If there is *no* base opposite the 5'-terminal target base, then return false.	
	
	const unsigned int terminal_query_3 = query_range.second + target_range.first;
	
	// Is there a target base alignable to the 5'-terminal query base (assuming a
	// gappless extrapolation)?
	return ( terminal_query_3 >= query.size() ) ? false :
		is_complemetary_base(target.front(), query[terminal_query_3]);
}

// Does the 3'-most base of the target have an exact match to the
// query?
bool NucCruc::match_terminal3_target() const
{
	// Compute the alignment range of both the query and target sequences
	pair<unsigned int, unsigned int> query_range;
	pair<unsigned int, unsigned int> target_range;
	
	alignment_range(query_range, target_range);
	
	// Is the query base opposite the 3'-terminal target base an exact match?
	// If there is *no* base opposite the 3'-terminal target base, then return false.	
	const int terminal_query_5 = int(query_range.first) - int(target.size() 
		- target_range.second) + 1;
		
	// Is there a query base alignable to the 5'-terminal target base (assuming a
	// gappless extrapolation)?
	return (terminal_query_5 < 0) ? false :
		is_complemetary_base(target.back(), query[terminal_query_5]);
}
