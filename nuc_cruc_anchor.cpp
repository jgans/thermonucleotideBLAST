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

#include "nuc_cruc.h"

using namespace std;
using namespace BASE;

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
		
		switch(query[query_index]){
			case A:
				if(target[target_index] != T){
					return anchor;
				}

				break;
			case T:
				if(target[target_index] != A){
					return anchor;
				}
				break;
			case G:
				if(target[target_index] != C){
					return anchor;
				}
				break;
			case C:
				if(target[target_index] != G){
					return anchor;
				}
				break;
			case I:
				if(target[target_index] == GAP){
					return anchor;
				}
				break;
			case E:
			case GAP:
				return anchor;
				break;
		};

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
		
		switch(query[query_index]){
			case A:
				if(target[target_index] != T){
					return anchor;
				}

				break;
			case T:
				if(target[target_index] != A){
					return anchor;
				}
				break;
			case G:
				if(target[target_index] != C){
					return anchor;
				}
				break;
			case C:
				if(target[target_index] != G){
					return anchor;
				}
				break;
			case I:
				if(target[target_index] == GAP){
					return anchor;
				}
				break;
			case E:
			case GAP:
				return anchor;
				break;
		};

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
		
		switch(query[query_index]){
			case A:
				if(target[target_index] != T){
					return anchor;
				}

				break;
			case T:
				if(target[target_index] != A){
					return anchor;
				}
				break;
			case G:
				if(target[target_index] != C){
					return anchor;
				}
				break;
			case C:
				if(target[target_index] != G){
					return anchor;
				}
				break;
			case I:
				if(target[target_index] == GAP){
					return anchor;
				}
				break;
			case E:
			case GAP:
				return anchor;
				break;
		};

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
		
		switch(query[query_index]){
			case A:
				if(target[target_index] != T){
					return anchor;
				}

				break;
			case T:
				if(target[target_index] != A){
					return anchor;
				}
				break;
			case G:
				if(target[target_index] != C){
					return anchor;
				}
				break;
			case C:
				if(target[target_index] != G){
					return anchor;
				}
				break;
			case I:
				if(target[target_index] == GAP){
					return anchor;
				}
				break;
			case E:
			case GAP:
				return anchor;
				break;
		};

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
			
			switch(*q_iter){
				case A:
					if(*t_iter != T){
						return false;
					}
					
					break;
				case T:
					if(*t_iter != A){
						return false;
					}
					break;
				case G:
					if(*t_iter != C){
						return false;
					}
					break;
				case C:
					if(*t_iter != G){
						return false;
					}
					break;
				case I:
					if(*t_iter == GAP){
						return false;
					}
					break;
				case E:
				case GAP:
					return false;
					break;
			};
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
		( int( query.front() ) == ( int(T) - int(target[terminal_target_3]) ) );
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
		( int( query.back() ) == ( int(T) - int(target[terminal_target_5]) ) );
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
		( int( target.front() ) == ( int(T) - int(query[terminal_query_3]) ) );
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
		( int( target.back() ) == ( int(T) - int(query[terminal_query_5]) ) );
}
