#include "nuc_cruc.h"
#include <iostream>
#include <list>
#include <sstream>

using namespace std;

ostream& operator << (ostream &s, const NucCruc &m_melt)
{
	// Order given by enum {A = 0, C, G, T, I, E}
	const char* base_map = "ACGTI$-MRSVWYHKDBN";
	
	if( m_melt.curr_align.query_align.size() != m_melt.curr_align.target_align.size() ){
		throw __FILE__ ":operator <<: alignment size mismatch!";
	}
	
	if(m_melt.tm_mode == NucCruc::HAIRPIN){
		
		s << "5' ";
	
		deque<BASE::nucleic_acid>::const_reverse_iterator q_iter, t_iter;

		for(t_iter = m_melt.curr_align.target_align.rbegin();t_iter != m_melt.curr_align.target_align.rend();t_iter++){

			s << base_map[*t_iter];
		}

		s << endl;

		t_iter = m_melt.curr_align.target_align.rbegin();

		s << "   ";

		for(q_iter = m_melt.curr_align.query_align.rbegin();q_iter != m_melt.curr_align.query_align.rend();q_iter++,t_iter++){

			enum {NO_MATCH, MATCH, INOSINE_MATCH};

			unsigned int match = NO_MATCH;

			if( is_complemetary_base(*q_iter, *t_iter) ){

				if( (*q_iter == BASE::I) || (*t_iter == BASE::I) ){
					match = INOSINE_MATCH;
				}
				else{
					match = MATCH;
				}
			}
			else{
				match = NO_MATCH;
			}

			switch(match){
				case NO_MATCH:
					s << ' ';
					break;
				case MATCH:
					s << '|';
					break;
				case INOSINE_MATCH:
					s << '!';
					break;
			};
		}

		s << "\n3' " ;

		for(q_iter = m_melt.curr_align.query_align.rbegin();q_iter != m_melt.curr_align.query_align.rend();q_iter++){
			s << base_map[*q_iter];
		}

		s << "\nhairpin";
	}
	else{ // m_melt.tm_mode != NucCruc::HAIRPIN

		deque<BASE::nucleic_acid>::const_iterator q_iter, t_iter;

		// Modified March 28, 2020 to include unaligned 5' and 3' bases. This is useful to illustrate
		// the presence or absence of 3' primer base clamps (which strongly effect PCR efficiency)

		//cerr << "first_match.first = " << m_melt.curr_align.first_match.first << endl;
		//cerr << "first_match.second = " << m_melt.curr_align.first_match.second << endl;

		//cerr << "last_match.first = " << m_melt.curr_align.last_match.first << endl;
		//cerr << "last_match.second = " << m_melt.curr_align.last_match.second << endl;

		const int query_len = m_melt.query.size();
		const int target_len = m_melt.target.size();

		const int prefix_len = max( 0, min(m_melt.curr_align.first_match.first, 
			target_len - 1 - m_melt.curr_align.first_match.second) );

		const int suffix_len = max( 0, min(query_len - 1 - m_melt.curr_align.last_match.first, 
			m_melt.curr_align.last_match.second) );

		s << "5' ";

		for(int i = 0;i < prefix_len;++i){
			s << base_map[ m_melt.query[m_melt.curr_align.first_match.first - prefix_len + i] ];
		}

		for(q_iter = m_melt.curr_align.query_align.begin();q_iter != m_melt.curr_align.query_align.end();q_iter++){
			s << base_map[*q_iter];
		}

		for(int i = 0;i < suffix_len;++i){
			s << base_map[ m_melt.query[m_melt.curr_align.last_match.first + 1 + i] ];
		}

		s << " 3'" << endl;

		s << "   ";

		// The prefix is formally unaligned. However, there may be complementary
		// bases between the query and target (that were not thermodynamically favorable, and hence
		// not aligned). Indicate complementary, but unaigned, bases with a ':'.

		for(int i = 0;i < prefix_len;++i){

			//if( is_complemetary_base(m_melt.query[i], m_melt.target[m_melt.curr_align.first_match.second + prefix_len - i]) ){
			if( is_complemetary_base(m_melt.query[m_melt.curr_align.first_match.first - prefix_len + i], 
									 m_melt.target[m_melt.curr_align.first_match.second + prefix_len - i]) ){
				s << ':';
			}
			else{
				s << ' ';
			}
		}

		q_iter = m_melt.curr_align.query_align.begin();
		
		for(t_iter = m_melt.curr_align.target_align.begin();t_iter != m_melt.curr_align.target_align.end();t_iter++,q_iter++){

			enum {NO_MATCH, MATCH, INOSINE_MATCH, DEGENERATE_MATCH};

			unsigned int match = NO_MATCH;

			if( is_complemetary_base(*t_iter, *q_iter) ){

				if( (*t_iter == BASE::I) || (*q_iter == BASE::I) ){
					match = INOSINE_MATCH;
				}
				else{
					if( IS_DEGENERATE_BASE(*t_iter) | IS_DEGENERATE_BASE(*q_iter) ){
						match = DEGENERATE_MATCH;
					}
					else{
						match = MATCH;
					}
				}
			}
			else{
				match = NO_MATCH;
			}

			switch(match){
				case NO_MATCH:
					s << ' ';
					break;
				case MATCH:
					s << '|';
					break;
				case INOSINE_MATCH:
					s << '!';
					break;
				case DEGENERATE_MATCH:
					s << '*';
					break;
			};
		}

		// The suffix is formally unaligned. However, there may be complementary
		// bases between the query and target (that were not thermodynamically favorable, and hence
		// not aligned). Indicate complementary, but unaigned, bases with a ':'.
		for(int i = 0;i < suffix_len;++i){

			if( is_complemetary_base(m_melt.query[m_melt.curr_align.last_match.first + 1 + i], 
				m_melt.target[m_melt.curr_align.last_match.second - i - 1]) ){

				s << ':';
			}
			else{
				s << ' ';
			}
		}

		s << endl;

		s << "3' " ;

		// The target is stored in 5' -> 3' orientation, so we need to reverse it here
		for(int i = prefix_len;i > 0;--i){
			s << base_map[ m_melt.target[m_melt.curr_align.first_match.second + i] ];
		}

		for(t_iter = m_melt.curr_align.target_align.begin();t_iter != m_melt.curr_align.target_align.end();t_iter++){
			s << base_map[*t_iter];
		}

		for(int i = 1;i <= suffix_len;++i){
			s << base_map[ m_melt.target[m_melt.curr_align.last_match.second - i] ];
		}

		s << " 5'" << endl;

		s << "dimer";
	}
	
	s << " alignment size = " << m_melt.curr_align.query_align.size();
	
	return s;
}

string NucCruc::query_seq() const
{
	// Order given by enum {A = 0, C, G, T, I, E}
	const char* base_map = "ACGTI$-MRSVWYHKDBN";
	CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH>::const_iterator iter;
	
	stringstream sout;
	
	for(iter = query.begin();iter != query.end();iter++){
		sout << base_map[*iter];
	}
	
	return sout.str();
}

string NucCruc::target_seq() const
{

	// Order given by enum {A = 0, C, G, T, I, E}
	const char* base_map = "ACGTI$-MRSVWYHKDBN";
	CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH>::const_iterator iter;
	
	stringstream sout;
	
	for(iter = target.begin();iter != target.end();iter++){
		sout << base_map[*iter];
	}
	
	return sout.str();
}


