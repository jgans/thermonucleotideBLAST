#include "nuc_cruc.h"
#include <iostream>
#include <list>
#include <sstream>

using namespace std;

ostream& operator << (ostream &s, const NucCruc &m_melt)
{
	// Order given by enum {A = 0, C, G, T, I, E}
	const char* base_map = "ACGTI$-";
	
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

			switch(*q_iter){
				case BASE::A:
					match = (*t_iter == BASE::T) ? MATCH : NO_MATCH;
					break;
				case BASE::T:
					match = (*t_iter == BASE::A) ? MATCH : NO_MATCH;
					break;
				case BASE::G:
					match = (*t_iter == BASE::C) ? MATCH : NO_MATCH;
					break;
				case BASE::C:
					match = (*t_iter == BASE::G) ? MATCH : NO_MATCH;
					break;
				case BASE::I:
					match = (*t_iter == BASE::GAP) ? NO_MATCH : INOSINE_MATCH;
					break;
				case BASE::E:
				case BASE::GAP:
					match = NO_MATCH;
					break;
			};

			if( (*t_iter == BASE::I) && (*q_iter != BASE::GAP) ){
				match = INOSINE_MATCH;
			}

			switch(match){
				case NO_MATCH:
					s << ' ';
					break;
				case MATCH:
					s << '|';
					break;
				case INOSINE_MATCH:
					s << '*';
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
	
		s << "5'-";

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

		for(int i = 0;i < prefix_len;++i){
			s << base_map[ m_melt.query[i] ];
		}

		for(q_iter = m_melt.curr_align.query_align.begin();q_iter != m_melt.curr_align.query_align.end();q_iter++){

			s << base_map[*q_iter];
		}

		for(int i = 0;i < suffix_len;++i){
			s << base_map[ m_melt.query[m_melt.curr_align.last_match.first + 1 + i] ];
		}

		s << "-3'" << endl;

		q_iter = m_melt.curr_align.query_align.begin();

		s << "   ";

		// The the prefix is unaligned (and gets a ' ', not a '|')
		for(int i = 0;i < prefix_len;++i){
			s << ' ';
		}

		for(t_iter = m_melt.curr_align.target_align.begin();t_iter != m_melt.curr_align.target_align.end();t_iter++,q_iter++){

			enum {NO_MATCH, MATCH, INOSINE_MATCH};

			unsigned int match = NO_MATCH;

			switch(*t_iter){
				case BASE::A:
					match = (*q_iter == BASE::T) ? MATCH : NO_MATCH;
					break;
				case BASE::T:
					match = (*q_iter == BASE::A) ? MATCH : NO_MATCH;
					break;
				case BASE::G:
					match = (*q_iter == BASE::C) ? MATCH : NO_MATCH;
					break;
				case BASE::C:
					match = (*q_iter == BASE::G) ? MATCH : NO_MATCH;
					break;
				case BASE::I:
					match = (*q_iter == BASE::GAP) ? NO_MATCH : INOSINE_MATCH;
					break;
				case BASE::E:
				case BASE::GAP:
					match = NO_MATCH;
					break;
			};

			if( (*q_iter == BASE::I) && (*t_iter != BASE::GAP) ){
				match = INOSINE_MATCH;
			}

			switch(match){
				case NO_MATCH:
					s << ' ';
					break;
				case MATCH:
					s << '|';
					break;
				case INOSINE_MATCH:
					s << '*';
					break;
			};
		}

		s << endl;

		s << "3'-" ;

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

		s << "-5'" << endl;

		s << "dimer";
	}
	
	s << " alignment size = " << m_melt.curr_align.query_align.size();
	
	return s;
}

string NucCruc::query_seq() const
{
	// Order given by enum {A = 0, C, G, T, I, E}
	const char* base_map = "ACGTI$";
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
	const char* base_map = "ACGTI$";
	CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH>::const_iterator iter;
	
	stringstream sout;
	
	for(iter = target.begin();iter != target.end();iter++){
		sout << base_map[*iter];
	}
	
	return sout.str();
}

string print_alignment(const deque<BASE::nucleic_acid> &m_query, const deque<BASE::nucleic_acid> &m_target)
{
        stringstream s;
        const char* base_map = "ACGTI$-";

        s << "5' ";

        deque<BASE::nucleic_acid>::const_iterator q_iter, t_iter;

        for(q_iter = m_query.begin();q_iter != m_query.end();q_iter++){

                s << base_map[*q_iter];
        }

        s << " 3'" << endl;

        q_iter = m_query.begin();

        s << "   ";

        for(t_iter = m_target.begin();t_iter != m_target.end();t_iter++,q_iter++){

                enum {NO_MATCH, MATCH, INOSINE_MATCH};

                unsigned int match = NO_MATCH;

                switch(*t_iter){
                        case BASE::A:
                                match = (*q_iter == BASE::T) ? MATCH : NO_MATCH;
                                break;
                        case BASE::T:
                                match = (*q_iter == BASE::A) ? MATCH : NO_MATCH;
                                break;
                        case BASE::G:
                                match = (*q_iter == BASE::C) ? MATCH : NO_MATCH;
                                break;
                        case BASE::C:
                                match = (*q_iter == BASE::G) ? MATCH : NO_MATCH;
                                break;
                        case BASE::I:
                                match = (*q_iter == BASE::GAP) ? NO_MATCH : INOSINE_MATCH;
                                break;
                        case BASE::E:
                                match = (*q_iter == BASE::E) ? NO_MATCH : MATCH;
                                break;
                        case BASE::GAP:
                                // Do nothing
                                break;
                };

                if( (*q_iter == BASE::I) && (*t_iter != BASE::GAP) ){
                        match = INOSINE_MATCH;
                }

                switch(match){
                        case NO_MATCH:
                                s << ' ';
                                break;
                        case MATCH:
                                s << '|';
                                break;
                        case INOSINE_MATCH:
                                s << '*';
                                break;
                };
        }

        s << endl;

        s << "3' " ;

        for(t_iter = m_target.begin();t_iter != m_target.end();t_iter++){
                s << base_map[*t_iter];
        }

        s << " 5'" << endl;

        s << " alignment size = " << m_query.size();

        return s.str();
}
