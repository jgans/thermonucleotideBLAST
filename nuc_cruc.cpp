#include "nuc_cruc.h"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <sstream>

using namespace std;
using namespace BASE;

int find_loop_index(const CircleBuffer<nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q,
	const unsigned int &m_start, const unsigned int &m_len);

inline nucleic_acid resolve_degenerate(const nucleic_acid &m_target_base, const nucleic_acid &m_query_base)
{
	nucleic_acid best_target = m_target_base;

	// When a degenerate target is not complementary to the query, 
	// choose a "real" base to represent the target in the following, arbitrary, order:
	//		A, T, G, C
	switch(m_target_base){
		// Real bases
		case A: case C: case G: case T: case I: 
		// "Virtual bases"
		case E: case GAP:
			break;
		// IUPAC degenerate bases
		case M: // A or C
			switch(m_query_base){
				case T:
					best_target = A;
					break;
				case G:
					best_target = C;
					break;
				default:
					best_target = A;
					break;
			};

			break;
		case R: // G or A
			switch(m_query_base){
				case T:
					best_target = A;
					break;
				case C:
					best_target = G;
					break;
				default:
					best_target = A;
					break;
			};

			break;
		case S: // G or C
			switch(m_query_base){
				case G:
					best_target = C;
					break;
				case C:
					best_target = G;
					break;
				default:
					best_target = G;
					break;
			};

			break;
		case V: // G or C or A
			switch(m_query_base){
				case G:
					best_target = C;
					break;
				case C:
					best_target = G;
					break;
				case T:
					best_target = A;
					break;
				default:
					best_target = A;
					break;
			};

			break;
		case W: // A or T
			switch(m_query_base){
				case A:
					best_target = T;
					break;
				case T:
					best_target = A;
					break;
				default:
					best_target = A;
					break;
			};

			break;
		case Y: // T or C
			switch(m_query_base){
				case G:
					best_target = C;
					break;
				case A:
					best_target = T;
					break;
				default:
					best_target = T;
					break;
			};

			break;
		case H: // A or C or T
			switch(m_query_base){
				case T:
					best_target = A;
					break;
				case G:
					best_target = C;
					break;
				case A:
					best_target = T;
					break;
				default:
					best_target = A;
					break;
			};

			break;
		case K: // G or T
			switch(m_query_base){
				case C:
					best_target = G;
					break;
				case A:
					best_target = T;
					break;
				default:
					best_target = T;
					break;
			};

			break;
		case D: // G or A or T
			switch(m_query_base){
				case C:
					best_target = G;
					break;
				case T:
					best_target = A;
					break;
				case A:
					best_target = T;
					break;
				default:
					best_target = A;
					break;
			};

			break;
		case B: // G or T or C
			switch(m_query_base){
				case C:
					best_target = G;
					break;
				case A:
					best_target = T;
					break;
				case G:
					best_target = C;
					break;
				default:
					best_target = T;
					break;
			};
		case N: // A or T or G or C

			switch(m_query_base){
				case A:
					best_target = T;
					break;
				case T:
					best_target = A;
					break;
				case G:
					best_target = C;
					break;
				case C:
					best_target = G;
					break;
				default:
					best_target = A;
					break;
			};
			break;
	};

	return best_target;
}

inline int best_base_pair(const nucleic_acid &m_target_base, const nucleic_acid &m_query_base)
{	
	// Currently, degenerate bases are only allowed in either the target base or
	// the query base, but not both. When we have a degenerate nucleotide in either the
	// target or query position, return the most optimistic base pair -- that is, pick the
	// complementary pair if allowed by the degeneracy.
	return BASE_PAIR(
		resolve_degenerate(m_target_base, m_query_base),
		resolve_degenerate(m_query_base, m_target_base) 
		);
}

#ifdef __DEBUG
const char *bp[] = {
			"AA", "AC", "AG", "AT", "AI", "AE", "A_", 
			"CA", "CC", "CG", "CT", "CI", "CE", "C_", 
			"GA", "GC", "GG", "GT", "GI", "GE", "G_", 
			"TA", "TC", "TG", "TT", "TI", "TE", "T_", 
			"IA", "IC", "IG", "IT", "II", "IE", "I_",
			"EA", "EC", "EG", "ET", "EI", "EE", "E_",
			"_A", "_C", "_G", "_T", "_I", "_E", "__"};
#endif // __DEBUG

NucCruc::NucCruc(const unsigned int &m_param_set /*= SANTA_LUCIA*/, const float &m_tm /* = 310.15 */)
{	
	// Are a pair of bases a watson and crick match (i.e. A & T, G & C, I & any base)
	memset( watson_and_crick, 0, NUM_BASE_PAIR*sizeof(bool) );
	
	watson_and_crick[AT] = watson_and_crick[TA] = true;
	watson_and_crick[CG] = watson_and_crick[GC] = true;
	
	watson_and_crick[AI] = watson_and_crick[IA] = true;
	watson_and_crick[TI] = watson_and_crick[IT] = true;
	watson_and_crick[GI] = watson_and_crick[IG] = true;
	watson_and_crick[CI] = watson_and_crick[IC] = true;
	watson_and_crick[II] = true;
	
	// The supplementary parameters from the TM program of Leber and Kaderali.
	// The associated paper is Leber and Kaderali et. al. Bioiformatics, 2005
	// (but these parameters only show up in the program source code).
	
	// delta G@37C = 0.3
	//param_supp[LOOP_H] = 0.0f;
	//param_supp[LOOP_S] = -0.97e-3f;
	
	// delta G@37C = 0.4
	//param_supp[BULGE_H] = 0.0f;
	//param_supp[BULGE_S] = -1.3e-3f;
	
	// delta G@37C = 1.75
	//param_supp[TERMINAL_MATCH_AT_H] = -2.66f;
	//param_supp[TERMINAL_MATCH_AT_S] = -14.22e-3f;
	
	// delta G@37C = 1.75
	//param_supp[TERMINAL_MATCH_GC_H] = -2.66f;
	//param_supp[TERMINAL_MATCH_GC_S] = -14.22e-3f;
	
	// delta G@37C = 1.75
	//param_supp[TERMINAL_MATCH_I_H] = -2.66f;
	//param_supp[TERMINAL_MATCH_I_S] = -14.22e-3f;
	
	// delta G@37C = 2.0
	//param_supp[TERMINAL_MISMATCH_H] = 0.0f;
	//param_supp[TERMINAL_MISMATCH_S] = -6.45e-3f;
	
	//////////////////////////////////////////////////////////////////////////
	// Supplementary parameters determined by minimization of delta G
	// LOOP //////////////////////////////////////////////////////////////////
	param_supp[LOOP_H] = -5.779f;		
	param_supp[LOOP_S] = -2.330e-2f;
	
	// BULGE //////////////////////////////////////////////////////////////////
	param_supp[BULGE_H] = 5.247e-1f;	
	param_supp[BULGE_S] = 3.318e-4f;
	
	// TERMINAL MATCH AT ///////////////////////////////////////////////////////
	param_supp[TERMINAL_MATCH_AT_H] = -4.474f;
	param_supp[TERMINAL_MATCH_AT_S] = -2.091e-2f;

	// TERMINAL MATCH GC ///////////////////////////////////////////////////////
	param_supp[TERMINAL_MATCH_GC_H] = -3.000f;
	param_supp[TERMINAL_MATCH_GC_S] =  -1.318e-2f;
		
	// TERMINAL MATCH I ////////////////////////////////////////////////////////
	// For now, consider matches with inosine to be equivalent to
	// A-T matches
	param_supp[TERMINAL_MATCH_I_H] = param_supp[TERMINAL_MATCH_AT_H];
	param_supp[TERMINAL_MATCH_I_S] = param_supp[TERMINAL_MATCH_AT_S];
	
	// TERMINAL MISMATCH ///////////////////////////////////////////////////////
	param_supp[TERMINAL_MISMATCH_H] = -2.421f;
	param_supp[TERMINAL_MISMATCH_S] = -1.180e-2f;
	
	// Initialize the salt correction mutlipliers
	param_supp_salt[LOOP_SALT] = 3.08f;
	param_supp_salt[BULGE_SALT] = 0.69f;
	param_supp_salt[TERMINAL_MATCH_SALT] = 0.56f;
	param_supp_salt[TERMINAL_MISMATCH_SALT] = 1.31f;
	
	// Initialize the parameter set
	switch(m_param_set){
		case SANTA_LUCIA:
			init_param_Santa_Lucia();
			break;
		default:
			throw "NucCruc: Unknowm parameter set!";
	};
	
	target_T = -1.0f;
	
	use_dinkelbach = false;
	
	// By default, always perform full alignments
	diagonal_alignment = false;
	
	// The dynamic programming parameters will be initialized when the user sets the salt
	// concentraion
	target_T = m_tm;
	
	// Initialize the concentrations with non-physical parameters
	// to force the user to specify them!
	na_concentration = -1.0f;
	strand_concentration = -1.0f;
	
	// By default, allow dangling bases on both ends of the alignment
	enable_dangle = make_pair(true, true);
	
	// Explore at most max_dp_path_enum high scoring paths through
	// the dynamic programming matrix
	max_dp_path_enum = 16;
	
	dp_matrix = new NC_Elem [(MAX_SEQUENCE_LENGTH + 1)*(MAX_SEQUENCE_LENGTH + 1)];
	
	// There is *no* current valid alignment
	tm_mode = INVALID;	
	
	//////////////////////////////////////////////////////////////////////////////
	// DEBUG -- make sure that all thermodynamic parameters are symmetric
	#ifdef __DEBUG
	for(char i = A;i <= GAP;i++){
		for(char j = A;j <= GAP;j++){
			for(char k = A;k <= GAP;k++){
				
				for(char l = A;l <= GAP;l++){
				
					const int curr = BASE_PAIR(i, j);
					const int prev = BASE_PAIR(k, l);
					
					const int _curr = BASE_PAIR(l, k);
					const int _prev = BASE_PAIR(j, i);
					
					if(param_H[prev][curr] != param_H[_prev][_curr]){
						cerr << "i = " << int(i) << endl;
						cerr << "j = " << int(j) << endl;
						
						cerr << "k = " << int(k) << endl;
						cerr << "l = " << int(l) << endl;
						
						cerr << "param_H[" << bp[prev] << "][" << bp[curr] << "] = " 
							<< param_H[prev][curr] << endl;
						cerr << "param_H[" << bp[_prev] << "][" << bp[_curr] << "] = " 
							<< param_H[_prev][_curr] << endl;
						
						cerr << "curr = " << bp[curr] << endl;
						cerr << "_curr = " << bp[_curr] << endl;
						
						cerr << "prev = " << bp[prev] << endl;
						cerr << "_prev = " << bp[_prev] << endl;
						
						throw "H error";
					}
					
					if(param_S[prev][curr] != param_S[_prev][_curr]){
						throw "S error";
					}
					
					//if(param_loop_terminal_H[prev][curr] != param_loop_terminal_H[_prev][_curr]){
					//	throw "loop terminal H";
					//}
					
					//if(param_loop_terminal_S[prev][curr] != param_loop_terminal_S[_prev][_curr]){
					//	throw "loop terminal S";
					//}
				}
			}
		}
	}
	#endif // __DEBUG
};

void NucCruc::update_dp_param()
{

	if(target_T < 0.0f){
		throw "update_dp_param: target_T < 0";
	}
	
	if(na_concentration <= 0.0f){
		throw "update_dp_param: na_concentration <= 0";
	}
	
	// The total salt correction has the form:
	// dS_total = param_SALT * (0.5f*num_base - 1) * log(na_concentration)
	// Where num_base is the total number of aligned bases in the *duplex*. The salt correction
	// per pair of duplex bases is approximately:
	// dS_per_base = param_SALT * log(na_concentration)
	
	const float salt_correction = param_SALT*log(na_concentration);
	
	const float loop_sc = salt_correction*param_supp_salt[LOOP_SALT];
	const float bulge_sc = salt_correction*param_supp_salt[BULGE_SALT];
	const float term_match_sc = salt_correction*param_supp_salt[TERMINAL_MATCH_SALT];
	const float term_mismatch_sc = salt_correction*param_supp_salt[TERMINAL_MISMATCH_SALT];
	
	for(unsigned int i = 0;i < NUM_BASE_PAIR;i++){
		for(unsigned int j = 0;j < NUM_BASE_PAIR;j++){
			
			// For these pairs, assume that all four bases receive a salt
			// correction (the pairs for which this is not true are handled below).
			// Since each pair of bound bases is counted twice in the NN model,
			// only apply a single pair salt correction to each set of NN terms.
			delta_g[i][j] = NC_SCORE_SCALE(param_H[i][j] - 
				target_T*(param_S[i][j] + salt_correction) );
		}
	}
	
	///////////////////////////////////////////////////////////////////////////////////
	// Double mismatches and pairs with gaps. These parameters are not specified
	// by SantaLucia et. al. The score for these pairs is *always* greater than
	// or equal to zero (i.e. always energetically unfavored).
	///////////////////////////////////////////////////////////////////////////////////
    
	for(char i = BASE::A;i <= BASE::I;i++){
		for(char j = BASE::A;j <= BASE::I;j++){
			
			const int curr = BASE_PAIR(i, j);
			
			for(char k = BASE::A;k <= BASE::I;k++){

				const int prev1 = BASE_PAIR(k, BASE::GAP);
				const int prev2 = BASE_PAIR(BASE::GAP, k);
				
				if (watson_and_crick[curr]){
					
					// bulge terminal match
					if( (curr == AT) || (curr == TA) ){
						
						NC_Score tmp_dg = 
							NC_SCORE_SCALE(param_supp[TERMINAL_MATCH_AT_H] - 
								target_T*(param_supp[TERMINAL_MATCH_AT_S] + 
								term_match_sc) );
								
						delta_g[curr][prev1] = delta_g[prev1][curr] = 
						delta_g[curr][prev2] = delta_g[prev2][curr] = 
							max(0, tmp_dg);
					}
					else{ 
						
						if( (curr == GC) || (curr == CG) ){
						
							NC_Score tmp_dg = 
								NC_SCORE_SCALE(param_supp[TERMINAL_MATCH_GC_H] - 
									target_T*(param_supp[TERMINAL_MATCH_GC_S] +
									term_match_sc) );
									
							delta_g[curr][prev1] = delta_g[prev1][curr] = 
							delta_g[curr][prev2] = delta_g[prev2][curr] = 
								max(0, tmp_dg);
						}
						else{ // curr contains inosine
							
							NC_Score tmp_dg = 
								NC_SCORE_SCALE(param_supp[TERMINAL_MATCH_I_H] - 
									target_T*(param_supp[TERMINAL_MATCH_I_S] + 
									term_match_sc) );
									
							delta_g[curr][prev1] = delta_g[prev1][curr] = 
							delta_g[curr][prev2] = delta_g[prev2][curr] = 
								max(0, tmp_dg);
						}
					}
				}
				else{ // watson_and_crick[curr] == false
					
					// mismatch bulge opening -- no salt correction
					// since this pair is internal to a loop
					NC_Score tmp_dg = 
						NC_SCORE_SCALE(param_supp[TERMINAL_MISMATCH_H] - 
							target_T*(param_supp[TERMINAL_MISMATCH_S] + 
							term_mismatch_sc) );
							
					delta_g[curr][prev1] = delta_g[prev1][curr] = 
					delta_g[curr][prev2] = delta_g[prev2][curr] = 
						max(0, tmp_dg);
				}
			}
			
			for(char k = BASE::A;k <= BASE::I;k++){
				for(char l = BASE::A;l <= BASE::I;l++){
					
					const int prev = BASE_PAIR(k, l);
					
					// Is this pair a double mismatch?
					if(!watson_and_crick[curr] && !watson_and_crick[prev]){
						
						// No salt correction for pairs of mismatches
						// (by definition, this NN-pair is internal to a loop).
						NC_Score tmp_dg = 
							NC_SCORE_SCALE(param_supp[LOOP_H] - 
								target_T*(param_supp[LOOP_S] + loop_sc) );

						delta_g[curr][prev] =
							max(0, tmp_dg);
					}
					
				}
			}
		}
	}
	
	for(char i = BASE::A;i <= BASE::I; i++){
		for(char j = BASE::A;j <= BASE::I;j++){
		
			int curr = BASE_PAIR(i, BASE::GAP);
			int prev = BASE_PAIR(j, BASE::GAP);
			
			NC_Score tmp_dg = 
				NC_SCORE_SCALE(param_supp[BULGE_H] - 
					target_T*(param_supp[BULGE_S] + bulge_sc) );
			
			delta_g[curr][prev] = max(0, tmp_dg);
										
			curr = BASE_PAIR(BASE::GAP, i);
			prev = BASE_PAIR(BASE::GAP, j);
			
			tmp_dg = NC_SCORE_SCALE(param_supp[BULGE_H] - 
				target_T*(param_supp[BULGE_S] + bulge_sc) );
			
			delta_g[curr][prev] = max(0, tmp_dg);
		}
	}
}

// Compute the Smith-Waterman local alignment between two sequences; the query 
// and the target. Both the query and target sequences are assumed to be
// in 5'-3' orientation.
NC_Score NucCruc::align_dimer(const CircleBuffer<nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q,
			const CircleBuffer<nucleic_acid, MAX_SEQUENCE_LENGTH> &m_t)
{		
	// Reset the largest element pointer buffer
	max_ptr.clear();
	
	const unsigned int query_len = m_q.size();
	const unsigned int target_len = m_t.size();
	
	// Copy the target sequence into the target buffer
	for(unsigned int i = 0;i < target_len;i++){
		target_buffer[i] = m_t[i];
	}
	
	NC_Score max_score = -1;
	
	for(unsigned int i = 1;i <= query_len;i++){
	
		// Since both query and target sequences are in 5'-3' orientation, we must
		// reverse one of the sequences. Reverse the query sequence.
		const nucleic_acid query_base = m_q[query_len - i];
		const nucleic_acid prev_query_base = (i == 1) ? GAP : m_q[query_len - (i - 1)];
		
		// A B
		// C X <-- dp[i][j]
		
		NC_Elem *C_ptr = dp_matrix + i*(MAX_SEQUENCE_LENGTH + 1);
		NC_Elem *X_ptr = C_ptr + 1;
		NC_Elem *A_ptr = C_ptr - (MAX_SEQUENCE_LENGTH + 1);
		NC_Elem *B_ptr = A_ptr + 1;
		
		for(unsigned int j = 1;j <= target_len;j++, 
			A_ptr++, B_ptr++, C_ptr++, X_ptr++){
			
			const nucleic_acid target_base = target_buffer[j - 1];
			const nucleic_acid prev_target_base = ( (j == 1) ? GAP : target_buffer[j - 2] );
			
			int cur_base_pair = best_base_pair(target_base, query_base);
			
			// Match or mismatch
			int prev_base_pair = best_base_pair(prev_target_base, prev_query_base);
			
			const NC_Score dg1 = (NC_Score(0) < A_ptr->M) ? A_ptr->M - delta_g[prev_base_pair][cur_base_pair] :
				- delta_g[prev_base_pair][cur_base_pair];
			
			// gap the query
			prev_base_pair = best_base_pair(prev_target_base, GAP);
			
			const NC_Score dg2 = (NC_Score(0) < A_ptr->I_query) ? A_ptr->I_query - delta_g[prev_base_pair][cur_base_pair] :
				- delta_g[prev_base_pair][cur_base_pair];
			
			// gap the target
			prev_base_pair = best_base_pair(GAP, prev_query_base);
			
			const NC_Score dg3 = (NC_Score(0) < A_ptr->I_target) ? A_ptr->I_target - delta_g[prev_base_pair][cur_base_pair] :
				- delta_g[prev_base_pair][cur_base_pair];
			
			if(dg1 >= dg2){
				
				if(dg1 >= dg3){
				
					X_ptr->M = dg1;
					X_ptr->M_trace = im1_jm1;
					
					#ifdef ENUMERATE_PATH
					// Allow for enumeration of multiple, equal scoring paths through
					// the dynamic programing matrix	
					if(dg1 == dg2){
						X_ptr->M_trace |= i_jm1;
					}
					
					if(dg1 == dg3){
						X_ptr->M_trace |= im1_j;
					}
					#endif // ENUMERATE_PATH
				}
				else{ // dg1 < dg3
					
					X_ptr->M = dg3;
					X_ptr->M_trace = im1_j;
				}
			}
			else{ // dg1 < dg2
				
				if(dg2 >= dg3){
				
					X_ptr->M = dg2;
					X_ptr->M_trace = i_jm1;
					
					#ifdef ENUMERATE_PATH
					// Allow for enumeration of multiple, equal scoring paths through
					// the dynamic programing matrix	
					if(dg2 == dg3){
						X_ptr->M_trace |= im1_j;
					}
					#endif // ENUMERATE_PATH
				}
				else{ // dg2 < dg3
					
					X_ptr->M = dg3;
					X_ptr->M_trace = im1_j;
				}
			}
			
			cur_base_pair = best_base_pair(target_base, GAP);
			prev_base_pair = best_base_pair(prev_target_base, query_base);
			
			NC_Score insert_gap = (NC_Score(0) < C_ptr->M) ? C_ptr->M - delta_g[prev_base_pair][cur_base_pair] :
				- delta_g[prev_base_pair][cur_base_pair];
			
			prev_base_pair = best_base_pair(prev_target_base, GAP);
			
			NC_Score extend_gap = (NC_Score(0) < C_ptr->I_query) ? C_ptr->I_query - delta_g[prev_base_pair][cur_base_pair] :
				- delta_g[prev_base_pair][cur_base_pair];
						
			if(insert_gap >= extend_gap){
			
				X_ptr->I_query = insert_gap;
				X_ptr->I_query_trace = im1_jm1;
				
				#ifdef ENUMERATE_PATH
				// Allow for enumeration of multiple, equal scoring paths through
				// the dynamic programing matrix	
				if(insert_gap == extend_gap){
					
					X_ptr->I_query_trace |= i_jm1;
				}
				#endif // ENUMERATE_PATH
				
			}
			else{
				X_ptr->I_query = extend_gap;
				X_ptr->I_query_trace = i_jm1;
			}
			
			cur_base_pair = best_base_pair(GAP, query_base);
			prev_base_pair = best_base_pair(target_base, prev_query_base);
			
			insert_gap = (NC_Score(0) < B_ptr->M) ? B_ptr->M - delta_g[prev_base_pair][cur_base_pair] :
				- delta_g[prev_base_pair][cur_base_pair];
			
			prev_base_pair = best_base_pair(GAP, prev_query_base);
			
			extend_gap = (NC_Score(0) < B_ptr->I_target) ? B_ptr->I_target - delta_g[prev_base_pair][cur_base_pair] :
				- delta_g[prev_base_pair][cur_base_pair];
			
			if(insert_gap >= extend_gap){

				X_ptr->I_target = insert_gap;
				X_ptr->I_target_trace = im1_jm1;
				
				#ifdef ENUMERATE_PATH
				// Allow for enumeration of multiple, equal scoring paths through
				// the dynamic programing matrix	
				if(insert_gap == extend_gap){
					
					X_ptr->I_target_trace |= im1_j;
				}
				#endif // ENUMERATE_PATH
			}
			else{
			
				X_ptr->I_target = extend_gap;
				X_ptr->I_target_trace = im1_j;
			}
		
			// Since gaps can only decrease the score, only test the match state for new,
			// high scoring cells
			if(X_ptr->M >= max_score){
				
				#ifdef ENUMERATE_PATH
				if(X_ptr->M > max_score){
					
					max_score = X_ptr->M;
					
					max_ptr.clear();
					max_ptr.push_back(X_ptr);
				}
				else{
					if(X_ptr->M == max_score){
						max_ptr.push_back(X_ptr);
					}
				}
				#else
				max_score = X_ptr->M;
					
				max_ptr.clear();
				max_ptr.push_back(X_ptr);
				#endif // ENUMERATE_PATH
			}			
		}
	}

	return max_score;
}

// Fast, diagonal only alignment
// 1) Both the query and target sequences are assumed to be in 5'-3' orientation. 
// 2) ** No gaps are allowed ** (The only states are mismatch or perfect match)
NC_Score NucCruc::align_dimer_diagonal(const CircleBuffer<nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q,
			const CircleBuffer<nucleic_acid, MAX_SEQUENCE_LENGTH> &m_t)
{		
	// Reset the largest element pointer buffer
	max_ptr.clear();
	
	const unsigned int query_len = m_q.size();
	const unsigned int target_len = m_t.size();
	
	// We are only aligning along the diagonal. Use the minimum(query, target) length
	const unsigned int len = min(query_len, target_len);
	
	// Copy the target sequence into the target buffer
	for(unsigned int i = 0;i < target_len;i++){
		target_buffer[i] = m_t[i];
	}
	
	NC_Score max_score = -1;
	
	NC_Elem *A_ptr = dp_matrix;
	NC_Elem *X_ptr = dp_matrix + (MAX_SEQUENCE_LENGTH + 2);
	
	int cur_base_pair;
	int prev_base_pair = best_base_pair(GAP, GAP);
	
	for(unsigned int i = 1;i <= len;i++, 
		A_ptr += (MAX_SEQUENCE_LENGTH + 2), X_ptr += (MAX_SEQUENCE_LENGTH + 2),
		prev_base_pair = cur_base_pair){

		// A B
		// C X <-- dp[i][j]
		cur_base_pair = best_base_pair(target_buffer[i - 1], m_q[query_len - i]);
			
		// Match or mismatch
		X_ptr->M = (NC_Score(0) < A_ptr->M) ? A_ptr->M - delta_g[prev_base_pair][cur_base_pair] :
			- delta_g[prev_base_pair][cur_base_pair];
			
		X_ptr->M_trace = im1_jm1;
					
		// Since gaps can only decrease the score, only test the match state for new,
		// high scoring cells
		if(X_ptr->M >= max_score){

			#ifdef ENUMERATE_PATH
			if(X_ptr->M > max_score){

				max_score = X_ptr->M;

				max_ptr.clear();
				max_ptr.push_back(X_ptr);
			}
			else{
				if(X_ptr->M == max_score){
					max_ptr.push_back(X_ptr);
				}
			}
			#else
			max_score = X_ptr->M;

			max_ptr.clear();
			max_ptr.push_back(X_ptr);
			#endif // ENUMERATE_PATH
		}			
	}

	return max_score;
}

// Compute the Smith-Waterman local alignment between a query sequence and itself to
// check for hairpin secondary structure
NC_Score NucCruc::align_hairpin(const CircleBuffer<nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q)
{	
	// We are performing a hairpin alignment (this info is needed to print the alignment
	// correctly).
	tm_mode = HAIRPIN;
	
	// Reset the largest element pointer buffer
	max_ptr.clear();
	
	// The smallest sterically allowed hairpin loop + a 1 base anchor
	const unsigned int steric_limit = 3u + 1u;  
	const unsigned int query_len = m_q.size();
	
	if(query_len == 0){
		throw ":NucCruc::align_hairpin: Empty query sequence";
	}
	
	int max_stem_len = (int)query_len - (int)steric_limit;
	
	// If max_stem_len <= 0, no alignment is performed!
	
	NC_Score max_score = -1;
	
	for(int i = 1;i <= max_stem_len;i++){
	
		// Reverse the query sequence
		const nucleic_acid query_base = m_q[query_len - i];
		const nucleic_acid prev_query_base = (i == 1) ? GAP : m_q[query_len - (i - 1)];
		
		const int upper_j = max_stem_len - (i - 1);
		
		// A B
		// C X <-- dp[i][j]
		
		NC_Elem *C_ptr = dp_matrix + i*(MAX_SEQUENCE_LENGTH + 1);
		NC_Elem *X_ptr = C_ptr + 1;
		NC_Elem *A_ptr = C_ptr - (MAX_SEQUENCE_LENGTH + 1);
		NC_Elem *B_ptr = A_ptr + 1;
		
		for(int j = 0;j < upper_j;j++,
			A_ptr++, B_ptr++, C_ptr++, X_ptr++){
			
			const nucleic_acid target_base = m_q[j];
			const nucleic_acid prev_target_base = ( (j == 0) ? GAP : m_q[j - 1] );
			
			int cur_base_pair = best_base_pair(target_base, query_base);
			
			// Match or mismatch
			int prev_base_pair = best_base_pair(prev_target_base, prev_query_base);
			
			const NC_Score dg1 = (NC_Score(0) < A_ptr->M) ? A_ptr->M - delta_g[prev_base_pair][cur_base_pair] :
				- delta_g[prev_base_pair][cur_base_pair];
			
			// gap the query
			prev_base_pair = best_base_pair(prev_target_base, GAP);
			
			const NC_Score dg2 = (NC_Score(0) < A_ptr->I_query) ? A_ptr->I_query - delta_g[prev_base_pair][cur_base_pair] :
				- delta_g[prev_base_pair][cur_base_pair];
			
			// gap the target
			prev_base_pair = best_base_pair(GAP, prev_query_base);
			
			const NC_Score dg3 = (NC_Score(0) < A_ptr->I_target) ? A_ptr->I_target - delta_g[prev_base_pair][cur_base_pair] :
				- delta_g[prev_base_pair][cur_base_pair];
			
			if(dg1 >= dg2){
				
				if(dg1 >= dg3){
				
					X_ptr->M = dg1;
					X_ptr->M_trace = im1_jm1;
					
					#ifdef ENUMERATE_PATH
					// Allow for enumeration of multiple, equal scoring paths through
					// the dynamic programing matrix	
					if(dg1 == dg2){
						X_ptr->M_trace |= i_jm1;
					}
					
					if(dg1 == dg3){
						X_ptr->M_trace |= im1_j;
					}
					#endif // ENUMERATE_PATH
					
				}
				else{ // dg1 < dg3
					
					X_ptr->M = dg3;
					X_ptr->M_trace = im1_j;
				}
			}
			else{ // dg1 < dg2
				
				if(dg2 >= dg3){
				
					X_ptr->M = dg2;
					X_ptr->M_trace = i_jm1;
					
					#ifdef ENUMERATE_PATH
					// Allow for enumeration of multiple, equal scoring paths through
					// the dynamic programing matrix	
					if(dg2 == dg3){
						X_ptr->M_trace |= im1_j;
					}
					#endif // ENUMERATE_PATH
				}
				else{ // dg2 < dg3
					
					X_ptr->M = dg3;
					X_ptr->M_trace = im1_j;
				}
			}
			
			cur_base_pair = best_base_pair(target_base, GAP);
			prev_base_pair = best_base_pair(prev_target_base, query_base);
			
			NC_Score insert_gap = (NC_Score(0) < C_ptr->M) ? C_ptr->M - delta_g[prev_base_pair][cur_base_pair] :
				- delta_g[prev_base_pair][cur_base_pair];
			
			prev_base_pair = best_base_pair(prev_target_base, GAP);
			
			NC_Score extend_gap = (NC_Score(0) < C_ptr->I_query) ? C_ptr->I_query - delta_g[prev_base_pair][cur_base_pair] :
				- delta_g[prev_base_pair][cur_base_pair];
			
			if(insert_gap >= extend_gap){
			
				X_ptr->I_query = insert_gap;
				X_ptr->I_query_trace = im1_jm1;
				
				#ifdef ENUMERATE_PATH
				// Allow for enumeration of multiple, equal scoring paths through
				// the dynamic programing matrix	
				if(insert_gap == extend_gap){
					X_ptr->I_query_trace |= i_jm1;
				}
				#endif // ENUMERATE_PATH
			}
			else{
				X_ptr->I_query = extend_gap;
				X_ptr->I_query_trace = i_jm1;
			}
			
			cur_base_pair = best_base_pair(GAP, query_base);
			prev_base_pair = best_base_pair(target_base, prev_query_base);
			
			insert_gap = (NC_Score(0) < B_ptr->M) ? B_ptr->M - delta_g[prev_base_pair][cur_base_pair] :
				- delta_g[prev_base_pair][cur_base_pair];
			
			prev_base_pair = best_base_pair(GAP, prev_query_base);
			
			extend_gap = (NC_Score(0) < B_ptr->I_target) ? B_ptr->I_target - delta_g[prev_base_pair][cur_base_pair] :
				- delta_g[prev_base_pair][cur_base_pair];
			
			if(insert_gap >= extend_gap){

				X_ptr->I_target = insert_gap;
				X_ptr->I_target_trace = im1_jm1;
				
				#ifdef ENUMERATE_PATH
				// Allow for enumeration of multiple, equal scoring paths through
				// the dynamic programing matrix	
				if(insert_gap == extend_gap){
					X_ptr->I_target_trace |= im1_j;
				}
				#endif // ENUMERATE_PATH
			}
			else{
			
				X_ptr->I_target = extend_gap;
				X_ptr->I_target_trace = im1_j;
			}
			
			// Since gaps can only decrease the score, only test the match state for new,
			// high scoring cells
			if(X_ptr->M >= max_score){
								
				#ifdef ENUMERATE_PATH
				if(X_ptr->M > max_score){
					
					max_score = X_ptr->M;
					
					max_ptr.clear();
					max_ptr.push_back(X_ptr);
				}
				else{
					if(X_ptr->M == max_score){
						max_ptr.push_back(X_ptr);
					}
				}
				#else
				max_score = X_ptr->M;
					
				max_ptr.clear();
				max_ptr.push_back(X_ptr);
				#endif // ENUMERATE_PATH
			}
		}		
	}

	return max_score;
}

void NucCruc::enumerate_dimer_alignments(NC_Elem *m_dp_matrix, NC_Elem *m_max_ptr,
			alignment &m_best_align,
			const CircleBuffer<nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q,
			const CircleBuffer<nucleic_acid, MAX_SEQUENCE_LENGTH> &m_t,
			const mode &m_mode)
{
	bool first_time = true;
	
	deque<trace_branch> trace_stack;
	int zero_count = -1;
	
	// Keep track of the number of alternate paths through the DP matrix
	// that we have explored
	unsigned int trace_count = 0;
	
	float best_dg = m_best_align.dH - target_T*m_best_align.dS;;
	
	const int query_len = int( m_q.size() );
	const int target_len = int( m_t.size() );
		
	while(true){

		// Test for zero_count <= 0 to catch 
		// 1) single paths with any number of zeros (zero_count == 0 at the end)
		// 2) multiple paths that will have zero_count == -1 after the last path
		//    has been completed
		if( !first_time && trace_stack.empty() && (zero_count <= 0) ){
			break;
		}
		
		// Stop searching once we've exceeded the maximum numnber
		// of paths through the DP matrix
		if( (max_dp_path_enum != 0) && (max_dp_path_enum < trace_count) ){
			break;
		}
		
		trace_count ++;
		
		first_time = false;
		
		alignment local_align;
				
		trace_back(m_dp_matrix, m_max_ptr, trace_stack, zero_count, local_align, m_q, m_t);
				
		////////////////////////////////////////////////////////////////////////////////////////////////////
		// Trim the ends of the alignment to remove "frayed ends" -- i.e. those ends that don't start
		// with a Watson and Crick base pair
		//
		// Trim the back of the alignment
		while( (local_align.query_align.empty() == false) &&
			!watson_and_crick[best_base_pair( local_align.query_align.back(), local_align.target_align.back() )] ){
			
			// If this pair does not include a gap, we need to adjust the last match as well
			if( !IS_VIRTUAL_BASE(local_align.query_align.back() ) ){
				--local_align.last_match.first;
			}
			
			if( !IS_VIRTUAL_BASE(local_align.target_align.back() ) ){

				++local_align.last_match.second;
			}
						
			local_align.query_align.pop_back();
			local_align.target_align.pop_back();			
		}
				
		// Trim the front of the alignment
		while( (local_align.query_align.empty() == false) &&
			!watson_and_crick[best_base_pair(local_align.query_align.front(), local_align.target_align.front())] ){
			
			// If this pair does not include a gap, we need to adjust the first match as well
			if( !IS_VIRTUAL_BASE(local_align.query_align.front() ) ){
				++local_align.first_match.first;
			}
			
			if( !IS_VIRTUAL_BASE(local_align.target_align.front() ) ){
				--local_align.first_match.second;
			}
			
			local_align.query_align.pop_front();
			local_align.target_align.pop_front();			
		}
		
		//////////////////////////////////////////////////////////////////////////////////////////////
		// If we've enumerated all of the zeros in the current path, update the trace stack
		if( (zero_count == 0) && !trace_stack.empty() ){
			
			while( !trace_stack.empty() && (trace_stack.back().next_trace() == false) ){
				trace_stack.pop_back();
			}
			
			// Reset the zero count to -1 to force a recalculation of the number of zeros in this new path
			zero_count = -1;
		}
		
		//////////////////////////////////////////////////////////////////////////////////////////////
		// How to handle dangling end bases. There are four possibilities to consider for each end:
		// 1) A blunt end -> no modification
		//	5'-AG...
		//	   ||
		//	3'-TC...
		// 2) 5' over-hang -> add 3' dangling end virtual base
		//	5'-ACG...
		//	    ||
		//	3'- TC...
		// 3) 3' over-hang -> add 5' dangling end virtual base
		//	5'- AG...
		//	    ||
		//	3'-CTC...
		// 4) frayed end  -> add the mismatched bases
		//	5'-CAG...
		//	    ||
		//	3'-CTC...
		
		// Do we need to attach bases to the front (5-query) of the alignment?
		if( enable_dangle.first && 
			( (local_align.first_match.first != 0) || ( local_align.first_match.second != (target_len - 1) ) ) ){
			
			// We've ruled out a blunt end
			if(local_align.first_match.first == 0){
				// 3' over-hang
				local_align.query_align.push_front(E);
			}
			else{
				// 5' over-hang or frayed end
				local_align.first_match.first --;
				local_align.query_align.push_front( m_q[local_align.first_match.first] );
			}
			
			if( local_align.first_match.second == (target_len - 1) ){
				// 3' over-hang
				local_align.target_align.push_front(E);
			}
			else{
				// 5' over-hang or frayed end
				local_align.first_match.second ++;
				local_align.target_align.push_front( m_t[local_align.first_match.second] );
			}
		}
				
		// Do we need to attach bases to the back (3-query) of the alignment?
		if( enable_dangle.second && 
			( (local_align.last_match.first != (query_len - 1) ) || (local_align.last_match.second != 0) ) ){
			
			// We've ruled out a blunt end
			if( local_align.last_match.first == (query_len - 1) ){
				// 5' over-hang
				local_align.query_align.push_back(E);
			}
			else{
				// 3' over-hang or frayed end
				local_align.last_match.first ++;
				local_align.query_align.push_back( m_q[local_align.last_match.first] );
			}
			
			if(local_align.last_match.second == 0){
				// 3' over-hang
				local_align.target_align.push_back(E);
			}
			else{
				// 3' over-hang or frayed end
				local_align.last_match.second --;
				local_align.target_align.push_back( m_t[local_align.last_match.second] );
			}
		}
		
		const unsigned int align_size = local_align.query_align.size();

		#ifdef __DEBUG
		if( align_size != local_align.target_align.size() ){
			throw __FILE__ ":NucCruc::tm_dimer: alignment size mismatch!";
		}
		#endif //__DEBUG

		// Do we have enough bases to melt?
		if(align_size < 3){
			continue;
		}
			
		if(evaluate_alignment(local_align, m_mode) == true){
		
			// Save this alignment if it is better than the current
			// best alignment. 
			// *********************************
			// * Use delta G to rank alignments*
			// *********************************
			
			const float local_dg = local_align.dH - target_T*local_align.dS;
			
			if(!m_best_align.valid || (local_dg < best_dg) ){

				m_best_align = local_align;
				m_best_align.valid = true;
				best_dg = local_dg;				
			}
		}
	}
}

void NucCruc::enumerate_hairpin_alignments(NC_Elem *m_dp_matrix, NC_Elem *m_max_ptr, 
			alignment &m_best_align, const CircleBuffer<BASE::nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q)
{
	const unsigned int min_hairpin_size = 3;
	
	bool first_time = true;
	
	deque<trace_branch> trace_stack;
	int zero_count = -1;
	
	// Keep track of the number of alternate paths through the DP matrix
	// that we have explored
	unsigned int trace_count = 0;
	
	float best_dg = m_best_align.dH - target_T*m_best_align.dS;
	
	const int query_len = int( m_q.size() );
	
	while(true){
		
		// Test for zero_count <= 0 to catch 
		// 1) single paths with any number of zeros (zero_count == 0 at the end)
		// 2) multiple paths that will have zero_count == -1 after the last path
		//    has been completed
		if( !first_time && trace_stack.empty() && (zero_count <= 0) ){
			break;
		}
		
		// Stop searching once we've exceeded the maximum number
		// of paths through the DP matrix
		if( (max_dp_path_enum != 0) && (max_dp_path_enum < trace_count) ){
			break;
		}
		
		trace_count ++;
		
		first_time = false;
		
		alignment local_align;
				
		trace_back(m_dp_matrix, m_max_ptr, trace_stack, zero_count, local_align, m_q, m_q);
				
		////////////////////////////////////////////////////////////////////////////////////////////////////
		// Trim the ends of the alignment to remove "frayed ends" -- i.e. those ends that don't start
		// with a Watson and Crick base pair
		//
		// Trim the back of the alignment
		while( (local_align.query_align.empty() == false) &&
			!watson_and_crick[best_base_pair( local_align.query_align.back(), local_align.target_align.back() )] ){
			
			// If this pair does not include a gap, we need to adjust the last match as well
			if( !IS_VIRTUAL_BASE(local_align.query_align.back() ) ){
				local_align.last_match.first --;
			}
			
			if( !IS_VIRTUAL_BASE(local_align.target_align.back() ) ){
				local_align.last_match.second ++;
			}
						
			local_align.query_align.pop_back();
			local_align.target_align.pop_back();
		}
				
		// Trim the front of the alignment
		while( (local_align.query_align.empty() == false) &&
			!watson_and_crick[best_base_pair(local_align.query_align.front(), local_align.target_align.front())] ){
			
			// If this pair does not include a gap, we need to adjust the first match as well
			if( !IS_VIRTUAL_BASE( local_align.query_align.front() ) ){
				local_align.first_match.first ++;
			}
			
			if( !IS_VIRTUAL_BASE( local_align.target_align.front() ) ){
				local_align.first_match.second --;
			}
			
			local_align.query_align.pop_front();
			local_align.target_align.pop_front();
		}
		
		//////////////////////////////////////////////////////////////////////////////////////////////
		// If we've enumerated all of the zeros in the current path, update the trace stack
		if( (zero_count == 0) && !trace_stack.empty() ){
			
			while( !trace_stack.empty() && (trace_stack.back().next_trace() == false) ){
				trace_stack.pop_back();
			}
			
			// Reset the zero count to -1 to force a recalculation of the number of zeros in this new path
			zero_count = -1;
		}
		
		// Check the alignment before we handle dangling end bases
		if( (local_align.query_align.size() >= min_hairpin_size) && evaluate_hairpin_alignment(local_align) ){
		
			// Save this alignment if it is better than the current
			// best alignment. 
			// *********************************
			// * Use delta G to rank alignments*
			// *********************************
			
			const float local_dg = local_align.dH - target_T*local_align.dS;
			
			if(!m_best_align.valid || (local_dg < best_dg) ){

				m_best_align = local_align;
				m_best_align.valid = true;
				best_dg = local_dg;
				
				/// DEBUG
				//cerr << "best_dg(0) = " << best_dg << endl;
				//cerr << "best_dh(0) = " << local_align.dH << endl;
				//cerr << "best_ds(0) = " << local_align.dS << endl;
			}
		}
		
		//////////////////////////////////////////////////////////////////////////////////////////////
		// How to handle dangling end bases. There are four possibilities to consider for each end:
		// 1) A blunt end -> no modification
		//	5'-AG...
		//	   ||
		//	3'-TC...
		// 2) 5' over-hang -> add 3' dangling end virtual base
		//	5'-ACG...
		//	    ||
		//	3'- TC...
		// 3) 3' over-hang -> add 5' dangling end virtual base
		//	5'- AG...
		//	    ||
		//	3'-CTC...
		// 4) frayed end  -> add the mismatched bases
		//	5'-CAG...
		//	    ||
		//	3'-CTC...
		// Do we need to attach bases to the front (5-query) of the alignment?
		if( (local_align.last_match.second != 0) || (local_align.last_match.first != (query_len - 1) ) ){

			if(local_align.last_match.second == 0){
				local_align.target_align.push_back(E);
			}
			else{
				local_align.last_match.second --;
				local_align.target_align.push_back(m_q[local_align.last_match.second]);
			}
			
			if( local_align.last_match.first == (query_len - 1) ){
				local_align.query_align.push_back(E);
			}
			else{
				local_align.last_match.first ++;
				local_align.query_align.push_back(m_q[local_align.last_match.first]);
			}
		}
		
		const unsigned int align_size = local_align.query_align.size();

		// Do we have enough bases to melt?
		if(align_size < 3){
			continue;
		}
		
		if( (align_size >= min_hairpin_size) && evaluate_hairpin_alignment(local_align)){
		
			// Save this alignment if it is better than the current
			// best alignment. 
			// *********************************
			// * Use delta G to rank alignments*
			// *********************************
			
			const float local_dg = local_align.dH - target_T*local_align.dS;
			
			if(!m_best_align.valid || (local_dg < best_dg) ){

				m_best_align = local_align;
				m_best_align.valid = true;
				best_dg = local_dg;
				
				/// DEBUG
				//cerr << "best_dg(1) = " << best_dg << endl;
				//cerr << "best_dh(1) = " << local_align.dH << endl;
				//cerr << "best_ds(1) = " << local_align.dS << endl;
			}
		}
		
		// If the closing basepair is an "A-T" or "T-A", then it carries a penalty that we may be better off
		// without. Try removing this pair an see if the alignment score improves. Note that there still needs to be
		// 3 or more hairpin stem base pairs *after* we remove the closing base pair.
		
		if(align_size <= 3){
			continue;
		}
		
		const int last_3 = local_align.first_match.first;
		const int last_5 = local_align.first_match.second;
		const int last_base_pair = best_base_pair(query[last_5], query[last_3]);
		
		if( (last_base_pair == GC) || (last_base_pair == CG) ){
			continue;
		}
		
		// DEBUG
		//continue;
		//cerr << "** Closing AT special case ** " << endl;
		
		// If we get here, the last base pair is an AT or TA. Remove this base pair	
		local_align.first_match.first ++;
		local_align.first_match.second --;
		
		// Shorten the actual alignment
		local_align.query_align.pop_front();
		local_align.target_align.pop_front();
			
		if( evaluate_hairpin_alignment(local_align) ){

			// Save this alignment if it is better than the current
			// best alignment. 
			// *********************************
			// * Use delta G to rank alignments*
			// *********************************

			const float local_dg = local_align.dH - target_T*local_align.dS;

			if(!m_best_align.valid || (local_dg < best_dg) ){

				m_best_align = local_align;
				m_best_align.valid = true;
				best_dg = local_dg;
				
				/// DEBUG
				//cerr << "best_dg(2) = " << best_dg << endl;
				//cerr << "best_dh(2) = " << local_align.dH << endl;
				//cerr << "best_ds(2) = " << local_align.dS << endl;
			}
		}
	}
}

void NucCruc::trace_back(NC_Elem *m_dp_matrix, NC_Elem *m_cell_ptr, 
	deque<trace_branch> &m_trace_stack, int &m_zero_count, 
	alignment &m_local_align,
	const CircleBuffer<nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q,
	const CircleBuffer<nucleic_acid, MAX_SEQUENCE_LENGTH> &m_t)
{	
	const int query_len = int( m_q.size() );

	int last_i = (m_cell_ptr - m_dp_matrix)/(MAX_SEQUENCE_LENGTH + 1);
	int last_j = (m_cell_ptr - m_dp_matrix)%(MAX_SEQUENCE_LENGTH + 1);
		
	// Track the first (5'-query/3'-target) and last (3'-query/5'-target) matches
	// so we can attach dangling end bases if needed
	// Assume that the first match of the alignment is an exact match
	m_local_align.first_match = make_pair(query_len - last_i, last_j - 1);
	
	int truncate_at_zero = 0;
	
	#ifdef ENUMERATE_PATH
	bool count_zeros = false;
	
	///////////////////////////////////////////////////////////////////////////////////////////////////
	// In order to enumerate all of the equally high scoring paths thought the dynamic programming
	// matrix (i.e. the equally high scoring alignments) we need to do the following:
	// 1) generate all of the alignments that stop at the highest scoring cells
	// 2) for each of these alignments, test all of the sub alignments that produced by truncating
	// alignments at cells that have a score of exactly 0.
	///////////////////////////////////////////////////////////////////////////////////////////////////
	
	// If this function is run with zero_count < 0, then it simply counts the number of
	// zero score values found along the dynamic programming path for match or mismatch cells
	// (not gaps though).
	if(m_zero_count < 0){
	
		m_zero_count = 0;
		count_zeros = true;
	}
	else{
		
		// If we're not counting zeros, then we're truncating at the n^th zero we encounter
		truncate_at_zero = m_zero_count--;
	}
	#else
	bool count_zeros = true;
	#endif // ENUMERATE_PATH
	
	// The first match state in an alignment must be an match or mismatch (no gaps allowed)
	static unsigned char first_match = query_target;
	
	unsigned char *match_ptr = &first_match;
	
	while(true){

		bool valid_alignment = true;
		
		unsigned char local_match;
		
		#ifdef ENUMERATE_PATH
		if( PATH_SPLIT(*match_ptr) == true ){
			
			// A split occurs at this element. Have we already found this split?
			deque<trace_branch>::iterator trace_iter = 
				find(m_trace_stack.begin(), m_trace_stack.end(), *match_ptr);

			if( trace_iter == m_trace_stack.end() ){

				// This is a new split
				m_trace_stack.push_back(*match_ptr);
				local_match = m_trace_stack.back().trace();				
			}
			else{

				// This is an existing split
				local_match = trace_iter->trace();
			}
		}
		else{

			local_match = *match_ptr;
		}
		#else
		local_match = *match_ptr;
		#endif // ENUMERATE_PATH
		
		switch(local_match){
			
			case query_target:
			
				if( (last_i > query_len) || (last_j < 1) ){
					valid_alignment = false;
				}
				else{
					if(m_cell_ptr->M < 0){
						valid_alignment = false;
					}
					else{
						if(m_cell_ptr->M == 0){

							if(count_zeros){
								m_zero_count++;
							}
							else{
								truncate_at_zero --;

								if(truncate_at_zero == 0){
									valid_alignment = false;
								}
							}
						}
					}

					#ifdef __DEBUG
					if(last_i > query_len){
						throw "last_i > query_len";
					}

					if(last_j < 1){
						throw "last_j < 1";
					}
					#endif // __DEBUG

					m_local_align.query_align.push_back(m_q[query_len - last_i]);

					// Reverse the target
					m_local_align.target_align.push_back(m_t[last_j - 1]);

					m_local_align.last_match = make_pair(query_len - last_i, last_j - 1);

					match_ptr = &(m_cell_ptr->M_trace);
					
					last_i --;
					last_j --;
				}
				
				break;
			case gap_target:
				
				if(last_j < 1){
					valid_alignment = false;
				}
				else{

					if(m_cell_ptr->I_query < 0){
						valid_alignment = false;
					}

					// gap the query
					m_local_align.query_align.push_back(GAP);

					#ifdef __DEBUG
					if(last_j < 1){
						throw "last_j < 1";
					}
					#endif // __DEBUG

					// Reverse the target
					m_local_align.target_align.push_back(m_t[last_j - 1]);

					m_local_align.last_match = make_pair(query_len - last_i + 1, last_j - 1);
					
					match_ptr = &(m_cell_ptr->I_query_trace);
					
					last_j --;
				}
				break;
			case query_gap:
			
				if(last_i > query_len){
					valid_alignment = false;
				}
				else{
					if(m_cell_ptr->I_target < 0){
						valid_alignment = false;
					}

					#ifdef __DEBUG
					if(last_i > query_len){
						throw "last_i > query_len";
					}
					#endif // __DEBUG

					// gap the target
					m_local_align.query_align.push_back(m_q[query_len - last_i]);

					// Reverse the target
					m_local_align.target_align.push_back(GAP);
					
					m_local_align.last_match = make_pair(query_len - last_i, last_j);
					
					match_ptr = &(m_cell_ptr->I_target_trace);
					
					last_i --;
				}
				break;
			default:
				throw "invalid_match in trace back";
		};
		
		if(!valid_alignment){
			break;
		}
		
		m_cell_ptr = m_dp_matrix + last_i*(MAX_SEQUENCE_LENGTH + 1) + last_j;
	}
	
	#ifndef ENUMERATE_PATH
	// If we're not enumerating paths, always return a "0" for the m_zero_count
	m_zero_count = 0;
	#endif // ENUMERATE_PATH
}

bool NucCruc::evaluate_alignment(alignment &local_align, const mode &m_mode)
{
	int terminal_base_pair = __;
	int last_last_base_pair = __;
	int last_base_pair = __;
	int cur_base_pair = __;

	if(m_mode != HAIRPIN){
	
		local_align.dH = param_init_H;

		// param_symmetry_S is only added for homo-dimers
		local_align.dS = param_init_S + ( (m_mode == HOMO_DIMER) ? param_symmetry_S : 0.0f);
	}
	
	//#define TRACK_ENERGY
	//#define UNAFOLD_COMPATIBILITY
	
	#ifdef TRACK_ENERGY
	cout << "Init dG(" 
		<< local_align.dH - target_T*local_align.dS
		<< ") = dH(" 
		<< local_align.dH
		<< ") + dS(" 
		<< local_align.dS
		<< ")" << endl;
	#endif // TRACK_ENERGY
	
	// Track the number of gaps and mismatches (for internal loops and bulges)
	unsigned int num_query_gap = 0;
	unsigned int num_target_gap = 0;
	unsigned int num_mismatch = 0;

	unsigned int num_base = 0;

	// Have we computed the 5' terminal exact match bases?
	bool terminal_5 = false;

	// The target sequence string has been reversed:
	// 5'- query  -3'
	// 3'- target -5'
	deque<nucleic_acid>::const_iterator q_iter = local_align.query_align.begin();
	deque<nucleic_acid>::const_iterator t_iter = local_align.target_align.begin();

	cur_base_pair = best_base_pair(*q_iter, *t_iter);

	if(watson_and_crick[cur_base_pair] == true){

		terminal_5 = true;

		if( (cur_base_pair == AT) || (cur_base_pair == TA) ){

			local_align.dH += param_AT_closing_H;
			local_align.dS += param_AT_closing_S;			
		}
	}

	// Count the number of bases in the alignment (for the salt correction)
	num_base += IS_VIRTUAL_BASE(*q_iter) ? 0 : 1;
	num_base += IS_VIRTUAL_BASE(*t_iter) ? 0 : 1;

	++q_iter;
	++t_iter;

	// The alignment index is used to test for dangling ends at the begining and end of an
	// alignment
	unsigned int alignment_index = 1;
	
	const unsigned int align_size = local_align.query_align.size();
	
	for(;q_iter != local_align.query_align.end();++q_iter, ++t_iter, ++alignment_index){
		
		#ifdef TRACK_ENERGY
		float last_dH = local_align.dH;
		float last_dS = local_align.dS;
		#endif // TRACK_ENERGY
	
		last_last_base_pair = last_base_pair;
		last_base_pair = cur_base_pair;
		cur_base_pair = best_base_pair(*q_iter, *t_iter);
		
		#ifdef TRACK_ENERGY			
		cout << "--------------------------------------------------" << endl;
		cout << "last last bp = " << bp[last_last_base_pair] << endl;
		cout << "last bp = " << bp[last_base_pair] << endl;
		cout << "curr bp = " << bp[cur_base_pair] << endl;		
		#endif // TRACK_ENERGY
		
		// Are we at the beging or end of an alignment?
		const bool align_start = (alignment_index == 1);
		const bool align_stop = (alignment_index == (align_size - 1) );
		
		// Are we inside a loop or a bulge? If so, the current base pair does *not* get
		// added to the total number of bases for the salt correction.
		const bool in_loop_or_bulge = (*q_iter == GAP) || (*t_iter == GAP) || 
			(!watson_and_crick[last_base_pair] && !watson_and_crick[cur_base_pair]);

		#ifdef TRACK_ENERGY
		if(in_loop_or_bulge){
			cout << "\tin_loop_or_bulge == true" << endl;
		}
		#endif // TRACK_ENERGY
				
		if(in_loop_or_bulge == false){
			
			// The contribution of a frayed end is the sum of the two dangling end configurations
			if(align_start && !watson_and_crick[last_base_pair] && NON_VIRTUAL_BASE_PAIR(last_base_pair) ){
				
				// Frayed-end at the begining of the alignment
				const nucleic_acid tmp_q = nucleic_acid( QUERY_BASE(last_base_pair) );
				const nucleic_acid tmp_t = nucleic_acid( TARGET_BASE(last_base_pair) );
				
				int tmp_pair = best_base_pair(tmp_q, E);
				
				local_align.dH += param_H[tmp_pair][cur_base_pair];
				local_align.dS += param_S[tmp_pair][cur_base_pair];
				
				tmp_pair = best_base_pair(E, tmp_t);
				
				local_align.dH += param_H[tmp_pair][cur_base_pair];
				local_align.dS += param_S[tmp_pair][cur_base_pair];
				
				#ifdef TRACK_ENERGY
				cout << "\tfrayed start" << endl;
				#endif // TRACK_ENERGY
			}
			else{
				if(align_stop && !watson_and_crick[cur_base_pair] && NON_VIRTUAL_BASE_PAIR(cur_base_pair) ){
				
					// Frayed-end at the end of the alignment
					int tmp_pair = best_base_pair(*q_iter, E);
				
					local_align.dH += param_H[last_base_pair][tmp_pair];
					local_align.dS += param_S[last_base_pair][tmp_pair];

					tmp_pair = best_base_pair(E, *t_iter);

					local_align.dH += param_H[last_base_pair][tmp_pair];
					local_align.dS += param_S[last_base_pair][tmp_pair];
					
					#ifdef TRACK_ENERGY
					cout << "\tfrayed stop" << endl;
					#endif // TRACK_ENERGY
				}
				else{
					local_align.dH += param_H[last_base_pair][cur_base_pair];
					local_align.dS += param_S[last_base_pair][cur_base_pair];
					
					#ifdef TRACK_ENERGY
					cout << "\tnormal pair = " 
						<< param_H[last_base_pair][cur_base_pair] - target_T*param_S[last_base_pair][cur_base_pair]
						<< endl;
					#endif // TRACK_ENERGY
				}
			}
			
			// Count the number of bases in the alignment (for the salt correction).
			// Bases that are internal to a loop do not get a salt correction (but the
			// terminal bases do).
			// From "The thermodynamics of DNA structural motifs"
			// J. SantaLucia and D. Hicks:
			// The terminal mismatches in internal loops are
			// assumed to have the same salt dependence as base pairs 
			// (Equation 3), whereas the stability of the remainder of 
			// the internal loop nucleotides are assumed to be salt
			// independent.
			num_base += IS_VIRTUAL_BASE(*q_iter) ? 0 : 1;
			num_base += IS_VIRTUAL_BASE(*t_iter) ? 0 : 1;			
		}

		// If these bases are an exact watson and crick match,
		// pay any free energy penalties that arise from 
		// loop or bulge closing
		//if( watson_and_crick[cur_base_pair] ){
		if( watson_and_crick[cur_base_pair] || (cur_base_pair == EE) ){ // Modified on March 14, 2019

			// Update the terminal_base_pair (since we don't know the location
			// of the last watson and crick match
			terminal_base_pair = cur_base_pair;

			// Do we still need to compute the terminal_5 initiation?
			if(!terminal_5){

				terminal_5 = true;

				if( (cur_base_pair == AT) || (cur_base_pair == TA) ){

					local_align.dH += param_AT_closing_H;
					local_align.dS += param_AT_closing_S;
					
					#ifdef TRACK_ENERGY
					cout << "\tterminal 5" << endl;
					#endif // TRACK_ENERGY
				}
			}

			const unsigned int max_gap = max(num_query_gap, num_target_gap);
			
			#ifdef TRACK_ENERGY
			cout << "\tmax_gap = " << int(max_gap) << endl;
			cout << "\tnum_mismatch = " << int(num_mismatch) << endl;			
			#endif // TRACK_ENERGY
					
			// Are we closing an internal loop?
			if( (num_mismatch > 1) || ( (max_gap > 0) && (num_mismatch == 1) ) ){
				
				// How many nucleotides are in this internal loop?
				const unsigned int gap_difference = (num_query_gap > num_target_gap) ?  
					num_query_gap - num_target_gap : 
					num_target_gap - num_query_gap;
				
				const unsigned int loop_size = num_mismatch*2 + gap_difference;

				// There are experimental parameters for the following
				// base pair combinations:
				// 	gg/tt, tt/gg, gt/tg and tg/gt
				if( (loop_size == 2) &&
				    ( (last_base_pair == GT) || (last_base_pair == TG) ) &&
				    ( (last_last_base_pair == GT) || (last_last_base_pair == TG) ) ){

					local_align.dH += param_H[last_last_base_pair][last_base_pair];
					local_align.dS += param_S[last_last_base_pair][last_base_pair];

					// Increment the number of bases to account for the terminal TG or GT pair
					// that were
					num_base += 2;
					
					#ifdef TRACK_ENERGY
					cout << "\tinternal loop close 2" << endl;
					#endif // TRACK_ENERGY
				}
				else{
					
					#ifdef TRACK_ENERGY
					float loop_dH = 0.0f;
					float loop_dS = 0.0f;
					#endif // TRACK_ENERGY
					
					// The internal loop closing penalty is entropic
					local_align.dS += param_loop_S[loop_size];

					// If this is an asymmetric internal loop, apply the asymmetry penalty
					local_align.dS += gap_difference * param_asymmetric_loop_dS;
					
					#ifdef TRACK_ENERGY
					cout << "\tloop dG = " <<  -target_T*param_loop_S[loop_size] << endl;
					cout << "\tloop penalty dG = " <<  -target_T*(gap_difference * 
						param_asymmetric_loop_dS) << endl;
					cout << "\tloop_size = " << int(loop_size) << endl;
					cout << "\tgap_difference = " << int(gap_difference) << endl;
					
					loop_dS += param_loop_S[loop_size] + gap_difference * param_asymmetric_loop_dS;
					#endif // TRACK_ENERGY
					
					///////////////////////////////////////////////////////////////////////
					// The loop score includes both *left* and *right* terminal mismatches 
					// If the loop starts or stops with a GAP, then we will need to add the
					// terminal mismatch ourselves.
					///////////////////////////////////////////////////////////////////////

					deque<nucleic_acid>::const_iterator rhs_q_iter = q_iter;
					deque<nucleic_acid>::const_iterator rhs_t_iter = t_iter;

					rhs_q_iter --;
					rhs_t_iter --;

					///////////////////////////////////////////////////////////////////////////
					// There are separate parameters for terminal mismatches that are internal
					// to a loop
					
					//***********************************************************************
					// Apply the correction here by first removing the previous contribution
					// from the right hand terminal base that we added above
					local_align.dH -= param_H[last_base_pair][cur_base_pair];
					local_align.dS -= param_S[last_base_pair][cur_base_pair];
					//***********************************************************************
					
					//local_align.dH += param_loop_terminal_H[last_base_pair][cur_base_pair] - 
					//	param_H[last_base_pair][cur_base_pair];
					//local_align.dS += param_loop_terminal_S[last_base_pair][cur_base_pair] - 
					//	param_S[last_base_pair][cur_base_pair];
					
					#ifdef TRACK_ENERGY
					cout << "\tremoving dg = " << param_H[last_base_pair][cur_base_pair] - 
						target_T*param_S[last_base_pair][cur_base_pair]
						<< endl;
					cout << "\tadding dg = " << param_loop_terminal_H[last_base_pair][cur_base_pair] - 
						target_T*param_loop_terminal_S[last_base_pair][cur_base_pair]
						<< endl;
					
					cout << "\tparam_loop_terminal_H = " << param_loop_terminal_H[last_base_pair][cur_base_pair] << endl;
					cout << "\tparam_loop_terminal_S = " << param_loop_terminal_S[last_base_pair][cur_base_pair] << endl;
					
					cout << "\tparam_H = " << param_H[last_base_pair][cur_base_pair] << endl;
					cout << "\tparam_S = " << param_S[last_base_pair][cur_base_pair] << endl;
					
					loop_dH += param_loop_terminal_H[last_base_pair][cur_base_pair] - 
						param_H[last_base_pair][cur_base_pair];
					loop_dS += param_loop_terminal_S[last_base_pair][cur_base_pair] - 
						param_S[last_base_pair][cur_base_pair];
					#endif // TRACK_ENERGY
						
					// The right terminal mistmatch. 
					if(HAS_GAP(last_base_pair) == false){
						
						// We've already subtracted off the contribution that was added before
						// we know that this was a terminal pair (see comment surrounded by **** above)
						local_align.dH += param_loop_terminal_H[last_base_pair][cur_base_pair];
						local_align.dS += param_loop_terminal_S[last_base_pair][cur_base_pair];
						
						#ifdef TRACK_ENERGY
						cout << "\tHAS_GAP(last_base_pair) == false"<< endl;
						cout << "\tAdded dH = " << param_loop_terminal_H[last_base_pair][cur_base_pair] << endl;
						cout << "\tAdded dS = " << param_loop_terminal_S[last_base_pair][cur_base_pair] << endl;
						#endif // TRACK_ENERGY
					}
					else{ // There is a gap in the right hand terminal pair

						int mm_base_pair = __;

						// The right terminal mismatch
						if(QUERY_BASE(last_base_pair) == GAP){

							// Extract the query base to use by walking backwards on the
							// query strand until we find a non-gap base
							while(true){

								if( !IS_VIRTUAL_BASE(*rhs_q_iter) ){

									mm_base_pair = best_base_pair( *rhs_q_iter, nucleic_acid( TARGET_BASE(last_base_pair) ) );
									break;
								}

								if(rhs_q_iter == local_align.query_align.begin() ){
									break;
								}

								rhs_q_iter --;
							}

						}
						else{ // TARGET_BASE(last_base_pair) == GAP

							// Extract the target base to use by walking backwards on the
							// target strand until we find a non-gap base
							while(true){

								if( !IS_VIRTUAL_BASE(*rhs_t_iter) ){

									mm_base_pair = best_base_pair(nucleic_acid( QUERY_BASE(last_base_pair) ), *rhs_t_iter);
									break;
								}

								if(rhs_t_iter == local_align.target_align.begin() ){
									break;
								}

								rhs_t_iter --;
							}
						}

						local_align.dH += param_loop_terminal_H[mm_base_pair][cur_base_pair];
						local_align.dS += param_loop_terminal_S[mm_base_pair][cur_base_pair];
						
						#ifdef TRACK_ENERGY
						cout << "\tright terminal mismatch = " 
							<< param_H[mm_base_pair][cur_base_pair] - target_T*param_S[mm_base_pair][cur_base_pair]
							<< endl;
						loop_dH += param_loop_terminal_H[mm_base_pair][cur_base_pair];
						loop_dS += param_loop_terminal_S[mm_base_pair][cur_base_pair];
						#endif // TRACK_ENERGY
					
					}

					// The left terminal mismatch. I'm not sure the optimal way to identify
					// these bases, other that to walk backwards until we encounter a pair
					// bases that are exactly complementary *and* are flanked (on the right)
					// by one or more gaps (i.e. a gap in the query and/or target strand). If
					// there are no flanking gaps, then this terminal mismatch has *already*
					// been accounted for.
					deque<nucleic_acid>::const_iterator lhs_q_iter = q_iter;
					deque<nucleic_acid>::const_iterator lhs_t_iter = t_iter;

					lhs_q_iter --;
					lhs_t_iter --;

					while(true){

						const int pm_base_pair = best_base_pair(*lhs_q_iter, *lhs_t_iter);

						// If we've found a perfect match, stop back tracking and
						// start reading ahead (into the loop)
						if(watson_and_crick[pm_base_pair] == true){

							lhs_q_iter ++;
							lhs_t_iter ++;

							// If the terminal bases do not contain gaps, then
							// we've already accounted for this terminal mismatch above
							// (including the salt correction)
							if( (*lhs_q_iter != GAP) && (*lhs_t_iter != GAP) ){
								
								// Since we previously added a contribution for this base, remove that
								// contribution now
								const int mm_base_pair = best_base_pair(*lhs_q_iter, *lhs_t_iter);
								
								local_align.dH -= param_H[pm_base_pair][mm_base_pair];
								local_align.dS -= param_S[pm_base_pair][mm_base_pair];
								
								#ifdef TRACK_ENERGY
								loop_dH -= param_H[pm_base_pair][mm_base_pair];
								loop_dS -= param_S[pm_base_pair][mm_base_pair];
								
								cout << "\tRemoved left terminal mismatch dG = "
									<< param_H[pm_base_pair][mm_base_pair] + 
									target_T*param_S[pm_base_pair][mm_base_pair] << endl;
								#endif // TRACK_ENERGY
							}
							else{
								// Since one the left hand side bases was a gap,
								// these bases were *not* counted in the salt correction.
								num_base += 2;
							
								while(*lhs_q_iter == GAP){
									lhs_q_iter ++;
								}

								while(*lhs_t_iter == GAP){
									lhs_t_iter ++;
								}
							}
							
							const int mm_base_pair = best_base_pair(*lhs_q_iter, *lhs_t_iter);

							local_align.dH += param_loop_terminal_H[pm_base_pair][mm_base_pair];
							local_align.dS += param_loop_terminal_S[pm_base_pair][mm_base_pair];
							
							#ifdef TRACK_ENERGY
							cout << "\tleft terminal mismatch = " 
								<< param_loop_terminal_H[pm_base_pair][mm_base_pair] - 
									target_T*param_loop_terminal_S[pm_base_pair][mm_base_pair]
								<< endl;
							loop_dH += param_loop_terminal_H[pm_base_pair][mm_base_pair];
							loop_dS += param_loop_terminal_S[pm_base_pair][mm_base_pair];
							#endif // TRACK_ENERGY

							break;
						}

						// If we hit the begining of the alignment, then we can't compute
						// the left terminal mismatch
						if( lhs_q_iter == local_align.query_align.begin() ){
							break;
						}

						lhs_q_iter --;
						lhs_t_iter --;
					}

					// The two left hand side bases have already been accounted for
					// in the num_base salt correction. Check the status of the right 
					// hand side bases
					if(rhs_q_iter != lhs_q_iter){
						num_base ++;
					}

					if(rhs_t_iter != lhs_t_iter){
						num_base ++;
					}
					
					#ifdef TRACK_ENERGY
					cout << "\tinternal loop close" << endl;
					cout << "\t\tcost dG(" 
						<< loop_dH - target_T*loop_dS
						<< ") = dH(" 
						<< loop_dH
						<< ") + dS(" 
						<< loop_dS
						<< ")" << endl;
					#endif // TRACK_ENERGY
				}

			}
			else{

				// Is this a bulge?
				if(num_query_gap || num_target_gap){
					
					#ifdef __DEBUG
					if( (num_query_gap > 0) && (num_target_gap > 0) ){
						
						cout << endl;
						cout << "q = " << query_seq() << endl;
						cout << "t = " << target_seq() << endl;
						cout << print_alignment(local_align.query_align, local_align.target_align);
						
						throw __FILE__ ":NucCruc::evaluate_alignment: (num_query_gap > 0) && (num_target_gap > 0)";
					}
					#endif // __DEBUG

					const unsigned int bulge_size = (num_query_gap > num_target_gap) ?  
						num_query_gap : num_target_gap;
					
					// Bulges of size 1 must include the NN terms from the flanking bases
					if(bulge_size == 1){

						local_align.dH += param_H[last_last_base_pair][cur_base_pair];
						local_align.dS += param_S[last_last_base_pair][cur_base_pair];
						
						#ifdef TRACK_ENERGY
						cout << "\tbulge 1" << endl;
						cout << "\t\tflanking dg = " 
							<< param_H[last_last_base_pair][cur_base_pair] - 
							   target_T*param_S[last_last_base_pair][cur_base_pair] << endl;
						#endif // TRACK_ENERGY
					}

					// The bulge loop closing penalty is entropic
					local_align.dS += param_bulge_S[bulge_size];
					
					#ifdef TRACK_ENERGY
						cout << "\tbulge loop closing penalty = " 
							<< -target_T*param_bulge_S[bulge_size] << endl;
					#endif // TRACK_ENERGY
						
					// The DNA folding program UNAFold does not apply the AT_closing
					// penalty to single base bulges (personal commiunication from Nick Markham,
					// one of the authors of UNAFold)
					#ifdef UNAFOLD_COMPATIBILITY
					// If the bulge is terminated by an A/T pair, apply the A/T closing penalty
					if( (bulge_size != 1) && ( (*q_iter == A) || (*q_iter == T) ) ){
						local_align.dS += param_bulge_AT_closing_S;						
					}
					#else
					// If the bulge is terminated by an A/T pair, apply the A/T closing penalty
					if( (*q_iter == A) || (*q_iter == T) ){
						local_align.dS += param_bulge_AT_closing_S;						
					}
					#endif // UNAFOLD_COMPATIBILITY
					
					// The DNA folding program UNAFold does not apply the AT_closing
					// penalty to single base bulges (personal commiunication from Nick Markham,
					// one of the authors of UNAFold)
					#ifdef UNAFOLD_COMPATIBILITY
					// Search backwards in the alignment to test for A/T initiation, which
					// will also receive the param_bulge_AT_closing_S penalty.
					if( (bulge_size != 1) && 
						has_AT_initiation( q_iter, local_align.query_align.begin(), t_iter, 
						local_align.target_align.begin() ) ){

						local_align.dS += param_bulge_AT_closing_S;						
					}
					#else
					// Search backwards in the alignment to test for A/T initiation, which
					// will also receive the param_bulge_AT_closing_S penalty.
					if( has_AT_initiation( q_iter, local_align.query_align.begin(), t_iter, 
						local_align.target_align.begin() ) ){

						local_align.dS += param_bulge_AT_closing_S;						
					}
					#endif // UNAFOLD_COMPATIBILITY
					
					#ifdef TRACK_ENERGY
					cout << "\tbulge" << endl;
					#endif // TRACK_ENERGY
				}
			}

			// Reset the gap and mismatch counters
			num_query_gap = 0;
			num_target_gap = 0;

			num_mismatch = 0;
		}
		else{
			// If either base is a GAP or $, then this is *not* a mismatch
			num_mismatch += ( !IS_VIRTUAL_BASE(*q_iter) && !IS_VIRTUAL_BASE(*t_iter) ) ? 1 : 0;
		}

		num_query_gap += (*q_iter == GAP) ? 1 : 0;
		num_target_gap += (*t_iter == GAP) ? 1 : 0;
		
		#ifdef TRACK_ENERGY
		
		last_dH = local_align.dH - last_dH;
		last_dS = local_align.dS - last_dS;
		
		cout << "change dG(" 
			<< last_dH - target_T*last_dS
			<< ") = dH(" 
			<< last_dH
			<< ") + dS(" 
			<< last_dS
			<< ")" << endl;
			
		cout << "partial dG(" 
			<< local_align.dH - target_T*local_align.dS
			<< ") = dH(" 
			<< local_align.dH
			<< ") + dS(" 
			<< local_align.dS
			<< ")" << endl;
		#endif // TRACK_ENERGY
	}

	// Account for a possible 3' AT base pair
	if( (terminal_base_pair == AT) || (terminal_base_pair == TA) ){

		local_align.dH += param_AT_closing_H;
		local_align.dS += param_AT_closing_S;
		
		#ifdef TRACK_ENERGY
			cout << "terminal AT dG(" 
				<< param_AT_closing_H - target_T*param_AT_closing_S
				<< ") = dH(" 
				<< param_AT_closing_H
				<< ") + dS(" 
				<< param_AT_closing_S
				<< ")" << endl;
		#endif // TRACK_ENERGY
	}

	#ifdef TRACK_ENERGY
	cout << "Final dG(" 
		<< local_align.dH - target_T*local_align.dS
		<< ") = dH(" 
		<< local_align.dH
		<< ") + dS(" 
		<< local_align.dS
		<< ")" << endl;
	#endif // TRACK_ENERGY
	
	// Binding must be enthalpically driven
	if(local_align.dH >= 0.0f){
		//return 0.0f;
		return false;
	}

	// Note that the nucleation contributions has already been included
	// in dH_d and dS_d. **IMPORTANT** It is now up the calling function to 
	// properly set the strand concentration!
	// For A + B -> D,
	// Ct = C_excess - C_limit/2
	// If C_excess == C_limit, Ct = C_excess/2
	// If C_excess >> C_limit, Ct = C_excess
	//const float heterodimer_inv_alpha = (m_mode == HETERO_DIMER) ? 1.0f/4.0f : 1.0f;
	const float heterodimer_inv_alpha = 1.0f;
	//const float heterodimer_inv_alpha = 1.0f/19.0f; // For the Citrus MOL-PCR assays, use the Tm @ 95% bound (instead of 50%)

	// Apply the salt correction
	local_align.dS += param_SALT * (0.5f*num_base - 1) * log(na_concentration);

	#ifdef TRACK_ENERGY
	cout << "Salt correction dS term (" << int(num_base) << " bp) = " 
		<< param_SALT * (0.5f*num_base - 1) * log(na_concentration) << endl;
	cout << "Salt correction dG(" 
		<< local_align.dH - target_T*local_align.dS
		<< ") = dH(" 
		<< local_align.dH
		<< ") + dS(" 
		<< local_align.dS
		<< ")" << endl;
	#endif // TRACK_ENERGY
	
	float __tm;
	
	if(m_mode == HAIRPIN){
		// Note that the hairpin calculation does not use a strand concentration
		__tm = local_align.dH/local_align.dS - NC_ZERO_C;
	}
	else{
		__tm = local_align.dH/(NC_R*log(strand_concentration*heterodimer_inv_alpha) + 
			local_align.dS) - NC_ZERO_C;
	}
	
	// Return physical temperatures (i.e. no negative Tm's)
	local_align.tm = max(0.0f, __tm);

	return true;
}

bool NucCruc::evaluate_hairpin_alignment(alignment &local_align)
{
	
	// We need to save the locations of the bases that flank the hairpin loop
	// to compute base-dependent loop energy terms
	const int last_3 = local_align.first_match.first;
	const int last_5 = local_align.first_match.second;
	
	const unsigned int hairpin_loop_len = last_3 - last_5 - 1; // query_len - (last_i + last_j);
	
	// DEBUG
	//cerr << "------------------------------------------------------" << endl;
	//cerr << "hairpin_loop_len = " << hairpin_loop_len << endl;
	
	// Hairpins do not pay the initiation cost
	local_align.dH = 0.0f;
	local_align.dS = 0.0f;

	// Add the entropic contribution from the hairpin loop here
	local_align.dS += param_hairpin_S[hairpin_loop_len];
	
	// DEBUG
	//cerr << "param_hairpin_S[hairpin_loop_len] = " << param_hairpin_S[hairpin_loop_len] << endl;
	
	//////////////////////////////////////////////////////////////////
	// Hairpin loop -- length dependent terms
	//////////////////////////////////////////////////////////////////
	
	// The closing bases for the hairpin loop are requried for loops of all sizes
	const int last_base_pair = best_base_pair(query[last_5], query[last_3]);
	int cur_base_pair = __;
	
	switch(hairpin_loop_len){
		case 3:
			{
				// Add tri-loop bonus
				const int loop_index = find_loop_index(query, last_5, 5 /* search for a 5 base match */);
				
				// DEBUG
				//cerr << "loop_index = " << loop_index << endl;
	
				if(loop_index >= 0){
					
					local_align.dH += param_hairpin_special_H[loop_index];
					local_align.dS += param_hairpin_special_S[loop_index];
				}
				
				// Add closing AT penalty
				if( (last_base_pair == AT) || (last_base_pair == TA) ){
					
					// DEBUG
					//cerr << "param_bulge_AT_closing_S = " << param_bulge_AT_closing_S << endl;
	
					local_align.dS += param_bulge_AT_closing_S;
				}
			}
			break;
		case 4:
			{
				// Add tetra-loop bonus
				const int loop_index = find_loop_index(query, last_5, 6 /* search for a 6 base match */);
								
				if(loop_index >= 0){
					
					local_align.dH += param_hairpin_special_H[loop_index];
					local_align.dS += param_hairpin_special_S[loop_index];					
				}
			}
			
			// Fall through to compute the terminal mismatch
		default:
			// Add terminal mismatch
			// last_base_pair is computed before the switch statement
			cur_base_pair = best_base_pair(query[last_5 + 1], query[last_3 - 1]);
			
			//local_align.dH += param_H[last_base_pair][cur_base_pair];
			//local_align.dS += param_S[last_base_pair][cur_base_pair];
			
			// Use the unpublished, terminal hairpin parameters (instead of the
			// "normal" base stacking parameters that are commented out above).
			local_align.dH += param_hairpin_terminal_H[last_base_pair][cur_base_pair];
			local_align.dS += param_hairpin_terminal_S[last_base_pair][cur_base_pair];
			
			// DEBUG
			//cerr << "param_loop_terminal_H[" << bp[last_base_pair] << "][" << bp[cur_base_pair] 
			//	<< "] = " << param_loop_terminal_H[last_base_pair][cur_base_pair] << endl;
			//cerr << "param_loop_terminal_S[" << bp[last_base_pair] << "][" << bp[cur_base_pair] 
			//	<< "] = " << param_loop_terminal_S[last_base_pair][cur_base_pair] << endl;
			
			break;
	};
		
	return evaluate_alignment(local_align, HAIRPIN);
}

// Given a query and target sequence, compute the melting temperature
float NucCruc::approximate_tm_heterodimer()
{
	if(use_dinkelbach){

		// Make a copy of the initial target_T (as we will need to 
		// overwrite it
		const float init_T = target_T;
		float q = -999999.9f;
		float last_q = q;
		float local_tm = 0.0f;
		NC_Score max_score = 0;
		
		temperature(NC_ZERO_C);
		
		unsigned int iter = 0;
				
		do{
			
			// Reset the alignment parameters
			curr_align.clear();

			// Compute the alignment 
			max_score = align_heterodimer();
			
			local_tm = tm_dimer(query, target, HETERO_DIMER);
			
			last_q = q;
			q = delta_G();
			
			// Set target_T = local_tm to compute -(dh - local_tm*ds) on the next iteration
			temperature(NC_ZERO_C + local_tm);
						
			iter ++;
		}
		while( (q < 0.0) && (q > last_q) );
		
		// Restore the initial temperature (and recompute the DP parameters)
		temperature(init_T);
		
		// Save the delta G as computed by dynamic programming alone
		curr_align.dp_dg = -NC_INV_SCORE_SCALE(max_score);
		
		return local_tm;
	}
	else{
		// Reset the alignment parameters
		curr_align.clear();

		// Compute the alignment 
		const NC_Score max_score = align_heterodimer();
		
		const float local_tm = tm_dimer(query, target, HETERO_DIMER);
		
		// Save the delta G as computed by dynamic programming alone
		curr_align.dp_dg = -NC_INV_SCORE_SCALE(max_score);
		
		return local_tm;
	}
};

float NucCruc::approximate_tm_homodimer()
{
	if(use_dinkelbach){

		// Make a copy of the initial target_T (as we will need to 
		// overwrite it
		const float init_T = target_T;
		float q = -999999.9f;
		float last_q = q;
		float local_tm = 0.0f;
		NC_Score max_score = 0;
		
		temperature(NC_ZERO_C);
		
		unsigned int iter = 0;
				
		do{
			
			// Reset the alignment parameters
			curr_align.clear();

			// Compute the alignment 
			max_score = align_homodimer();
			
			local_tm = tm_dimer(query, query, HOMO_DIMER);
			
			last_q = q;
			q = delta_G();
			
			// Set target_T = local_tm to compute -(dh - local_tm*ds) on the next iteration
			temperature(NC_ZERO_C + local_tm);
						
			iter ++;
		}
		while( (q < 0.0) && (q > last_q) );
				
		// Restore the initial temperature (and recompute the DP parameters)
		temperature(init_T);
		
		// Save the delta G as computed by dynamic programming alone
		curr_align.dp_dg = -NC_INV_SCORE_SCALE(max_score);
		
		return local_tm;
	}
	else{
		// Reset the alignment parameters
		curr_align.clear();

		// Compute the alignment 
		const NC_Score max_score = align_homodimer();
		
		const float local_tm = tm_dimer(query, query, HOMO_DIMER);
		
		// Save the delta G as computed by dynamic programming alone
		curr_align.dp_dg = -NC_INV_SCORE_SCALE(max_score);
		
		return local_tm;
	}
}

float NucCruc::tm_dimer(const CircleBuffer<nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q,
			const CircleBuffer<nucleic_acid, MAX_SEQUENCE_LENGTH> &m_t,
			const mode &m_mode)
{
	// Assume that the alignment(s) have already been computed (but have not been enumerated)
	
	if(na_concentration <= 0.0f){
		throw ":NucCruc::tm_dimer: Invalid na_concentration";
	}
	
	if(strand_concentration <= 0.0f){
		throw ":NucCruc::tm_dimer: Invalid strand_concentration";
	}
	
	vector<NC_Elem*>::const_iterator max_iter;
	
	for(max_iter = max_ptr.begin();max_iter != max_ptr.end();max_iter++){
		
		enumerate_dimer_alignments( dp_matrix, *max_iter, curr_align, 
			m_q, m_t, m_mode);
	}
	
	return curr_align.tm;
}

float NucCruc::approximate_tm_hairpin()
{
	if(na_concentration <= 0.0f){
		throw ":NucCruc::approximate_tm_hairpin: Invalid na_concentration";
	}

	if(use_dinkelbach){

		// Make a copy of the initial target_T (as we will need to 
		// overwrite it
		const float init_T = target_T;
		float q = -999999.9f;
		float last_q = q;
		float local_tm = 0.0f;
		NC_Score max_score = 0;
		
		temperature(NC_ZERO_C);
		
		unsigned int iter = 0;
				
		do{
			
			// Reset the alignment parameters
			curr_align.clear();

			// Compute the hairpin alignment
			max_score = align_hairpin(query);
			
			vector<NC_Elem*>::const_iterator max_iter;

			for(max_iter = max_ptr.begin();max_iter != max_ptr.end();max_iter++){

				enumerate_hairpin_alignments(dp_matrix, *max_iter, curr_align, query);
			}
		
			local_tm = curr_align.tm;
			
			last_q = q;
			q = delta_G();
			
			// Set target_T = local_tm to compute -(dh - local_tm*ds) on the next iteration
			temperature(NC_ZERO_C + local_tm);
						
			iter ++;
		}
		while( (q < 0.0) && (q > last_q) );
				
		// Restore the initial temperature (and recompute the DP parameters)
		temperature(init_T);
		
		// Save the delta G as computed by dynamic programming alone
		curr_align.dp_dg = -NC_INV_SCORE_SCALE(max_score);
		
		return local_tm;
	}
	else{
	
		// Reset the alignment parameters
		curr_align.clear();

		// Compute the hairpin alignment
		const NC_Score max_score = align_hairpin(query);

		vector<NC_Elem*>::const_iterator max_iter;

		for(max_iter = max_ptr.begin();max_iter != max_ptr.end();max_iter++){
			enumerate_hairpin_alignments(dp_matrix, *max_iter, curr_align, query);
		}
		
		// Save the delta G as computed by dynamic programming alone
		curr_align.dp_dg = -NC_INV_SCORE_SCALE(max_score);
		
		return curr_align.tm;
	}
}

// m_len must be either 5 or 6
int NucCruc::find_loop_index(const CircleBuffer<nucleic_acid, MAX_SEQUENCE_LENGTH> &m_q,
	const unsigned int &m_start, const unsigned int &m_len)
{
	// This lookup table was produced using:
	// cat nuc_cruc_santa_lucia_data.txt | awk '{print "\""$1"\","}'
	static const char *hairpin_special_loop[] = {
		"AAAAAT",
		"AAAACT",
		"AAACAT",
		"ACTTGT",
		"AGAAAT",
		"AGAAT",
		"AGAGAT",
		"AGATAT",
		"AGCAAT",
		"AGCAT",
		"AGCGAT",
		"AGCTTT",
		"AGGAAT",
		"AGGAT",
		"AGGGAT",
		"AGGGGT",
		"AGTAAT",
		"AGTAT",
		"AGTGAT",
		"AGTTCT",
		"ATTCGT",
		"ATTTGT",
		"ATTTTT",
		"CAAAAG",
		"CAAACG",
		"CAACAG",
		"CAACCG",
		"CCTTGG",
		"CGAAAG",
		"CGAAG",
		"CGAGAG",
		"CGATAG",
		"CGCAAG",
		"CGCAG",
		"CGCGAG",
		"CGCTTG",
		"CGGAAG",
		"CGGAG",
		"CGGGAG",
		"CGGGGG",
		"CGTAAG",
		"CGTAG",
		"CGTGAG",
		"CGTTCG",
		"CTTCGG",
		"CTTTGG",
		"CTTTTG",
		"GAAAAC",
		"GAAAAT",
		"GAAACC",
		"GAAACT",
		"GAACAC",
		"GAACAT",
		"GCTTGC",
		"GCTTGT",
		"GGAAAC",
		"GGAAAT",
		"GGAAC",
		"GGAGAC",
		"GGAGAT",
		"GGATAC",
		"GGATAT",
		"GGCAAC",
		"GGCAAT",
		"GGCAC",
		"GGCGAC",
		"GGCGAT",
		"GGCTTC",
		"GGCTTT",
		"GGGAAC",
		"GGGAAT",
		"GGGAC",
		"GGGGAC",
		"GGGGAT",
		"GGGGGC",
		"GGGGGT",
		"GGTAAC",
		"GGTAAT",
		"GGTAC",
		"GGTGAC",
		"GGTGAT",
		"GGTTCC",
		"GTATAT",
		"GTTCGC",
		"GTTCGT",
		"GTTTGC",
		"GTTTGT",
		"GTTTTC",
		"GTTTTT",
		"TAAAAA",
		"TAAAAG",
		"TAAACA",
		"TAAACG",
		"TAACAA",
		"TAACAG",
		"TCTTGA",
		"TCTTGG",
		"TGAAA",
		"TGAAAA",
		"TGAAAG",
		"TGAGAA",
		"TGAGAG",
		"TGATAA",
		"TGATAG",
		"TGCAA",
		"TGCAAA",
		"TGCAAG",
		"TGCGAA",
		"TGCGAG",
		"TGCTTA",
		"TGCTTG",
		"TGGAA",
		"TGGAAA",
		"TGGAAG",
		"TGGGAA",
		"TGGGAG",
		"TGGGGA",
		"TGGGGG",
		"TGTAA",
		"TGTAAA",
		"TGTAAG",
		"TGTGAA",
		"TGTGAG",
		"TGTTCA",
		"TTTCGA",
		"TTTCGG",
		"TTTTAG",
		"TTTTGA",
		"TTTTGG",
		"TTTTTA",
		"TTTTTG"
	};
	
	const char *base_name = "ACGTE";
	
	unsigned int match = 0;
	
	#define	BASE_LEN	6
	
	char base[BASE_LEN]; // <- a maximum of BASE_LEN bases will be tested! 
	
	// Initialize to avoid compiler warnings
	memset(base, 0, BASE_LEN);
	
	base[0] = base_name[ m_q[m_start] ];
	 
	for(unsigned int i = 0;i < NUM_SPECIAL_HAIRPIN_LOOP;i++){
	
		const char* seq = hairpin_special_loop[i];
		
		if(seq[0] == base[0]){
			
			if(match == 0){
			
				match = 1;
				
				base[1] = base_name[ m_q[m_start + 1] ];
			}
			
			if(seq[1] == base[1]){
			
				if(match == 1){
				
					match = 2;

					base[2] = base_name[ m_q[m_start + 2] ];
				}

				if(seq[2] == base[2]){
			
					if(match == 2){
					
						match = 3;

						base[3] = base_name[ m_q[m_start + 3] ];
					}

					if(seq[3] == base[3]){
			
						if(match == 3){
						
							match = 4;

							base[4] = base_name[ m_q[m_start + 4] ];
						}

						if(seq[4] == base[4]){
			
							if(match == 4){

								match = 5;
								
								if(m_len == 5){
									
									// Did we match a five base pattern or five
									// bases of a longer pattern?	
									if(seq[5] == '\0'){
										return i;
									}
									else{
										continue;
									}
								}
								
								base[5] = base_name[ m_q[m_start + 5] ];
							}
							
							if(seq[5] == base[5]){
							
								return i;
							}
						}
						else{
							if(match > 4){
								return -1;
							}
						}
					}
					else{
						if(match > 3){
							return -1;
						}
					}
				}
				else{
					if(match > 2){
						return -1;
					}
				}
			}
			else{
				if(match > 1){
					return -1;
				}
			}
		}
		else{
			if(match > 0){
				return -1;
			}
		}
	}
	
	return -1;
}

bool NucCruc::is_internal_to_loop( deque<nucleic_acid>::const_iterator m_q, 
	const deque<nucleic_acid>::const_iterator &m_query_end, 
	deque<nucleic_acid>::const_iterator m_t, 
	const deque<nucleic_acid>::const_iterator &m_target_end )
{
	
	// Advance the iterator *copies*
	m_q ++;
	m_t ++;
	
	if( (m_q == m_query_end) || (m_t == m_target_end) ){
		return false;
	}

	return (watson_and_crick[best_base_pair(*m_q, *m_t)] == false);
}

bool NucCruc::has_AT_initiation( deque<nucleic_acid>::const_iterator m_q, 
	const deque<nucleic_acid>::const_iterator &m_query_begin, 
	deque<nucleic_acid>::const_iterator m_t, 
	const deque<nucleic_acid>::const_iterator &m_target_begin )
{
	
	// Advance the iterator *copies*
	do{
		
		m_q --;
		m_t --;
	}
	while( (m_q != m_query_begin) && (m_t != m_target_begin) && ( (*m_q == GAP) || (*m_t == GAP) ) );
	
	int bp = best_base_pair(*m_q, *m_t);
	
	return ( (bp == AT) || (bp == TA) );
}

#ifdef __DEBUG

void NucCruc::dump_tables()
{
	const char* pairs[] = {
		"AA", "AC", "AG", "AT", "AI", "AE", "A_", 
		"CA", "CC", "CG", "CT", "CI", "CE", "C_", 
		"GA", "GC", "GG", "GT", "GI", "GE", "G_", 
		"TA", "TC", "TG", "TT", "TI", "TE", "T_", 
		"IA", "IC", "IG", "IT", "II", "IE", "I_", 
		"EA", "EC", "EG", "ET", "EI", "EE", "E_",
		"_A", "_C", "_G", "_T", "_I", "_E", "__"
	};
	
	static const char *hairpin_special_loop[] = {
		"AAAAAT",
		"AAAACT",
		"AAACAT",
		"ACTTGT",
		"AGAAAT",
		"AGAAT",
		"AGAGAT",
		"AGATAT",
		"AGCAAT",
		"AGCAT",
		"AGCGAT",
		"AGCTTT",
		"AGGAAT",
		"AGGAT",
		"AGGGAT",
		"AGGGGT",
		"AGTAAT",
		"AGTAT",
		"AGTGAT",
		"AGTTCT",
		"ATTCGT",
		"ATTTGT",
		"ATTTTT",
		"CAAAAG",
		"CAAACG",
		"CAACAG",
		"CAACCG",
		"CCTTGG",
		"CGAAAG",
		"CGAAG",
		"CGAGAG",
		"CGATAG",
		"CGCAAG",
		"CGCAG",
		"CGCGAG",
		"CGCTTG",
		"CGGAAG",
		"CGGAG",
		"CGGGAG",
		"CGGGGG",
		"CGTAAG",
		"CGTAG",
		"CGTGAG",
		"CGTTCG",
		"CTTCGG",
		"CTTTGG",
		"CTTTTG",
		"GAAAAC",
		"GAAAAT",
		"GAAACC",
		"GAAACT",
		"GAACAC",
		"GAACAT",
		"GCTTGC",
		"GCTTGT",
		"GGAAAC",
		"GGAAAT",
		"GGAAC",
		"GGAGAC",
		"GGAGAT",
		"GGATAC",
		"GGATAT",
		"GGCAAC",
		"GGCAAT",
		"GGCAC",
		"GGCGAC",
		"GGCGAT",
		"GGCTTC",
		"GGCTTT",
		"GGGAAC",
		"GGGAAT",
		"GGGAC",
		"GGGGAC",
		"GGGGAT",
		"GGGGGC",
		"GGGGGT",
		"GGTAAC",
		"GGTAAT",
		"GGTAC",
		"GGTGAC",
		"GGTGAT",
		"GGTTCC",
		"GTATAT",
		"GTTCGC",
		"GTTCGT",
		"GTTTGC",
		"GTTTGT",
		"GTTTTC",
		"GTTTTT",
		"TAAAAA",
		"TAAAAG",
		"TAAACA",
		"TAAACG",
		"TAACAA",
		"TAACAG",
		"TCTTGA",
		"TCTTGG",
		"TGAAA",
		"TGAAAA",
		"TGAAAG",
		"TGAGAA",
		"TGAGAG",
		"TGATAA",
		"TGATAG",
		"TGCAA",
		"TGCAAA",
		"TGCAAG",
		"TGCGAA",
		"TGCGAG",
		"TGCTTA",
		"TGCTTG",
		"TGGAA",
		"TGGAAA",
		"TGGAAG",
		"TGGGAA",
		"TGGGAG",
		"TGGGGA",
		"TGGGGG",
		"TGTAA",
		"TGTAAA",
		"TGTAAG",
		"TGTGAA",
		"TGTGAG",
		"TGTTCA",
		"TTTCGA",
		"TTTCGG",
		"TTTTAG",
		"TTTTGA",
		"TTTTGG",
		"TTTTTA",
		"TTTTTG"
	};
		
	for(unsigned int i = AA;i <= __;i++){
		
		for(unsigned int j = AA;j <= __;j++){
		
			cout << "param_H[" 
				<< pairs[i]
				<< "][" 
				<< pairs[j]
				<< "] = " << param_H[i][j] << endl;
				
			cout << "param_S[" 
				<< pairs[i]
				<< "][" 
				<< pairs[j]
				<< "] = " << param_S[i][j] << endl;
				
			cout << endl;
		}
	}
	
	for(unsigned int i = 0;i < NUM_SPECIAL_HAIRPIN_LOOP;i++){
		
		cout << "param_hairpin_special_H[" << hairpin_special_loop[i]
			<< "] = " << param_hairpin_special_H[i] << endl;
			
		cout << "param_hairpin_special_S[" << hairpin_special_loop[i]
			<< "] = " << param_hairpin_special_S[i] << endl;
	}
	
	for(unsigned int i = 0;i < MAX_LOOP_LENGTH + 1;i++){
		
		cout << "param_loop_S[" << i
			<< "] = " << param_loop_S[i] << endl;
	}
	
	for(unsigned int i = 0;i < MAX_BULGE_LENGTH + 1;i++){
		
		cout << "param_bulge_S[" << i
			<< "] = " << param_bulge_S[i] << endl;
	}
}

#endif // __DEBUG

string print_alignment(deque<nucleic_acid> &m_query, deque<nucleic_acid> &m_target)
{
	stringstream s;

	// This base map must match the order of the enumerated base values
	const char* base_map = "ACGTI$-MRSVWYHKDBN";
	
	s << "5' ";
	
	deque<nucleic_acid>::const_iterator q_iter, t_iter;
	
	for(q_iter = m_query.begin();q_iter != m_query.end();q_iter++){
	
		s << base_map[*q_iter];
	}
	
	s << " 3'" << endl;
	
	q_iter = m_query.begin();
	
	s << "   ";
	
	for(t_iter = m_target.begin();t_iter != m_target.end();t_iter++,q_iter++){
	
		enum {NO_MATCH, MATCH, INOSINE_MATCH, DEGENERATE_MATCH};
		
		unsigned int match = NO_MATCH;
		
		switch(*t_iter){
			case A:
				match = (*q_iter == T) ? MATCH : NO_MATCH;
				break;
			case T:
				match = (*q_iter == A) ? MATCH : NO_MATCH;
				break;
			case G:
				match = (*q_iter == C) ? MATCH : NO_MATCH;
				break;
			case C:
				match = (*q_iter == G) ? MATCH : NO_MATCH;
				break;
			case I:
				match = (*q_iter == GAP) ? NO_MATCH : INOSINE_MATCH;
				break;
			case E:
				match = (*q_iter == E) ? NO_MATCH : MATCH;
				break;
			case M:
				match = ( (*q_iter == A) || (*q_iter == C) ) ? DEGENERATE_MATCH : NO_MATCH;
				break;
			case R:
				match = ( (*q_iter == G) || (*q_iter == A) ) ? DEGENERATE_MATCH : NO_MATCH;
				break;
			case S:
				match = ( (*q_iter == G) || (*q_iter == C) ) ? DEGENERATE_MATCH : NO_MATCH;
				break;
			case V:
				match = ( (*q_iter == G) || (*q_iter == C) || (*q_iter == A) ) ? DEGENERATE_MATCH : NO_MATCH;
				break;
			case W:
				match = ( (*q_iter == A) || (*q_iter == T) ) ? DEGENERATE_MATCH : NO_MATCH;
				break;
			case Y:
				match = ( (*q_iter == T) || (*q_iter == C) ) ? DEGENERATE_MATCH : NO_MATCH;
				break;
			case H:
				match = ( (*q_iter == A) || (*q_iter == C) || (*q_iter == T) ) ? DEGENERATE_MATCH : NO_MATCH;
				break;
			case K:
				match = ( (*q_iter == G) || (*q_iter == T) ) ? DEGENERATE_MATCH : NO_MATCH;
				break;
			case D:
				match = ( (*q_iter == G) || (*q_iter == A) || (*q_iter == T) ) ? DEGENERATE_MATCH : NO_MATCH;
				break;
			case B:
				match = ( (*q_iter == G) || (*q_iter == T) || (*q_iter == C) ) ? DEGENERATE_MATCH : NO_MATCH;
				break;
			case N:
				// N matches everything!
				match = DEGENERATE_MATCH;
				break;
			case GAP:
				// Do nothing
				break;
		};
		
		if( (*q_iter == I) && (*t_iter != GAP) ){
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
				s << '!';
				break;
			case DEGENERATE_MATCH:
				s << '*';
				break;
		};
	}
	
	s << endl;
	
	s << "3' " ;
	
	for(t_iter = m_target.begin();t_iter != m_target.end();t_iter++){

		// dEBUG
		cerr << int(*t_iter) << '\t' << base_map[*t_iter] << endl;

		s << base_map[*t_iter];
	}
	
	s << " 5'" << endl;
	
	s << " alignment size = " << m_query.size();
	
	return s.str();
}

void NucCruc::set_supp_param(double m_param[NUM_SUPP_PARAM])
{

	// Set the linear models to the specified parameters
	for(unsigned int i = 0;i < NUM_SUPP_PARAM;i++){
		param_supp[i] = m_param[i];
	}
		
	// Update the parameters
	update_dp_param();
}

void NucCruc::set_supp_salt_param(double m_param[NUM_SALT_PARAM])
{

	// Set the linear models to the specified parameters
	for(unsigned int i = 0;i < NUM_SALT_PARAM;i++){
		param_supp_salt[i] = m_param[i];
	}
		
	// Update the parameters
	update_dp_param();
}

// Return the change in free energy on binding for the most recent melting temperature
// calculation as computing by the best scoring cell found in the dynamic programming
// matrix
float NucCruc::delta_G_dp() const
{
	//return curr_align.dp_dg;

	return curr_align.dp_dg + param_init_H - target_T*param_init_S /* initiation energy = 1.96 @ 37C */;
};
