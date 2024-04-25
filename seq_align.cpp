#include "seq_align.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace std;
using namespace SA;

SeqAlign::SeqAlign(const AlignmentMode &m_mode, const bool &m_is_na):
	is_na(m_is_na)
{
	query_len = SIMD_SET1(0);
	target_len = SIMD_SET1(0);
	
	dp_matrix = NULL;
	query = NULL;
	target = NULL;
	
	max_loc.simd = SIMD_SET1(-1);

	mode = m_mode;
	
	for(unsigned int i = 0;i < SA_LEN;i++){
		curr_align[i].set_na(m_is_na);
	}

	// By default, do not use N as mask for na sequences
	disable_mask_N_na();
	
	if(is_na){
	
		// These are blastn defaults, note that other nucleotide alignment algorithms (like FASTA, megablast, etc.) 
		// use different paramters
		match.simd = SIMD_SET1(DEFAULT_NA_ALIGN_MATCH);
		mask.simd = SIMD_SET1(DEFAULT_NA_ALIGN_MASK);
		mismatch.simd = SIMD_SET1(DEFAULT_NA_ALIGN_MISMATCH);
		gap_existance.simd = SIMD_SET1(DEFAULT_NA_ALIGN_GAP_OPEN);
		gap_extension.simd = SIMD_SET1(DEFAULT_NA_ALIGN_GAP_EXTEND);
	}
	else{
		// Gaps costs from BLASTP default parameters
		gap_existance.simd = SIMD_SET1(DEFAULT_AA_ALIGN_GAP_OPEN);
		gap_extension.simd = SIMD_SET1(DEFAULT_AA_ALIGN_GAP_EXTEND);
		
		// There are number of possible protein score matricies: 
		// PAM30, PAM70,  BLOSUM45, BLOSUM62, BLOSUM80
		// By default, BLASTP used BLOSUM62
		init_BLOSUM62(aa_score_matrix);
	}
};

// Compute the alignment between two sequences; the query 
// and the target. Both the query and target sequences are assumed to be
// in 5'-3' orientation.
void SeqAlign::align_overlap()
{
	const simd_elem m_minus_mm( SIMD_SUB(match.simd, mismatch.simd) );
	const simd_elem zero(0);
	const simd_elem all_N(NA::N);
	const simd_elem query_length_minus_1( SIMD_SUB( query_len.simd, SIMD_SET1(1) ) );
	const simd_elem target_length_minus_1( SIMD_SUB( target_len.simd, SIMD_SET1(1) ) );
	
	const simd_elem _im1_jm1(im1_jm1);
	const simd_elem _i_jm1(i_jm1);
	const simd_elem _im1_j(im1_j);
	
	simd_elem max_score(-100); // A very small score!	
	simd_elem tmp_A;
	simd_elem tmp_B;
	simd_elem tmp_C;
		
	if(dp_matrix != NULL){
		
		SIMD_FREE(dp_matrix);
		dp_matrix = NULL;
	}
	
	const SA_Score max_query_len = query_len.max();
	const SA_Score max_target_len = target_len.max();
		
	dp_matrix = (SA_Elem*)SIMD_MALLOC( (max_query_len + 1)*(max_target_len + 1)*sizeof(SA_Elem), 
		MEMORY_ALIGNMENT);
	
	if(dp_matrix == NULL){
		throw __FILE__ ":SeqAlign::align_overlap: Unable to allocate dp_matrix";
	}
	
	// Initialize the dynamic programming matrix to have zeros along the first row and first column
	for(SA_Score i = 0;i <= max_target_len;i++){

		SA_Elem &elem_ref = dp_matrix[i];

		elem_ref.M = zero;
		elem_ref.I_query.simd = elem_ref.I_target.simd = gap_existance.simd;
		elem_ref.M_trace.simd = SIMD_SET1(query_target);
	}

	for(SA_Score i = 0;i <= max_query_len;i++){

		SA_Elem &elem_ref = dp_matrix[i*(max_target_len + 1)];

		elem_ref.M = zero;
		elem_ref.I_query.simd = elem_ref.I_target.simd = gap_existance.simd;
		elem_ref.M_trace.simd = SIMD_SET1(query_target);		
	}
	
	// Reset the max loc to -1
	max_loc.simd = SIMD_SET1(-1);
		
	for(SA_Score i = 0;i < max_query_len;i++){
			
		// The dp matrix has query_size rows and target_size columns
		// A B
		// C X <-- dp[i][j]
		SA_Elem *A_ptr = dp_matrix + i*(max_target_len + 1);
		SA_Elem *B_ptr = A_ptr + 1;
		SA_Elem *C_ptr = A_ptr + (max_target_len + 1);
		SA_Elem *X_ptr = C_ptr + 1;

		const simd_elem all_i(i);
		const simd_elem valid_query( SIMD_CMPGT(query_len.simd, all_i.simd) );
		const simd_elem query_is_N( SIMD_CMPEQ(query[i].simd, all_N.simd) );
		
		for(SA_Score j = 0;j < max_target_len;j++, A_ptr++, B_ptr++, C_ptr++, X_ptr++){
			
			tmp_A.simd = SIMD_ADD(
					SIMD_BITWISE_AND(
						SIMD_CMPGT(SIMD_BITWISE_AND(query[i].simd, target[j].simd), zero.simd),
						m_minus_mm.simd),
					mismatch.simd);

			if(mask_N_na){
				
				tmp_B.simd = SIMD_BITWISE_OR(query_is_N.simd, SIMD_CMPEQ(target[j].simd, all_N.simd) );

				tmp_A.simd = SIMD_BITWISE_OR( 
						SIMD_BITWISE_AND(tmp_B.simd, mask.simd), 
						SIMD_BITWISE_AND_NOT(tmp_B.simd, tmp_A.simd) );
			}

			X_ptr->M.simd = SIMD_ADD(
				SIMD_MAX( SIMD_MAX(A_ptr->M.simd, A_ptr->I_query.simd), A_ptr->I_target.simd ),
				 tmp_A.simd);
					
			X_ptr->M_trace.simd =  SIMD_BITWISE_AND_NOT(
							SIMD_BITWISE_OR(
								SIMD_CMPGT(A_ptr->I_query.simd, A_ptr->M.simd), 
								SIMD_CMPGT(A_ptr->I_target.simd, A_ptr->M.simd) ),
							_im1_jm1.simd);
							
			X_ptr->M_trace.simd =  SIMD_BITWISE_OR( X_ptr->M_trace.simd,
							SIMD_BITWISE_AND(
								SIMD_BITWISE_AND_NOT(
									SIMD_CMPGT(A_ptr->I_target.simd, A_ptr->I_query.simd),
									SIMD_CMPGT(A_ptr->I_query.simd, A_ptr->M.simd) ),
								_i_jm1.simd) );
							
			X_ptr->M_trace.simd =  SIMD_BITWISE_OR( X_ptr->M_trace.simd,
							SIMD_BITWISE_AND(
								SIMD_BITWISE_AND(
									SIMD_CMPGT(A_ptr->I_target.simd, A_ptr->M.simd),
									SIMD_CMPGT(A_ptr->I_target.simd, A_ptr->I_query.simd) ),
								_im1_j.simd) );


			// Unlike the Smith-Waterman alignments, don't clamp to zero
			tmp_A.simd = SIMD_ADD(C_ptr->M.simd, gap_existance.simd);
			
			tmp_B.simd = SIMD_ADD(C_ptr->I_query.simd, gap_extension.simd);
			
			X_ptr->I_query.simd = SIMD_MAX(tmp_A.simd, tmp_B.simd );

			tmp_C.simd = SIMD_CMPGT(tmp_B.simd, tmp_A.simd);
			
			X_ptr->I_query_trace.simd = SIMD_BITWISE_OR(
							SIMD_BITWISE_AND_NOT(
								tmp_C.simd, _im1_jm1.simd),
							SIMD_BITWISE_AND(
								tmp_C.simd, _i_jm1.simd) );
								
			// Unlike the Smith-Waterman alignments, don't clamp to zero
			tmp_A.simd =  SIMD_ADD(B_ptr->M.simd, gap_existance.simd);
			
			tmp_B.simd = SIMD_ADD(B_ptr->I_target.simd, gap_extension.simd);
			
			X_ptr->I_target.simd = SIMD_MAX(tmp_A.simd, tmp_B.simd);

			tmp_C.simd = SIMD_CMPGT(tmp_B.simd, tmp_A.simd);
			
			X_ptr->I_target_trace.simd = SIMD_BITWISE_OR(
							SIMD_BITWISE_AND_NOT(tmp_C.simd, _im1_jm1.simd),
							SIMD_BITWISE_AND(tmp_C.simd, _im1_j.simd) );
			
			// Force the alignment to go to the end of the query or the end of the target
			const simd_elem all_j( SIMD_SET1(j) );
			const simd_elem valid_query_and_target( 
				SIMD_BITWISE_AND(valid_query.simd, SIMD_CMPGT(target_len.simd, all_j.simd) ) );
							
			tmp_A.simd = SIMD_BITWISE_AND(valid_query_and_target.simd,
					SIMD_BITWISE_AND_NOT(
						SIMD_CMPGT(max_score.simd, X_ptr->M.simd),
						SIMD_BITWISE_OR(
							SIMD_CMPEQ(all_i.simd, query_length_minus_1.simd),
							SIMD_CMPEQ(all_j.simd, target_length_minus_1.simd) ) ) );
			
			max_score.simd = SIMD_BITWISE_OR(
						SIMD_BITWISE_AND_NOT(tmp_A.simd, max_score.simd),
						SIMD_BITWISE_AND(tmp_A.simd, X_ptr->M.simd) );
			
			max_loc.simd = SIMD_BITWISE_OR(
						SIMD_BITWISE_AND_NOT( tmp_A.simd, max_loc.simd),
						SIMD_BITWISE_AND(tmp_A.simd, SIMD_SET1(i*(max_target_len + 1) + j) ) );
		}
	}
}

// Compute the alignment between two sequences; the query 
// and the target. Both the query and target sequences are assumed to be
// in 5'-3' orientation.
void SeqAlign::align_global()
{
	const simd_elem m_minus_mm( SIMD_SUB(match.simd, mismatch.simd) );
	const simd_elem zero(0);
	const simd_elem all_N(NA::N);
	const simd_elem query_length_minus_1( SIMD_SUB( query_len.simd, SIMD_SET1(1) ) );
	const simd_elem target_length_minus_1( SIMD_SUB( target_len.simd, SIMD_SET1(1) ) );
	
	const simd_elem _im1_jm1(im1_jm1);
	const simd_elem _i_jm1(i_jm1);
	const simd_elem _im1_j(im1_j);
	
	simd_elem tmp_A;
	simd_elem tmp_B;
	simd_elem tmp_C;
		
	if(dp_matrix != NULL){
		
		SIMD_FREE(dp_matrix);
		dp_matrix = NULL;
	}
	
	const SA_Score max_query_len = query_len.max();
	const SA_Score max_target_len = target_len.max();
		
	dp_matrix = (SA_Elem*)SIMD_MALLOC( (max_query_len + 1)*(max_target_len + 1)*sizeof(SA_Elem), 
		MEMORY_ALIGNMENT);
	
	if(dp_matrix == NULL){
		throw __FILE__ ":SeqAlign::align_global: Unable to allocate dp_matrix";
	}
	
	// Initialize the dynamic programming matrix for a global alignment
	simd_elem boundary_penalty = gap_existance;
	
	for(SA_Score i = 1;i <= max_target_len;i++){

		// The first row ...
		SA_Elem &elem_ref = dp_matrix[i];

		// Add an additional penalty (of gap_existance) to the match score along the edge
		elem_ref.M = SIMD_ADD(boundary_penalty.simd, gap_existance.simd);
		
		elem_ref.I_query = boundary_penalty;
		elem_ref.I_target = boundary_penalty;

		elem_ref.M_trace.simd = SIMD_SET1(gap_target);
		elem_ref.I_target_trace.simd = SIMD_SET1(gap_target);
		elem_ref.I_query_trace.simd = SIMD_SET1(gap_target);

		boundary_penalty = SIMD_ADD(boundary_penalty.simd, gap_extension.simd);
	}

	boundary_penalty = gap_existance;
	
	for(SA_Score i = 1;i <= max_query_len;i++){

		// The first column ...
		SA_Elem &elem_ref = dp_matrix[i*(max_target_len + 1)];

		// Add an additional penalty (of gap_existance) to the match score along the edge
		elem_ref.M = SIMD_ADD(boundary_penalty.simd, gap_existance.simd);

		elem_ref.I_query = boundary_penalty;
		elem_ref.I_target = boundary_penalty;

		elem_ref.M_trace.simd = SIMD_SET1(query_gap);
		elem_ref.I_target_trace.simd = SIMD_SET1(query_gap);
		elem_ref.I_query_trace.simd = SIMD_SET1(query_gap);

		boundary_penalty = SIMD_ADD(boundary_penalty.simd, gap_extension.simd);
	}
	
	// The upper-left hand cell has a zero match score
	dp_matrix[0].M = dp_matrix[0].I_query = dp_matrix[0].I_target = zero;
	dp_matrix[0].M_trace = dp_matrix[0].I_target_trace = dp_matrix[0].I_query_trace = query_target;

	for(SA_Score i = 0;i < max_query_len;i++){
			
		// The dp matrix has query_size rows and target_size columns
		// A B
		// C X <-- dp[i][j]
		SA_Elem *A_ptr = dp_matrix + i*(max_target_len + 1);
		SA_Elem *B_ptr = A_ptr + 1;
		SA_Elem *C_ptr = A_ptr + (max_target_len + 1);
		SA_Elem *X_ptr = C_ptr + 1;

		const simd_elem query_is_N( SIMD_CMPEQ(query[i].simd, all_N.simd) );
		
		for(SA_Score j = 0;j < max_target_len;j++, A_ptr++, B_ptr++, C_ptr++, X_ptr++){
			
			tmp_A.simd = SIMD_ADD(
					SIMD_BITWISE_AND(
						SIMD_CMPGT(SIMD_BITWISE_AND(query[i].simd, target[j].simd), zero.simd),
						m_minus_mm.simd),
					mismatch.simd);

			if(mask_N_na){
				
				tmp_B.simd = SIMD_BITWISE_OR(query_is_N.simd, SIMD_CMPEQ(target[j].simd, all_N.simd) );

				tmp_A.simd = SIMD_BITWISE_OR( 
						SIMD_BITWISE_AND(tmp_B.simd, mask.simd), 
						SIMD_BITWISE_AND_NOT(tmp_B.simd, tmp_A.simd) );
			}

			X_ptr->M.simd = SIMD_ADD(
				SIMD_MAX( SIMD_MAX(A_ptr->M.simd, A_ptr->I_query.simd), A_ptr->I_target.simd ),
				 tmp_A.simd);
					
			X_ptr->M_trace.simd =  SIMD_BITWISE_AND_NOT(
							SIMD_BITWISE_OR(
								SIMD_CMPGT(A_ptr->I_query.simd, A_ptr->M.simd), 
								SIMD_CMPGT(A_ptr->I_target.simd, A_ptr->M.simd) ),
							_im1_jm1.simd);
							
			X_ptr->M_trace.simd =  SIMD_BITWISE_OR( X_ptr->M_trace.simd,
							SIMD_BITWISE_AND(
								SIMD_BITWISE_AND_NOT(
									SIMD_CMPGT(A_ptr->I_target.simd, A_ptr->I_query.simd),
									SIMD_CMPGT(A_ptr->I_query.simd, A_ptr->M.simd) ),
								_i_jm1.simd) );
							
			X_ptr->M_trace.simd =  SIMD_BITWISE_OR( X_ptr->M_trace.simd,
							SIMD_BITWISE_AND(
								SIMD_BITWISE_AND(
									SIMD_CMPGT(A_ptr->I_target.simd, A_ptr->M.simd),
									SIMD_CMPGT(A_ptr->I_target.simd, A_ptr->I_query.simd) ),
								_im1_j.simd) );


			// Don't clamp to zero in a global alignment
			tmp_A.simd = SIMD_ADD(C_ptr->M.simd, gap_existance.simd);
			
			tmp_B.simd = SIMD_ADD(C_ptr->I_query.simd, gap_extension.simd);
			
			X_ptr->I_query.simd = SIMD_MAX(tmp_A.simd, tmp_B.simd );

			tmp_C.simd = SIMD_CMPGT(tmp_B.simd, tmp_A.simd);
			
			X_ptr->I_query_trace.simd = SIMD_BITWISE_OR(
							SIMD_BITWISE_AND_NOT(
								tmp_C.simd, _im1_jm1.simd),
							SIMD_BITWISE_AND(
								tmp_C.simd, _i_jm1.simd) );
								
			// Don't clamp to zero in a global alignment
			tmp_A.simd = SIMD_ADD(B_ptr->M.simd, gap_existance.simd);
			
			tmp_B.simd = SIMD_ADD(B_ptr->I_target.simd, gap_extension.simd);
			
			X_ptr->I_target.simd = SIMD_MAX(tmp_A.simd, tmp_B.simd);

			tmp_C.simd = SIMD_CMPGT(tmp_B.simd, tmp_A.simd);
			
			X_ptr->I_target_trace.simd = SIMD_BITWISE_OR(
							SIMD_BITWISE_AND_NOT(tmp_C.simd, _im1_jm1.simd),
							SIMD_BITWISE_AND(tmp_C.simd, _im1_j.simd) );
		}
	}
	
	// Force the alignment to go to the end of the query *and* the end of the target
	// To do this, we set the max_loc value that corresponds to the 
	// end of each query paired with the end of each matching target.
	for(unsigned int i = 0;i < SA_LEN;++i){
		
		// Index from query len - 1 and target len - 1, since we will be adding
		// one to each of the matrix coordinates in the trace_back function.
		// Check for zero length sequences to make sure we don't generate an
		// invalid index.
		const size_t index = ( (query_len.v[i] > 0) ? query_len.v[i] - 1 : 0)*(max_target_len + 1) 
			+ ( (target_len.v[i] > 0) ? target_len.v[i] - 1 : 0 );
		
		max_loc.v[i] = index;
	}
}

void SeqAlign::align_smith_waterman()
{
	simd_elem max_score(-100); // A very small score!
	simd_elem tmp_A;
	simd_elem tmp_B;
	simd_elem tmp_C;
	const simd_elem m_minus_mm( SIMD_SUB(match.simd, mismatch.simd) );
	const simd_elem zero(0);
	const simd_elem all_N(NA::N);
	
	const simd_elem _im1_jm1(im1_jm1);
	const simd_elem _i_jm1(i_jm1);
	const simd_elem _im1_j(im1_j);
		
	if(dp_matrix != NULL){
		
		SIMD_FREE(dp_matrix);
		dp_matrix = NULL;
	}
	
	const SA_Score max_query_len = query_len.max();
	const SA_Score max_target_len = target_len.max();
		
	dp_matrix = (SA_Elem*)SIMD_MALLOC( (max_query_len + 1)*(max_target_len + 1)*sizeof(SA_Elem), 
		MEMORY_ALIGNMENT);
	
	if(dp_matrix == NULL){
		throw __FILE__ ":SeqAlign::align_smith_waterman: Unable to allocate dp_matrix";
	}
	
	// Initialize the dynamic programming matrix to have zeros along the first row and first column
	for(SA_Score i = 0;i <= max_target_len;i++){

		// The first row ...
		SA_Elem &elem_ref = dp_matrix[i];

		elem_ref.M = zero;
		elem_ref.I_query.simd = elem_ref.I_target.simd = gap_existance.simd;
		elem_ref.M_trace.simd = SIMD_SET1(query_target);
	}

	for(SA_Score i = 0;i <= max_query_len;i++){

		// The first column ...
		SA_Elem &elem_ref = dp_matrix[i*(max_target_len + 1)];

		elem_ref.M = zero;
		elem_ref.I_query.simd = elem_ref.I_target.simd = gap_existance.simd;
		elem_ref.M_trace.simd = SIMD_SET1(query_target);		
	}
	
	// Reset the max loc to -1
	max_loc.simd = SIMD_SET1(-1);
		
	for(SA_Score i = 0;i < max_query_len;i++){
			
		// The dp matrix has query_size rows and target_size columns
		// A B
		// C X <-- dp[i][j]
		SA_Elem *A_ptr = dp_matrix + i*(max_target_len + 1);
		SA_Elem *B_ptr = A_ptr + 1;
		SA_Elem *C_ptr = A_ptr + (max_target_len + 1);
		SA_Elem *X_ptr = C_ptr + 1;

		const simd_elem query_is_N( SIMD_CMPEQ(query[i].simd, all_N.simd) );
		
		for(SA_Score j = 0;j < max_target_len;j++, A_ptr++, B_ptr++, C_ptr++, X_ptr++){
			
			tmp_A.simd = SIMD_ADD(
					SIMD_BITWISE_AND(
						SIMD_CMPGT(SIMD_BITWISE_AND(query[i].simd, target[j].simd), zero.simd),
						m_minus_mm.simd),
					mismatch.simd);

			if(mask_N_na){
				
				tmp_B.simd = SIMD_BITWISE_OR(query_is_N.simd, SIMD_CMPEQ(target[j].simd, all_N.simd) );

				tmp_A.simd = SIMD_BITWISE_OR( 
						SIMD_BITWISE_AND(tmp_B.simd, mask.simd), 
						SIMD_BITWISE_AND_NOT(tmp_B.simd, tmp_A.simd) );
			}

			X_ptr->M.simd = SIMD_ADD(
				// Clamp to zero for a smith-waterman alignment
				SIMD_MAX(SIMD_MAX( SIMD_MAX(A_ptr->M.simd, A_ptr->I_query.simd), A_ptr->I_target.simd ), zero.simd),
				 tmp_A.simd);
					
			X_ptr->M_trace.simd =  SIMD_BITWISE_AND_NOT(
							SIMD_BITWISE_OR(
								SIMD_CMPGT(A_ptr->I_query.simd, A_ptr->M.simd), 
								SIMD_CMPGT(A_ptr->I_target.simd, A_ptr->M.simd) ),
							_im1_jm1.simd);
							
			X_ptr->M_trace.simd =  SIMD_BITWISE_OR( X_ptr->M_trace.simd,
							SIMD_BITWISE_AND(
								SIMD_BITWISE_AND_NOT(
									SIMD_CMPGT(A_ptr->I_target.simd, A_ptr->I_query.simd),
									SIMD_CMPGT(A_ptr->I_query.simd, A_ptr->M.simd) ),
								_i_jm1.simd) );
							
			X_ptr->M_trace.simd =  SIMD_BITWISE_OR( X_ptr->M_trace.simd,
							SIMD_BITWISE_AND(
								SIMD_BITWISE_AND(
									SIMD_CMPGT(A_ptr->I_target.simd, A_ptr->M.simd),
									SIMD_CMPGT(A_ptr->I_target.simd, A_ptr->I_query.simd) ),
								_im1_j.simd) );


			tmp_A.simd = SIMD_ADD(SIMD_MAX(C_ptr->M.simd, zero.simd), gap_existance.simd);
			
			tmp_B.simd = SIMD_ADD(SIMD_MAX(C_ptr->I_query.simd, zero.simd), gap_extension.simd);
			
			X_ptr->I_query.simd = SIMD_MAX(tmp_A.simd, tmp_B.simd );

			tmp_C.simd = SIMD_CMPGT(tmp_B.simd, tmp_A.simd);
			
			X_ptr->I_query_trace.simd = SIMD_BITWISE_OR(
							SIMD_BITWISE_AND_NOT(
								tmp_C.simd, _im1_jm1.simd),
							SIMD_BITWISE_AND(
								tmp_C.simd, _i_jm1.simd) );
			
			tmp_A.simd =  SIMD_ADD(SIMD_MAX(B_ptr->M.simd, zero.simd),gap_existance.simd);
			
			tmp_B.simd = SIMD_ADD(SIMD_MAX(B_ptr->I_target.simd, zero.simd), gap_extension.simd);
			
			X_ptr->I_target.simd = SIMD_MAX(tmp_A.simd, tmp_B.simd);

			tmp_C.simd = SIMD_CMPGT(tmp_B.simd, tmp_A.simd);
			
			X_ptr->I_target_trace.simd = SIMD_BITWISE_OR(
							SIMD_BITWISE_AND_NOT(tmp_C.simd, _im1_jm1.simd),
							SIMD_BITWISE_AND(tmp_C.simd, _im1_j.simd) );

			tmp_A.simd = SIMD_CMPGT(X_ptr->M.simd, max_score.simd);
			
			max_score.simd = SIMD_BITWISE_OR(
						SIMD_BITWISE_AND(tmp_A.simd, X_ptr->M.simd),
						SIMD_BITWISE_AND_NOT(tmp_A.simd, max_score.simd) );
			
			max_loc.simd = SIMD_BITWISE_OR(
						SIMD_BITWISE_AND( tmp_A.simd, SIMD_SET1(i*(max_target_len + 1) + j) ),
						SIMD_BITWISE_AND_NOT(tmp_A.simd, max_loc.simd) );
		}
	}
}

void SeqAlign::compute_score(const unsigned int &m_index)
{
	if( (max_loc.v[m_index] < 0) ){
		//throw __FILE__ ":SeqAlign::compute_score: Invalid pointer";
		return;
	}

	Alignment &align_ref = curr_align[m_index];
	
	// Remove any existing alignment information
	align_ref.clear();
	
	const size_t max_target_len = target_len.max();
	
	// Add 1 to last_i and last_j since we need the location of X (not A) in the
	// DP programming scheme.
	
	const int last_i = max_loc.v[m_index]/(max_target_len + 1) + 1;
	const int last_j = max_loc.v[m_index]%(max_target_len + 1) + 1;
			
	SA_Elem *max_ptr = dp_matrix + last_i*(max_target_len + 1) + last_j;

	if(max_ptr->M.v[m_index] >= max_ptr->I_query.v[m_index]){
		
		if(max_ptr->M.v[m_index] >= max_ptr->I_target.v[m_index]){
			align_ref.set_score(max_ptr->M.v[m_index]);
			
		}
		else{ // max_ptr->M.v[m_index] < max_ptr->I_target.v[m_index]
			align_ref.set_score(max_ptr->I_target.v[m_index]);
		}
	}
	else{ // max_ptr->M.v[m_index] < max_ptr->I_query.v[m_index]
		
		if(max_ptr->I_query.v[m_index] >= max_ptr->I_target.v[m_index]){
			align_ref.set_score(max_ptr->I_query.v[m_index]);
		}
		else{ // max_ptr->I_query.v[m_index] < max_ptr->I_target.v[m_index]
			align_ref.set_score(max_ptr->I_target.v[m_index]);
		}
	}
	
	// Compute the relative score if we have non-zero length sequences
	if( (query_len.v[m_index] > 0) && (target_len.v[m_index] > 0) ){
	
		SA_Score best_possible_score = 0;
		SA_Score worst_possible_score = 0;

		const int delta_len = max(query_len.v[m_index], target_len.v[m_index]) - 
					min(query_len.v[m_index], target_len.v[m_index]);

		switch(mode){
			case Overlap:
			case SmithWaterman:

				worst_possible_score = 0.0;
				best_possible_score = min(query_len.v[m_index], target_len.v[m_index])*
					get_match_score();
				break;
			case Global:

				//	---------- ACGATGCATGCATCG
				//	ATCGATAGCTG---------------
				worst_possible_score = (query_len.v[m_index] - 1)*get_gap_extension_score() + 
					(target_len.v[m_index] - 1)*get_gap_extension_score() +
					2.0*get_gap_open_score();

				// 	ACATAGCATGCATGCA----------------
				//	||||||||||||||||
				//	ACATAGCATGCATGCAAACGCGATAGCGCTAG
				best_possible_score = min(query_len.v[m_index], target_len.v[m_index])*
					get_match_score();

				if(delta_len > 0){
					best_possible_score += get_gap_open_score() + 
						(delta_len - 1)*get_gap_extension_score();
				}

				break;
			default:
				throw __FILE__ ":SeqAlign::compute_score: Unknown mode";
		};

		// Can we normalize the relative score?
		if(best_possible_score <= worst_possible_score){
			throw __FILE__ ":SeqAlign::compute_score: best_possible_score <= worst_possible_score";	
		}

		if(align_ref.get_score() < worst_possible_score){
			throw __FILE__ ":SeqAlign::compute_score: best_possible_score <= worst_possible_score";	
		}

		align_ref.set_relative_score(
			float(align_ref.get_score() - worst_possible_score)/
			(best_possible_score - worst_possible_score) );
	}
	else{
		align_ref.set_relative_score(0.0f);
	}
}

void SeqAlign::trace_back(const unsigned int &m_index)
{		
	if( (max_loc.v[m_index] < 0) ){
		throw __FILE__ ":SeqAlign::trace_back: Invalid pointer";
	}

	Alignment &align_ref = curr_align[m_index];
	
	// After computing the trace back, this alignment will be valid
	align_ref.set_valid();
	
	const size_t max_target_len = target_len.max();
	
	// Add 1 to last_i and last_j since we need the location of X (not A) in the
	// DP programming scheme.
	
	int last_i = max_loc.v[m_index]/(max_target_len + 1) + 1;
	int last_j = max_loc.v[m_index]%(max_target_len + 1) + 1;

	// DEBUG
	//cerr << "init last_i = " << last_i << endl;
	//cerr << "init last_j = " << last_j << endl << endl;
	
	SA_Elem *max_ptr = dp_matrix + last_i*(max_target_len + 1) + last_j;
	SA_Elem *curr_ptr = max_ptr;
	
	bool found_last_match = false;

	unsigned char curr_trace = max_ptr->M_trace.v[m_index];
	
	if(max_ptr->M.v[m_index] >= max_ptr->I_query.v[m_index]){
		
		if(max_ptr->M.v[m_index] >= max_ptr->I_target.v[m_index]){
			curr_trace = query_target;
		}
		else{ // max_ptr->M.v[m_index] < max_ptr->I_target.v[m_index]
			curr_trace = query_gap;
		}
	}
	else{ // max_ptr->M.v[m_index] < max_ptr->I_query.v[m_index]
		
		if(max_ptr->I_query.v[m_index] >= max_ptr->I_target.v[m_index]){
			curr_trace = gap_target;
		}
		else{ // max_ptr->I_query.v[m_index] < max_ptr->I_target.v[m_index]
			curr_trace = query_gap;
		}
	}
	
	while(true){
		
		switch(curr_trace){
			
			case query_target:

				switch(mode){
					case Overlap:
						// The alignment stops when we exhaust either the target *or* the query sequenecs
						if( (last_i < 1) || (last_j < 1) ){
							return;
						}
						break;
					case Global:
						// The alignment stops when we exhaust both the target *and* the query sequenecs
						if( (last_i < 1) && (last_j < 1) ){
							return;
						}
						break;
					case SmithWaterman:
						if( (last_i < 1) || (last_j < 1) || (curr_ptr->M.v[m_index] <= 0) ){
							return;
						}
						break;
				};
				
				--last_i;
				--last_j;
				
				if(!found_last_match){

					found_last_match = true;

					align_ref.set_query_stop(last_i);
					align_ref.set_target_stop(last_j);
				}

				align_ref.push_front(query[last_i].v[m_index], target[last_j].v[m_index]);
			
				curr_trace = curr_ptr->M_trace.v[m_index];
								
				break;
			case gap_target:
				
				switch(mode){
					case Overlap:
					
						// The alignment stops when we exhaust either the target or the query sequenecs
						if(last_j < 1){
							return;
						}
						break;
					case Global:
						// The alignment stops when we exhaust both the target and the query sequenecs
						if( (last_i < 1) && (last_j < 1) ){
							return;
						}
						break;
					case SmithWaterman:
					
						if( (last_j < 1) || (curr_ptr->I_query.v[m_index] <= 0)){
							return;
						}
						break;
				};

				if(last_j > 0){

					--last_j;
					
					// gap the query
					align_ref.push_front(NA::GAP, target[last_j].v[m_index]);
				}

				curr_trace = curr_ptr->I_query_trace.v[m_index];
				
				break;
			case query_gap:
			
				switch(mode){
					case Overlap:

						// The alignment stops when we exhaust either the target or the query sequenecs
						if(last_i < 1){
							return;
						}
						break;
					case Global:
						// The alignment stops when we exhaust both the target and the query sequenecs
						if( (last_i < 1) && (last_j < 1) ){
							return;
						}
						break;
					case SmithWaterman:
						if( (last_i < 1) || (curr_ptr->I_target.v[m_index] <= 0) ){
							return;
						}
						break;
				};

				if(last_i > 0){
				
					--last_i;
					
					// gap the target
					align_ref.push_front(query[last_i].v[m_index], NA::GAP);
				}

				curr_trace = curr_ptr->I_target_trace.v[m_index];

				break;
			default:
				
				// DEBUG
				//cerr << "last_i = " << last_i << endl;
				//cerr << "last_j = " << last_j << endl;
				
				throw __FILE__ ":SeqAlign::trace_back: Invalid_match in trace back";
		};
		
		align_ref.set_query_start(last_i);
		align_ref.set_target_start(last_j);
				
		if( (last_i < 0) || (last_j < 0) ){

			//cerr << "last_i = " << last_i << endl;
			//cerr << "last_j = " << last_j << endl;

			throw __FILE__ ":SeqAlign::trace_back: Index underflow in trace back";
		}
		
		// DEBUG
		//cerr << "curr_i = " << last_i << endl;
		//cerr << "curr_j = " << last_j << endl;
		//cerr << "curr_index = " << last_i*(max_target_len + 1) + last_j << endl;

		//switch(curr_trace){
		//	case query_target:
		//		cerr << "query_target" << endl << endl;
		//		break;
		//	case query_gap:
		//		cerr << "query_gap" << endl << endl;
		//		break;
		//	case gap_target:
		//		cerr << "gap_target" << endl << endl;
		//		break;
		//	case invalid_match:
		//		cerr << "invalid_match" << endl << endl;
		//		break;
		//	default:
		//		cerr << "Unknown trace!" << endl << endl;
		//		break;
		//};
				
		curr_ptr = dp_matrix + last_i*(max_target_len + 1) + last_j;
	}

	// It's an error if we get here
	throw __FILE__ ":trace_back: Error in trace_back";
}

string Alignment::print() const
{
	stringstream s;
	
	if(is_na){
		s << "5' ";
	}
	
	deque<base_type>::const_iterator q_iter, t_iter;
	
	for(q_iter = query_align.begin();q_iter != query_align.end();q_iter++){
	
		s << ( is_na ? bits_to_na(*q_iter) : bits_to_aa(*q_iter) );
	}
	
	if(is_na){
		s << " 3'" << endl;
	}
	else{
		s << endl;
	}
	
	q_iter = query_align.begin();
	
	if(is_na){
		s << "   ";
	}
	
	for(t_iter = target_align.begin();t_iter != target_align.end();t_iter++,q_iter++){
	
		enum {NO_MATCH, MATCH, DEGENERATE_MATCH};
		
		unsigned int match = NO_MATCH;
		
		if(*q_iter == *t_iter){
			match = MATCH;
		}
		else{
			if( is_na && (*q_iter & *t_iter) ){
				match = DEGENERATE_MATCH;
			}
		}
				
		switch(match){
			case NO_MATCH:
				s << ' ';
				break;
			case MATCH:
				s << '|';
				break;
			case DEGENERATE_MATCH:
				s << '*';
				break;
		};
	}
	
	s << endl;
	
	if(is_na){
		s << "5' " ;
	}
	
	for(t_iter = target_align.begin();t_iter != target_align.end();t_iter++){
		s << ( is_na ? bits_to_na(*t_iter) : bits_to_aa(*t_iter) );
	}
	
	if(is_na){
		s << " 3'" << endl;
	}
	else{
		s << endl;
	}
	
	s << "Alignment size = " << query_align.size();
	
	return s.str();
}

void SeqAlign::init_BLOSUM62(vector<SA_Score> &m_matrix)
{
	// From the NCBI blast source distribution: ~/ncbi/data/BLOSUM62
	const SA_Score matrix[AA::MATRIX_SIZE] = {
		4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, -1,
		-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1,
		-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1,
		-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1,
		0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -1,
		-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1,
		-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1,
		0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1,
		-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1,
		-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1,
		-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1,
		-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1,
		-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1,
		-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1,
		-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -1,
		1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, -1,
		0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, -1,
		-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -1,
		-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1,
		0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1,
		-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1,
		-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
	};

	m_matrix = vector<SA_Score>(matrix, matrix + AA::MATRIX_SIZE);
}

ostream& operator << (ostream &m_s, const Alignment &m_align)
{
	m_s << m_align.print();
	
	return m_s;
}


char SA::bits_to_na(base_type m_bits)
{
	switch(m_bits){
		case NA::A:
			return 'A';
		case NA::C:
			return 'C';
		case NA::G:
			return 'G';
		case NA::T:
			return 'T';
		case NA::M:
			return 'M';
		case NA::R:
			return 'R';
		case NA::S:
			return 'S';
		case NA::V:
			return 'V';
		case NA::W:
			return 'W';
		case NA::Y:
			return 'Y';
		case NA::H:
			return 'H';
		case NA::K:
			return 'K';
		case NA::D:
			return 'D';
		case NA::B:
			return 'B';
		case NA::N:
			return 'N';
		case NA::GAP:
			return '-';
	};

	throw __FILE__ ":bit_to_na: Unknown base!";
	return 'X'; // Keep the compiler happy
};

base_type SA::na_to_bits(char m_base)
{
	switch(m_base){
		case 'A': case 'a':
			return NA::A;
		case 'C': case 'c':
			return NA::C;
		case 'G': case 'g':
			return NA::G;
		case 'T': case 't':
			return NA::T;
		case 'M': case 'm':
			return NA::M;
		case 'R': case 'r':
			return NA::R;
		case 'S': case 's':
			return NA::S;
		case 'V': case 'v':
			return NA::V;
		case 'W': case 'w': 
			return NA::W;
		case 'Y': case  'y':
			return NA::Y;
		case 'H': case 'h':
			return NA::H;
		case 'K': case 'k':
			return NA::K;
		case 'D': case 'd':
			return NA::D;
		case 'B': case 'b':
			return NA::B;
		case 'N': case 'n':
			return NA::N;
		case '-':
			return NA::GAP;
	};

	throw __FILE__ ":na_to_bits: Unknown base!";
	return NA::GAP; // Keep the compiler happy
};

char SA::bits_to_aa(base_type m_bits)
{
	switch(m_bits){
		case AA::A:
			return 'A';
		case AA::R:
			return 'R';
		case AA::N:
			return 'N';
		case AA::D:
			return 'D';
		case AA::C:
			return 'C';
		case AA::Q:
			return 'Q';
		case AA::E:
			return 'E';
		case AA::G:
			return 'G';
		case AA::H:
			return 'H';
		case AA::I:
			return 'I';
		case AA::L:
			return 'L';
		case AA::K:
			return 'K';
		case AA::M:
			return 'M';
		case AA::F:
			return 'F';
		case AA::P:
			return 'P';
		case AA::S:
			return 'S';
		case AA::T:
			return 'T';
		case AA::W:
			return 'W';
		case AA::Y:
			return 'Y';
		case AA::V:
			return 'V';
		case AA::B:
			return 'B';
		case AA::Z:
			return 'Z';
		case AA::X:
			return 'X';
		case AA::GAP:
			return '-';
	};

	throw __FILE__ ":bit_to_aa: Unknown base!";
	return '?'; // Keep the compiler happy
};

base_type SA::aa_to_bits(char m_base)
{
	switch(m_base){
		case 'A': case 'a':
			return AA::A;
		case 'R': case 'r':
			return AA::R;
		case 'N': case 'n':
			return AA::N;
		case 'D': case 'd':
			return AA::D;
		case 'C': case 'c':
			return AA::C;
		case 'Q': case 'q':
			return AA::Q;
		case 'E': case 'e':
			return AA::E;
		case 'G': case 'g':
			return AA::G;
		case 'H': case 'h':
			return AA::H;
		case 'I': case 'i':
			return AA::I;
		case 'L': case 'l':
			return AA::L;
		case 'K': case 'k':
			return AA::K;
		case 'M': case 'm':
			return AA::M;
		case 'F': case 'f':
			return AA::F;
		case 'P': case 'p':
			return AA::P;
		case 'S': case 's':
			return AA::S;
		case 'T': case 't':
			return AA::T;
		case 'W': case 'w':
			return AA::W;
		case 'Y': case 'y':
			return AA::Y;
		case 'V': case 'v':
			return AA::V;
		case 'B': case 'b':
			return AA::B;
		case 'Z': case 'z':
			return AA::Z;
		case 'X': case 'x':
			return AA::X;
		case '-':
			return AA::GAP;
	};

	throw __FILE__ ":aa_to_bit: Unknown base!";
	return AA::GAP; // Keep the compiler happy
};

base_type SA::complement(base_type m_base)
{
	#ifdef _DEBUG
	if(!is_na){
		throw __FILE__ ":complement: Only valid for NA sequences";
	}
	#endif // _DEBUG

	switch(m_base){
		case NA::A:
			return NA::T;
		case NA::C:
			return NA::G;
		case NA::G:
			return NA::C;
		case NA::T:
			return NA::A;
		case NA::M:
			return NA::K;
		case NA::R:
			return NA::Y;
		case NA::S:
			return NA::S;
		case NA::V:
			return NA::B;
		case NA::W:
			return NA::W;
		case NA::Y:
			return NA::R;
		case NA::H:
			return NA::D;
		case NA::K:
			return NA::M;
		case NA::D:
			return NA::H;
		case NA::B:
			return NA::V;
		case NA::N:
			return NA::N;
		case NA::GAP:
			return NA::GAP;
	};

	throw __FILE__ ":complement: Unknown base";
	return NA::GAP; // Keep the compiler happy
};

deque<base_type> SA::reverse_complement(const deque<base_type> &m_seq)
{
	deque<base_type> ret(m_seq);
	
	reverse( ret.begin(), ret.end() );
	
	for(deque<base_type>::iterator i = ret.begin();i != ret.end();++i){
		*i = complement(*i);
	}
	
	return ret;
}
