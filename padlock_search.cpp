#include "tntblast.h"

using namespace std;

// Search for binding and ligation sites for Padlock probes / MOL-PCR (both
// have the same assay geometry).
//
// Here are all of the possible binding configurations for primers P1 and P2
// against a target sequence (for now, we assume that all sequences have both a 5'-3'
// and a 3'-5' strand -- obviously not true for single stranded viruses):
//
// ************************************
// The P1 primer is 5' of the P2 primer
// That is, 5'-P1-3' 5'-P2-3'
// ************************************
//
//                     3'-P1-5'3'-P2-5' (not allowed)
//	5'- =============================== -3'
//	3'- =============================== -5'
//         
//
//                     3'-P2-5'3'-P1-5'
//	5'- =============================== -3'
//	3'- =============================== -5'
//         
//                             
//	5'- =============================== -3'
//	3'- =============================== -5'
//         5'-P1-3'5'-P2-3'
//
//
//	5'- =============================== -3'
//	3'- =============================== -5'
//         5'-P2-3'5'-P1-3' (not allowed)
//
/////////////////////////////////////////////////////////////////////////////////
// Note that single primer reactions have been removed for padlock probes
// (they are not allowed due to assay contraints).
//
//                     3'-P1-5'3'-P1-5'
//	5'- =============================== -3'		* single primer amplification
//	3'- =============================== -5'
// 
//
//                     3'-P2-5'3'-P2-5'
//	5'- =============================== -3'		* single primer amplification
//	3'- =============================== -5'
// 
//
//
//	5'- =============================== -3'		* single primer amplification
//	3'- =============================== -5'
//         5'-P1-3'5'-P1-3
//
//
//	5'- =============================== -3'		* single primer amplification
//	3'- =============================== -5'
//         5'-P2-3'5'-P2-3
//

list<hybrid_sig> padlock(DNAHash &m_hash, const pair<string, SEQPTR> &m_seq, 
	const hybrid_sig &m_sig, NucCruc &m_melt,
	std::unordered_map<BindCacheKey, BindCacheValue> &m_plus_strand_melt_cache, 
	std::unordered_map<BindCacheKey, BindCacheValue> &m_minus_strand_melt_cache,
	const float &m_forward_primer_strand, 
	const float &m_reverse_primer_strand, 
	const float &m_min_probe_tm, const float &m_max_probe_tm,
	const float &m_min_probe_dg, const float &m_max_probe_dg,
	const unsigned int &m_probe_clamp_5, const unsigned int &m_probe_clamp_3,
	const unsigned int &m_max_gap, const unsigned int &m_max_mismatch,
	const int &m_target_strand)
{	
	const float forward_primer_strand = m_forward_primer_strand/m_sig.forward_degen;
	const float reverse_primer_strand = m_reverse_primer_strand/m_sig.reverse_degen;

	list<hybrid_sig> sig_list;
	
	list<oligo_info> upstream_bind;
	list<oligo_info> downstream_bind;
	
	list<oligo_info>::const_iterator up_iter;
	list<oligo_info>::const_iterator down_iter;
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// primer 2 upstream (5') and primer 1 downstream (3'): bind to minus strand
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// The upstream (5') primer has a 3' dangling base
	m_melt.set_query(m_sig.reverse_oligo);
	
	// Assume that the primer oligos are in vast excess to the target strands
	m_melt.strand(reverse_primer_strand, 0.0f);
	
	if(m_target_strand & Seq_strand_minus){
	
		bind_oligo_to_minus_strand(upstream_bind, 
			m_hash, m_seq.second, 
			m_sig.reverse_oligo,
			m_melt, m_minus_strand_melt_cache,
			m_min_probe_tm, m_max_probe_tm,
			m_min_probe_dg, m_max_probe_dg,
			m_probe_clamp_5, 	// Upstream primer has a 5' clamp,
			0,			// but no 3' clamp
			m_max_gap,
			m_max_mismatch); 	
	}
		
	// Test primer 1 in the downstream position
	m_melt.set_query(m_sig.forward_oligo);
	
	// Assume that the primer oligos are in vast excess to the target strands
	m_melt.strand(forward_primer_strand, 0.0f);
	
	if(m_target_strand & Seq_strand_minus){
	
		bind_oligo_to_minus_strand(downstream_bind, 
			m_hash, m_seq.second, 
			m_sig.forward_oligo,
			m_melt, m_minus_strand_melt_cache,
			m_min_probe_tm, m_max_probe_tm,
			m_min_probe_dg, m_max_probe_dg,
			0, 				// Downstream primer has a 3' clamp,
			m_probe_clamp_3,
			m_max_gap,
			m_max_mismatch); 	// but no 5' clamp
	}
	
	//   <---downstream                 upstream--->
	//
	//		         5'-P1-3'5'-P2-3'
	//	3'-======================================-5'  minus strand
	//	5'-======================================-3'  plus strand
	//
	for(up_iter = upstream_bind.begin();up_iter != upstream_bind.end();up_iter++){
		
		for(down_iter = downstream_bind.begin();down_iter != downstream_bind.end();down_iter++){
			
			// DEBUG
			//cerr << "\tdown = [" << down_iter->loc_5 << ", " << down_iter->loc_3 << "]" << endl;
			//cerr << "\tup = [" << up_iter->loc_5 << ", " << up_iter->loc_3 << "]" << endl;
			
			// loc coordinates are measured in the target plus strand
			if( (down_iter->loc_3 + 1) == up_iter->loc_5){
				
				// DEBUG
				//cerr << "\t21-" << endl;
				
				// These primers are next to each other!
				const int start = down_iter->loc_5;
				const int stop = up_iter->loc_3;
				
				if(start > stop){
					throw __FILE__ ":padlock: start > stop (4)";
				}
				
				const unsigned int len = stop - start + 1;
				
				hybrid_sig tmp;

				// Make a copy of the signature (including the id and the name)
				tmp = m_sig;

				// Downstream primer
				tmp.forward_oligo = m_sig.forward_oligo; 

				// Upstream primer
				tmp.reverse_oligo = m_sig.reverse_oligo;
				
				tmp.primer_strand = hybrid_sig::MINUS;
				tmp.amplicon_def = m_seq.first;
				
				tmp.amplicon_range.first = start;
				tmp.amplicon_range.second = stop;
				
				tmp.forward_tm = down_iter->tm;
				tmp.reverse_tm = up_iter->tm;
				
				tmp.forward_dH = down_iter->dH;
				tmp.reverse_dH = up_iter->dH;
				
				tmp.forward_dS = down_iter->dS;
				tmp.reverse_dS = up_iter->dS;
				
				tmp.forward_mm = down_iter->num_mm;
				tmp.reverse_mm = up_iter->num_mm;
				
				tmp.forward_gap = down_iter->num_gap;
				tmp.reverse_gap = up_iter->num_gap;
				
				tmp.forward_align = down_iter->alignment;
				tmp.reverse_align = up_iter->alignment;
				
				tmp.forward_primer_clamp = down_iter->anchor_3;
				tmp.reverse_primer_clamp = up_iter->anchor_5;
				
				// Extract the match sequence in the same orientation as the primers
				tmp.amplicon = string(len, '-');
				
				SEQPTR ptr = SEQ_START(m_seq.second) + max(0, start);
				
				for(unsigned int i = max(0, 1 - start);i < len;i++, ptr++){
					
					// Don't run past the end of the sequence
					if( (ptr - SEQ_START(m_seq.second) ) >= (int)SEQ_SIZE(m_seq.second) ){
						break;
					}
					
					tmp.amplicon[i] = hash_base_to_ascii(*ptr);
				}

				sig_list.push_back(tmp);
			}
		}
	}
	
	downstream_bind.clear();
	upstream_bind.clear();
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// primer 2 upstream (5') and primer 1 downstream (3'): bind to plus strand
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// The upstream (5') primer has a 3' dangling base
	m_melt.set_query(m_sig.reverse_oligo);
	
	// Assume that the primer oligos are in vast excess to the target strands
	m_melt.strand(reverse_primer_strand, 0.0f);
	
	if(m_target_strand & Seq_strand_plus){
	
		bind_oligo_to_plus_strand(upstream_bind, 
			m_hash, m_seq.second, 
			m_sig.reverse_oligo,
			m_melt, m_plus_strand_melt_cache,
			m_min_probe_tm, m_max_probe_tm,
			m_min_probe_dg, m_max_probe_dg,
			m_probe_clamp_5, 	// Upstream primer has a 5' clamp,
			0,			// but no 3' clamp
			m_max_gap,
			m_max_mismatch); 
	}
	
	// Test primer 1 in the downstream position
	m_melt.set_query(m_sig.forward_oligo);
	
	// Assume that the primer oligos are in vast excess to the target strands
	m_melt.strand(forward_primer_strand, 0.0f);
	
	if(m_target_strand & Seq_strand_plus){
	
		bind_oligo_to_plus_strand(downstream_bind, 
			m_hash, m_seq.second, 
			m_sig.forward_oligo,
			m_melt, m_plus_strand_melt_cache,
			m_min_probe_tm, m_max_probe_tm,
			m_min_probe_dg, m_max_probe_dg,
			0, 				// Downstream primer has a 3' clamp,
			m_probe_clamp_3, 		// but no 5' clamp
			m_max_gap,
			m_max_mismatch); 
	}
	
	//   <---upstream                 downstream--->
	//
	//		         3'-P2-5'3'-P1-5'
	//	5'-======================================-3'  plus strand
	//
	for(up_iter = upstream_bind.begin();up_iter != upstream_bind.end();up_iter++){
		
		for(down_iter = downstream_bind.begin();down_iter != downstream_bind.end();down_iter++){
			
			// loc coordinates are measured in the target plus strand
			if( (up_iter->loc_3 + 1) == down_iter->loc_5){
				
				// These primers are next to each other!
				const int start = up_iter->loc_5;
				const int stop = down_iter->loc_3;
				
				if(start > stop){
					throw __FILE__ ":padlock: start > stop (6)";
				}
				
				const unsigned int len = stop - start + 1;
				
				hybrid_sig tmp;

				// Make a copy of the signature (including the id and the name)
				tmp = m_sig;

				// Downstream primer
				tmp.forward_oligo = m_sig.forward_oligo; 

				// Upstream primer
				tmp.reverse_oligo = m_sig.reverse_oligo;
				
				tmp.primer_strand = hybrid_sig::PLUS;
				tmp.amplicon_def = m_seq.first;
				
				tmp.amplicon_range.first = start;
				tmp.amplicon_range.second = stop;
				
				tmp.forward_tm = down_iter->tm;
				tmp.reverse_tm = up_iter->tm;
				
				tmp.forward_dH = down_iter->dH;
				tmp.reverse_dH = up_iter->dH;
				
				tmp.forward_dS = down_iter->dS;
				tmp.reverse_dS = up_iter->dS;
				
				tmp.forward_mm = down_iter->num_mm;
				tmp.reverse_mm = up_iter->num_mm;
				
				tmp.forward_gap = down_iter->num_gap;
				tmp.reverse_gap = up_iter->num_gap;
				
				tmp.forward_align = down_iter->alignment;
				tmp.reverse_align = up_iter->alignment;
				
				tmp.forward_primer_clamp = down_iter->anchor_3;
				tmp.reverse_primer_clamp = up_iter->anchor_5;
				
				// Extract the match sequence in the same orientation as the primers
				tmp.amplicon = string(len, '-');
								
				SEQPTR ptr = SEQ_START(m_seq.second) + min( stop, (int)SEQ_SIZE(m_seq.second) - 1);
				
				for(unsigned int i = max(0, stop - (int)SEQ_SIZE(m_seq.second) - 1);i < len;i++, ptr--){
			
					// Don't run past the end of the sequence
					if( ptr < SEQ_START(m_seq.second) ){
						break;
					}
					
					tmp.amplicon[i] = hash_base_to_ascii_complement(*ptr);
				}

				sig_list.push_back(tmp);
			}
		}
	}
		
	return sig_list;
}
