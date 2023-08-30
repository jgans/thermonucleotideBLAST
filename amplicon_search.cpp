#include "tntblast.h"

using namespace std;

// A bit mask to pull out the strand information
#define STRAND_INFO (oligo_info::PLUS_STRAND | oligo_info::MINUS_STRAND)

pair<unsigned int, unsigned int> cull_oligo_match(list<oligo_info> &m_match_list, const unsigned int &m_max_amplicon_len, 
	const bool &m_has_probe, const bool &m_single_primer_pcr);

struct sort_by_oligo_loc
{
	inline bool operator()(const oligo_info &m_a, const oligo_info &m_b) const
	{
		if( !(m_a.loc_5 + m_a.loc_3) || !(m_b.loc_5 + m_b.loc_3) ){
			return m_a.target_loc < m_b.target_loc;
		}
		
		if(m_a.loc_5 == m_b.loc_5){
			return m_a.loc_3 < m_b.loc_3;
		}
		
		return m_a.loc_5 < m_b.loc_5;
	};
};

// Terminology:
//	P1 = Forward primer
//	P2 = Reverse primer
//
// Here are all of the possible binding configurations for primers P1 and P2
// against a target sequence:
//
//                             3'-P2-5'
//	5'- =============================== -3'	= plus strand
//	3'- =============================== -5' = minus strand
//         5'-P1-3'
//
//                             3'-P1-5'
//	5'- =============================== -3'
//	3'- =============================== -5'
//         5'-P2-3'
//
//////////////////////////////////////////////////////////////////////////////////
// The following cases are only tested if m_single_primer_pcr == true
//////////////////////////////////////////////////////////////////////////////////
//                             3'-P1-5'
//	5'- =============================== -3'		* single primer amplification
//	3'- =============================== -5'
//         5'-P1-3'
//
//                             3'-P2-5'
//	5'- =============================== -3'		* single primer amplification
//	3'- =============================== -5'
//         5'-P2-3'
//
list<hybrid_sig> amplicon(DNAHash &m_hash, 
	const pair<string, SEQPTR> &m_seq, 
	const hybrid_sig &m_sig, NucCruc &m_melt,
	unordered_map<BindCacheKey, BindCacheValue> &m_plus_strand_melt_cache,
	unordered_map<BindCacheKey, BindCacheValue> &m_minus_strand_melt_cache,
	const float &m_forward_primer_strand, 
	const float &m_reverse_primer_strand, const float &m_probe_strand, 
	const float &m_min_primer_tm, const float &m_max_primer_tm,
	const float &m_min_primer_dg, const float &m_max_primer_dg,
	const float &m_min_probe_tm, const float &m_max_probe_tm, 
	const float &m_min_probe_dg, const float &m_max_probe_dg, 
	const unsigned int &m_primer_clamp, 
	const int &m_min_max_primer_clamp, 
	const unsigned int &m_probe_clamp_5,
	const unsigned int &m_probe_clamp_3,
	const unsigned int &m_max_gap,
	const unsigned int &m_max_mismatch,
	const unsigned int &m_max_amplicon_len,
	const bool &m_single_primer_pcr)
{	
	// Only apply the min primer clamp test when we have a sensible (>= 0) value
	const bool apply_min_max_primer_clamp = (m_min_max_primer_clamp >= 0);
	
	const float forward_primer_strand = m_forward_primer_strand/m_sig.forward_degen;
	const float reverse_primer_strand = m_reverse_primer_strand/m_sig.reverse_degen;
	const float probe_strand = m_probe_strand/m_sig.probe_degen;

	unsigned int min_max_primer_clamp = apply_min_max_primer_clamp ? (unsigned int)m_min_max_primer_clamp : 0;
	
	list<hybrid_sig> sig_list;
	
	// Assemble a list of hash matches so we can cull *before* computing sequence alignments
	list<oligo_info> match_list;
	
	match_oligo_to_minus_strand(match_list, m_hash, m_sig.forward_oligo, oligo_info::F);
	match_oligo_to_minus_strand(match_list, m_hash, m_sig.reverse_oligo, oligo_info::R);
	
	// Did we find any oligo locations for either primer to bind to the
	// minus strand? If not, then we can stop looking right now
	const unsigned int num_minus_match = match_list.size();
	
	if(num_minus_match == 0){
		return sig_list;
	}	
	
	match_oligo_to_plus_strand(match_list, m_hash, m_sig.forward_oligo, oligo_info::F);
	match_oligo_to_plus_strand(match_list, m_hash, m_sig.reverse_oligo, oligo_info::R);
	
	const unsigned int num_plus_match = match_list.size();

	// If the number of matches has not increased, then we did not match any oligos to the plus strand
	if(num_plus_match == num_minus_match){
		return sig_list;
	}	
	
	if( m_sig.has_probe() ){
	
		match_oligo_to_minus_strand(match_list, m_hash, m_sig.probe_oligo, oligo_info::P);
		match_oligo_to_plus_strand(match_list, m_hash, m_sig.probe_oligo, oligo_info::P);
		
		// If the number of matches has not increased, then we did not match any probe oligos
		if(match_list.size() == num_plus_match){
			return sig_list;
		}
	}

	const pair<unsigned int, unsigned int> strand_count = 
		cull_oligo_match(match_list, m_max_amplicon_len, m_sig.has_probe(), m_single_primer_pcr);

	if(strand_count.first < strand_count.second){ // num minus < num plus
		
		m_melt.set_query(m_sig.forward_oligo);

		// Assume that the primer oligos are in vast excess to the target strands
		m_melt.strand(forward_primer_strand, 0.0f);
	
		// Attempt to bind forward primers to both the minus *and* plus strands
		bind_oligo_to_minus_strand(match_list, 
			oligo_info::F,
			m_seq.second, 
			m_sig.forward_oligo,
			m_melt, m_minus_strand_melt_cache,
			m_min_primer_tm, m_max_primer_tm,
			m_min_primer_dg, m_max_primer_dg,
			0, // no 5' clamp for primers
			m_primer_clamp,
			m_max_gap, m_max_mismatch);

		// Cull orphaned primers and probes
		cull_oligo_match(match_list, m_max_amplicon_len, m_sig.has_probe(), m_single_primer_pcr);
		
		if( match_list.empty() ){
			return sig_list;
		}
		
		m_melt.set_query(m_sig.reverse_oligo);
		
		// Assume that the primer oligos are in vast excess to the target strands
		m_melt.strand(reverse_primer_strand, 0.0f);
	
		// Attempt to bind forward primers to both the minus *and* plus strands
		bind_oligo_to_minus_strand(match_list, 
			oligo_info::R,
			m_seq.second, 
			m_sig.reverse_oligo,
			m_melt, m_minus_strand_melt_cache,
			m_min_primer_tm, m_max_primer_tm,
			m_min_primer_dg, m_max_primer_dg,
			0, // no 5' clamp for primers
			m_primer_clamp,
			m_max_gap, m_max_mismatch);
		
		// Cull orphaned primers and probes
		cull_oligo_match(match_list, m_max_amplicon_len, m_sig.has_probe(), m_single_primer_pcr);

		if( match_list.empty() ){
			return sig_list;
		}
	
		m_melt.set_query(m_sig.forward_oligo);
		
		// Assume that the primer oligos are in vast excess to the target strands
		m_melt.strand(forward_primer_strand, 0.0f);
	
		// Attempt to bind forward primers to both the minus *and* plus strands
		bind_oligo_to_plus_strand(match_list, 
			oligo_info::F,
			m_seq.second, 
			m_sig.forward_oligo,
			m_melt, m_plus_strand_melt_cache,
			m_min_primer_tm, m_max_primer_tm,
			m_min_primer_dg, m_max_primer_dg,
			0, // no 5' clamp for primers
			m_primer_clamp,
			m_max_gap, m_max_mismatch);

		// Cull orphaned primers and probes
		cull_oligo_match(match_list, m_max_amplicon_len, m_sig.has_probe(), m_single_primer_pcr);
		
		m_melt.set_query(m_sig.reverse_oligo);
		
		// Assume that the primer oligos are in vast excess to the target strands
		m_melt.strand(reverse_primer_strand, 0.0f);
	
		// Attempt to bind forward primers to both the minus *and* plus strands
		bind_oligo_to_plus_strand(match_list, 
			oligo_info::R,
			m_seq.second, 
			m_sig.reverse_oligo,
			m_melt, m_plus_strand_melt_cache,
			m_min_primer_tm, m_max_primer_tm,
			m_min_primer_dg, m_max_primer_dg,
			0, // no 5' clamp for primers
			m_primer_clamp,
			m_max_gap, m_max_mismatch);
	}
	else{ // num plus >= num minus

		m_melt.set_query(m_sig.forward_oligo);
		
		// Assume that the primer oligos are in vast excess to the target strands
		m_melt.strand(forward_primer_strand, 0.0f);
	
		// Attempt to bind forward primers to both the minus *and* plus strands
		bind_oligo_to_plus_strand(match_list, 
			oligo_info::F,
			m_seq.second, 
			m_sig.forward_oligo,
			m_melt, m_plus_strand_melt_cache,
			m_min_primer_tm, m_max_primer_tm,
			m_min_primer_dg, m_max_primer_dg,
			0, // no 5' clamp for primers
			m_primer_clamp,
			m_max_gap, m_max_mismatch);

		// Cull orphaned primers and probes
		cull_oligo_match(match_list, m_max_amplicon_len, m_sig.has_probe(), m_single_primer_pcr);
		
		if( match_list.empty() ){
			return sig_list;
		}
		
		m_melt.set_query(m_sig.reverse_oligo);
		
		// Assume that the primer oligos are in vast excess to the target strands
		m_melt.strand(reverse_primer_strand, 0.0f);
	
		// Attempt to bind forward primers to both the minus *and* plus strands
		bind_oligo_to_plus_strand(match_list, 
			oligo_info::R,
			m_seq.second, 
			m_sig.reverse_oligo,
			m_melt, m_plus_strand_melt_cache,
			m_min_primer_tm, m_max_primer_tm,
			m_min_primer_dg, m_max_primer_dg,
			0, // no 5' clamp for primers
			m_primer_clamp,
			m_max_gap, m_max_mismatch);
		
		// Cull orphaned primers and probes
		cull_oligo_match(match_list, m_max_amplicon_len, m_sig.has_probe(), m_single_primer_pcr);
		
		if( match_list.empty() ){
			return sig_list;
		}
		
		m_melt.set_query(m_sig.forward_oligo);
		
		// Assume that the primer oligos are in vast excess to the target strands
		m_melt.strand(forward_primer_strand, 0.0f);
	
		// Attempt to bind forward primers to both the minus *and* plus strands
		bind_oligo_to_minus_strand(match_list, 
			oligo_info::F,
			m_seq.second, 
			m_sig.forward_oligo,
			m_melt, m_minus_strand_melt_cache,
			m_min_primer_tm, m_max_primer_tm,
			m_min_primer_dg, m_max_primer_dg,
			0, // no 5' clamp for primers
			m_primer_clamp,
			m_max_gap, m_max_mismatch);

		// Cull orphaned primers and probes
		cull_oligo_match(match_list, m_max_amplicon_len, m_sig.has_probe(), m_single_primer_pcr);
		
		if( match_list.empty() ){
			return sig_list;
		}
		
		m_melt.set_query(m_sig.reverse_oligo);
		
		// Assume that the primer oligos are in vast excess to the target strands
		m_melt.strand(reverse_primer_strand, 0.0f);
	
		// Attempt to bind forward primers to both the minus *and* plus strands
		bind_oligo_to_minus_strand(match_list, 
			oligo_info::R,
			m_seq.second, 
			m_sig.reverse_oligo,
			m_melt, m_minus_strand_melt_cache,
			m_min_primer_tm, m_max_primer_tm,
			m_min_primer_dg, m_max_primer_dg,
			0, // no 5' clamp for primers
			m_primer_clamp,
			m_max_gap, m_max_mismatch);
	}

	// Reuse the melting engine for binding a probe (if present)
	if( m_sig.has_probe() ){
		
		// Cull orphaned primers and probes *before* attempting to bind probes
		cull_oligo_match(match_list, m_max_amplicon_len, m_sig.has_probe(), m_single_primer_pcr);

		if( match_list.empty() ){
			return sig_list;
		}
	
		m_melt.set_query(m_sig.probe_oligo);
		
		// Assume that the probes are in vast excess to the amplicons
		m_melt.strand(probe_strand, 0.0f);
		
		/////////////////////////////////////////////////////////////////////////////
		// Compute all probe binding locations to the plus and minus
		// strands. Originally, this calculation was delayed until
		// an amplicon was known to be produced. However, for certain sequences
		// (Bordetella pertusiss, I'm looking at you) that contain many repeats of
		// a primer motif, this lead to an explosion of amplicons to compute probes
		// for (performing the same calculation over and over again!)
		/////////////////////////////////////////////////////////////////////////////
				
		// Does the probe bind to the minus strand?
		bind_oligo_to_minus_strand(match_list, 
			oligo_info::P, m_seq.second, 
			m_sig.probe_oligo,
			m_melt, m_minus_strand_melt_cache,
			m_min_probe_tm, m_max_probe_tm,
			m_min_probe_dg, m_max_probe_dg,
			m_probe_clamp_5, m_probe_clamp_3,
			m_max_gap, m_max_mismatch);
		
		// Does the probe bind to the plus strand?
		bind_oligo_to_plus_strand(match_list,
			oligo_info::P, m_seq.second, 
			m_sig.probe_oligo,
			m_melt, m_plus_strand_melt_cache,
			m_min_probe_tm, m_max_probe_tm,
			m_min_probe_dg, m_max_probe_dg,
			m_probe_clamp_5, m_probe_clamp_3,
			m_max_gap, m_max_mismatch);
	}
	
	// We need one final sort before testing for assay matches
	match_list.sort( sort_by_oligo_loc() );

	// For every primer that can serve as a forward primer
	for(list<oligo_info>::iterator f = match_list.begin();f != match_list.end();f++){
		
		// Exclude probes and primers that bind to the plus strand
		if( (f->mask & (oligo_info::PLUS_STRAND | oligo_info::P) ) ){
			continue;
		}
		
		list<oligo_info>::iterator r = f;
		r++;
		
		for(;r != match_list.end();r++){

			// Exclude probes and primers that bind to the minus strand
			if( r->mask & (oligo_info::MINUS_STRAND | oligo_info::P) ){
				continue;
			}
			
			if( !m_single_primer_pcr && 
				( ( f->mask & (oligo_info::R | oligo_info::F) ) == 
				  ( r->mask & (oligo_info::R | oligo_info::F) ) ) ){
				  continue;
			}
			
			// Check the primer orientation
			if(f->loc_3 >= r->loc_5){
				continue;
			}

			// Check the amplicon length against the allowed length
			if( (r->loc_3 - f->loc_5 + 1) > (int)m_max_amplicon_len){
				continue;
			}

			if(apply_min_max_primer_clamp && 
				(max(f->anchor_3, r->anchor_3) <= min_max_primer_clamp) ){

				// Only exclude this assay if *both* primers have a small clamp
				continue;
			}

			if( m_sig.has_probe() ){
				
				list<oligo_info>::iterator p = f;
				p++;
				
				for(;p != r;p++){

					if(p->mask & oligo_info::P){

						// If we get here, we have a primer pair and probe!
						const int amp_start = f->loc_5;
						const int amp_stop = r->loc_3;

						if(amp_start > amp_stop){
							throw __FILE__ ":amplicon: amp_start > amp_stop";
						}
						
						// The original version only checked that the probe
						// bound to the amplicon
						if( !( (p->loc_5 >= amp_start) && (p->loc_3 <= amp_stop) ) ){
							continue;
						}
						
						// Now, we also check to make sure that the probe does not overlap
						// the binding site of the primer that binds to the same strand as
						// the probe (which would prevent probe hydrolysis and therfore
						// prevent the release of dye).
						if( (p->mask & STRAND_INFO) == (f->mask & STRAND_INFO) ){
							
							// The probe and *forward* primer are on the same strand
							if(p->loc_5 <= f->loc_3){							
								// The probe overlaps the forward primer
								continue;
							}
						}
						else{ 
						
							// The probe and *reverse* primer are on the same strand
							if(p->loc_3 >= r->loc_5){
								// The probe overlaps the reverse primer
								continue;
							}
						}
						
						///////////////////////////////////////////////////////////////
						// If we get here, we have a valid pair of primers and a probe
						const unsigned int amp_len = amp_stop - amp_start + 1;

						hybrid_sig tmp;

						// Make a copy of the signature (including the id and the name)
						tmp = m_sig;
												
						// Is this amplicon produced by a single primer?
						if( (f->mask & oligo_info::R) && (r->mask & oligo_info::R) ){
							
							// Two reverse oligos
							tmp.forward_oligo = m_sig.reverse_oligo;
						}
						
						if( (f->mask & oligo_info::F) && (r->mask & oligo_info::F) ){
							
							// Two forward oligos
							tmp.reverse_oligo = m_sig.forward_oligo;
						}
						
						tmp.primer_strand = (f->mask & oligo_info::F) ? hybrid_sig::PLUS : hybrid_sig::MINUS;
						
						tmp.amplicon_def = m_seq.first;

						tmp.amplicon_range.first = amp_start;
						tmp.amplicon_range.second = amp_stop;
						
						// Always print match information for the forward primer in the upstream position and
						// the reverse primer in the downstream position. If the primers bind in the reverse
						// orientation, then we need to swap them (for output only).
						list<oligo_info>::iterator f_output = f;
						list<oligo_info>::iterator r_output = r;
						
						if( (f->mask & oligo_info::R) && (r->mask & oligo_info::F) ){
							
							// Swap. The match order is *not* the output order
							swap(f_output, r_output);
							
						}
						
						tmp.forward_tm = f_output->tm;
						tmp.reverse_tm = r_output->tm;

						tmp.forward_dH = f_output->dH;
						tmp.reverse_dH = r_output->dH;

						tmp.forward_dS = f_output->dS;
						tmp.reverse_dS = r_output->dS;
						
						tmp.forward_mm = f_output->num_mm;
						tmp.reverse_mm = r_output->num_mm;
						
						tmp.forward_gap = f_output->num_gap;
						tmp.reverse_gap = r_output->num_gap;
						
						tmp.forward_primer_clamp = f_output->anchor_3;
						tmp.reverse_primer_clamp = r_output->anchor_3;

						tmp.forward_align = f_output->alignment;
						tmp.reverse_align = r_output->alignment;

						// Copy the amplicon bases in the orientation of primer 1
						// (starting from the first valid base)	
						if(tmp.primer_strand == hybrid_sig::PLUS){
						
							SEQPTR ptr = SEQ_START(m_seq.second) + max(0, amp_start);

							tmp.amplicon = string(amp_len, '-');

							for(unsigned int i = max(0, -amp_start);i < amp_len;i++, ptr++){

								// Don't run past the end of the sequence
								if( (ptr - SEQ_START(m_seq.second) ) >= (int)SEQ_SIZE(m_seq.second) ){
									break;
								}

								tmp.amplicon[i] = hash_base_to_ascii(*ptr);
							}
						}
						else{
							// Taking the complement of the amplicon
							SEQPTR ptr = SEQ_START(m_seq.second) + min(amp_stop, (int)SEQ_SIZE(m_seq.second) - 1);

							tmp.amplicon = string(amp_len, '-');

							for(unsigned int i = max(0, int(amp_stop) - (int)SEQ_SIZE(m_seq.second) + 1);i < amp_len;i++, ptr--){

								// Don't run past the end of the sequence
								if( ptr < SEQ_START(m_seq.second) ){
									break;
								}

								tmp.amplicon[i] = hash_base_to_ascii_complement(*ptr);
							}
						}
						
						tmp.probe_range = make_pair(p->loc_5, p->loc_3);
						tmp.probe_tm = p->tm;
						tmp.probe_dH = p->dH;
						tmp.probe_dS = p->dS;
						tmp.probe_mm = p->num_mm;
						tmp.probe_gap = p->num_gap;
						tmp.probe_strand = (p->mask & oligo_info::PLUS_STRAND) ? hybrid_sig::PLUS : hybrid_sig::MINUS;
						tmp.probe_align = p->alignment;
						
						// This is a valid solution
						sig_list.push_back(tmp);
					}
				}
			}
			else{

				// If we get here, we have a valid primer pair!
				const int amp_start = f->loc_5;
				const int amp_stop = r->loc_3;

				if(amp_start > amp_stop){
					throw __FILE__ ":amplicon: amp_start > amp_stop";
				}

				const unsigned int amp_len = amp_stop - amp_start + 1;

				hybrid_sig tmp;

				// Make a copy of the signature (including the id and the name)
				tmp = m_sig;

				// Is this amplicon produced by a single primer?
				if( (f->mask & oligo_info::R) && (r->mask & oligo_info::R) ){

					// Two reverse oligos
					tmp.forward_oligo = m_sig.reverse_oligo;
				}

				if( (f->mask & oligo_info::F) && (r->mask & oligo_info::F) ){

					// Two forward oligos
					tmp.reverse_oligo = m_sig.forward_oligo;
				}
				
				tmp.primer_strand = (f->mask & oligo_info::F) ? hybrid_sig::PLUS : hybrid_sig::MINUS;
				
				tmp.amplicon_def = m_seq.first;

				tmp.amplicon_range.first = amp_start;
				tmp.amplicon_range.second = amp_stop;
				
				// Always print match information for the forward primer in the upstream position and
				// the reverse primer in the downstream position. If the primers bind in the reverse
				// orientation, then we need to swap them (for output only).
				list<oligo_info>::iterator f_output = f;
				list<oligo_info>::iterator r_output = r;

				if( (f->mask & oligo_info::R) && (r->mask & oligo_info::F) ){

					// Swap. The match order is *not* the output order
					swap(f_output, r_output);
				}

				tmp.forward_tm = f_output->tm;
				tmp.reverse_tm = r_output->tm;

				tmp.forward_dH = f_output->dH;
				tmp.reverse_dH = r_output->dH;

				tmp.forward_dS = f_output->dS;
				tmp.reverse_dS = r_output->dS;
				
				tmp.forward_mm = f_output->num_mm;
				tmp.reverse_mm = r_output->num_mm;
								
				tmp.forward_gap = f_output->num_gap;
				tmp.reverse_gap = r_output->num_gap;

				tmp.forward_primer_clamp = f_output->anchor_3;
				tmp.reverse_primer_clamp = r_output->anchor_3;

				tmp.forward_align = f_output->alignment;
				tmp.reverse_align = r_output->alignment;

				// Copy the amplicon bases in the orientation of primer 1
				// (starting from the first valid base)	
				if(tmp.primer_strand == hybrid_sig::PLUS){

					SEQPTR ptr = SEQ_START(m_seq.second) + max(0, amp_start);

					tmp.amplicon = string(amp_len, '-');

					for(unsigned int i = max(0, -amp_start);i < amp_len;i++, ptr++){

						// Don't run past the end of the sequence
						if( (ptr - SEQ_START(m_seq.second) ) >= (int)SEQ_SIZE(m_seq.second) ){
							break;
						}

						tmp.amplicon[i] = hash_base_to_ascii(*ptr);
					}
				}
				else{

					// Taking the complement of the amplicon
					SEQPTR ptr = SEQ_START(m_seq.second) + min(amp_stop, (int)SEQ_SIZE(m_seq.second) - 1);

					tmp.amplicon = string(amp_len, '-');

					for(unsigned int i = max(0, int(amp_stop) - (int)SEQ_SIZE(m_seq.second) + 1);i < amp_len;i++, ptr--){

						// Don't run past the end of the sequence
						if( ptr < SEQ_START(m_seq.second) ){
							break;
						}

						tmp.amplicon[i] = hash_base_to_ascii_complement(*ptr);
					}
				}
				
				// This is a valid solution
				sig_list.push_back(tmp);
			}
		}
	}
	
	return sig_list;
}

pair<unsigned int, unsigned int> cull_oligo_match(list<oligo_info> &m_match_list, const unsigned int &m_max_amplicon_len, 
	const bool &m_has_probe, const bool &m_single_primer_pcr)
{	
	const unsigned int threshold = m_max_amplicon_len + 50;
	pair<unsigned int, unsigned int> ret = make_pair(0, 0);
	
	// ** It may be better to move this sort out of this function and into bind oligo **
	
	// Sort the oligo matches by the target_loc in ascending order
	m_match_list.sort( sort_by_oligo_loc() );
	
	// Unmask the valid bit
	for(list<oligo_info>::iterator i = m_match_list.begin();i != m_match_list.end();i++){
		i->mask &= ~oligo_info::VALID;
	}
	
	// For every primer that can serve as a forward primer
	for(list<oligo_info>::iterator f = m_match_list.begin();f != m_match_list.end();f++){
		
		if( (f->mask & (oligo_info::PLUS_STRAND | oligo_info::P) ) ){
			continue;
		}
		
		list<oligo_info>::iterator r = f;
		r++;
		
		for(;r != m_match_list.end();r++){
			
			if( (r->target_loc - f->target_loc) > threshold){
				break;
			}
			
			if( r->mask & (oligo_info::MINUS_STRAND | oligo_info::P) ){
				continue;
			}
			
			if( !m_single_primer_pcr && 
				( ( f->mask & (oligo_info::R | oligo_info::F) ) == 
				  ( r->mask & (oligo_info::R | oligo_info::F) ) ) ){
				  continue;
			}

			if(m_has_probe){
				
				list<oligo_info>::iterator p = f;
				p++;
				
				for(;p != r;p++){
				
					if(p->mask & oligo_info::P){
					
						p->mask |= oligo_info::VALID;
						f->mask |= oligo_info::VALID;
						r->mask |= oligo_info::VALID;
					}
				}
			}
			else{
				f->mask |= oligo_info::VALID;
				r->mask |= oligo_info::VALID;
			}
		}
	}
	
	// Remove all list elements that are not valid
	list<oligo_info>::iterator i = m_match_list.begin();
	
	while( i != m_match_list.end() ){
	
		if(i->mask & oligo_info::VALID){
		
			i++;
			
			ret.first += (i->mask & oligo_info::MINUS_STRAND) ? 1 : 0;
			ret.second += (i->mask & oligo_info::PLUS_STRAND) ? 1 : 0;
			
			continue;
		}
		
		list<oligo_info>::iterator del = i;
		i++;
		
		m_match_list.erase(del);
	}
	
	return ret;
}
