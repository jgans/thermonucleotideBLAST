#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include "tntblast.h"
#include <sstream>

using namespace std;

#ifdef PROFILE
extern unsigned int num_plus_tm_eval;
extern unsigned int num_minus_tm_eval;
#endif // PROFILE

struct sort_by_delta
{
	inline bool operator()(const pair<unsigned int, unsigned int> &m_a, 
		const pair<unsigned int, unsigned int> &m_b) const
	{
		return ( int(m_a.first) - int(m_a.second) ) < ( int(m_b.first) - int(m_b.second) );
	};
};

struct unique_by_delta
{
	inline bool operator()(const pair<unsigned int, unsigned int> &m_a, 
		const pair<unsigned int, unsigned int> &m_b) const
	{
		return ( int(m_a.first) - int(m_a.second) ) == ( int(m_b.first) - int(m_b.second) );
	};
};

struct sort_by_hash_match
{
	inline bool operator()(const oligo_info &m_a, const oligo_info &m_b) const
	{
		return ( int(m_a.query_loc) - int(m_a.target_loc) ) < ( int(m_b.query_loc) - int(m_b.target_loc) );
	};
};

struct unique_by_hash_match
{
	inline bool operator()(const oligo_info &m_a, const oligo_info &m_b) const
	{
		return ( ( int(m_a.query_loc) - int(m_a.target_loc) ) == ( int(m_b.query_loc) - int(m_b.target_loc) ) );
	};
};

struct sort_by_bound_match
{
	inline bool operator()(const oligo_info &m_a, const oligo_info &m_b) const
	{
		if(m_a.loc_5 != m_b.loc_5){
			return (m_a.loc_5 < m_b.loc_5);
		}
		
		if(m_a.loc_3 != m_b.loc_3){
			return (m_a.loc_3 < m_b.loc_3);
		}
		
		// If we get here, then the same target range was found
		if(m_a.tm == m_b.tm){
		
			// Due to boundary effects (i.e. how much sequence did we push into the 
			// comparison buffer) it is possible to get the same alignments but with differing
			// numbers of mismatches. Always take the alignment with the *highest* number of mismatches
			// (since lower mismatch counts are an artifact!). I wish there were a more elegant way
			// to handle this!
			if(m_a.num_mm == m_b.num_mm){

				// Choose the *longest* (most informative) alignment to preferentially to make sure that we show
				// any dangling ends to the user
				return ( m_a.alignment.size() > m_b.alignment.size() );
			}

			return (m_a.num_mm > m_b.num_mm);
		}
		
		// Reverse the tm comparison to sort by descending tm
		return (m_a.tm > m_b.tm);
	};
};

void match_oligo_to_minus_strand(list<oligo_info> &info_list, 
		DNAHash &m_hash, const string &m_oligo, const unsigned char &m_mask)
{
	list<oligo_info> hit_list;
	
	const unsigned char local_mask = m_mask | oligo_info::MINUS_STRAND;
	
	// Iterate over all possible word arrangements within this fragment of window bases. Note that the prefix iterator
	// increment (i.e. ++hash_iter) is more efficient than the postfix interator increment (i.e. hash_iter++) since the
	// latter requires an addition iterator copy operation.
	for(DNAHash::iterator hash_iter = m_hash.find(m_oligo);hash_iter != m_hash.end();++hash_iter){
		hit_list.push_back( oligo_info(hash_iter.offset(), *hash_iter, local_mask) );
	}
	
	hit_list.sort( sort_by_hash_match() );
	hit_list.unique( unique_by_hash_match() );

	info_list.merge(hit_list);
}

void match_oligo_to_plus_strand(list<oligo_info> &info_list, 
		DNAHash &m_hash, const string &m_oligo, const unsigned char &m_mask)
{
	list<oligo_info> hit_list;
		
	const unsigned char local_mask = m_mask | oligo_info::PLUS_STRAND;
		
	// Iterate over all possible word arrangements within this fragment of window bases. Note that the prefix iterator
	// increment (i.e. ++hash_iter) is more efficient than the postfix interator increment (i.e. hash_iter++) since the
	// latter requires an addition iterator copy operation.
	for(DNAHash::iterator hash_iter = m_hash.find_complement(m_oligo);hash_iter != m_hash.end();++hash_iter){
		hit_list.push_back( oligo_info(hash_iter.offset(), *hash_iter, local_mask) );
	}
	
	hit_list.sort( sort_by_hash_match() );
	hit_list.unique( unique_by_hash_match() );

	info_list.merge(hit_list);
}
	
void bind_oligo_to_minus_strand(list<oligo_info> &info_list, 
		DNAHash &m_hash, SEQPTR m_seq, 
		const string &m_oligo,
		NucCruc &m_melt, unordered_map<BindCacheKey, BindCacheValue> &m_melt_cache,
		const float &m_min_oligo_tm, const float &m_max_oligo_tm,
		const float &m_min_oligo_dg, const float &m_max_oligo_dg,
		const unsigned int &m_clamp_5,
		const unsigned int &m_clamp_3,
		const unsigned int &m_max_gap,
		const unsigned int &m_max_mismatch)
{		
	const unsigned int window = m_melt.size_query();

	// How many target bases should we add to both the 5' and 3' tails of the target sequence? 
	const unsigned int flanking_bases_5 = NUM_FLANK_BASE;
	const unsigned int flanking_bases_3 = NUM_FLANK_BASE;
	
	const unsigned int target_length = window + flanking_bases_5 + flanking_bases_3;	
	
	list<oligo_info> hit_list;

	// Find a unique list of all valid match seeds between the query (i.e. m_oligo_seq) and 
	// target (i.e. m_background_seq) sequences.
	list< pair<unsigned int, unsigned int> > match_seed; // query, target
		
	// Iterate over all possible word arrangements within this fragment of window bases. Note that the prefix iterator
	// increment (i.e. ++hash_iter) is more efficient than the postfix interator increment (i.e. hash_iter++) since the
	// latter requires an addition iterator copy operation.
	for(DNAHash::iterator hash_iter = m_hash.find(m_oligo);hash_iter != m_hash.end();++hash_iter){
		match_seed.push_back( make_pair(hash_iter.offset(), *hash_iter) );
	}
		
	match_seed.sort( sort_by_delta() );
	match_seed.unique( unique_by_delta() );
		
	// Iterate over *unique* word arrangements within this fragment of window bases
	for(list< pair<unsigned int, unsigned int> >::const_iterator match_iter = match_seed.begin();match_iter != match_seed.end();
		++match_iter){
		
		// The desired target offset (if we can get it)
		unsigned int target_start = max( int(match_iter->second) - int(match_iter->first /* query offset */ + flanking_bases_5),
			0);
		
		unsigned int target_stop = min( target_start + target_length, SEQ_SIZE(m_seq) );
		
		// Is there a result already in the cache?
		BindCacheKey key(m_oligo, target_start, target_stop);

		unordered_map<BindCacheKey, BindCacheValue>::const_iterator cache_iter = m_melt_cache.find(key);

		if( cache_iter == m_melt_cache.end() ){ // Cache miss
			
			// Clear the target buffer to make room for the bases we're about to add
			m_melt.clear_target();
		
			// The first valid base of the target
			SEQPTR seq_ptr = SEQ_START(m_seq) + target_start;
	
			// Attempt to extract target_length bases from the start
			for(unsigned int i = target_start;i < target_stop;i++, seq_ptr++){
			
				// Bind to the plus strand
				switch(*seq_ptr){
					case DB_A:
						m_melt.push_front_target(BASE::T);
						break;
					case DB_T:
						m_melt.push_front_target(BASE::A);
						break;
					case DB_C:
						m_melt.push_front_target(BASE::G);
						break;
					case DB_G:
						m_melt.push_front_target(BASE::C);
						break;
					case DB_I:
						m_melt.push_front_target(BASE::I);
						break;
					// IUPAC degenerate bases
					case DB_M: // A or C
						m_melt.push_front_target(BASE::K);
						break;
					case DB_R: // G or A
						m_melt.push_front_target(BASE::Y);
						break;
					case DB_S: // G or C
						m_melt.push_front_target(BASE::S);
						break;
					case DB_V: // G or C or A
						m_melt.push_front_target(BASE::B);
						break;
					case DB_W: // A or T
						m_melt.push_front_target(BASE::W);
						break;
					case DB_Y: // T or C
						m_melt.push_front_target(BASE::R);
						break;
					case DB_H: // A or C or T
						m_melt.push_front_target(BASE::D);
						break;
					case DB_K: // G or T
						m_melt.push_front_target(BASE::M);
						break;
					case DB_D: // G or A or T
						m_melt.push_front_target(BASE::H);
						break;
					case DB_B: // G or T or C
						m_melt.push_front_target(BASE::V);
						break;
					case DB_N: // A or T or G or C
						m_melt.push_front_target(BASE::N);
						break;
					default:
						
						// Are we closer to the start or the end?
						if( (i - target_start) > (target_stop - (i + 1) ) ){
						
							// Closer to the end
							target_stop = i;
						}
						else{
							// Closer to the start
							m_melt.clear_target();
							target_start = i + 1;
						}
						
						break;
				};
			}

			#ifdef PROFILE
			num_minus_tm_eval ++;
			#endif // PROFILE
		
			const float tm = m_melt.approximate_tm_heterodimer();
		
			if( (tm < m_min_oligo_tm) || (tm > m_max_oligo_tm) ){

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, 0.0 /*dg*/, 0.0 /*dH*/, 0.0 /*dS*/,
					0 /*anchor_5*/, 0 /*anchor_3*/, 0 /*target_5*/, 0 /*target_3*/,
					0 /*num_mismatch*/, 0 /*num_gap*/, 0.0 /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}

			const float dg = m_melt.delta_G();
		
			if( (dg < m_min_oligo_dg) || (dg > m_max_oligo_dg) ){

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					0 /*anchor_5*/, 0 /*anchor_3*/, 0 /*target_5*/, 0 /*target_3*/,
					0 /*num_mismatch*/, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);
					
				continue;
			}
		
			const unsigned int anchor_5 = m_melt.anchor5_query();
		
			if(anchor_5 < m_clamp_5){

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, 0 /*anchor_3*/, 0 /*target_5*/, 0 /*target_3*/,
					0 /*num_mismatch*/, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
		
			const unsigned int anchor_3 = m_melt.anchor3_query();
		
			if(anchor_3 < m_clamp_3){

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, anchor_3, 0 /*target_5*/, 0 /*target_3*/,
					0 /*num_mismatch*/, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
		
			const unsigned int num_mismatch = m_melt.num_mismatch();
		
			if(num_mismatch > m_max_mismatch){
				
				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, anchor_3, 0 /*target_5*/, 0 /*target_3*/,
					num_mismatch, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
		
			const unsigned int num_gap = m_melt.num_gap();
		
			if(num_gap > m_max_gap){
				
				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, anchor_3, 0 /*target_5*/, 0 /*target_3*/,
					num_mismatch, num_gap, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
		
			// How many "real" (non-degenerate, non-virtual) bases make up the target sequence?
			// Since degenerate bases get turned into 'N' by the hash function, we can
			// end up with spurious matches to poly-N sequences. Require that at least
			// half of the query can be aligned to real (ATGCI) bases
			const float fraction_real_base_pairs = m_melt.fraction_aligned_real_base_pairs();

			if(fraction_real_base_pairs < 0.5f){

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, anchor_3, 0 /*target_5*/, 0 /*target_3*/,
					num_mismatch, num_gap, fraction_real_base_pairs, string() /*align*/);

				continue;
			}

			pair<unsigned int, unsigned int> query_range;
			pair<unsigned int, unsigned int> target_range;
			
			m_melt.alignment_range(query_range, target_range);
		
			unsigned int target_5 = target_start;
			unsigned int target_3 = target_5;

			target_5 += target_stop - target_start - 1 - (int)target_range.second;
			target_3 += target_stop - target_start - 1 - (int)target_range.first;
		
			// Adjust the extent of the binding region to incorporate
			// unbound, non-virtual query bases
			target_5 -= (int)query_range.first;

			target_3 += (int)(window - 1) - (int)query_range.second;

			// Save the alignment
			stringstream ss_align;
			ss_align << m_melt;

			const string align = ss_align.str();

			hit_list.push_back(oligo_info( target_5, target_3, tm, 
				m_melt.delta_H(), m_melt.delta_S(), anchor_5, anchor_3, 
				num_mismatch, num_gap, align) );
			
			// Add this result to the cache
			m_melt_cache[key] = BindCacheValue(tm, dg, m_melt.delta_H(), m_melt.delta_S(),
				anchor_5, anchor_3, target_5, target_3,
				num_mismatch, num_gap, fraction_real_base_pairs, align);
		}
		else{ // Cache hit

			if( (cache_iter->second.tm < m_min_oligo_tm) || (cache_iter->second.tm > m_max_oligo_tm) ){
				continue;
			}

			if( (cache_iter->second.dg < m_min_oligo_dg) || (cache_iter->second.dg > m_max_oligo_dg) ){
				continue;
			}

			if(cache_iter->second.anchor_5 < m_clamp_5){
				continue;
			}

			if(cache_iter->second.num_mismatch > m_max_mismatch){
				continue;
			}

			if(cache_iter->second.num_gap > m_max_gap){
				continue;
			}

			if(cache_iter->second.fraction_real_base_pairs < 0.5f){
				continue;
			}

			hit_list.push_back(oligo_info( cache_iter->second.target_5, cache_iter->second.target_3, 
				cache_iter->second.tm, 
				cache_iter->second.dH, cache_iter->second.dS, 
				cache_iter->second.anchor_5, cache_iter->second.anchor_3, 
				cache_iter->second.num_mismatch, cache_iter->second.num_gap, cache_iter->second.seq_align) );
		}

	}
		
	info_list.clear();
	
	if(hit_list.empty() == true){
		return;
	}

	// Oligos can appear multiple times in the info_list list (but with the same 
	// target start and stop range). Remove these duplicate entries.
	hit_list.sort();
	
	// Only return unique oligos. Take the highest tm oligo out of every set of oligos that
	// matches a particular target range.
	info_list.push_back( hit_list.front() );
	hit_list.pop_front();
	
	while(hit_list.empty() == false){
	
		if( info_list.back() != hit_list.front() ){
			info_list.push_back( hit_list.front() );
		}
			
		hit_list.pop_front();
	}
}

void bind_oligo_to_minus_strand(list<oligo_info> &info_list, 
		const unsigned char &m_oligo_mask, SEQPTR m_seq, 
		const string &m_oligo,
		NucCruc &m_melt, unordered_map<BindCacheKey, BindCacheValue> &m_melt_cache,
		const float &m_min_oligo_tm, const float &m_max_oligo_tm,
		const float &m_min_oligo_dg, const float &m_max_oligo_dg,
		const unsigned int &m_clamp_5,
		const unsigned int &m_clamp_3,
		const unsigned int &m_max_gap,
		const unsigned int &m_max_mismatch)
{
	const unsigned int window = m_melt.size_query();
	
	// How many target bases should we add to both the 5' and 3' tails of the target sequence? 
	const unsigned int flanking_bases_5 = NUM_FLANK_BASE;
	const unsigned int flanking_bases_3 = NUM_FLANK_BASE;
	
	const unsigned int target_length = window + flanking_bases_5 + flanking_bases_3;	

	// Extract the target oligos into a separate list
	list<oligo_info> curr_oligo;
	
	list<oligo_info>::iterator info_iter = info_list.begin();
	
	const char strand_and_oligo_mask = m_oligo_mask | oligo_info::MINUS_STRAND;
	
	while( info_iter != info_list.end() ){
		
		//if( !(info_iter->mask & m_oligo_mask) || !(info_iter->mask & oligo_info::MINUS_STRAND) ){
		if( (info_iter->mask & strand_and_oligo_mask) != strand_and_oligo_mask ){

			++info_iter;
			continue;
		}
				
		// Add this element to the curr_oligo list and erase it from the info list 
		list<oligo_info>::iterator tmp_iter = info_iter;
		++info_iter;
		
		curr_oligo.push_front(*tmp_iter);
		info_list.erase(tmp_iter);
		
		list<oligo_info>::iterator match_iter = curr_oligo.begin();
				
		// The desired target offset (if we can get it)
		unsigned int target_start = max( int(match_iter->target_loc) - int(match_iter->query_loc /* query offset */ + flanking_bases_5),
			0);
		
		unsigned int target_stop = min( target_start + target_length, SEQ_SIZE(m_seq) );

		// Is there a result already in the cache?
		BindCacheKey key(m_oligo, target_start, target_stop);
		
		unordered_map<BindCacheKey, BindCacheValue>::const_iterator cache_iter = m_melt_cache.find(key);

		if( cache_iter == m_melt_cache.end() ){ // Cache miss

			// Clear the target buffer to make room for the bases we're about to add
			m_melt.clear_target();
			
			// The first valid base of the target
			SEQPTR seq_ptr = SEQ_START(m_seq) + target_start;

			// Attempt to extract target_length bases from the start
			for(unsigned int i = target_start;i < target_stop;i++, seq_ptr++){
				
				// Bind to the minus strand
				switch(*seq_ptr){
					case DB_A:
						m_melt.push_front_target(BASE::T);
						break;
					case DB_T:
						m_melt.push_front_target(BASE::A);
						break;
					case DB_C:
						m_melt.push_front_target(BASE::G);
						break;
					case DB_G:
						m_melt.push_front_target(BASE::C);
						break;
					case DB_I:
						m_melt.push_front_target(BASE::I);
						break;
					// IUPAC degenerate bases
					case DB_M: // A or C
						m_melt.push_front_target(BASE::K);
						break;
					case DB_R: // G or A
						m_melt.push_front_target(BASE::Y);
						break;
					case DB_S: // G or C
						m_melt.push_front_target(BASE::S);
						break;
					case DB_V: // G or C or A
						m_melt.push_front_target(BASE::B);
						break;
					case DB_W: // A or T
						m_melt.push_front_target(BASE::W);
						break;
					case DB_Y: // T or C
						m_melt.push_front_target(BASE::R);
						break;
					case DB_H: // A or C or T
						m_melt.push_front_target(BASE::D);
						break;
					case DB_K: // G or T
						m_melt.push_front_target(BASE::M);
						break;
					case DB_D: // G or A or T
						m_melt.push_front_target(BASE::H);
						break;
					case DB_B: // G or T or C
						m_melt.push_front_target(BASE::V);
						break;
					case DB_N: // A or T or G or C
						m_melt.push_front_target(BASE::N);
						break;
					default:
						
						// Are we closer to the start or the end?
						if( (i - target_start) > (target_stop - (i + 1) ) ){
						
							// Closer to the end
							target_stop = i;
						}
						else{
							// Closer to the start
							m_melt.clear_target();
							target_start = i + 1;
						}
						
						break;
				};
			}
			
			#ifdef PROFILE
			num_minus_tm_eval ++;
			#endif // PROFILE
			
			const float tm = m_melt.approximate_tm_heterodimer();

			if( (tm < m_min_oligo_tm) || (tm > m_max_oligo_tm) ){

				curr_oligo.pop_front();

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, 0.0f /*dg*/, 0.0f /*dH*/, 0.0f /*dS*/,
					0 /*anchor_5*/, 0 /*anchor_3*/, 0 /*target_5*/, 0 /*target_3*/,
					0 /*num_mismatch*/, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			const float dg = m_melt.delta_G();
			
			if( (dg < m_min_oligo_dg) || (dg > m_max_oligo_dg) ){

				curr_oligo.pop_front();

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					0 /*anchor_5*/, 0 /*anchor_3*/, 0 /*target_5*/, 0 /*target_3*/,
					0 /*num_mismatch*/, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			const unsigned int anchor_5 = m_melt.anchor5_query();
			
			if(anchor_5 < m_clamp_5){

				curr_oligo.pop_front();

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, 0 /*anchor_3*/, 0 /*target_5*/, 0 /*target_3*/,
					0 /*num_mismatch*/, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			const unsigned int anchor_3 = m_melt.anchor3_query();
			
			if(anchor_3 < m_clamp_3){

				curr_oligo.pop_front();

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, anchor_3, 0 /*target_5*/, 0 /*target_3*/,
					0 /*num_mismatch*/, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			const unsigned int num_mismatch = m_melt.num_mismatch();
			
			if(num_mismatch > m_max_mismatch){

				curr_oligo.pop_front();

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, anchor_3, 0 /*target_5*/, 0 /*target_3*/,
					num_mismatch, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			const unsigned int num_gap = m_melt.num_gap();
			
			if(num_gap > m_max_gap){

				curr_oligo.pop_front();

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, anchor_3, 0 /*target_5*/, 0 /*target_3*/,
					num_mismatch, num_gap, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			// How many "real" (non-degenerate, non-virtual) bases make up the target sequence?
			// Since degenerate bases get turned into 'N' by the hash function, we can
			// end up with spurious matches to poly-N sequences. Require that at least
			// half of the query can be aligned to real (ATGCI) bases
			const float fraction_real_base_pairs = m_melt.fraction_aligned_real_base_pairs();

			if(fraction_real_base_pairs < 0.5f){

				curr_oligo.pop_front();

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, anchor_3, 0 /*target_5*/, 0 /*target_3*/,
					num_mismatch, num_gap, fraction_real_base_pairs, string() /*align*/);

				continue;
			}

			pair<unsigned int, unsigned int> query_range;
			pair<unsigned int, unsigned int> target_range;
			
			m_melt.alignment_range(query_range, target_range);
			
			int target_5 = target_start;
			int target_3 = target_5;

			target_5 += target_stop - target_start - 1 - (int)target_range.second;
			target_3 += target_stop - target_start - 1 - (int)target_range.first;
			
			// Adjust the extent of the binding region to incorporate
			// unbound, non-virtual query bases
			target_5 -= (int)query_range.first;

			target_3 += (int)(window - 1) - (int)query_range.second;
			
			// Save the alignment
			stringstream ss_align;
			ss_align << m_melt;

			const string align = ss_align.str();

			match_iter->loc_5 = target_5;
			match_iter->loc_3 = target_3;
			match_iter->tm = tm;
			match_iter->dH = m_melt.delta_H();
			match_iter->dS = m_melt.delta_S();
			match_iter->anchor_5 = anchor_5;
			match_iter->anchor_3 = anchor_3;
			match_iter->num_mm = num_mismatch;
			match_iter->num_gap = num_gap;
			match_iter->alignment = align;

			// Add this result to the cache
			m_melt_cache[key] = BindCacheValue(tm, dg, m_melt.delta_H(), m_melt.delta_S(),
				anchor_5, anchor_3, target_5, target_3,
				num_mismatch, num_gap, fraction_real_base_pairs, align);
		}
		else{ // Cache hit

			if( (cache_iter->second.tm < m_min_oligo_tm) || (cache_iter->second.tm > m_max_oligo_tm) ){

				curr_oligo.pop_front();
				continue;
			}

			if( (cache_iter->second.dg < m_min_oligo_dg) || (cache_iter->second.dg > m_max_oligo_dg) ){

				curr_oligo.pop_front();
				continue;
			}

			if(cache_iter->second.anchor_5 < m_clamp_5){

				curr_oligo.pop_front();
				continue;
			}

			if(cache_iter->second.num_mismatch > m_max_mismatch){

				curr_oligo.pop_front();
				continue;
			}

			if(cache_iter->second.num_gap > m_max_gap){

				curr_oligo.pop_front();
				continue;
			}

			if(cache_iter->second.fraction_real_base_pairs < 0.5f){

				curr_oligo.pop_front();
				continue;
			}

			match_iter->loc_5 = cache_iter->second.target_5;
			match_iter->loc_3 = cache_iter->second.target_3;
			match_iter->tm = cache_iter->second.tm;
			match_iter->dH =cache_iter->second.dH;
			match_iter->dS = cache_iter->second.dS;
			match_iter->anchor_5 = cache_iter->second.anchor_5;
			match_iter->anchor_3 = cache_iter->second.anchor_3;
			match_iter->num_mm = cache_iter->second.num_mismatch;
			match_iter->num_gap = cache_iter->second.num_gap;
			match_iter->alignment = cache_iter->second.seq_align;
		}
	}
	
	if( curr_oligo.empty() ){
		return;
	}
	
	// Oligos can appear multiple times in the info_list list (but with the same 
	// target start and stop range). Remove these duplicate entries.
	curr_oligo.sort( sort_by_bound_match() );
	
	// Only return unique oligos. Take the highest tm oligo out of every set of oligos that
	// matches a particular target range.
	info_list.push_back( curr_oligo.front() );
	curr_oligo.pop_front();
	
	while( !curr_oligo.empty() ){
	
		if( info_list.back() != curr_oligo.front() ){
			info_list.push_back( curr_oligo.front() );
		}
			
		curr_oligo.pop_front();
	}
}

void bind_oligo_to_plus_strand(list<oligo_info> &info_list, 
		DNAHash &m_hash, SEQPTR m_seq, 
		const string &m_oligo,
		NucCruc &m_melt, unordered_map<BindCacheKey, BindCacheValue> &m_melt_cache,
		const float &m_min_oligo_tm, const float &m_max_oligo_tm,
		const float &m_min_oligo_dg, const float &m_max_oligo_dg,
		const unsigned int &m_clamp_5,
		const unsigned int &m_clamp_3,
		const unsigned int &m_max_gap,
		const unsigned int &m_max_mismatch)
{
	const unsigned int window = m_melt.size_query();
	
	// How many target bases should we add to both the 5' and 3' tails of the target sequence? 
	const unsigned int flanking_bases_5 = NUM_FLANK_BASE;
	const unsigned int flanking_bases_3 = NUM_FLANK_BASE;
	
	const unsigned int target_length = window + flanking_bases_5 + flanking_bases_3;	
	
	// Find a unique list of all valid match seeds between the query (i.e. m_oligo_seq) and 
	// target (i.e. m_background_seq) sequences.
	list< pair<unsigned int, unsigned int> > match_seed; // query, target
		
	// Iterate over all possible word arrangements within this fragment of window bases. Note that the prefix iterator
	// increment (i.e. ++hash_iter) is more efficient than the postfix interator increment (i.e. hash_iter++) since the
	// latter requires an addition iterator copy operation.
	for(DNAHash::iterator hash_iter = m_hash.find_complement(m_oligo);hash_iter != m_hash.end();++hash_iter){
		match_seed.push_back( make_pair(hash_iter.offset(), *hash_iter) );
	}
	
	match_seed.sort( sort_by_delta() );
	match_seed.unique( unique_by_delta() );
	
	list<oligo_info> hit_list;
	
	// Iterate over *unique* word arrangements within this fragment of window bases
	for(list< pair<unsigned int, unsigned int> >::const_iterator match_iter = match_seed.begin();match_iter != match_seed.end();
		++match_iter){
		
		// The desired target offset (if we can get it)
		unsigned int target_start = max( int(match_iter->second) - int(match_iter->first /* query offset */ + flanking_bases_5),
			0);
		
		unsigned int target_stop = min( target_start + target_length, SEQ_SIZE(m_seq) );
		
		// Is there a result already in the cache?
		BindCacheKey key(m_oligo, target_start, target_stop);

		unordered_map<BindCacheKey, BindCacheValue>::const_iterator cache_iter = m_melt_cache.find(key);

		if( cache_iter == m_melt_cache.end() ){ // Cache miss
			
			// Clear the target buffer to make room for the bases we're about to add
			m_melt.clear_target();
			
			// The first valid base of the target
			SEQPTR seq_ptr = SEQ_START(m_seq) + target_start;
		
			// Attempt to extract target_length bases from the start
			for(unsigned int i = target_start;i < target_stop;i++, seq_ptr++){
				
				// Bind to the plus strand
				switch(*seq_ptr){
					case DB_A:
						m_melt.push_back_target(BASE::A);
						break;
					case DB_T:
						m_melt.push_back_target(BASE::T);
						break;
					case DB_C:
						m_melt.push_back_target(BASE::C);
						break;
					case DB_G:
						m_melt.push_back_target(BASE::G);
						break;
					case DB_I:
						m_melt.push_back_target(BASE::I);
						break;
					// IUPAC degenerate bases
					case DB_M: // A or C
						m_melt.push_back_target(BASE::M);
						break;
					case DB_R: // G or A
						m_melt.push_back_target(BASE::R);
						break;
					case DB_S: // G or C
						m_melt.push_back_target(BASE::S);
						break;
					case DB_V: // G or C or A
						m_melt.push_back_target(BASE::V);
						break;
					case DB_W: // A or T
						m_melt.push_back_target(BASE::W);
						break;
					case DB_Y: // T or C
						m_melt.push_back_target(BASE::Y);
						break;
					case DB_H: // A or C or T
						m_melt.push_back_target(BASE::H);
						break;
					case DB_K: // G or T
						m_melt.push_back_target(BASE::K);
						break;
					case DB_D: // G or A or T
						m_melt.push_back_target(BASE::D);
						break;
					case DB_B: // G or T or C
						m_melt.push_back_target(BASE::B);
						break;
					case DB_N: // A or T or G or C
						m_melt.push_back_target(BASE::N);
						break;
					default:
						
						// Are we closer to the start or the end?
						if( (i - target_start) > (target_stop - (i + 1) ) ){
						
							// Closer to the end
							target_stop = i;
						}
						else{
							// Closer to the start
							m_melt.clear_target();
							target_start = i + 1;
						}
						
						break;
				};
			}

			#ifdef PROFILE
			num_plus_tm_eval ++;
			#endif // PROFILE
			
			const float tm = m_melt.approximate_tm_heterodimer();
			
			if( (tm < m_min_oligo_tm) || (tm > m_max_oligo_tm) ){

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, 0.0f /*dg*/, 0.0f /*dH*/, 0.0f /*dS*/,
					0 /*anchor_5*/, 0 /*anchor_3*/, 0 /*target_5*/, 0 /*target_3*/,
					0 /*num_mismatch*/, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			const float dg = m_melt.delta_G();
			
			if( (dg < m_min_oligo_dg) || (dg > m_max_oligo_dg) ){

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					0 /*anchor_5*/, 0 /*anchor_3*/, 0 /*target_5*/, 0 /*target_3*/,
					0 /*num_mismatch*/, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			const unsigned int anchor_5 = m_melt.anchor5_query();
			
			if(anchor_5 < m_clamp_5){

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, 0 /*anchor_3*/, 0 /*target_5*/, 0 /*target_3*/,
					0 /*num_mismatch*/, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			const unsigned int anchor_3 = m_melt.anchor3_query();
			
			if(anchor_3 < m_clamp_3){

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, anchor_3, 0 /*target_5*/, 0 /*target_3*/,
					0 /*num_mismatch*/, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			const unsigned int num_mismatch = m_melt.num_mismatch();
			
			if(num_mismatch > m_max_mismatch){

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, anchor_3, 0 /*target_5*/, 0 /*target_3*/,
					num_mismatch, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			const unsigned int num_gap = m_melt.num_gap();
			
			if(num_gap > m_max_gap){

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0 /*dH*/, 0.0 /*dS*/,
					anchor_5, anchor_3, 0 /*target_5*/, 0 /*target_3*/,
					num_mismatch, num_gap, 0.0 /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			// How many "real" (non-degenerate, non-virtual) bases make up the target sequence?
			// Since degenerate bases get turned into 'N' by the hash function, we can
			// end up with spurious matches to poly-N sequences. Require that at least
			// half of the query can be aligned to real (ATGCI) bases
			const float fraction_real_base_pairs = m_melt.fraction_aligned_real_base_pairs();

			if(fraction_real_base_pairs < 0.5f){

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0 /*dH*/, 0.0 /*dS*/,
					anchor_5, anchor_3, 0 /*target_5*/, 0 /*target_3*/,
					num_mismatch, num_gap, fraction_real_base_pairs, string() /*align*/);

				continue;
			}

			pair<unsigned int, unsigned int> query_range;
			pair<unsigned int, unsigned int> target_range;
			
			m_melt.alignment_range(query_range, target_range);
			
			int target_5 = target_start;
			int target_3 = target_5;

			target_5 += (int)target_range.first;
			target_3 += (int)target_range.second;

			// Adjust the extent of the binding region to incorporate
			// unbound, non-virtual query bases
			target_3 += (int)query_range.first;

			target_5 -= (int)(window - 1) - (int)query_range.second;
			
			// Save the alignment
			stringstream ss_align;
			ss_align << m_melt;
			
			const string align = ss_align.str();

			hit_list.push_back( oligo_info(target_5, target_3, tm, 
				m_melt.delta_H(), m_melt.delta_S(), anchor_5, anchor_3, 
				num_mismatch, num_gap, align) );
			
			// Add this result to the cache
			m_melt_cache[key] = BindCacheValue(tm, dg, m_melt.delta_H(), m_melt.delta_S(),
				anchor_5, anchor_3, target_5, target_3,
				num_mismatch, num_gap, fraction_real_base_pairs, align);
		}
		else{ // Cache hit

			if( (cache_iter->second.tm < m_min_oligo_tm) || (cache_iter->second.tm > m_max_oligo_tm) ){
				continue;
			}

			if( (cache_iter->second.dg < m_min_oligo_dg) || (cache_iter->second.dg > m_max_oligo_dg) ){
				continue;
			}

			if(cache_iter->second.anchor_5 < m_clamp_5){
				continue;
			}

			if(cache_iter->second.num_mismatch > m_max_mismatch){
				continue;
			}

			if(cache_iter->second.num_gap > m_max_gap){
				continue;
			}

			if(cache_iter->second.fraction_real_base_pairs < 0.5f){
				continue;
			}

			hit_list.push_back(oligo_info( cache_iter->second.target_5, cache_iter->second.target_3, 
				cache_iter->second.tm, 
				cache_iter->second.dH, cache_iter->second.dS, 
				cache_iter->second.anchor_5, cache_iter->second.anchor_3, 
				cache_iter->second.num_mismatch, cache_iter->second.num_gap, cache_iter->second.seq_align) );
		}
	}
	
	info_list.clear();
	
	if(hit_list.empty() == true){
		return;
	}

	// Oligos can appear multiple times in the info_list list (but with the same 
	// target start and stop range). Remove these duplicate entries.
	hit_list.sort();
	
	// Only return unique oligos. Take the highest tm oligo out of every set of oligos that
	// matches a particular target range.
	info_list.push_back( hit_list.front() );
	hit_list.pop_front();
	
	while(hit_list.empty() == false){
	
		if( info_list.back() != hit_list.front() ){
			info_list.push_back( hit_list.front() );			
		}
			
		hit_list.pop_front();
	}
}		

void bind_oligo_to_plus_strand(list<oligo_info> &info_list, 
		const unsigned char &m_oligo_mask, SEQPTR m_seq, 
		const string &m_oligo,
		NucCruc &m_melt, unordered_map<BindCacheKey, BindCacheValue> &m_melt_cache,
		const float &m_min_oligo_tm, const float &m_max_oligo_tm,
		const float &m_min_oligo_dg, const float &m_max_oligo_dg,
		const unsigned int &m_clamp_5,
		const unsigned int &m_clamp_3,
		const unsigned int &m_max_gap,
		const unsigned int &m_max_mismatch)
{
	const unsigned int window = m_melt.size_query();
	
	// How many target bases should we add to both the 5' and 3' tails of the target sequence? 
	const unsigned int flanking_bases_5 = NUM_FLANK_BASE;
	const unsigned int flanking_bases_3 = NUM_FLANK_BASE;
	
	const unsigned int target_length = window + flanking_bases_5 + flanking_bases_3;	
	
	// Extract the target oligos into a separate list
	list<oligo_info> curr_oligo;
	
	list<oligo_info>::iterator info_iter = info_list.begin();
	
	const char strand_and_oligo_mask = m_oligo_mask | oligo_info::PLUS_STRAND;
	
	while( info_iter != info_list.end() ){
		
		//if( !(info_iter->mask & m_oligo_mask) || !(info_iter->mask & oligo_info::PLUS_STRAND) ){
		if( (info_iter->mask & strand_and_oligo_mask) != strand_and_oligo_mask ){

			++info_iter;
			continue;
		}
				
		// Add this element to the curr_oligo list and erase it from the info list 
		list<oligo_info>::iterator tmp_iter = info_iter;
		++info_iter;
		
		curr_oligo.push_front(*tmp_iter);
		info_list.erase(tmp_iter);
		
		list<oligo_info>::iterator match_iter = curr_oligo.begin();
		
		// The desired target offset (if we can get it)
		unsigned int target_start = max( int(match_iter->target_loc) - int(match_iter->query_loc /* query offset */ + flanking_bases_5),
			0);
		
		unsigned int target_stop = min( target_start + target_length, SEQ_SIZE(m_seq) );
		
		// Is there a result already in the cache?
		BindCacheKey key(m_oligo, target_start, target_stop);

		unordered_map<BindCacheKey, BindCacheValue>::const_iterator cache_iter = m_melt_cache.find(key);

		if( cache_iter == m_melt_cache.end() ){ // Cache miss
			
			// Clear the target buffer to make room for the bases we're about to add
			m_melt.clear_target();
			
			// The first valid base of the target
			SEQPTR seq_ptr = SEQ_START(m_seq) + target_start;
		
			// Attempt to extract target_length bases from the start
			for(unsigned int i = target_start;i < target_stop;i++, seq_ptr++){
				
				// Bind to the plus strand
				switch(*seq_ptr){
					case DB_A:
						m_melt.push_back_target(BASE::A);
						break;
					case DB_T:
						m_melt.push_back_target(BASE::T);
						break;
					case DB_C:
						m_melt.push_back_target(BASE::C);
						break;
					case DB_G:
						m_melt.push_back_target(BASE::G);
						break;
					case DB_I:
						m_melt.push_back_target(BASE::I);
						break;
					// IUPAC degenerate bases
					case DB_M: // A or C
						m_melt.push_back_target(BASE::M);
						break;
					case DB_R: // G or A
						m_melt.push_back_target(BASE::R);
						break;
					case DB_S: // G or C
						m_melt.push_back_target(BASE::S);
						break;
					case DB_V: // G or C or A
						m_melt.push_back_target(BASE::V);
						break;
					case DB_W: // A or T
						m_melt.push_back_target(BASE::W);
						break;
					case DB_Y: // T or C
						m_melt.push_back_target(BASE::Y);
						break;
					case DB_H: // A or C or T
						m_melt.push_back_target(BASE::H);
						break;
					case DB_K: // G or T
						m_melt.push_back_target(BASE::K);
						break;
					case DB_D: // G or A or T
						m_melt.push_back_target(BASE::D);
						break;
					case DB_B: // G or T or C
						m_melt.push_back_target(BASE::B);
						break;
					case DB_N: // A or T or G or C
						m_melt.push_back_target(BASE::N);
						break;
					default:
										
						// Are we closer to the start or the end?
						if( (i - target_start) > (target_stop - (i + 1) ) ){
						
							// Closer to the end
							target_stop = i;
						}
						else{
							// Closer to the start
							m_melt.clear_target();
							target_start = i + 1;
						}

						break;
				};
			}

			#ifdef PROFILE
			num_plus_tm_eval ++;
			#endif // PROFILE
			
			const float tm = m_melt.approximate_tm_heterodimer();
			
			if( (tm < m_min_oligo_tm) || (tm > m_max_oligo_tm) ){

				curr_oligo.pop_front();

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, 0.0f /*dg*/, 0.0f /*dH*/, 0.0f /*dS*/,
					0 /*anchor_5*/, 0 /*anchor_3*/, 0 /*target_5*/, 0 /*target_3*/,
					0 /*num_mismatch*/, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			const float dg = m_melt.delta_G();
			
			if( (dg < m_min_oligo_dg) || (dg > m_max_oligo_dg) ){

				curr_oligo.pop_front();

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					0 /*anchor_5*/, 0 /*anchor_3*/, 0 /*target_5*/, 0 /*target_3*/,
					0 /*num_mismatch*/, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			const unsigned int anchor_5 = m_melt.anchor5_query();
			
			if(anchor_5 < m_clamp_5){

				curr_oligo.pop_front();

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, 0 /*anchor_3*/, 0 /*target_5*/, 0 /*target_3*/,
					0 /*num_mismatch*/, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			const unsigned int anchor_3 = m_melt.anchor3_query();
			
			if(anchor_3 < m_clamp_3){

				curr_oligo.pop_front();

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, anchor_3, 0 /*target_5*/, 0 /*target_3*/,
					0 /*num_mismatch*/, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			const unsigned int num_mismatch = m_melt.num_mismatch();
			
			if(num_mismatch > m_max_mismatch){

				curr_oligo.pop_front();

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, anchor_3, 0 /*target_5*/, 0 /*target_3*/,
					num_mismatch, 0 /*num_gap*/, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			const unsigned int num_gap = m_melt.num_gap();
			
			if(num_gap > m_max_gap){

				curr_oligo.pop_front();

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0f /*dH*/, 0.0f /*dS*/,
					anchor_5, anchor_3, 0 /*target_5*/, 0 /*target_3*/,
					num_mismatch, num_gap, 0.0f /*fraction_real_base_pairs*/, string() /*align*/);

				continue;
			}
			
			// How many "real" (non-degenerate, non-virtual) bases make up the target sequence?
			// Since degenerate bases get turned into 'N' by the hash function, we can
			// end up with spurious matches to poly-N sequences. Require that at least
			// half of the query can be aligned to real (ATGCI) bases
			const float fraction_real_base_pairs = m_melt.fraction_aligned_real_base_pairs();

			if(fraction_real_base_pairs < 0.5f){

				curr_oligo.pop_front();

				// Add this partial result to the cache
				m_melt_cache[key] = BindCacheValue(tm, dg, 0.0 /*dH*/, 0.0 /*dS*/,
					anchor_5, anchor_3, 0 /*target_5*/, 0 /*target_3*/,
					num_mismatch, num_gap, fraction_real_base_pairs, string() /*align*/);

				continue;
			}

			pair<unsigned int, unsigned int> query_range;
			pair<unsigned int, unsigned int> target_range;
			
			m_melt.alignment_range(query_range, target_range);
			
			int target_5 = target_start;
			int target_3 = target_5;

			target_5 += (int)target_range.first;
			target_3 += (int)target_range.second;

			// Adjust the extent of the binding region to incorporate
			// unbound, non-virtual query bases
			target_3 += (int)query_range.first;

			target_5 -= (int)(window - 1) - (int)query_range.second;

			// Save the alignment
			stringstream ss_align;
			ss_align << m_melt;
			
			const string align = ss_align.str();

			match_iter->loc_5 = target_5;
			match_iter->loc_3 = target_3;
			match_iter->tm = tm;
			match_iter->dH = m_melt.delta_H();
			match_iter->dS = m_melt.delta_S();
			match_iter->anchor_5 = anchor_5;
			match_iter->anchor_3 = anchor_3;
			match_iter->num_mm = num_mismatch;
			match_iter->num_gap = num_gap;
			match_iter->alignment = align;

			// Add this result to the cache
			m_melt_cache[key] = BindCacheValue(tm, dg, m_melt.delta_H(), m_melt.delta_S(),
				anchor_5, anchor_3, target_5, target_3,
				num_mismatch, num_gap, fraction_real_base_pairs, align);
		}
		else{ // Cache miss

			if( (cache_iter->second.tm < m_min_oligo_tm) || (cache_iter->second.tm > m_max_oligo_tm) ){

				curr_oligo.pop_front();
				continue;
			}

			if( (cache_iter->second.dg < m_min_oligo_dg) || (cache_iter->second.dg > m_max_oligo_dg) ){

				curr_oligo.pop_front();
				continue;
			}

			if(cache_iter->second.anchor_5 < m_clamp_5){

				curr_oligo.pop_front();
				continue;
			}

			if(cache_iter->second.num_mismatch > m_max_mismatch){

				curr_oligo.pop_front();
				continue;
			}

			if(cache_iter->second.num_gap > m_max_gap){

				curr_oligo.pop_front();
				continue;
			}

			if(cache_iter->second.fraction_real_base_pairs < 0.5f){

				curr_oligo.pop_front();
				continue;
			}

			match_iter->loc_5 = cache_iter->second.target_5;
			match_iter->loc_3 = cache_iter->second.target_3;
			match_iter->tm = cache_iter->second.tm;
			match_iter->dH =cache_iter->second.dH;
			match_iter->dS = cache_iter->second.dS;
			match_iter->anchor_5 = cache_iter->second.anchor_5;
			match_iter->anchor_3 = cache_iter->second.anchor_3;
			match_iter->num_mm = cache_iter->second.num_mismatch;
			match_iter->num_gap = cache_iter->second.num_gap;
			match_iter->alignment = cache_iter->second.seq_align;
		}
	}
	
	if( curr_oligo.empty() ){
		return;
	}
	
	// Oligos can appear multiple times in the info_list list (but with the same 
	// target start and stop range). Remove these duplicate entries.
	curr_oligo.sort( sort_by_bound_match() );
	
	// Only return unique oligos. Take the highest tm oligo out of every set of oligos that
	// matches a particular target range.
	info_list.push_back( curr_oligo.front() );
	curr_oligo.pop_front();
	
	while( !curr_oligo.empty() ){
	
		if( info_list.back() != curr_oligo.front() ){
			info_list.push_back( curr_oligo.front() );
		}
			
		curr_oligo.pop_front();
	}
}
