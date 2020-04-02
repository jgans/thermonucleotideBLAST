// ThermonucleotideBLAST
// 
// Copyright (c) 2008, Los Alamos National Security, LLC
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

// DNA Sequence hashing code (for fast sequence searching)
// J. D. Gans
// Los Alamos National Laboratory
//
// Version 1
//	- Note this code is *not* thread safe due to the use of global
// 	variables in hash_db.cpp. In addition, the current use of these global
//	variables only allow for a *single* hash per process!
//
// Version 1.1 (10/25/06)
//	- Almost a complete rewrite. Functionality has been moved into a class (DNAHash) which
//	is accessed via an iterator. In addition to the double word hashing that was initially
//	used, I have added single word hashing (from words of length 2 to 12).
//
// Version 1.2 (2/7/07)
//	- Use a bitmask to track memory access for faster performance
//
// Version 2.0 (8/20/08)
// 	- A complete rewrite! "Better, faster, stronger, we can rebuild it ..."
//	- This version employs a repacking scheme to allow efficient access (and avoid cache thrashing).

#ifndef __HASH_DBASE
#define __HASH_DBASE

#include <string>
#include <string.h> // for memset
#include <vector>
#include "seq.h"

#define	SEQ_HASH_END	0xffffffff

const unsigned short int word_mask_lookup[] = 
{
	0x0, // 0 bits -- a place holder to align array elements starting at 0
	0x3, // 2 bits -- a place holder to align array elements starting at 0
	0xf, // 2^4 - 1
	0x3f, // 2^6 - 1
	0xff, // 2^8 - 1
	0x3ff, // 2^10 - 1
	0xfff, // 2^12 - 1
	0x3fff, // 2^14 - 1
	0xffff, // 2^16 - 1
};

class DNAHash_iterator;

class DNAHash {

	public:
		friend class DNAHash_iterator;
		typedef DNAHash_iterator iterator;
	private:
	
		std::pair<size_t, size_t> *hash_table;
		size_t *index_table;
		size_t index_table_size;
		unsigned char word_length;
		
		template <class SEQ>
		iterator _find(const SEQ &m_seq, const size_t &m_len, const bool &m_complement) const;
	public:
		
		DNAHash(const unsigned char &m_word_length = 8)
		{
			word_length = m_word_length;

			if( (word_length < 2) || (word_length > 8) ){
				throw __FILE__ ":DNAHash: Unsupported word length";
			}

			hash_table = NULL;
			index_table = NULL;
			index_table_size = 0;
		};
		
		~DNAHash()
		{
			clear();
		};
		
		inline void clear()
		{
			if(hash_table){
			
				delete [] hash_table;
				hash_table = NULL;
			}
						
			if(index_table){
			
				delete [] index_table;
				index_table = NULL;
				index_table_size = 0;
			}			
		};
		
		inline unsigned char word_size() const
		{
			return word_length;
		};
		
		// The smallest sequence that can be hashed
		inline size_t min_sequence_size() const
		{
			return size_t(word_length);
		};
		
		inline unsigned short int word_mask() const
		{
			return word_mask_lookup[word_length];
		};
		
		// The range of sequence data hashed is [m_start, m_stop).
		template <class SEQ>
		inline void hash(const SEQ &m_seq, const size_t &m_len, 
			const size_t &m_start, const size_t &m_stop);
		
		inline void hash(const std::string &m_seq, const size_t &m_start, const size_t &m_stop)
		{
			hash( m_seq, m_seq.size(), m_start, m_stop);
		};
		
		inline bool empty() const
		{
			return (hash_table == NULL);
		};
		
		template <class SEQ>
		inline iterator find(const SEQ &m_seq, const size_t &m_len) const;
		
		template <class SEQ>
		inline iterator find_complement(const SEQ &m_seq, const size_t &m_len) const;
		
		inline iterator find(const std::string &m_seq) const;
		inline iterator find_complement(const std::string &m_seq) const;
		
		inline iterator find(const SEQPTR &m_seq) const;
		inline iterator find_complement(const SEQPTR &m_seq) const;
		
		inline iterator end() const;
};

class DNAHash_iterator
{
	private:
		friend class DNAHash;
		
		const std::pair<size_t, size_t> *hash_table;
		const size_t *index_table;
		unsigned char word_length;
		std::vector<unsigned short int> word_list;
		
		size_t word_index;
		size_t hash_index;
		std::pair<size_t, size_t> hash_range;
		
		inline void clear()
		{
			hash_table = NULL;
			index_table = NULL;
			word_list.clear();
			word_index = 0;
			hash_index = SEQ_HASH_END;
		};
		
	public:
		DNAHash_iterator()
		{
			hash_table = NULL;
                        index_table = NULL;
                        word_list.clear();
                        word_index = 0;
                        hash_index = SEQ_HASH_END;
		};
		
		DNAHash_iterator(const DNAHash &m_hash) : hash_table(m_hash.hash_table),
			index_table(m_hash.index_table), word_length(m_hash.word_length)
		{
			word_index = 0;
			hash_index = SEQ_HASH_END;
		};
		
		DNAHash_iterator(const DNAHash_iterator &m_copy) : hash_table(m_copy.hash_table),
			index_table(m_copy.index_table), word_length(m_copy.word_length), 
			word_list(m_copy.word_list), word_index(m_copy.word_index), 
			hash_index(m_copy.hash_index), hash_range(m_copy.hash_range)
		{
			
		};
		
		~DNAHash_iterator()
		{
			// Do nothing
		};
		
		// offset() returns the index of the first matching base in
		// the query sequence
		inline size_t offset() const
		{
			return word_index;
		};
		
		// operator * () returns the index to the first matching base in
		// the target sequence
		inline size_t operator * ()
		{
			#ifdef _DEBUG
			// Is the user trying to end()++ ?
			if(hash_table == NULL){
				throw __FILE__ ":operator*: hash_ptr == NULL";
			}
			#endif // _DEBUG

			return index_table[hash_index];
		};
		
		inline bool operator == (const DNAHash_iterator &m_rhs) const
		{
			return ( (hash_table == m_rhs.hash_table) && (hash_index == m_rhs.hash_index) );
		};
		
		inline bool operator != (const DNAHash_iterator &m_rhs) const
		{
			return ( (hash_table != m_rhs.hash_table) || (hash_index != m_rhs.hash_index) );
		};
		
		// Prefix
		inline DNAHash_iterator& operator ++ ()
		{
			#ifdef _DEBUG
			// Is the user trying to end()++ ?
			if(hash_table == NULL){
				throw __FILE__ ":++DNAHash_iterator: hash_table == NULL";
			}
			#endif // _DEBUG
			
			if(++hash_index < hash_range.second){
				return (*this);
			}
			
			do{
				word_index ++;

				if( word_index >= word_list.size() ){

					// No matches were found!
					clear();
					return (*this);
				}

				hash_range = hash_table[ word_list[word_index] ];
			}
			while(hash_range.first == hash_range.second);
			
			hash_index = hash_range.first;
			
			return (*this);
		};

		// Postfix -- not as efficient as the prefix operator (has the additional overhead of 
		// the iterator copy).
		inline DNAHash_iterator operator ++ (int)
		{
			DNAHash_iterator ret(*this);

			++(*this);
			
			return ret;
		};
		
		template <class SEQ>
		inline void build_word_list(const SEQ &m_seq, const size_t &m_len, const bool &m_complement)
		{
			#ifdef _DEBUG
			// Is the user trying to end()++ ?
			if(hash_table == NULL){
				throw __FILE__ ":build_word_list: hash_table == NULL";
			}
			#endif // _DEBUG

			unsigned short int word = 0;
			size_t j = 1;
			
			word_list.clear();
						
			// Return now if there are not enough bases to build a single word
			if(word_length > m_len){
				return;
			}
			
			word_list.reserve( m_len - word_length + 1);
			
			if(m_complement == true){
				
				for(long int i = (long int)(m_len) - 1;i >= 0;i--, j++){

					switch(m_seq[i]){
						case 'A':
						case 'a':
							word = (word << 2) | DB_T;
							break;
						case 'T':
						case 't':
							word = (word << 2) | DB_A;
							break;
						case 'C':
						case 'c':
							word = (word << 2) | DB_G;
							break;
						case 'G':
						case 'g':
							word = (word << 2) | DB_C;
							break;
						default:
							// By default, treat non-ATGC bases as an 'A'
							word = (word << 2) | DB_A;
							break;
					};
					
					// Is this a valid word?
					if( j >= word_length ){
						word_list.push_back( word & word_mask_lookup[word_length] );
					}
				}
			}
			else{ // m_complement == false
				
				for(size_t i = 0;i < m_len;i++, j++){

					switch(m_seq[i]){
						case 'A':
						case 'a':
							word = (word << 2) | DB_A;
							break;
						case 'T':
						case 't':
							word = (word << 2) | DB_T;
							break;
						case 'C':
						case 'c':
							word = (word << 2) | DB_C;
							break;
						case 'G':
						case 'g':
							word = (word << 2) | DB_G;
							break;
						default:
							// By default, treat non-ATGC bases as an 'A'
							word = (word << 2) | DB_A;
							break;
					};
					
					// Is this a valid word?
					if( j >= word_length ){
						word_list.push_back( word & word_mask_lookup[word_length] );
					}
				}
			}
		};
};

// The use of a template function requires the following code to go in a 
// header file. The range of sequence data hashed is [m_start, m_stop).
template <class SEQ>
inline void DNAHash::hash(const SEQ &m_seq, const size_t &m_len, 
	const size_t &m_start, const size_t &m_stop)
{
	#ifdef _DEBUG
	if(m_start > m_stop){
		throw __FILE__ ":DNAHash::hash: m_start > m_stop";
	}
	
	if(m_stop > m_len){
		throw __FILE__ ":DNAHash::hash: m_stop > m_len";
	}
	#endif // _DEBUG
	
	const size_t table_size = word_mask() + 1;
	
	// Is this the first time we've been called? If so, allocate memory for
	// the hash table
	if(hash_table == NULL){
		
		hash_table = new std::pair<size_t, size_t> [table_size];
		
		if(hash_table == NULL){
			throw __FILE__ ":DNAHash::hash: Unable to allocate hash_table";
		}
		
	}
	
	memset( hash_table, 0, table_size*sizeof(std::pair<size_t, size_t>) );
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// The current hash scheme makes *two* passes through the sequence to be hashed! The first
	// pass computes the size of the contiguous match list at every position in the hash table. The second
	// pass populates those lists. On modern processors (i.e. Intel Xeon CPUs cache size ~ 2 Mb),
	// this scheme provides faster table look ups that a single pass hash table construction that
	// use discontinuous match lists (at the expense of slower hash table construction). It appears that
	// accessing the data in a discontinuous match list (where neighboring elements may be separated
	// by 10s of Mb) thrashes the CPU cache and kills performance!
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Compute the sizes of the hash table list required to store this sequence
	unsigned short int word = 0;
	size_t word_index = 0;
		
	for(size_t i = m_start;i < m_stop;i++){
		
		switch(m_seq[i]){
			case 'A':
			case 'a':
				word = (word << 2) | DB_A;
				word_index ++;
				break;
			case 'T':
			case 't':
				word = (word << 2) | DB_T;
				word_index ++;
				break;
			case 'C':
			case 'c':
				word = (word << 2) | DB_C;
				word_index ++;
				break;
			case 'G':
			case 'g':
				word = (word << 2) | DB_G;
				word_index ++;
				break;
			default:
				// By default, treat non-ATGC bases as an 'A'
				word = (word << 2) | DB_A;
				word_index ++;
				break;
		};
				
		// Is this a valid word?
		if(word_index >= word_length){
			
			hash_table[ word & word_mask() ].first++;
		}
	}
	
	// Partition the index_table into separate, continguous hash element lists
	size_t hash_start = 0;
	
	for(size_t i = 0;i < table_size;i++){
		
		const size_t tmp = hash_table[i].first;
		
		hash_table[i].first = hash_table[i].second = hash_start;
		hash_start += tmp;
	}
	
	// hash_start now contains the total number of hashed elements
	
	if(index_table_size < hash_start){
		
		if(index_table != NULL){
			delete [] index_table;
		}
		
		index_table_size = hash_start;
		
		index_table = new size_t [index_table_size];
		
		if(index_table == NULL){
			throw __FILE__ ":DNAHash::hash: Unable to allocate index_table";
		}
	}
	
	word = 0;
	word_index = 0;
	
	// Build the hash table using scratch space. This resulting array of indicies (stored in scratch)
	// will be repacked to allow sequential memory access. This step is crucial for efficient iteration
	// through the hash table (a factor of 5 speed up was observed on an Intel Xeon CPU due to a reduction
	// in cache misses).
	for(size_t i = m_start;i < m_stop;i++){
		
		switch(m_seq[i]){
			case 'A':
			case 'a':
				word = (word << 2) | DB_A;
				word_index ++;
				break;
			case 'T':
			case 't':
				word = (word << 2) | DB_T;
				word_index ++;
				break;
			case 'C':
			case 'c':
				word = (word << 2) | DB_C;
				word_index ++;
				break;
			case 'G':
			case 'g':
				word = (word << 2) | DB_G;
				word_index ++;
				break;
			default:
				// By default, treat non-ATGC bases as an 'A'
				word = (word << 2) | DB_A;
				word_index ++;
		};
				
		// Is this a valid word?
		if(word_index >= word_length){
			
			// Save the corresponding starting index value in the hash table
			//const size_t start_index = (i + 1) - word_length;
			//const size_t tmp_index = word & word_mask();
			
			index_table[ hash_table[ word & word_mask() ].second++ ] = (i + 1) - word_length;
		}
	}	
}

// The use of a template function requires the following code to go in a 
// header file. The range of sequence data hashed is [m_start, m_stop).
template <>
inline void DNAHash::hash(const SEQPTR &m_seq, const size_t &m_len, 
	const size_t &m_start, const size_t &m_stop)
{
	#ifdef _DEBUG
	if(m_start > m_stop){
		throw __FILE__ ":DNAHash::hash: m_start > m_stop";
	}
	
	if(m_stop > m_len){
		throw __FILE__ ":DNAHash::hash: m_stop > m_len";
	}
	#endif // _DEBUG
	
	const size_t table_size = word_mask() + 1;
	
	// Is this the first time we've been called? If so, allocate memory for
	// the hash table
	if(hash_table == NULL){
		
		hash_table = new std::pair<size_t, size_t> [table_size];
		
		if(hash_table == NULL){
			throw __FILE__ ":DNAHash::hash: Unable to allocate hash_table";
		}
		
	}
	
	memset( hash_table, 0, table_size*sizeof(std::pair<size_t, size_t>) );
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// The current hash scheme makes *two* passes through the sequence to be hashed! The first
	// pass computes the size of the contiguous match list at every position in the hash table. The second
	// pass populates those lists. On modern processors (i.e. Intel Xeon CPUs cache size ~ 2 Mb),
	// this scheme provides faster table look ups that a single pass hash table construction that
	// use discontinuous match lists (at the expense of slower hash table construction). It appears that
	// accessing the data in a discontinuous match list (where neighboring elements may be separated
	// by 10s of Mb) thrashes the CPU cache and kills performance!
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Compute the sizes of the hash table list required to store this sequence
	unsigned short int word = 0;
	size_t word_index = 0;
	
	SEQPTR ptr = SEQ_START(m_seq) + m_start;
	
	for(size_t i = m_start;i < m_stop;i++, ptr++){
		
		word = (word << 2) | (*ptr & DB_MASK);
		word_index ++;
				
		// Is this a valid word?
		if(word_index >= word_length){
			
			hash_table[ word & word_mask() ].first++;
		}
	}
	
	// Partition the index_table into separate, continguous hash element lists
	size_t hash_start = 0;
	
	for(size_t i = 0;i < table_size;i++){
		
		const size_t tmp = hash_table[i].first;
		
		hash_table[i].first = hash_table[i].second = hash_start;
		hash_start += tmp;
	}
	
	// hash_start now contains the total number of hashed elements
	
	if(index_table_size < hash_start){
		
		if(index_table != NULL){
			delete [] index_table;
		}
		
		index_table_size = hash_start;
		
		index_table = new size_t [index_table_size];
		
		if(index_table == NULL){
			throw __FILE__ ":DNAHash::hash: Unable to allocate index_table";
		}
	}
	
	word = 0;
	word_index = 0;
	
	ptr = SEQ_START(m_seq) + m_start;
	
	// Build the hash table using scratch space. This resulting array of indicies (stored in scratch)
	// will be repacked to allow sequential memory access. This step is crucial for efficient iteration
	// through the hash table (a factor of 5 speed up was observed on an Intel Xeon CPU due to a reduction
	// in cache misses).
	for(size_t i = m_start;i < m_stop;i++, ptr++){
		
		word = (word << 2) | (*ptr & DB_MASK);
		word_index ++;
				
		// Is this a valid word?
		if(word_index >= word_length){
			
			// Save the corresponding starting index value in the hash table
			//const size_t start_index = (i + 1) - word_length;
			//const size_t tmp_index = word & word_mask();
			
			index_table[ hash_table[ word & word_mask() ].second++ ] = (i + 1) - word_length;
		}
	}	
}

template <>
inline void DNAHash_iterator::build_word_list(const SEQPTR &m_seq, const size_t &m_len, const bool &m_complement)
{
	#ifdef _DEBUG
	// Is the user trying to end()++ ?
	if(hash_table == NULL){
		throw __FILE__ ":build_word_list: hash_table == NULL";
	}
	#endif // _DEBUG

	unsigned short int word = 0;
	size_t j = 1;

	word_list.clear();

	// Return now if there are not enough bases to build a single word
	if(word_length > m_len){
		return;
	}

	word_list.reserve( m_len - word_length + 1);

	SEQPTR ptr = SEQ_START(m_seq);
	
	if(m_complement == true){

		for(long int i = (long int)(m_len) - 1;i >= 0;i--, j++){

			word = (word << 2) | (DB_T - ( *ptr & DB_MASK) );
			
			// Is this a valid word?
			if( j >= word_length ){
				word_list.push_back( word & word_mask_lookup[word_length] );
			}
		}
	}
	else{ // m_complement == false

		for(size_t i = 0;i < m_len;i++, j++){

			word = (word << 2) | (*ptr & DB_MASK);

			// Is this a valid word?
			if( j >= word_length ){
				word_list.push_back( word & word_mask_lookup[word_length] );
			}
		}
	}
};


inline DNAHash::iterator DNAHash::end() const
{
	return DNAHash_iterator();
}
	
template <class SEQ>
inline DNAHash::iterator DNAHash::find(const SEQ &m_seq, const size_t &m_len) const
{
	return _find(m_seq, m_len, false);
};

template <class SEQ>
inline DNAHash::iterator DNAHash::find_complement(const SEQ &m_seq, const size_t &m_len) const
{
	return _find(m_seq, m_len, true);
};

inline DNAHash::iterator DNAHash::find(const std::string &m_seq) const
{
	return find( m_seq, m_seq.size() );
};

inline DNAHash::iterator DNAHash::find_complement(const std::string &m_seq) const
{
	return find_complement( m_seq, m_seq.size() );
};

inline DNAHash::iterator DNAHash::find(const SEQPTR &m_seq) const
{
	return find( m_seq, SEQ_SIZE(m_seq) );
};

inline DNAHash::iterator DNAHash::find_complement(const SEQPTR &m_seq) const
{
	return find_complement( m_seq, SEQ_SIZE(m_seq) );
};


template <class SEQ>
DNAHash::iterator DNAHash::_find(const SEQ &m_seq, const size_t &m_len, 
	const bool &m_complement) const
{
	DNAHash::iterator iter(*this);
	
	iter.build_word_list(m_seq, m_len, m_complement);
	
	if(iter.word_list.empty() == true){
		return end();
	}
	
	iter.hash_range = hash_table[ iter.word_list[iter.word_index] ];
	
	while(iter.hash_range.first == iter.hash_range.second){
		
		iter.word_index ++;
		
		if( iter.word_index >= iter.word_list.size() ){
		
			// No matches were found!
			return end();
		}
		
		iter.hash_range = hash_table[ iter.word_list[iter.word_index] ];
	}
	
	iter.hash_index = iter.hash_range.first;
	
	return iter;
}

#endif // __HASH_DBASE
