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

// bitmask v. 0.1
// J.D. Gans 9/14/05
//
// A class, very similar to bitset, that provides:
// 1) dynamic sizing
// 2) Comperable speed to bitset
// 3) The serialization funtions:
//	mpi_pack/mpi_unpack, read/write
//
// 12/14/05	Fixed bug in operator=(). Bits were not properly set.
// 12/19/05	Added mpi_pack and mpi_pack member functions
// 01/10/06	Added Set operations and helper functions:
//			operator=(const set<unsigned int> &m_set);
// 02/01/06	Fixed bug in clear() -- bits was not zero to zero.
//			Added addition converstions between bitmask and std::set<unsigned int>
// 08/29/06	Added read() and write() for binary I/O
// 09/06/06	Added operator=(const __bitelem &m_bit). This operator is required
//			for chaining bit copies; i.e. A[0] = B[0] = 1.
// 01/09/07	Added a toggle() member function to __bitelem to allow efficient bit
//			flipping.

#ifndef 	__BIT_MASK
#define	__BIT_MASK

#include <string.h> // For memset
#include <set>
#include <iostream>
#include <fstream>
#include <string>

class __bitelem{
	private:
		unsigned char *byte;
		unsigned char offset;
		
	public:
		inline __bitelem(unsigned char *m_byte, const unsigned char &m_offset)
		{
			byte = m_byte;
			offset = m_offset;
		};
		
		// Wow -- performance is dramatically improved by an explicit 
		// destructor.
		~__bitelem() 
		{ 
			// Do nothing!
		};
		
		inline __bitelem& operator=(const bool &m_bool)
		{
			*byte = (m_bool) ? ( *byte | (1 << offset) ) : ( *byte & ~(1 << offset) );
			
			return *this;
		};
		
		inline __bitelem& operator=(const __bitelem &m_bit)
		{
			*byte = bool(m_bit) ? ( *byte | (1 << offset) ) : ( *byte & ~(1 << offset) );
			
			return *this;
		};
		
		inline operator bool() const
		{
			return ( (*byte >> offset) & 1 );
		};
		
		inline void toggle()
		{
			// Similar code to operator=(const __bitelem &m_bit) but:
			// 1) explicity inlines the the bool() and
			// 2) swaps the arguments to ? (i.e. bool ? A : B --> bool ? B : A)
			*byte = ( (*byte >> offset) & 1 ) ? 
				( *byte & ~(1 << offset) ) : ( *byte | (1 << offset) );
		};
};

class bitmask {
	private:
		enum {BITS_PER_BYTE = 8};
		
		unsigned char* mask;
		size_t len;
		size_t bits;
		
		// Make sure that the left-over bits are always zero
		// (this will allow for efficient BYTE-wise comparisons)
		inline void zero_remainder()
		{
			// Set any unused bits to zero
			const unsigned char last_bit = (unsigned char)(bits%BITS_PER_BYTE);
			
			const unsigned char remainder_mask[] = {
				0xFF, // 11111111
				0x01, // 00000001,
				0x03, // 00000011
				0x07, // 00000111
				0x0F, // 00001111
				0x1F, // 00011111
				0x3F, // 00111111
				0x7F  // 01111111
				};
			
			if(last_bit != 0){
				
				// len == 0 only if last_bit == 0
				mask[len - 1] &= remainder_mask[last_bit];
			}
			
		};
		
	public:
		bitmask()
		{
			mask = NULL;
			len = 0;
			bits = 0;
		};
		
		bitmask(const unsigned int &m_num_bits, const bool &m_bool = false)
		{
			bits = m_num_bits;
			
			if(bits == 0){
				mask = NULL;
				len = 0;
				return;
			}
			
			len = (bits/BITS_PER_BYTE) + 
				(bits%BITS_PER_BYTE != 0);
			
			mask = new unsigned char [len];
			
			if(mask == NULL){
				throw "bitmask: Unable to allocate memory!";
			}
			
			// Flip bits to either all 1's or all 0's
			memset(mask, m_bool ? 0xff : 0x00, len);
			
			// Set all remainder bits (if any) to zero
			zero_remainder();
		};
		
		bitmask(const bitmask &m_mask)
		{
			if(m_mask.mask == NULL){
				mask = NULL;
				len = 0;
				bits = 0;
				
				return;
			}
			
			len = m_mask.len;
			bits = m_mask.bits;
			
			mask = new unsigned char [len];
			
			if(mask == NULL){
				throw "bitmask: Unable to allocate memory!";
			}
			
			// Copy the memory
			memcpy(mask, m_mask.mask, len);
		};
		
		~bitmask()
		{
			if(mask){
				delete [] mask;
				mask = NULL;
			}
		};
		
		// No bounds checking for efficient access
		inline __bitelem operator[](const unsigned int &m_index)
		{
			const unsigned int core = m_index/BITS_PER_BYTE;
			const unsigned char remainder = m_index%BITS_PER_BYTE;
			
			return __bitelem(mask + core, remainder);
		};
		
		inline bool operator[](const unsigned int &m_index) const
		{
			const unsigned int core = m_index/BITS_PER_BYTE;
			const unsigned char remainder = m_index%BITS_PER_BYTE;
			
			return ( (mask[core] >> remainder) & 1 );
		};
		
		// Check the bounds
		inline __bitelem at(const unsigned int &m_index)
		{
			if(m_index >= bits){
				throw "bitmask:at: index out of bounds!";
			}
			
			const unsigned int core = m_index/BITS_PER_BYTE;
			const unsigned char remainder = m_index%BITS_PER_BYTE;
			
			return __bitelem(mask + core, remainder);
		};
		
		inline bool at(const unsigned int &m_index) const
		{
			if(m_index >= bits){
				throw "bitmask:at: index out of bounds!";
			}
			
			const unsigned int core = m_index/BITS_PER_BYTE;
			const unsigned char remainder = m_index%BITS_PER_BYTE;
			
			return ( (mask[core] >> remainder) & 1 );
		};
		
		// Reset -> 0
		inline void reset()
		{
			memset(mask, 0x00, len);
			
			// We don't need to zero the remainder since all bits
			// are already 0 .
			// zero_remainder();
		};
		
		// Overwrite this bitmask with the contents of the m_set. Note that 
		// m_set is assumed to be a zero based set! Sets whose elements are 
		// not clustered around zero may result in a bit set that can be very large
		// (i.e inefficient)!
		void reset(const unsigned int &m_num_bits, const std::set<unsigned int> &m_set);
		
		inline void clear()
		{
			if(mask){
				delete [] mask;
				mask = NULL;
			}
			
			len = 0;
			bits = 0;
		};
		
		inline size_t size() const
		{
			return bits;
		};
		
		// Similar to std::includes()
		bool includes(const bitmask &m_set) const;
		
		bool operator==(const bitmask &m_mask) const;
		bool operator==(const std::set<unsigned int> &m_mask) const;
		bool operator!=(const bitmask &m_mask) const;
		bool operator<(const bitmask &m_mask) const;
		bool operator>(const bitmask &m_mask) const;
		
		// Set difference
		bitmask operator-(const bitmask &m_mask) const;
		
		// Set union
		bitmask operator+(const bitmask &m_mask) const;
		
		// Make (*this) the union of (*this) and m_set
		bitmask& operator += (const bitmask &m_set);
		
		// Set intersection
		bitmask operator&(const bitmask &m_mask) const;
		
		// Make (*this) the intersection of (*this) and m_set
		bitmask& operator &= (const bitmask &m_set);
		
		// True if all bits are zero, false otherwise
		inline bool empty() const
		{
			unsigned char *ptr = mask;
			unsigned char offset = 0;
			
			for(unsigned int i = 0;i < bits;i++){
				
				// If this bit is 1 (i.e. true) return false
				// right away.
				if( ( (*ptr >> offset) & 1) == 1 ){
					return false;
				}
				
				offset++;
				
				if(offset == BITS_PER_BYTE){
					ptr++;
					offset = 0;
				}
			}
			
			return true;
		};
		
		// Count the number of bits == 1
		inline unsigned int count() const
		{
			unsigned int tmp = 0;
			unsigned char *ptr = mask;
			unsigned char offset = 0;
			
			for(unsigned int i = 0;i < bits;i++){
				
				tmp += ( ( (*ptr >> offset) & 1) == 1 );
				
				offset++;
				
				if(offset == BITS_PER_BYTE){
					ptr++;
					offset = 0;
				}
			}
			
			return tmp;
		};
		
		inline void resize(const unsigned int &m_num_bits, const bool &m_bool)
		{
			// Deallocate existing memory
			if(mask){
				delete [] mask;
			}
			
			bits = m_num_bits;
			
			if(bits == 0){
				mask = NULL;
				len = 0;
				return;
			}
			
			len = (bits/BITS_PER_BYTE) + 
				(bits%BITS_PER_BYTE != 0);
			
			mask = new unsigned char [len];
			
			if(mask == NULL){
				throw "bitmask: Unable to allocate memory!";
			}
			
			// Flip bits to either all 1's or all 0's
			memset(mask, m_bool ? 0xff : 0x00, len);
			
			// Set all remainder bits (if any) to zero
			zero_remainder();
		};
		
		inline bitmask& operator=(const bitmask &m_mask)
		{
			if(mask){
				delete [] mask;
			}
			
			if(m_mask.mask == NULL){
				mask = NULL;
				len = 0;
				bits = 0;
				
				return *this;
			}
			
			len = m_mask.len;
			bits = m_mask.bits;
			
			mask = new unsigned char [len];
			
			if(mask == NULL){
				throw "bitmask: Unable to allocate memory!";
			}
			
			// Copy the memory
			memcpy(mask, m_mask.mask, len);
			
			return *this;
		};
		
		inline bitmask& operator=(const int &m_int)
		{
			// Flip bits to either all 0's or all 1's
			memset(mask, (m_int == 0) ? 0x00 : 0xff, len);
			
			// Set all remainder bits (if any) to zero
			zero_remainder();
			
			return *this;
		};
		
		inline bitmask& operator=(const bool &m_bool)
		{
			// Flip bits to either all 1's or all 0's
			memset(mask, m_bool ? 0xff : 0x00, len);
			
			// Set all remainder bits (if any) to zero
			zero_remainder();
				
			return *this;
		};
		
		bitmask& operator=(const std::string &m_mask);
		
		bitmask& operator=(const char* m_mask);
		
		inline size_t mpi_size() const
		{
			return sizeof(unsigned int) + // bits
				  sizeof(unsigned int) + // len
				  len;
		};
		
		inline unsigned char *mpi_pack(unsigned char *m_ptr) const
		{
			// The number of bits
			memcpy(m_ptr, &bits, sizeof(unsigned int) );
			m_ptr += sizeof(unsigned int);

			// The length of the bit array (in bytes)
			memcpy(m_ptr, &len, sizeof(unsigned int) );
			m_ptr += sizeof(unsigned int);
			
			// The bit array
			memcpy(m_ptr, mask, len);
			m_ptr += len;
			
			return m_ptr;
		};
		
		inline unsigned char *mpi_unpack(unsigned char *m_ptr)
		{
			// Clean up any existing memory
			if(mask){
				delete [] mask;
				mask = NULL;
			}
			
			// The number of bits
			memcpy(&bits, m_ptr, sizeof(unsigned int) );
			m_ptr += sizeof(unsigned int);

			// The length of the bit array (in bytes)
			memcpy(&len, m_ptr, sizeof(unsigned int) );
			m_ptr += sizeof(unsigned int);
			
			mask = new unsigned char [len];

			if(mask == NULL){
				throw "bitmask:mpi_unpack: Unable to allocate memory!";
			}

			// Copy the bit array
			memcpy(mask, m_ptr, len);
			m_ptr += len*sizeof(unsigned char);
			
			return m_ptr;
		};
		
		void write(std::ofstream &m_fout) const
		{
			// The number of bits
			m_fout.write( (const char*)&bits, sizeof(bits) );

			// The length of the bit array (in bytes)
			m_fout.write( (const char*)&len, sizeof(len) );
			
			// The bit array
			m_fout.write( (const char*)mask, len );
		};
		
		void read(std::ifstream &m_fin)
		{
			// Clean up any existing memory
			if(mask){
				delete [] mask;
				mask = NULL;
			}
			
			// The number of bits
			m_fin.read( (char*)&bits, sizeof(bits) );

			// The length of the bit array (in bytes)
			m_fin.read( (char*)&len, sizeof(len) );
			
			mask = new unsigned char [len];

			if(mask == NULL){
				throw "bitmask:read: Unable to allocate memory!";
			}
			
			// The bit array
			m_fin.read( (char*)mask, len );
		
		}
	
		// Convert this bitmask into a std::set<unsigned int>
		operator std::set<unsigned int>() const;
};

std::ostream& operator<<(std::ostream &s, const bitmask &m_set);

#endif // __BIT_MASK
