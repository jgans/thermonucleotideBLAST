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

#include "bitmask.h"

using namespace std;

void bitmask::reset(const unsigned int &m_num_bits, const std::set<unsigned int> &m_set)
{
	// Only delete memory if we can not rescue the current memory
	// associated with this class
	if( mask && (m_num_bits != bits) ){
		delete [] mask;
		mask = NULL;
	}
	
	if(m_num_bits == 0){	
		mask = NULL;
		len = 0;
		bits = 0;
		return;
	}
	
	if(mask == NULL){
		bits = m_num_bits;

		len = (bits/BITS_PER_BYTE) + 
			(bits%BITS_PER_BYTE != 0);

		mask = new unsigned char [len];

		if(mask == NULL){
			throw "bitmask::reset: Unable to allocate memory!";
		}
	}
	
	// Flip bits to all 0's to initialize the bitmask
	memset(mask, 0x00, len);
	
	// We don't need to zero the remainder since all bits
	// are already 0 .
	// zero_remainder();
			
	set<unsigned int>::const_iterator iter;
	
	for(iter = m_set.begin();iter != m_set.end();iter++){
		// Each set element gets its own bit. For safety, check 
		// the bounds on each assignment
		this->at(*iter) = true;
	}
}

// Similar to std::includes()
bool bitmask::includes(const bitmask &m_set) const
{
	// Make sure we are comparing bit strings of the same length
	// (not formally required, but in practice a good reality check and
	// will allow use to efficiently search both sets with a single for-
	// loop).
	if(m_set.bits != bits){
		throw "bitmask::includes: bit length mis-match";
	}
	
	unsigned char *ptr_super = mask;
	unsigned char *ptr_sub = m_set.mask;
	unsigned char offset = 0;
	
	// Assume that m_set is a subset of (*this) until proven otherwise  
	for(unsigned int i = 0;i < bits;i++){
		
		// Every 1 bit in the subset must be a 1 in the superset as well
		// (to maintain the superset/subset relationship).
		if( ( ( (*ptr_sub >> offset) & 1 ) == 1 ) && 
		    ( ( (*ptr_super >> offset) & 1 ) == 0 ) ){
			
			// m_set is *not* a subset of (*this)
			return false;
		}
		
		offset++;

		if(offset == BITS_PER_BYTE){
			ptr_super++;
			ptr_sub++;
			offset = 0;
		}
	}
	
	// If we get here, m_set is included in (*this)
	return true;
}

// Make (*this) the union of (*this) and m_set
bitmask& bitmask::operator += (const bitmask &m_set)
{
	// Make sure we are comparing bit strings of the same length
	// (not formally required, but in practice a good reality check and
	// will allow use to efficiently search both sets with a single for-
	// loop).
	if(m_set.bits != bits){
		throw "bitmask::operator |=: bit length mis-match";
	}
	
	for(unsigned int i = 0;i < len;i++){
		mask[i] = mask[i] | m_set.mask[i];
	}
	
	return (*this);
}

// Make (*this) the intersection of (*this) and m_set
bitmask& bitmask::operator &= (const bitmask &m_set)
{
	// Make sure we are comparing bit strings of the same length
	// (not formally required, but in practice a good reality check and
	// will allow use to efficiently search both sets with a single for-
	// loop).
	if(m_set.bits != bits){
		throw "bitmask::operator &=: bit length mis-match";
	}
	
	for(unsigned int i = 0;i < len;i++){
		mask[i] = mask[i] & m_set.mask[i];
	}
	
	return (*this);
}

bool bitmask::operator==(const bitmask &m_mask) const
{
	// Make sure we are comparing bit strings of the same length
	// (not formally required, but in practice a good reality check and
	// will allow use to efficiently compare both sets with memcmp);
	if(m_mask.bits != bits){
		throw "bitmask::operator==: bit length mis-match";
	}
	
	return (memcmp(mask,  m_mask.mask, len) == 0);	
}

bool bitmask::operator==(const set<unsigned int> &m_mask) const
{
	// Compare this bitstring to the std::set.
	// First, their must be the same number of elements in both
	// set representations
	if( count() != m_mask.size() ){
		return false;
	}
	
	set<unsigned int>::const_iterator iter;
	
	for(iter = m_mask.begin();iter != m_mask.end();iter++){
		if(*iter >= bits){
			// Instead of simply returning false, throw an error 
			// to let the user know that set dimensions are not equal!
			// return false;
			throw "bitmask::operator==: set element out of bounds";
		}
		
		// Is the element (*iter) present in the bitstring?
		if( (*this)[*iter] == false){
			// The element *iter was not found, these sets can
			// not be equal!
			return false;
		}
	}
	
	return true;
}

bool bitmask::operator!=(const bitmask &m_mask) const
{
	return !(*this == m_mask);
}

bool bitmask::operator<(const bitmask &m_mask) const
{
	// Make sure we are comparing bit strings of the same length
	// (not formally required, but in practice a good reality check and
	// will allow use to efficiently compare both sets with memcmp);
	if(m_mask.bits != bits){
		throw "bitmask::operator<: bit length mis-match";
	}
	
	return (memcmp(mask,  m_mask.mask, len) < 0);	
}

bool bitmask::operator>(const bitmask &m_mask) const
{
	// Make sure we are comparing bit strings of the same length
	// (not formally required, but in practice a good reality check and
	// will allow use to efficiently compare both sets with memcmp);
	if(m_mask.bits != bits){
		throw "bitmask::operator>: bit length mis-match";
	}
	
	return (memcmp(mask,  m_mask.mask, len) > 0);	
}

// Return a bitmask containing those bits found in (*this) but not m_mask
// (i.e the set difference).
bitmask bitmask::operator-(const bitmask &m_mask) const
{
	// Make sure we are comparing bit strings of the same length
	// (not formally required, but in practice a good reality check and
	// will allow use to efficiently compare both sets with memcmp);
	if(bits != m_mask.bits){
		throw "bitmask::operator-: bit length mis-match";
	}
	
	bitmask tmp(bits, false);
	
	for(unsigned int i = 0;i < bits;i++){
		tmp[i] = ( (*this)[i] == true ) && (m_mask[i] == false);
	}
	
	return tmp;
}

// Return a bitmask containing the union of bits found in (*this) and m_mask.
bitmask bitmask::operator+(const bitmask &m_mask) const
{
	// Make sure we are comparing bit strings of the same length
	// (not formally required, but in practice a good reality check and
	// will allow use to efficiently compare both sets with memcmp);
	if(bits != m_mask.bits){
		throw "bitmask::operator+: bit length mis-match";
	}
	
	bitmask tmp(bits, false);
	
	for(unsigned int i = 0;i < bits;i++){
		tmp[i] = (*this)[i] || m_mask[i];
	}
	
	return tmp;
}

// Return a bitmask containing the intersection of bits found in (*this) and m_mask.
bitmask bitmask::operator&(const bitmask &m_mask) const
{
	// Make sure we are comparing bit strings of the same length
	// (not formally required, but in practice a good reality check and
	// will allow use to efficiently compare both sets with memcmp);
	if(bits != m_mask.bits){
		throw "bitmask::operator&: bit length mis-match";
	}
	
	bitmask tmp(bits, false);
	
	for(unsigned int i = 0;i < bits;i++){
		tmp[i] = (*this)[i] && m_mask[i];
	}
	
	return tmp;
}

bitmask& bitmask::operator=(const string &m_mask)
{
	if(mask){
		delete [] mask;
	}

	bits = m_mask.size();
	
	cerr << "bits = " << bits << endl;
	
	if(bits == 0){
		mask = NULL;
		len = 0;
		return (*this);
	}

	len = (bits/BITS_PER_BYTE) + 
		(bits%BITS_PER_BYTE != 0);

	mask = new unsigned char [len];

	if(mask == NULL){
		throw "bitmask: Unable to allocate memory!";
	}
	
	// Flip bits to all 0's
	memset(mask, 0x00, len);
			
	// Set all remainder bits (if any) to zero
	zero_remainder();
	
	for(unsigned int i = 0;i < bits;i++){
		switch(m_mask[i]){
			case '1':
				(*this)[i] = true;
				break;
			case '0':
				(*this)[i] = false;
				break;
			default:
				throw "bitmask::operator=: Illegal string character!";
		};
	}
	
	return (*this);
}

bitmask& bitmask::operator=(const char* m_mask)
{
	if(mask){
		delete [] mask;
	}

	bits = strlen(m_mask);
	
	if(bits == 0){
		mask = NULL;
		len = 0;
		return (*this);
	}

	len = (bits/BITS_PER_BYTE) + 
		(bits%BITS_PER_BYTE != 0);

	mask = new unsigned char [len];

	if(mask == NULL){
		throw "bitmask: Unable to allocate memory!";
	}
	
	// Flip bits to all 0's
	memset(mask, 0x00, len);
			
	// Set all remainder bits (if any) to zero
	zero_remainder();
	
	for(unsigned int i = 0;i < bits;i++){
		switch(m_mask[i]){
			case '1':
				(*this)[i] = true;
				break;
			case '0':
				(*this)[i] = false;
				break;
			default:
				throw "bitmask::operator=: Illegal string character!";
		};
	}
	
	return (*this);
}

// Convert this bitmask into a std::set<unsigned int>
bitmask::operator set<unsigned int>() const
{
	set<unsigned int> ret;
	
	for(unsigned int i = 0;i < bits;i++){
		if( (*this)[i] == true ){
			ret.insert(i);
		}
	}
	
	return  ret;
}

ostream& operator<<(ostream &s, const bitmask &m_set)
{
	// Print the set as human readable binary. LSB first.
	const unsigned int len = m_set.size();
	
	for(unsigned int i = 0;i < len;i++){
		s << (int)(m_set[i]);
	}
	
	return s;
}
