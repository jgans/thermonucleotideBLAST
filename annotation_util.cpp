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

#include "annotation.h"

#include <sstream>
#include <string.h>

using namespace std;

unsigned long int line_number = 1;

// Note that genbank entries are 1's based (not zero based like asn.1
// entires).
bool read_range(ifstream &m_fin, pair<unsigned int, unsigned int> &m_range, 
				list< pair<unsigned int, unsigned int> > &m_seg_list)
{
	unsigned int i = 0;
	string buffer;
	
	getline(m_fin, buffer);
	line_number ++;
	
	unsigned int len = buffer.size();

	// Make the range parsing robust to Mac/PC?Unix EOL differences
	if(buffer[len - 1] == '\r'){
		len--;
	}
	
	// Skip to the start of the range entry.
	while(isspace(buffer[i])){
		i++;
	}
	
	// For now, ignore "<" and ">"
	if((buffer[i] == '<') || (buffer[i] == '>')){
		i++;
	}

	// Is this a basic range (i.e. xxx..yyy) ?
	if( isdigit(buffer[i]) ){
	
		// Read a basic range: The first entry ...
		list<char> number;

		while( isdigit(buffer[i]) ){
		
			number.push_back(buffer[i] - '0');
			i++;
		}
	
		// Convert this list to an zero based base location
		m_range.first = list_to_int(number) - 1;

		if(i == len){

			// We have found a single number entry for the
			// range (i.e. NC_005816.gbk)
			m_range.second = m_range.first;

			// Not a complement
			return false;
		}

		// The second entry. Read the spacers ".."
		while((i != len) && !isdigit(buffer[i])){
			i++;
		}

		while((i != len) && isdigit(buffer[i])){
			number.push_back(buffer[i] - '0');
			i++;
		}
		
		// Convert this list to an zero based base location
		m_range.second = list_to_int(number) - 1;

		// Not a complement
		return false;
	}

	// Is this range a simple complement?
	if(buffer[i] == 'c'){
	
		unsigned int j = i + 11; // + strlen( "complement(" )
		
		// For now, ignore "<" and ">"
		if((buffer[j] == '<') || (buffer[j] == '>')){
			j++;
		}

		// Look -- not checking the bounds on buffer for speed reasons
		if( isdigit(buffer[j]) ){
		
			// Read a basic range: The first entry ...
			list<char> number;

			while(isdigit(buffer[j])){
				number.push_back(buffer[j] - '0');
				j++;
			}
			
			// Convert this list to an zero based base location
			m_range.first = list_to_int(number) - 1;

			if(j == len - 1 /* skip the closing ')' of complement(...) */){

				// We have found a single number entry for the
				// range (i.e. NC_005816.gbk)
				m_range.second = m_range.first;

				// Is a complement
				return true;
			}

			// The second entry. Read the spacers ".."
			while( !isdigit(buffer[j]) ){
				j++;
			}

			while( isdigit(buffer[j]) ){
			
				number.push_back(buffer[j] - '0');
				j++;
			}
			
			// Convert this list to an zero based base location
			m_range.second = list_to_int(number) - 1;

			// Is a complement
			return true;
		}
	}
	
	// We are now into the complicated ranges! i.e. join and complement(join( .. ) )
	// First, make sure that we have the entire range!
	unsigned int left_paren = 0;
	unsigned int right_paren = 0;
	unsigned int j;

	for(j = i;buffer[j] != '\0';j++){
		if(buffer[j] == '('){
			left_paren ++;
		}

		if(buffer[j] == ')'){
			right_paren ++;
		}
	}

	while(left_paren != right_paren){
		string tmp;

		getline(m_fin, tmp);
		line_number ++;

		buffer += tmp;
		
		left_paren = right_paren = 0;

		for(unsigned int j = i;buffer[j] != '\0';j++){
			if(buffer[j] == '('){
				left_paren ++;
			}

			if(buffer[j] == ')'){
				right_paren ++;
			}
		}
	}

	bool is_complement = false;
	bool is_join = false;

	// Now test for complement
	if(buffer[i] == 'c'){
		is_complement = true;

		i += 11; // + strlen( "complement(" )
	}

	// Read the join info
	if(buffer[i] == 'j'){
		is_join = true;

		i += 5; // + strlen( "join(" )
	}
	else{
		if(buffer[i] == 'o'){
			is_join = true;

			i += 6; // + strlen( "order(" )
		}
		else{
			throw error_msg(":read_range: Did not find expected \"join\"");
		}
	}
	
	// Test for "complement" *again* to catch the extreamly rare "join(complement" -- see 
	// NC_005213.gbk for an example.
	if(buffer[i] == 'c'){
		is_complement = true;

		i += 11; // + strlen( "complement(" )
	}
	
	// For now, ignore "<" and ">"
	if((buffer[i] == '<') || (buffer[i] == '>')){
		i++;
	}

	j = buffer.size();

	while(i < j){
		// Read a range: The first entry ...
		list<char> number;
		
		pair<unsigned int, unsigned int> tmp_range;

		while((i < j) && isdigit(buffer[i])){
			number.push_back(buffer[i] - '0');
			i++;
		}
		
		// Convert this list to an zero based base location
		tmp_range.first = list_to_int(number) - 1;
		
		// Is there a second value in this range, or is it a single
		// base entry?
		bool is_single_base = false;

		while((i < j) && !isdigit(buffer[i])){
			if((buffer[i] == ',') || (buffer[i] == ')')){
				is_single_base = true;
			}

			i++;
		}

		// Is there a second entry to read?
		if(is_single_base){
			tmp_range.second = tmp_range.first;
		}
		else{
			// The second entry
			while((i < j) && isdigit(buffer[i])){
				number.push_back(buffer[i] - '0');
				i++;
			}
			
			// Convert this list to an zero based base location
			tmp_range.second = list_to_int(number) - 1;
		}

		m_seg_list.push_back(tmp_range);

		// Skip all non-digit characters
		while((i < j) && !isdigit(buffer[i])){
			i++;
		}
	}

	if(m_seg_list.empty()){
		throw error_msg("read_range: Unable to read join(...)");
	}

	// Sort the seg list (important for annotations that span the origin)
	m_seg_list.sort();

	m_range.first = m_seg_list.front().first;
	m_range.second = m_seg_list.back().second;

	return is_complement;

	
}

// Convert the input list of char's into an unsinged integer.
// The list is exhausted in the process
unsigned int list_to_int(list<char> &m_number)
{
#ifdef _DEBUG
	if(m_number.empty()){
		throw error_msg("Empty number list!");
	}
#endif // _DEBUG

	unsigned int num = m_number.back();
			
	m_number.pop_back();

	int power = 10;

	while(m_number.empty() == false){
		num += power*m_number.back();

		power *= 10;
		m_number.pop_back();
	}
			
	return num;
}

char *error_msg(const char *m_error)
{
	static char error[1024];
	
	stringstream ssout;
	
	ssout << m_error << " : line #" << line_number;
	
	strcpy( error, ssout.str().c_str() );
	
	return error;
}

