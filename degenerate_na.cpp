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

#include "degenerate_na.h"

#include <ctype.h>
#include <vector>
#include <iostream>

using namespace std;

list< pair<string, string> > 
	expand_nucleic_acid(const list< pair<string, string> > &m_seq)
{
	list< pair<string, string> > seq_list;
	
	list< pair<string, string> >::const_iterator iter;
	
	for(iter = m_seq.begin();iter != m_seq.end();iter++){
	
		list<string> primer1 = expand_nucleic_acid(iter->first);
		list<string> primer2 = expand_nucleic_acid(iter->second);
		
		list<string>::const_iterator iter1;
		list<string>::const_iterator iter2;
				
		for(iter1 = primer1.begin();iter1 != primer1.end();iter1++){
			for(iter2 = primer2.begin();iter2 != primer2.end();iter2++){
				seq_list.push_back( make_pair(*iter1, *iter2) );
			}
		}
		
	}
		
	return seq_list;
}

list<string> expand_nucleic_acid(const string &m_seq)
{
	const unsigned int len = (unsigned int)( m_seq.size() );
	 
	list<string> seq_list;
	
	vector<string> degen_vec( len );
	unsigned int i;

	for(i = 0;i < len;i++){
	
		switch( toupper(m_seq[i]) ){
			case 'A':
				degen_vec[i] = "A";
				break;
			case 'T':
				degen_vec[i] = "T";
				break;
			case 'G':
				degen_vec[i] = "G";
				break;
			case 'C':
				degen_vec[i] = "C";
				break;
			case 'I': // Inosine
				degen_vec[i] = "I";
				break;
			case 'M':
				degen_vec[i] = "AC";
				break;
			case 'R':
				degen_vec[i] = "GA";
				break;
			case 'S':
				degen_vec[i] = "GC";
				break;
			case 'V':
				degen_vec[i] = "GCA";
				break;
			case 'W':
				degen_vec[i] = "AT";
				break;
			case 'Y':
				degen_vec[i] = "TC";
				break;
			case 'H':
				degen_vec[i] = "ACT";
				break;
			case 'K':
				degen_vec[i] = "GT";
				break;
			case 'D':
				degen_vec[i] = "GAT";
				break;
			case 'B':
				degen_vec[i] = "GTC";
				break;
			case 'N':
				degen_vec[i] = "ATGC";
				break;
			default:
			
				cerr << "Found the unknown letter \"" << m_seq[i] << "\"" << endl;
				throw "Unknown base!";
		};
	}
	
	size_t num_seq = 1;
	
	vector<string::const_iterator> iter_vec(len);
	
	for(i = 0;i < len;i++){

		iter_vec[i] = degen_vec[i].begin();
		
		num_seq *= degen_vec[i].size();
	}
	
	for(i = 0;i < num_seq;i++){

		string tmp(len, 'A');
		unsigned int j;

		for(j = 0;j < len;j++){
			tmp[j] = *iter_vec[j];
		}
		
		seq_list.push_back(tmp);
		
		for(j = 0;j < len;j++){
			iter_vec[j]++;
			
			if( iter_vec[j] == degen_vec[j].end() ){
				// Wrap around
				iter_vec[j] = degen_vec[j].begin();
			}
			else{
				break;
			}
		}
	}
	
	return seq_list;
}

unsigned int degeneracy(const std::string &m_oligo)
{
	unsigned int d = 1;

	for(string::const_iterator i = m_oligo.begin();i != m_oligo.end();i++){

		switch( toupper(*i) ){
			case 'A':
			case 'T':
			case 'G':
			case 'C':
			case 'I': // Inosine

				// d *= 1
				break;
			case 'M':
			case 'R':
			case 'S':
			case 'W':
			case 'Y':
			case 'K':
				d *= 2;
				break;
			case 'V':
			case 'H':
			case 'D':
			case 'B':
				d *= 3;
				break;
			case 'N':
				d *= 4;
				break;
			default:
			
				cerr << "Found the unknown letter '" << *i << "'" << endl;
				throw "Unknown base!";
		};
	}

	return d;
}
