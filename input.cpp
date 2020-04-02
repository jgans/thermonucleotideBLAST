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

#include "hybrid_sig.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

size_t read_input_file(const string &m_file, vector<hybrid_sig> &m_sig_list, 
	const bool &m_ignore_probe, const bool &m_force_probe)
{
	if(m_ignore_probe && m_force_probe){
		throw "Error: Can not both ignore and force probes at the same time!";
	}
	
	// Append assays to the existing assays in m_sig_list
	size_t sig_count = m_sig_list.size();
	
	ifstream fin( m_file.c_str() );
	
	if(!fin){
		cerr << "Unable to open " << m_file << " for reading" << endl;
		throw __FILE__ ":read_input_file: I/O Error";
	}
	
	string line;
	
	while( getline(fin, line) ){
	
		vector<string> tmp_sig;
		string entry;
		
		// Trim any characters that appear after a comment symbol on a line
		const string::size_type pos = line.find(COMMENT_SYMBOL);
		
		if(pos != string::npos ){
			line = line.substr(0, pos);
		}
		
		stringstream ss(line);
		
		while(ss >> entry){
			tmp_sig.push_back(entry);
		}
		
		// Skip empty lines
		if(tmp_sig.empty() == true){
			continue;
		}
		
		switch( tmp_sig.size() ){
			case 2: // Name + probe
				{
					if(m_ignore_probe){
						throw "ignore_probe is true but only probes have been provided!";
					}

					m_sig_list.push_back( hybrid_sig(tmp_sig[0], tmp_sig[1], sig_count++) );
				}
				break;
			case 3: // Name + forward + reverse
				{
					if(m_force_probe){
					
						// When we turn a pair of oligos into probes, append "_F" and "_R" to
						// make the names unique
						m_sig_list.push_back( hybrid_sig(tmp_sig[0] + "_F", tmp_sig[1], sig_count++) );
						m_sig_list.push_back( hybrid_sig(tmp_sig[0] + "_R", tmp_sig[2], sig_count++) );
					}
					else{
					
						m_sig_list.push_back( hybrid_sig(tmp_sig[0], tmp_sig[1], tmp_sig[2], sig_count++) );
					}
				}
				break;
			case 4: // Name + forward + reverse + probe
				{
					if(m_ignore_probe){
					
						// As the user has requested, do not include the probe sequence
						m_sig_list.push_back( hybrid_sig(tmp_sig[0], tmp_sig[1], tmp_sig[2], sig_count++) );
					}
					else{
						if(m_force_probe){
						
							// When we turn a pair of oligos into probes, append "_F", "_R" and "_P" to
							// make the names unique
							m_sig_list.push_back( hybrid_sig(tmp_sig[0] + "_F", tmp_sig[1], 
								sig_count++) );
							m_sig_list.push_back( hybrid_sig(tmp_sig[0] + "_R", tmp_sig[2],
								sig_count++) );
							m_sig_list.push_back( hybrid_sig(tmp_sig[0] + "_P", tmp_sig[3], 
								sig_count++) );
						}
						else{
						
							m_sig_list.push_back( hybrid_sig(tmp_sig[0], tmp_sig[1], tmp_sig[2], tmp_sig[3], 
								sig_count++) );
						}
					}
				}
				break;
			default:
				cerr << "Error reading " << m_file << endl;
				throw "Error: Invalid number of arguments input file";
		};
	}
	
	fin.close();
	
	return sig_count;
}
