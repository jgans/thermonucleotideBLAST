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

#include <iostream>
#include <string>
#include <sstream>
#include "primer.h"

using namespace std;
using namespace ASSY_HEURISTIC;

int PCRPrimer::operator() (const CircleBuffer<unsigned char, MAX_PCR_LENGTH> &m_primer,
	const unsigned int &m_mask) const
{
	int status = PCR_VALID;
	
	if( m_primer.empty() ){
		throw "PCRPrimer::operator(): Empty primer buffer";
	}

	// Is this the last base?
	if( (m_mask & NO_3_T) && (m_primer.back() == T) ){
		if(report_verbose){
			status |= NO_3_T;
		}
		else{
			return NO_3_T;
		}
	}
	
	// Is this the first base?
	if( (m_mask & NO_5_G) && (m_primer.front() == G) ){
		if(report_verbose){
			status |= NO_5_G;
		}
		else{
			return NO_5_G;
		}
	}
	
	// Do we need to check to first and last base?
	if( (m_mask & NO_5_G_3_C) && (m_primer.front() == G) && (m_primer.back() == C)){
		if(report_verbose){
			status |= NO_5_G_3_C;
		}
		else{
			return NO_5_G_3_C;
		}
	}
	
	const unsigned int len = m_primer.size();

	if(m_mask & NO_3_POLY_G){
		
		if(len >= 3){
			
			if( (m_primer[len - 1] == G) && (m_primer[len - 2] == G) && (m_primer[len - 3] == G) ){
			
				if(report_verbose){
					status |= NO_3_POLY_G;
				}
				else{
					return NO_3_POLY_G;
				}
			}
		}
		
		if(len >= 4){
			
			if( (m_primer[len - 1] == G) && (m_primer[len - 2] == A) && (m_primer[len - 3] == G) && 
				(m_primer[len - 4] == G) ){
				
				if(report_verbose){
					status |= NO_3_POLY_G;
				}
				else{
					return NO_3_POLY_G;
				}
			}
		}
		
	}
	
	CircleBuffer<unsigned char, MAX_PCR_LENGTH>::const_iterator iter;
	
	unsigned int TC_run = 0;
	unsigned int AG_run = 0;
	unsigned int G_run = 0;
	unsigned int C_run = 0;
	unsigned int A_run = 0;
	unsigned int T_run = 0;
	unsigned int max_TC_run = 0;
	unsigned int max_AG_run = 0;
	unsigned int num_GC = 0;

	unsigned int index = 0;

	// Use a heuristic definition of the "middle of the probe"
	const unsigned int lower_middle = (unsigned int)(len*MIDDLE_LOWER_BOUND);
	const unsigned int upper_middle = (unsigned int)(len*MIDDLE_UPPER_BOUND);
	
	// The list must be packed from 5' -> 3'
	for(iter = m_primer.begin();iter != m_primer.end();iter++, index++){
		switch(*iter){
			case A:
				AG_run ++;
				A_run++;
				TC_run = 0;
				T_run = 0;
				G_run = 0;
				C_run = 0;
				break;
			case T:
				TC_run++;
				T_run ++;
				AG_run = 0;
				G_run = 0;
				C_run = 0;
				A_run = 0;
				break;
			case G:
				AG_run ++;
				G_run ++;
				TC_run = 0;
				C_run = 0;
				A_run = 0;
				T_run = 0;
				num_GC ++;
				break;
			case C:
				TC_run ++;
				C_run ++;
				AG_run = 0;
				G_run = 0;
				A_run = 0;
				T_run = 0;
				num_GC ++;
				break;
		};

		if(max_TC_run < TC_run){
			max_TC_run = TC_run;
		}

		if(max_AG_run < AG_run){
			max_AG_run = AG_run;
		}
		
		if(m_mask & NO_IDENTICAL_RUNS){
			if(A_run >= run_len){
				if(report_verbose){
					status |= NO_IDENTICAL_RUNS;
				}
				else{
					return NO_IDENTICAL_RUNS;
				}
			}
			
			if(T_run >= run_len){
				if(report_verbose){
					status |= NO_IDENTICAL_RUNS;
				}
				else{
					return NO_IDENTICAL_RUNS;
				}
			}
			
			if(G_run >= run_len){
				if(report_verbose){
					status |= NO_IDENTICAL_RUNS;
				}
				else{
					return NO_IDENTICAL_RUNS;
				}
			}
			
			if(C_run >= run_len){
				if(report_verbose){
					status |= NO_IDENTICAL_RUNS;
				}
				else{
					return NO_IDENTICAL_RUNS;
				}
			}
		}
		
		if( (index == 4) && (m_mask & MULTI_5_GC) ){
			if(max_TC_run >= 2){
				if(num_GC > 2){
					if(report_verbose){
						status |= MULTI_5_GC;
					}
					else{
						return MULTI_5_GC;
					}
				}
			}
			else{
				if(num_GC > 3){
					if(report_verbose){
						status |= MULTI_5_GC;
					}
					else{
						return MULTI_5_GC;
					}
				}
			}
		}
		
		if( (index == 1) && (m_mask & NO_5_PENULTIMATE_G) && (*iter == G) ){
			
			if(report_verbose){
				status |= NO_5_PENULTIMATE_G;
			}
			else{
				return NO_5_PENULTIMATE_G;
			}
		}
		
		// Start at index == lower_middle + 1
		if( (m_mask & NO_MIDDLE_CC) && (index > lower_middle) && (index <= upper_middle) && (C_run >= 2) ){
			
			if(report_verbose){
				status |= NO_MIDDLE_CC;
			}
			else{
				return NO_MIDDLE_CC;
			}
		}
	}

	if(m_mask & NO_POLY_RUNS){
		if(max_TC_run >= run_len){
			if(report_verbose){
				status |= NO_POLY_RUNS;
			}
			else{
				return NO_POLY_RUNS;
			}
		}

		if(max_AG_run >= run_len){
			if(report_verbose){
				status |= NO_POLY_RUNS;
			}
			else{
				return NO_POLY_RUNS;
			}
		}
	}

	if(m_mask & POLY_3_GC){
		if(C_run >= 3){
			if(report_verbose){
				status |= POLY_3_GC;
			}
			else{
				return POLY_3_GC;
			}
		}

		if(G_run >= 3){
			if(report_verbose){
				status |= POLY_3_GC;
			}
			else{
				return POLY_3_GC;
			}
		}
	}

	if(m_mask & GC_CONTENT){
		const float gc_content = (float)(num_GC)/len;

		if(gc_content < gc_min){
			if(report_verbose){
				status |= GC_CONTENT;
			}
			else{
				return GC_CONTENT;
			}
		}

		if(gc_content > gc_max){
			if(report_verbose){
				status |= GC_CONTENT;
			}
			else{
				return GC_CONTENT;
			}
		}
	}

	// If we get here, this primer is valid!
	return status;
}

int PCRPrimer::operator() (const string &m_primer, const unsigned int &m_mask) const
{
	CircleBuffer<unsigned char, MAX_PCR_LENGTH> primer_list;

	const size_t len = m_primer.size();

	for(size_t i = 0;i < len;i++){
		switch(m_primer[i]){
			case 'A':
			case 'a':
				primer_list.push_back(A);
				break;
			case 'T':
			case 't':
				primer_list.push_back(T);
				break;
			case 'G':
			case 'g':
				primer_list.push_back(G);
				break;
			case 'C':
			case 'c':
				primer_list.push_back(C);
				break;
			default:
				return BAD_BASE;
		};
	}

	// If we get here, this primer is valid!
	return (*this)(primer_list, m_mask);
}


int PCRPrimer::operator() () const
{
	return (*this)(primer_buffer, mask);
}

int PCRPrimer::operator() (const unsigned int &m_mask) const
{
	return (*this)(primer_buffer, m_mask);
}

string PCRPrimer::str() const
{
	string primer_str;

	primer_str.resize( primer_buffer.size() );

	CircleBuffer<unsigned char, MAX_PCR_LENGTH>::const_iterator iter;

	unsigned int index = 0;

	// The list must be packed from 5' -> 3'
	for(iter = primer_buffer.begin();iter != primer_buffer.end();iter++, index++){
		switch(*iter){
			case A:
				primer_str[index] = 'A';
				break;
			case T:
				primer_str[index] = 'T';
				break;
			case G:
				primer_str[index] = 'G';
				break;
			case C:
				primer_str[index] = 'C';
				break;
		};
	}
	
	return primer_str;
}


string  PCRPrimer::error(const int &m_code) const
{
	stringstream sout;
	bool delim = false;
	
	if(m_code == PCR_VALID){
		return "PCR_VALID";
	}
	
	if(m_code == BAD_BASE){
		return "BAD_BASE";
	}
	
	if(m_code & POLY_3_GC){
		sout << "POLY_3_GC";
		delim = true;
	}
		
	if(m_code & MULTI_5_GC){
		if(delim){
			sout << ", ";
		}
		
		sout << "MULTI_5_GC";
		delim = true;
	}
	
	if(m_code & NO_POLY_RUNS){
		if(delim){
			sout << ", ";
		}
		
		sout << "NO_POLY_RUNS";
		delim = true;
	}
	
	if(m_code & NO_3_T){
		if(delim){
			sout << ", ";
		}
		
		sout << "NO_3_T";
		delim = true;
	}
	
	if(m_code & GC_CONTENT){
		if(delim){
			sout << ", ";
		}
		
		sout << "GC_CONTENT";
		delim = true;
	}
	
	if(m_code & NO_5_G){
		if(delim){
			sout << ", ";
		}
		
		sout << "NO_5_G";
		delim = true;
	}
	
	if(m_code & NO_5_G_3_C){
		if(delim){
			sout << ", ";
		}
		
		sout << "NO_5_G_3_C";
		delim = true;
	}
	
	if(m_code & NO_IDENTICAL_RUNS){
		if(delim){
			sout << ", ";
		}
		
		sout << "NO_IDENTICAL_RUNS";
		delim = true;
	}
	
	if(m_code & NO_5_PENULTIMATE_G){
		if(delim){
			sout << ", ";
		}
		
		sout << "NO_5_PENULTIMATE_G";
		delim = true;
	}
	
	if(m_code & NO_3_POLY_G){
		if(delim){
			sout << ", ";
		}
		
		sout << "NO_3_POLY_G";
		delim = true;
	}
	
	if(m_code & NO_MIDDLE_CC){
		if(delim){
			sout << ", ";
		}
		
		sout << "NO_MIDDLE_CC";
		delim = true;
	}
	
	return sout.str();
}

pair<unsigned int, unsigned int> PCRPrimer::num_g_num_c(const string &m_seq) const
{
	unsigned int num_g = 0;
	unsigned int num_c = 0;
	
	string::const_iterator iter;
	
	for(iter = m_seq.begin();iter != m_seq.end();iter++){
		switch(*iter){
			case 'G':
			case 'g':
				num_g ++;
				break;
			case 'C':
			case 'c':
				num_c++;
				break;
		};
	}
	
	return make_pair(num_g, num_c);
}

// Helper function
unsigned int parse_pcr_filter_options(const vector<string> &m_opt)
{
	unsigned int ret = 0;

	vector<string>::const_iterator iter;
	
	for(iter = m_opt.begin();iter != m_opt.end();iter++){
		if(*iter == "POLY_3_GC"){
			ret |= POLY_3_GC;
			continue;
		}
		
		if(*iter == "MULTI_5_GC"){
			ret |= MULTI_5_GC;
			continue;
		}
		
		if(*iter == "NO_POLY_RUNS"){
			ret |= NO_POLY_RUNS;
			continue;
		}
		
		if(*iter == "NO_3_T"){
			ret |= NO_3_T;
			continue;
		}
		
		if(*iter == "GC_CONTENT"){
			ret |= GC_CONTENT;
			continue;
		}
		
		if(*iter == "NO_5_G"){
			ret |= NO_5_G;
			continue;
		}
		
		if(*iter == "NO_5_G_3_C"){
			ret |= NO_5_G_3_C;
			continue;
		}
		
		if(*iter == "NO_IDENTICAL_RUNS"){
			ret |= NO_IDENTICAL_RUNS;
			continue;
		}
		
		if(*iter == "NO_5_PENULTIMATE_G"){
			ret |= NO_5_PENULTIMATE_G;
			continue;
		}
		
		if(*iter == "NO_3_POLY_G"){
			ret |= NO_3_POLY_G;
			continue;
		}
		
		if(*iter == "NO_MIDDLE_CC"){
			ret |= NO_MIDDLE_CC;
			continue;
		}
		
		// If we get here, we've encountered an unknown option
		cerr << "Unknown PCR filter option \"" << *iter << "\"" << endl;
		
		throw "Unknown PCR filter option";
	}
	
	return ret;
}

