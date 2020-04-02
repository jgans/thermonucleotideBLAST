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

#include "sequence_data.h"

// Needed for sort()
#include <algorithm>

using namespace std;

#ifdef USE_NCBI
void sequence_data::load_asn(const std::string &m_filename, const bool &m_sort_by_len)
{
	format = ASN_1;
	
	// Clear any existing annotation data
	mol.clear();
	streampos pos;
	
	while(true){
		
		mol.push_back( DNAMol() );
		
		if(mol.back().load(m_filename, DNAMol::ASN_1, pos) == false){
			break;
		}
	}
	
	if(mol.back().empty() == true){
		mol.pop_back();
	}
		
	list<DNAMol>::const_iterator iter = mol.begin();
	const size_t num_seq = size();

	seq_length.resize(num_seq);

	for(size_t i = 0;i < num_seq;i++, iter++){

		// Don't use readdb_get_sequence_length -- it's too slow on large databases
		const unsigned int seq_len = iter->num_bases();

		seq_length[i] = make_pair(seq_len, i);
	}

	if(m_sort_by_len == true){
		sort( seq_length.begin(), seq_length.end(), sequence_order() );
	}
}
#endif // USE_NCBI

void sequence_data::load_gbk(const std::string &m_filename, const bool &m_sort_by_len)
{
	format = GBK;
	
	// Clear any existing annotation data
	mol.clear();
	streampos pos;
	
	while(true){
		
		mol.push_back( DNAMol() );
		
		if(mol.back().load(m_filename, DNAMol::GBK, pos) == false){
			break;
		}
	}
	
	if(mol.back().empty() == true){
		mol.pop_back();
	}
		
	list<DNAMol>::const_iterator iter = mol.begin();
	const size_t num_seq = size();

	seq_length.resize(num_seq);

	for(size_t i = 0;i < num_seq;i++, iter++){

		// Don't use readdb_get_sequence_length -- it's too slow on large databases
		const unsigned int seq_len = iter->num_bases();

		seq_length[i] = make_pair(seq_len, i);
	}

	if(m_sort_by_len == true){
		sort( seq_length.begin(), seq_length.end(), sequence_order() );
	}
}

void sequence_data::load_embl(const std::string &m_filename, const bool &m_sort_by_len)
{
	format = EMBL;
	
	// Clear any existing annotation data
	mol.clear();
	streampos pos;
	
	while(true){
		
		mol.push_back( DNAMol() );
		
		if(mol.back().load(m_filename, DNAMol::EMBL, pos) == false){
			break;
		}
	}
	
	if(mol.back().empty() == true){
		mol.pop_back();
	}
		
	list<DNAMol>::const_iterator iter = mol.begin();
	const size_t num_seq = size();

	seq_length.resize(num_seq);

	for(size_t i = 0;i < num_seq;i++, iter++){

		// Don't use readdb_get_sequence_length -- it's too slow on large databases
		const unsigned int seq_len = iter->num_bases();

		seq_length[i] = make_pair(seq_len, i);
	}

	if(m_sort_by_len == true){
		sort( seq_length.begin(), seq_length.end(), sequence_order() );
	}
}

void sequence_data::load_gff3(const std::string &m_filename, const bool &m_sort_by_len)
{
	format = GFF3;
	
	// Clear any existing annotation data
	mol.clear();
	
	GFF3File fin(m_filename);
			
	if(!fin){
		throw "Unable to parse GFF3File";
	}

	vector<string> source = fin.feature_source();			
	
	for(vector<string>::const_iterator i = source.begin();i != source.end();i++){

		mol.push_back( DNAMol() );
		
		mol.back().load(fin, *i);
	}
		
	list<DNAMol>::const_iterator iter = mol.begin();
	const size_t num_seq = size();

	seq_length.resize(num_seq);

	for(size_t i = 0;i < num_seq;i++, iter++){

		// Don't use readdb_get_sequence_length -- it's too slow on large databases
		const unsigned int seq_len = iter->num_bases();

		seq_length[i] = make_pair(seq_len, i);
	}

	if(m_sort_by_len == true){
		sort( seq_length.begin(), seq_length.end(), sequence_order() );
	}
}

void sequence_data::load_ptt(const std::string &m_filename, const bool &m_sort_by_len)
{
	format = PTT;
	
	// Clear any existing annotation data
	mol.clear();
	streampos pos;
	
	while(true){
		
		mol.push_back( DNAMol() );
		
		if(mol.back().load(m_filename, DNAMol::PTT, pos) == false){
			break;
		}
	}
	
	if(mol.back().empty() == true){
		mol.pop_back();
	}
		
	list<DNAMol>::const_iterator iter = mol.begin();
	const size_t num_seq = size();

	seq_length.resize(num_seq);

	for(size_t i = 0;i < num_seq;i++, iter++){

		// Don't use readdb_get_sequence_length -- it's too slow on large databases
		const unsigned int seq_len = iter->num_bases();

		seq_length[i] = make_pair(seq_len, i);
	}

	if(m_sort_by_len == true){
		sort( seq_length.begin(), seq_length.end(), sequence_order() );
	}
}
