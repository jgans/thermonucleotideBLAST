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

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <string.h>

using namespace std;

// Turn a range or segment string into a start/stop and segment list
// Return false on failure and true on success
void GeneAnnotation::range(const string &m_range)
{
	// Clear any existing segments
	seg_list.clear();
	
	const unsigned int len = m_range.size();
	unsigned int i = 0;
	
	// Skip to the start of the range entry.
	while( (i < len) && isspace(m_range[i]) ){
		i++;
	}
	
	if(i == len){
		throw "Range: Empty range string";
	}
	
	while(i < len){
		
		list<char> number;
		
		if( !isdigit(m_range[i]) ){
			
			seg_list.clear();
			throw "Range: Could not read start";
		}
		
		while( (i < len) && isdigit(m_range[i]) ){
		
			number.push_back(m_range[i] - '0');
			i++;
		}
		
		const int start = list_to_int(number);
		
		while( (i < len) && isspace(m_range[i]) ){
			i++;
		}

		if( (i == len) || (m_range[i] != '.') ){
		
			seg_list.clear();
			throw "Range: Could not read \"..\"";
		}
		
		i++;
		
		if( (i == len) || (m_range[i] != '.') ){
		
			seg_list.clear();
			throw "Range: Could not read \"..\"";
		}
		
		i++;
				
		while( (i < len) && isspace(m_range[i]) ){
			i++;
		}
		
		if(i == len){
		
			seg_list.clear();
			throw "Range: Could not read stop";
		}
		
		
		if( !isdigit(m_range[i]) ){
		
			seg_list.clear();
			throw "Range: Could not read stop";
		}
		
		while( (i < len) && isdigit(m_range[i]) ){
		
			number.push_back(m_range[i] - '0');
			i++;
		}
		
		const int stop = list_to_int(number);
				
		while( (i < len) && ( isspace(m_range[i]) || (m_range[i] == ',') ) ){
			i++;
		}
		
		seg_list.push_back( make_pair(start, stop) );
	}
	
	seg_list.sort();
	
	start( seg_list.front().first );
	stop( seg_list.back().second );
	
	// Only keep the segment list around if it stores more than 1 element
	if(seg_list.size() == 1){
		seg_list.clear();
	}
}

// Using the start/stop and segment list, return a valid range string
string GeneAnnotation::range() const
{
	stringstream ss;
	
	if(is_segmented() == true){
		list< pair<unsigned int, unsigned int> >::const_iterator iter = 
			seg_list.begin();
		
		while( iter != seg_list.end() ){
			ss << iter->first << ".." << iter->second;
			
			iter++;
			
			if( iter != seg_list.end() ){
				ss << ", ";
			}
		}
	}
	else{
		ss << start() << ".." << stop();
	}
	
	return ss.str();
}

bool DNAMol::load(const std::string &m_filename, const unsigned int &m_type, streampos &m_pos)
{
	// Return true if there are multiple genomes to read in this file
	bool ret = false;
	
	switch(m_type){
		case PTT:
			loadPTT(m_filename);
			break;
		case GBK:
			ret = loadGBK(m_filename, m_pos);
			break;
		case EMBL:
			ret = loadEMBL(m_filename, m_pos);
			break;
		default:
			throw "Unknown file format!";
	};

	return ret;
}

void DNAMol::processGeneList(const bool &m_loading /* = false */)
{	
	list<GeneAnnotation>::iterator gene_iter;

	bool annot_overlaps_origin = false;
		
	for(gene_iter = gene_list.begin();gene_iter != gene_list.end();gene_iter++){
		
		// Handle the special case of a gene that overlaps the genome start site
		if( gene_iter->handle_gene_start_overlap(seq_len) == true){
			annot_overlaps_origin = true;
		}		
	}

	// Sort genes
	gene_list.sort();


	if( m_loading && gene_list.empty() ){
	
		// No annotations were found! Is there any sequence to be read?
		// If so, then it will be classified as intergenic space
		if(seq_len != 0){
		
			GeneAnnotation inter_genic_space;
				
			inter_genic_space.start(0);
			inter_genic_space.stop(seq_len - 1);

			gene_list.push_back(inter_genic_space);			
		}
	}
	else if(m_loading == true){

		unsigned int last_stop_plus_1;
		unsigned int first_start;
		bool annot_is_entire_mol = false;

		gene_iter = gene_list.begin();

		// Assume a linear genome
		last_stop_plus_1 = 0;
		first_start = 0;
		
		// Only add intergenic space if there is room -- i.e. no
		// molecule spanning single annotation!
		if(annot_is_entire_mol == false){
		
			// Add intergenic spaces to the gene list.
			for(/* Already initialized */;gene_iter != gene_list.end();gene_iter++){
				
				// Only consider annotations that do not overlap the origin
				// (i.e. only consider cases where start <= stop).
				if( ( gene_iter->start() <= gene_iter->stop() ) && 
				     (gene_iter->start() > last_stop_plus_1) ){
					
					GeneAnnotation inter_genic_space;
					
					inter_genic_space.start(last_stop_plus_1);
					inter_genic_space.stop(gene_iter->start() - 1);
					
					gene_list.insert(gene_iter, inter_genic_space);
				}
				
				last_stop_plus_1 = max(last_stop_plus_1, gene_iter->stop() + 1);
			}
			
			// Do we need to add intergenic space to the tail?
			last_stop_plus_1 = seq_len;

			first_start = seq_len - 1;
			
			// Only add intergenic space if no annotation overlaps the origin.
			// Also, when last_stop_plus_1 == 0 and first_start == seq_len - 1, there is
			// no intergenic space (these points are right next to each other on the
			// circle).
			if( !( (last_stop_plus_1 == 0) && (first_start == seq_len - 1) ) &&
				(last_stop_plus_1 - 1 != first_start) && (annot_overlaps_origin == false)){
				GeneAnnotation inter_genic_space;
				
				if(last_stop_plus_1 <= first_start){
					inter_genic_space.start(last_stop_plus_1);
					inter_genic_space.stop(first_start);
				}
				else{
					inter_genic_space.start(last_stop_plus_1);
					inter_genic_space.stop(seq_len + first_start);
				}
								
				gene_list.push_back(inter_genic_space);
			}
		}
	}
}

bool DNAMol::import(const DNAMol &m_mol, const unsigned int &m_type)
{
	// Import all annotations from m_mol that are within the sequence bounds of 
	// mol.
	
	// Delete all intergenic space prior to the addition of new annotation
	list<GeneAnnotation>::iterator iter;

	list<list<GeneAnnotation>::iterator> reaper;
	
	for(iter = gene_list.begin();iter != gene_list.end();iter++){
		if(iter->is_intergenic() == true){
			reaper.push_back(iter);
		}
	}
	
	while(reaper.empty() == false){
		gene_list.erase( reaper.back() );
		reaper.pop_back();
	}
	
	list<GeneAnnotation>::const_iterator m_iter;
	
	for(m_iter = m_mol.gene_list.begin();m_iter != m_mol.gene_list.end();m_iter++){
		// If this annotation *starts* before the sequence end, then we will include
		// it (to allow for annotations that overlap the oriC).
		if(m_iter->start() < seq_len){
		
			gene_list.push_back(*m_iter);
			
			if(m_iter->stop() > seq_len){
			
				// This annotation needs to wrap around the oriC
				gene_list.back().stop(m_iter->stop() - seq_len);
			}
		}
	}
	
	processGeneList(true /* recomputing intergenic space */);
	
	return true;
}

int file_type(const std::string &m_filename)
{
	// Open the file and read until we can determine the file type
	ifstream fin(m_filename.c_str(), ios::binary);

	int ret = -1; // return -1 on error

	if(!fin){
		return ret;
	}
	
	const unsigned int buffer_size = 512;
	char buffer[buffer_size];

	// Clear the array of bytes
	memset(buffer, 0, buffer_size);

	fin.read(buffer, buffer_size);

	if(fin){
		fin.close();
	}

	bool GBK_hint = false;
	unsigned int ptt_hint = 0;
	
	int first_non_space_char = -1;

	for(unsigned int i = 0;i < buffer_size;i++){
	
		if( isspace(buffer[i]) ){
			continue;
		}

		if(first_non_space_char < 0){
			first_non_space_char = i;
		}

		if(buffer[i] == '>'){
		
			ret = DNAMol::FASTA;
			break;
		}
		
		if(isupper(buffer[i])){
		
			// Genbank flat files typically start with upper case string,
			// i.e. "LOCUS".
			GBK_hint = true;
			
			// This is only a hint, so don't break -- keep reading!
		}
		
		if(buffer[i] == '.'){
			if( (i != 0) && (buffer[i - 1] == '.') ){
				ptt_hint++;
			}
		}
	}
	
	// No valid characters detected! This is an error!
	if(first_non_space_char == -1){
		return ret;
	}
	
	if(buffer[first_non_space_char] == '@'){
		return DNAMol::FASTQ;
	}
	
	// Quick test for GBK file
	if( strstr(buffer, "LOCUS") && strstr(buffer, "DEFINITION") ){
		return DNAMol::GBK;
	}
	
	// Test for GFF3 file
	if( strstr(buffer, "##") && strstr(buffer, "gff-version") ){
		return DNAMol::GFF3;
	}
	
	// Use the arbitrary threshold of two ".." occurances
	// to indicate a PTT file.
	if(ptt_hint >= 2){
		return DNAMol::PTT;
	}
	
	if((ret == -1) && GBK_hint){

		// Does this file start with the string "ID"?
		if(first_non_space_char < (int)buffer_size - 1){
			if( ( (buffer[first_non_space_char] == 'I') && 
				(buffer[first_non_space_char + 1] == 'D') ) || 
				strstr(buffer, "\nFT") ){

				return DNAMol::EMBL;
			}
		}

		ret = DNAMol::GBK;
	}

	return ret;
}
