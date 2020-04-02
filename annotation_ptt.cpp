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

#include <ctype.h>

using namespace std;


// Local functions
bool parse_ptt(GeneAnnotation &m_gene, const string &m_line);
bool is_annotation(const string &m_line);
SeqIdPtr write_gi_PTT(SeqIdPtr &m_sip, const unsigned int &m_gi);

bool DNAMol::loadPTT(const string &m_filename)
{
	// Read in file as ASCII
	ifstream fin( m_filename.c_str() );
	
	if(!fin){
		throw error_msg("Unable to open Protein Translation Table (ptt)");
	}
	
	string line;
	GeneAnnotation gene;
	
	// Skip the header
	while( getline(fin, line) ){
		// Annotation lines start with \d+\s*..\s*\d+
		if(is_annotation(line) == true){
		
			// Parse this line
			if(parse_ptt(gene, line) == true){
				gene_list.push_back(gene);
			}
			
			break;
		}
	}
	
	while( getline(fin, line) ){

		// Parse this line
		if(parse_ptt(gene, line) == true){
			gene_list.push_back(gene);
		}
	}
	
	fin.close();
	
	processGeneList(false /* Don't add intergenic space to PTT records */);
	
	return true;
}

bool is_annotation(const string &m_line)
{
	const unsigned int len = m_line.size();
	
	unsigned int pos = 0;
	
	// Skip any spaces
	while( (pos < len) && isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Is this a blank line?
	if(pos > len){
		return false;
	}
	
	unsigned int count = 0;
	
	while( (pos < len) && isdigit(m_line[pos]) ){
		pos ++;
		count ++;
	}
	
	// Did we find a least one digit?
	if(count == 0){
		return false;
	}
	
	// Skip any spaces
	while( (pos < len) && isspace(m_line[pos]) ){
		pos ++;
	}
	
	count = 0;
	
	while( (pos < len) && (m_line[pos] == '.') ){
		pos ++;
		count ++;
	}
	
	// Did we find two dots (i.e. "..")?
	if(count != 2){
		return false;
	}
	
	// Skip any spaces
	while( (pos < len) && isspace(m_line[pos]) ){
		pos ++;
	}
	
	count = 0;
	
	while( (pos < len) && isdigit(m_line[pos]) ){
		pos ++;
		count ++;
	}
	
	// Did we find a least one digit?
	if(count == 0){
		return false;
	}
	
	// We found a valid annotation line!
	return true;
}

bool parse_ptt(GeneAnnotation &m_gene, const string &m_line)
{
	// Clear any existing annotations
	m_gene.clear();
	
	SeqIdPtr sip = NULL;
	
	const unsigned int len = m_line.size();
	
	unsigned int pos = 0;
	unsigned int start = 0;
	
	// The annotation type is user until we find out otherwise
	unsigned int annot_type = GeneAnnotation::USER;
	
	// Skip any spaces
	while( (pos < len) && isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Is this a blank line?
	if(pos >= len){
		return false;
	}
	
	unsigned int count = 0;
	
	start = pos;
	
	while( (pos < len) && isdigit(m_line[pos]) ){
		pos ++;
		count ++;
	}
	
	// Did we find a least one digit?
	if(count == 0){
		throw "PTT: Could not read the start";
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record";
	}
	
	//////////////////////////////////////////////////////////////////
	// PTT files are 1-based, while our internal representation is
	// 0-based.
	m_gene.start(atoi( m_line.substr(start, pos - start).c_str() ) - 1);

	// Skip any spaces
	while( (pos < len) && isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record after start";
	}
	
	count = 0;
	
	while( (pos < len) && (m_line[pos] == '.') ){
		pos ++;
		count ++;
	}
	
	// Did we find two dots (i.e. "..")?
	if(count != 2){
		throw "PTT: could not read ..";
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record after ..";
	}
	
	// Skip any spaces
	while( (pos < len) && isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record before stop";
	}
	
	count = 0;
	
	start = pos;
	
	while( (pos < len) && isdigit(m_line[pos]) ){
		pos ++;
		count ++;
	}
	
	// Did we find a least one digit?
	if(count == 0){
		throw "PTT: Could not read the stop";
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record no stop";
	}
	
	//////////////////////////////////////////////////////////////////
	// PTT files are 1-based, while our internal representation is
	// 0-based.
	m_gene.stop(atoi( m_line.substr(start, pos - start).c_str() ) - 1);
	
	// Make sure that start <= stop
	if( m_gene.start() > m_gene.stop() ){
		const unsigned int swap_value = m_gene.start();
		
		// Swap the start and stop values
		m_gene.start( m_gene.stop() );
		m_gene.stop(swap_value);
	}
	
	// Skip any spaces
	while( (pos < len) && isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record after stop";
	}
	
	start = pos;
	
	// Read the strand
	while( (pos < len) && !isspace(m_line[pos]) ){
		pos ++;
	}
	
	if(m_line[start] == '+'){
		m_gene.is_complement(false);
	}
	else{
		if(m_line[start] == '-'){
			m_gene.is_complement(true);
		}
		else{
			// We did not read a valid strand symbol!
			throw "PTT: Could not read the strand [+|-]";
		}
	}
	
	// Skip any spaces
	while( (pos < len) && isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record after strand";
	}
	
	// Read the length
	while( (pos < len) && !isspace(m_line[pos]) ){
		pos ++;
	}
		
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record after length";
	}
		
	// Length is (a) not used and (b) the number of amino acids in the protein, 
	// not the number of nucleic acids in the gene.
	
	// Skip any spaces
	while( (pos < len) && isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record before pid";
	}
	
	start = pos;
	
	// Read the PID
	while( (pos < len) && !isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record no pid";
	}
	
	if(m_line[start] != '-'){
		sip = write_gi_PTT(sip, 
			atoi( m_line.substr(start, pos - start).c_str() ) );
		
		// Promote this annotation to a CDS
		annot_type = GeneAnnotation::CDS;
	}
	
	// Skip any spaces
	while( (pos < len) && isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record after pid";
	}
	
	start = pos;
	
	// Read the Gene
	while( (pos < len) && !isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record before gene";
	}
	
	if(m_line[start] != '-'){
		m_gene.info( GeneAnnotation::LOCUS, 
			m_line.substr(start, pos - start).c_str() );
		
		// Promote this annotation to a CDS
		annot_type = GeneAnnotation::CDS;
	}
	
	// Skip any spaces
	while( (pos < len) && isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record after gene";
	}
	
	start = pos;
	
	// Read the Synonym
	while( (pos < len) && !isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record before synonym";
	}
	
	if(m_line[start] != '-'){
		m_gene.info( GeneAnnotation::LOCUS_TAG, 
			m_line.substr(start, pos - start).c_str() );
		
		// Promote this annotation to a GENE (unless
		// it is already a CDS --  CDS < GENE).
		annot_type = min(annot_type, (unsigned int)GeneAnnotation::GENE);
	}
	
	// Skip any spaces
	while( (pos < len) && isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record after synonym";
	}
	
	start = pos;
	
	// Read the COG code
	while( (pos < len) && !isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record before COD code";
	}
	
	if(m_line[start] != '-'){
		m_gene.info( GeneAnnotation::COG_CODE, 
			m_line.substr(start, pos - start).c_str() );
	}
	
	// Skip any spaces
	while( (pos < len) && isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record after COG code";
	}
	
	start = pos;
	
	// Read the COG id
	while( (pos < len) && !isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record before COD id";
	}
	
	if(m_line[start] != '-'){
		m_gene.info( GeneAnnotation::COG_ID, 
			m_line.substr(start, pos - start).c_str() );
	}
	
	// Skip any spaces
	while( (pos < len) && isspace(m_line[pos]) ){
		pos ++;
	}
	
	// Did we run out of characters?
	if(pos >= len){
		throw "PTT: Truncated record after COG id";
	}
	
	// Read the Product -- i.e. the remainder of the line
	m_gene.info(GeneAnnotation::PRODUCT, 
		m_line.substr(pos, len - pos).c_str() );
	
	// Set the annotation type
	m_gene.type(annot_type);
	
	// Set the SeqId
	if(sip){
		m_gene.seqid(sip);

		sip = SeqIdSetFree(sip);
	}
	
	return true;
}

// Write a gi to the given SeqIdPtr. If a gi entry
// already exists, this code will overwrite it! The updated SedIdPtr is
// returned and the input pointer is invalid.
SeqIdPtr write_gi_PTT(SeqIdPtr &m_sip, const unsigned int &m_gi)
{
	SeqIdPtr sip = NULL;
	SeqIdPtr tmp_new = NULL;
	SeqIdPtr tmp_old = NULL;

	sip = ValNodeNew(NULL);
	sip->choice = SEQID_GI;
	sip->data.intvalue = m_gi;

	if(m_sip == NULL){
		return sip;
	}

	tmp_new = sip;
	tmp_old = m_sip;

	while(tmp_old != NULL){
		if(tmp_old->choice != SEQID_GI){
			tmp_new->next = SeqIdDup(tmp_old);
			tmp_new = tmp_new->next;
		}

		tmp_old = tmp_old->next;
	}

	// Free the old ptr
	m_sip = SeqIdSetFree(m_sip);

	return sip;
}
