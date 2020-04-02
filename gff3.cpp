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

#include "gff3.h"

#include <stdlib.h> // For atoi

#include <sstream>

// For debugging
#include <iostream>

using namespace std;

// These are static variables in GFF3Record
string GFF3Record::sofa_version;
map<string, int> GFF3Record::sofa;
		
// Local functions
bool ishex(const char &m_c);
unsigned int parse_type(const string &m_type);
size_t replace_escapes(string &m_str);
string extract_defline(const string &m_line);

void GFF3File::parse(ifstream &m_fin)
{
	unsigned int line_number = 0;
	
	try{
		while( getline(m_fin, line_buffer) ){

			parse_line(m_fin);

			line_number++;
		}
	}
	catch(const char* error){
		stringstream ssin;
		
		ssin << "(line #" << line_number << "): " << error;
		
		throw ssin.str();
	}
	catch(const string &error){
		stringstream ssin;
		
		ssin << "(line #" << line_number << "): " << error;
		
		throw ssin.str();
	}
	
	catch(...){
		stringstream ssin;
		
		ssin << "(line #" << line_number << "): Unhandled error";
		
		throw ssin.str();
	}
}

void GFF3File::parse_pragma(std::ifstream &m_fin)
{
	stringstream ss_pragma(line_buffer);
	
	string tmp;
	
	ss_pragma >> tmp;
	
	if(tmp == "##FASTA"){
	
		parse_fasta(m_fin);
		return;
	}
	
	if(tmp == "##gff-version"){
	
		if( !(ss_pragma >> pragma_version) ){
			throw "Unable to read gff-version";
		}
		
		return;
	}
	
	if(tmp == "##feature-ontology"){
	
		if( !(ss_pragma >> pragma_feature_ontology) ){
			throw "Unable to read feature-ontology";
		}
		
		return;
	}
	
	if(tmp == "##sequence-region"){
	
		string seq_id;
		int start;
		int stop;
		
		if( !(ss_pragma >> seq_id) ){
			throw "Unable to read sequence-region->seq_id";
		}
		
		if( !(ss_pragma >> start) ){
			throw "Unable to read sequence-region->start";
		}
		
		if( !(ss_pragma >> stop) ){
			throw "Unable to read sequence-region->stop";
		}
		
		pragma_sequence_region.insert( make_pair(seq_id, make_pair(start, stop) ) );

		return;
	}
	
	if(tmp == "##attribute-ontology"){
	
		if( !(ss_pragma >> pragma_attribute_ontology) ){
			throw "Unable to read attribute-ontology";
		}
		
		return;
	}
	
	if(tmp == "##source-ontology"){
	
		if( !(ss_pragma >> pragma_source_ontology) ){
			throw "Unable to read source-ontology";
		}
		
		return;
	}
	
	if(tmp == "##species"){
	
		if( !(ss_pragma >> pragma_species) ){
			throw "Unable to read species";
		}
		
		return;
	}
	
	if(tmp == "##genome-build"){
	
		if( !(ss_pragma >> pragma_genome_build_source) ){
			throw "Unable to read genome-build->source";
		}
		
		if( !(ss_pragma >> pragma_genome_build_name) ){
			throw "Unable to read genome-build->name";
		}
		
		return;
	}
	
	if(tmp == "###"){
		return;
	}
}

void GFF3File::parse_fasta(std::ifstream &m_fin)
{
	// Read the rest of the file as fasta
	string defline;
	string::size_type loc = line_buffer.find('>');
	
	if(loc == 0){
		defline = line_buffer;
	}
	else{
		if( !getline(m_fin, defline) ){
			throw "Unable to read fasta defline";
		}
		
		loc = defline.find('>');
	}
	
	do{
		// Is this is a valid defline?
		if( (loc != 0) || (defline.size() == 1) ){
			throw "Malformed fasta defline";
		}
		
		// Strip the '>' from the fasta defline
		defline = extract_defline(defline);
		
		// Truncate defline at the first white space character
		string &seq_ref = seq[defline];
		
		while( getline(m_fin, line_buffer) ){
			
			if(line_buffer.empty() == true){
				continue;
			}
			
			if(line_buffer[0] == '>'){
				defline = extract_defline(line_buffer);
				break;
			}
			
			// This is a line of sequence (either nucleotide or amino acid)
			seq_ref += line_buffer;
		}
	}
	while(m_fin);
}

void GFF3File::parse_line(ifstream &m_fin)
{
	size_t line_len = line_buffer.size();
	
	// Ignore blank lines
	if(line_len == 0){
		return;
	}

	// Is this a comment or a pragma line?
	if(line_buffer[0] == '#'){
		
		if(line_len == 1){
			return; // This is a comment
		}
		
		if(line_buffer[1] == '#'){
		
			// This is a pragma
			parse_pragma(m_fin);			
			return;
		}
		
		// This is a comment
		return;
	}
	
	// Backwards compatibility requires treating '>' in the first position as
	// an implicit FASTA pramga
	if(line_buffer[0] == '>'){
		parse_fasta(m_fin);
	}
	
	try{
		string seqid;
		string id;
		
		GFF3Record tmp(line_buffer, seqid, id);
		
		map<string, GFF3Record> &seqid_map = features[seqid];
		
		// According to the GFF3 specification, records can be split into multiple lines
		// (each with the same ID). Since using a multimap complicates the user interface,
		// I have opted to collapse mulitple records in to a single monolithic record (and
		// use a map)
		//seqid_map.insert( make_pair(id, tmp) );
		seqid_map[id] += tmp; 
	}
	catch(const char *error){
		throw error;
	}
	catch(const string &error){
		throw error;
	}
	catch(...){
		throw "Error parsing GFF3 file";
	}
}

void GFF3Record::parse(const string &m_line_buffer, string &m_seqid, string &m_id)
{
	init();

	/////////////////////////////////////////////////////////////////////////
	// Read the seqid
	string::size_type last_delim = 0;
	string::size_type next_delim = m_line_buffer.find(GFF3_DELIM, last_delim);
	
	if(next_delim == string::npos){
		throw "Unable to read the seqid from GFF3 file";
	}
	
	m_seqid = m_line_buffer.substr(last_delim, next_delim - last_delim);
	
	replace_escapes(m_seqid);
	
	last_delim = next_delim + 1;
	
	/////////////////////////////////////////////////////////////////////////
	// Read the source
	next_delim = m_line_buffer.find(GFF3_DELIM, last_delim);
	
	if(next_delim == string::npos){
		throw "Unable to read the source from GFF3 file";
	}
	
	rec_source = m_line_buffer.substr(last_delim, next_delim - last_delim);
	
	replace_escapes(rec_source);
	
	last_delim = next_delim + 1;
	
	/////////////////////////////////////////////////////////////////////////
	// Read the type
	next_delim = m_line_buffer.find(GFF3_DELIM, last_delim);
	
	if(next_delim == string::npos){
		throw "Unable to read the type from GFF3 file";
	}
	
	rec_type = parse_SOFA( m_line_buffer.substr(last_delim, next_delim - last_delim) );
	
	last_delim = next_delim + 1;
	
	/////////////////////////////////////////////////////////////////////////
	// Read the start
	next_delim = m_line_buffer.find(GFF3_DELIM, last_delim);
	
	if(next_delim == string::npos){
		throw "Unable to read the start from GFF3 file";
	}

	GFFSegment seg;

	// The start is 1's based
	seg.range.first = atoi( m_line_buffer.substr(last_delim, next_delim - last_delim).c_str() ) - 1;
	
	if(seg.range.first < 0){
		throw "Start is out of bounds (< 0) in GFF file";
	}
	
	last_delim = next_delim + 1;
	
	/////////////////////////////////////////////////////////////////////////
	// Read the stop
	next_delim = m_line_buffer.find(GFF3_DELIM, last_delim);
	
	if(next_delim == string::npos){
		throw "Unable to read the stop from GFF3 file";
	}
	
	// The stop is 1's based
	seg.range.second = atoi( m_line_buffer.substr(last_delim, next_delim - last_delim).c_str() ) - 1;
	
	if(seg.range.second < 0){
		throw "Stop is out of bounds (< 0) in GFF file";
	}
	
	if(seg.range.first > seg.range.second){
		throw "Bounds error: Start > Stop in GFF file";
	}
	
	last_delim = next_delim + 1;
	
	/////////////////////////////////////////////////////////////////////////
	// Read the score
	next_delim = m_line_buffer.find(GFF3_DELIM, last_delim);
	
	if(next_delim == string::npos){
		throw "Unable to read the score from GFF3 file";
	}
	
	const string feature_score = m_line_buffer.substr(last_delim, next_delim - last_delim);

	if(feature_score == "."){
		rec_score = 0.0;
	}
	else{
		rec_score = atof( m_line_buffer.substr(last_delim, next_delim - last_delim).c_str() );
	}
	
	last_delim = next_delim + 1;
	
	/////////////////////////////////////////////////////////////////////////
	// Read the strand
	next_delim = m_line_buffer.find(GFF3_DELIM, last_delim);
	
	if(next_delim == string::npos){
		throw "Unable to read the strand from GFF3 file";
	}
	
	const string feature_strand = m_line_buffer.substr(last_delim, next_delim - last_delim);
	
	if(feature_strand.size() != 1){
		throw "Error reading the strand from GFF3 file";
	}
	
	// Validate and assign the strand character
	switch(feature_strand[0]){
		case '-': // Minus strand
			rec_strand = GFF3_MINUS_STRAND;
			break;
		case '+': // Plus strand
			rec_strand = GFF3_PLUS_STRAND;
			break;
		case '.': // Not a stranded feature
			rec_strand = GFF3_NO_STRAND;
			break;
		case '?': // Strand information is not available
			rec_strand = GFF3_UKNOWN_STRAND;
			break;
		default:
			throw "Illegal strand character in GFF3 file";
	};
	
	last_delim = next_delim + 1;
	
	/////////////////////////////////////////////////////////////////////////
	// Read the phase
	next_delim = m_line_buffer.find(GFF3_DELIM, last_delim);
	
	if(next_delim == string::npos){
		throw "Unable to read the phase from GFF3 file";
	}
	
	const string feature_phase = m_line_buffer.substr(last_delim, next_delim - last_delim);

	if(feature_phase == "."){
		seg.phase = 0;
	}
	else{
		seg.phase = atoi( feature_phase.c_str() );
	}
	
	if( (seg.phase < 0) || (seg.phase > 2) ){
		throw "Phase is out of bounds in GFF3 file (phase must be one of 0, 1 or 2)";
	}
	
	last_delim = next_delim + 1;
	
	// Save the segment
	rec_seg.push_back(seg);

	/////////////////////////////////////////////////////////////////////////
	// Read the attributes. This is the last column.
	next_delim = m_line_buffer.size();
		
	// Parse the contents of the feature_attributes string and store the results in the attrib map and 
	// the m_id string
	m_id = parse_attrib( m_line_buffer.substr(last_delim, next_delim - last_delim) );
}

string GFF3Record::parse_attrib(const string &m_line_buffer)
{
	string id;
	
	// Split the input buffer on ';'
	string::size_type last_delim = 0;
	string::size_type next_delim = m_line_buffer.find(GFF3_ATTRIBUTE_DELIM, last_delim);
	
	do{
		if(next_delim == string::npos){
			
			parse_attrib_pair(m_line_buffer.substr(last_delim, m_line_buffer.size() - last_delim), id);
			last_delim = string::npos;
		}
		else{
		
			parse_attrib_pair( m_line_buffer.substr(last_delim, next_delim - last_delim), id );
			last_delim = next_delim + 1;
			next_delim = m_line_buffer.find(GFF3_ATTRIBUTE_DELIM, last_delim);
		}
	}
	while( (last_delim != string::npos) || (next_delim != string::npos) );
	
	return id;
}

// Attribute pairs are of the form: tag=value
void GFF3Record::parse_attrib_pair(const string &m_line_buffer, string &m_id)
{
	// If the input string is empty, do nothing
	if(m_line_buffer.empty() == true){
		return;
	}
	
	// Extract the tag
	string::size_type loc = m_line_buffer.find(GFF3_TAG_VALUE_DELIM);
	
	if(loc == string::npos){
		throw "Unable to find tag/value delimeter (i.e. =) in GFF3 file";
	}
	
	string tag = m_line_buffer.substr(0, loc);
	
	if(tag.empty() == true){
		throw "Empty tag in GFF3 file";
	}
	
	replace_escapes(tag);
	
	string value = m_line_buffer.substr(loc + 1, m_line_buffer.size() - loc);
	
	if(value.empty() == true){
		throw "Empty value in GFF3 file";
	}
	
	// Test for the special case of an ID tag
	if(tag == "ID"){
		
		replace_escapes(value);
		
		m_id = value;
		return;
	}
	
	// Split the value string on ','
	string::size_type last_loc = 0;
	
	do{
		loc = value.find(GFF3_VALUE_DELIM, last_loc);
		
		string local_value;
		
		if(loc == string::npos){
			local_value = value.substr(last_loc, value.size() - last_loc);
		}
		else{
			local_value = value.substr(last_loc, loc - last_loc);
		}
		
		replace_escapes(local_value);
		
		attrib.insert( make_pair( tag, local_value) );
		
		last_loc = loc + 1;
	}
	while(loc != string::npos);	
}

GFF3Record& GFF3Record::operator+=(const GFF3Record &m_rhs)
{
	// Has the left hand side (i.e. lhs += m_rhs) previously been initialized?
	if(rec_type == INVALID_TYPE){

		// No existing data -- just copy the right hand side and return
		*this = m_rhs;
		return *this;
	}
	
	// If we get here, there is existing data. Make sure that the two entries are compatible
	if(rec_source != m_rhs.rec_source){
		throw "GFF3Record: Lines with the same ID, differ in source";
	}

	if(rec_strand != m_rhs.rec_strand){
		throw "GFF3Record: Lines with the same ID, differ in strand";
	}

	if(rec_type != m_rhs.rec_type){
		throw "GFF3Record: Lines with the same ID, differ in type";
	}

	if(rec_type != m_rhs.rec_type){
		throw "GFF3Record: Lines with the same ID, differ in type";
	}

	// Merge the range and phase information
	for(list<GFFSegment>::const_iterator iter = m_rhs.rec_seg.begin();
		iter != m_rhs.rec_seg.end();iter++){

		rec_seg.push_back(*iter);
	}

	rec_seg.sort();
	rec_seg.unique();

	// Merge the attribute data. 
	multimap<string, string> new_attrib;

	// First, pool all of the map keys.
	set<string> attrib_keys;

	for(multimap<string, string>::const_iterator iter = attrib.begin();iter != attrib.end();iter++){
		attrib_keys.insert(iter->first);
	}

	for(multimap<string, string>::const_iterator iter = m_rhs.attrib.begin();
		iter != m_rhs.attrib.end();iter++){
		attrib_keys.insert(iter->first);
	}

	typedef multimap<string, string>::const_iterator I;

	// Now that the keys are pooled, make a unique list for the values belonging to each key
	for(set<string>::const_iterator iter = attrib_keys.begin();iter != attrib_keys.end();iter++){

		list<string> value;

		pair<I, I> range = attrib.equal_range(*iter);

		for(I i = range.first;i != range.second;i++){
			value.push_back(i->second);
		}
		range = m_rhs.attrib.equal_range(*iter);

		for(I i = range.first;i != range.second;i++){
			value.push_back(i->second);
		}

		value.sort();
		value.unique();

		for(list<string>::const_iterator i = value.begin();i != value.end();i++){
			new_attrib.insert( make_pair(*iter, *i) );
		}
	}

	// Replace the old attribute list with the merged attribute list
	attrib = new_attrib;

	return *this;
}

string extract_defline(const string &m_line)
{
	const size_t len = m_line.size();

	string::size_type start = m_line.find('>');

	// If we can't find the '>', assume that this is a defline that has
	// *already* been extracted
	if(start == string::npos){
		return m_line;
	}

	// Skip the '>' character
	start ++;

	// Skip any interveaning white space
	while( (start < len) && isspace(m_line[start]) ){
		start++;
	}

	string::size_type stop = start;

	// Find the end of string or first white space
	while( (stop < len) && !isspace(m_line[stop]) ){
		stop++;
	}

	return m_line.substr(start, stop - start);
}

size_t replace_escapes(string &m_str)
{
	string::size_type loc = 0;
	size_t len = m_str.size();
	
	while(true){
	
		loc = m_str.find('%', loc);
		
		if( (loc == string::npos) || (loc + 2 >= len) ){
			return len;
		}
		
		if( ishex(m_str[loc + 1]) && ishex(m_str[loc + 2]) ){
			
			// This is a valid escape sequence, replace it with the corresponding
			// charater
			string escape_code = m_str.substr(loc + 1, 2);
			
			escape_code[0] = toupper(escape_code[0]);
			escape_code[1] = toupper(escape_code[1]);
			
			string escape_str;
			
			switch(escape_code[0]){
				case '2':
					
					switch(escape_code[1]){
						case '0':
							escape_str = " ";
							break;

						case '2':
							escape_str = "\"";
							break;

						case '3':
							escape_str = "#";
							break;

						case '4':
							escape_str = "$";
							break;

						case '5':
							escape_str = "%";
							break;

						case '6':
							escape_str = "&";
							break;
						
						case 'B':
							escape_str = "+";
							break;
						
						case 'C':
							escape_str = ",";
							break;

						case 'F':
							escape_str = "/";
							break;
						default:
							throw "Unknown URL escape code";
					};

					break;
				case '3':
					
					switch(escape_code[1]){
						case 'A':
							escape_str = ":";
							break;

						case 'B':
							escape_str = ";";
							break;

						case 'C':
							escape_str = "<";
							break;
						
						case 'D':
							escape_str = "=";
							break;

						case 'E':
							escape_str = ">";
							break;

						case 'F':
							escape_str = "?";
							break;

						default:
							throw "Unknown URL escape code";
					};
					
					break;
				case '4':
					switch(escape_code[1]){
						case '0':
							escape_str = "@";
							break;
						default:
							throw "Unknown URL escape code";
					};
					
					break;
				case '5':
					switch(escape_code[1]){
						case 'B':
							escape_str = "[";
							break;

						case 'C':
							escape_str = "\\";
							break;

						case 'D':
							escape_str = "]";
							break;

						case 'E':
							escape_str = "^";
							break;

						default:
							throw "Unknown URL escape code";
					};
					
					break;
				case '6':
					switch(escape_code[1]){
						case '0':
							escape_str = "`";
							break;
						default:
							throw "Unknown URL escape code";
					};

					break;
				case '7':
					switch(escape_code[1]){
					
						case 'B':
							escape_str = "{";
							break;						

						case 'C':
							escape_str = "|";
							break;
					
						case 'D':
							escape_str = "}";
							break;

						case 'E':
							escape_str = "~";
							break;
						default:
							throw "Unknown URL escape code";
					};
					break;
				default:
					throw "Unknown URL escape code";
			};
			
			m_str.replace(loc, 3, escape_str);
			
			// Update the string length
			len = m_str.size();
		}
		
		loc++;
	}
	
	return len;
}

bool ishex(const char &m_c)
{
	switch(m_c){
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
		case 'A':
		case 'a':
		case 'B':
		case 'b':
		case 'C':
		case 'c':
		case 'D':
		case 'd':
		case 'E':
		case 'e':
		case 'F':
		case 'f':
			return true;
	};
	
	return false;
}

void GFF3Record::init()
{
	if(sofa_version.empty() == false){
		// Already initialized
		return;
	}

	sofa_version = "2.3";
	// 
	sofa["SEQUENCE_ONTOLOGY"] = sofa["SO:0000000"] = SEQUENCE_ONTOLOGY;

	// A sequence_feature with an extent greater than zero.
	sofa["REGION"] = sofa["SO:0000001"] = REGION;

	// 
	sofa["INTERIOR_CODING_EXON"] = sofa["SO:0000004"] = INTERIOR_CODING_EXON;

	// A region amplified by a PCR reaction.
	sofa["PCR_PRODUCT"] = sofa["SO:0000006"] = PCR_PRODUCT;

	// A pair of sequencing reads in which the two members of the pair are related by originating at either end of a clone insert.
	sofa["READ_PAIR"] = sofa["SO:0000007"] = READ_PAIR;

	// Any one of several small cytoplasmic RNA molecules present in the cytoplasm and sometimes nucleus of a eukaryote.
	sofa["SCRNA"] = sofa["SO:0000013"] = SCRNA;

	// A collection of match parts.
	sofa["MATCH_SET"] = sofa["SO:0000038"] = MATCH_SET;

	// A part of a match, for example an hsp from blast isa match_part.
	sofa["MATCH_PART"] = sofa["SO:0000039"] = MATCH_PART;

	// A part of a gene, that has no other route in the ontology back to region. This concept is necessary for logical inference as these parts must have the properties of region. It also allows us to associate all the parts of genes with a gene.
	sofa["GENE_PART"] = sofa["SO:0000050"] = GENE_PART;

	// A regulatory element of an operon to which activators or repressors bind thereby effecting translation of genes in that operon.
	sofa["OPERATOR"] = sofa["SO:0000057"] = OPERATOR;

	// A transposon or insertion sequence. An element that can insert in a variety of DNA sequences.
	sofa["TRANSPOSABLE_ELEMENT"] = sofa["SO:0000101"] = TRANSPOSABLE_ELEMENT;

	// A match to an EST or cDNA sequence.
	sofa["EXPRESSED_SEQUENCE_MATCH"] = sofa["SO:0000102"] = EXPRESSED_SEQUENCE_MATCH;

	// The end of the clone insert.
	sofa["CLONE_INSERT_END"] = sofa["SO:0000103"] = CLONE_INSERT_END;

	// A sequence of amino acids linked by peptide bonds which may lack appreciable tertiary structure and may not be liable to irreversible denaturation.
	sofa["POLYPEPTIDE"] = sofa["SO:0000104"] = POLYPEPTIDE;

	// A sequence_variant is a non exact copy of a sequence_feature or genome exhibiting one or more sequence_alteration.
	sofa["SEQUENCE_VARIANT_OBS"] = sofa["SO:0000109"] = SEQUENCE_VARIANT_OBS;

	// An extent of biological sequence.
	sofa["SEQUENCE_FEATURE"] = sofa["SO:0000110"] = SEQUENCE_FEATURE;

	// A short preexisting polynucleotide chain to which new deoxyribonucleotides can be added by DNA polymerase.
	sofa["PRIMER"] = sofa["SO:0000112"] = PRIMER;

	// A viral sequence which has integrated into a host genome.
	sofa["PROVIRAL_REGION"] = sofa["SO:0000113"] = PROVIRAL_REGION;

	// A methylated deoxy-cytosine.
	sofa["METHYLATED_C"] = sofa["SO:0000114"] = METHYLATED_C;

	// A primary transcript that, at least in part, encodes one or more proteins.
	sofa["PROTEIN_CODING_PRIMARY_TRANSCRIPT"] = sofa["SO:0000120"] = PROTEIN_CODING_PRIMARY_TRANSCRIPT;

	// Region in mRNA where ribosome assembles.
	sofa["RIBOSOME_ENTRY_SITE"] = sofa["SO:0000139"] = RIBOSOME_ENTRY_SITE;

	// A sequence segment located within the five prime end of an mRNA that causes premature termination of translation.
	sofa["ATTENUATOR"] = sofa["SO:0000140"] = ATTENUATOR;

	// The sequence of DNA located either at the end of the transcript that causes RNA polymerase to terminate transcription.
	sofa["TERMINATOR"] = sofa["SO:0000141"] = TERMINATOR;

	// A region of sequence which may be used to manufacture a longer assembled, sequence.
	sofa["ASSEMBLY_COMPONENT"] = sofa["SO:0000143"] = ASSEMBLY_COMPONENT;

	// A region of the transcript sequence within a gene which is not removed from the primary RNA transcript by RNA splicing.
	sofa["EXON"] = sofa["SO:0000147"] = EXON;

	// One or more contigs that have been ordered and oriented using end-read information. Contains gaps that are filled with N's.
	sofa["SUPERCONTIG"] = sofa["SO:0000148"] = SUPERCONTIG;

	// A contiguous sequence derived from sequence assembly. Has no gaps, but may contain N's from unvailable bases.
	sofa["CONTIG"] = sofa["SO:0000149"] = CONTIG;

	// A sequence obtained from a single sequencing experiment. Typically a read is produced when a base calling program interprets information from a chromatogram trace file produced from a sequencing machine.
	sofa["READ"] = sofa["SO:0000150"] = READ;

	// A piece of DNA that has been inserted in a vector so that it can be propagated in E. coli or some other organism.
	sofa["CLONE"] = sofa["SO:0000151"] = CLONE;

	// The point at which a deletion occured.
	sofa["DELETION"] = sofa["SO:0000159"] = DELETION;

	// A methylated adenine.
	sofa["METHYLATED_A"] = sofa["SO:0000161"] = METHYLATED_A;

	// The position where intron is excised.
	sofa["SPLICE_SITE"] = sofa["SO:0000162"] = SPLICE_SITE;

	// The junction between the 3 prime end of an exon and the following intron.
	sofa["FIVE_PRIME_SPLICE_SITE"] = sofa["SO:0000163"] = FIVE_PRIME_SPLICE_SITE;

	// The junction between the 3 prime end of an intron and the following exon.
	sofa["THREE_PRIME_SPLICE_SITE"] = sofa["SO:0000164"] = THREE_PRIME_SPLICE_SITE;

	// A cis-acting sequence that increases the utilization of (some) eukaryotic promoters, and can function in either orientation and in any location (upstream or downstream) relative to the promoter.
	sofa["ENHANCER"] = sofa["SO:0000165"] = ENHANCER;

	// A regulatory_region composed of the TSS(s) and binding sites for TF_complexes of the basal transcription machinery.
	sofa["PROMOTER"] = sofa["SO:0000167"] = PROMOTER;

	// A nucleotide match against a sequence from another organism.
	sofa["CROSS_GENOME_MATCH"] = sofa["SO:0000177"] = CROSS_GENOME_MATCH;

	// A group of contiguous genes transcribed as a single (polycistronic) mRNA from a single regulatory region.
	sofa["OPERON"] = sofa["SO:0000178"] = OPERON;

	// The start of the clone insert.
	sofa["CLONE_INSERT_START"] = sofa["SO:0000179"] = CLONE_INSERT_START;

	// A match against a translated sequence.
	sofa["TRANSLATED_NUCLEOTIDE_MATCH"] = sofa["SO:0000181"] = TRANSLATED_NUCLEOTIDE_MATCH;

	// A region of the gene which is not transcribed.
	sofa["NON_TRANSCRIBED_REGION"] = sofa["SO:0000183"] = NON_TRANSCRIBED_REGION;

	// A transcript that in its initial state requires modification to be functional.
	sofa["PRIMARY_TRANSCRIPT"] = sofa["SO:0000185"] = PRIMARY_TRANSCRIPT;

	// A group of characterized repeat sequences.
	sofa["REPEAT_FAMILY"] = sofa["SO:0000187"] = REPEAT_FAMILY;

	// A segment of DNA that is transcribed, but removed from within the transcript by splicing together the sequences (exons) on either side of it.
	sofa["INTRON"] = sofa["SO:0000188"] = INTRON;

	// A polymorphism detectable by the size differences in DNA fragments generated by a restriction enzyme.
	sofa["RFLP_FRAGMENT"] = sofa["SO:0000193"] = RFLP_FRAGMENT;

	// The sequence of the 5' exon that encodes for protein.
	sofa["FIVE_PRIME_EXON_CODING_REGION"] = sofa["SO:0000196"] = FIVE_PRIME_EXON_CODING_REGION;

	// The sequence of the 3' exon that encodes for protein.
	sofa["THREE_PRIME_EXON_CODING_REGION"] = sofa["SO:0000197"] = THREE_PRIME_EXON_CODING_REGION;

	// An exon that does not contain any codons.
	sofa["NONCODING_EXON"] = sofa["SO:0000198"] = NONCODING_EXON;

	// Messenger RNA sequences that are untranslated and lie five prime and three prime to sequences which are translated.
	sofa["UTR"] = sofa["SO:0000203"] = UTR;

	// A region at the 5' end of a mature transcript (preceding the initiation codon) that is not translated into a protein.
	sofa["FIVE_PRIME_UTR"] = sofa["SO:0000204"] = FIVE_PRIME_UTR;

	// A region at the 3' end of a mature transcript (following the stop codon) that is not translated into a protein.
	sofa["THREE_PRIME_UTR"] = sofa["SO:0000205"] = THREE_PRIME_UTR;

	// A transcript which has undergone the necessary modifications for its function. In eukaryotes this includes, for example, processing of introns, cleavage, base modification, and modifications to the 5' and/or the 3' ends, other than addition of bases. In bacteria functional mRNAs are usually not modified.
	sofa["PROCESSED_TRANSCRIPT"] = sofa["SO:0000233"] = PROCESSED_TRANSCRIPT;

	// Messenger RNA is the intermediate molecule between DNA and protein. It includes UTR and coding sequences. It does not contain introns.
	sofa["MRNA"] = sofa["SO:0000234"] = MRNA;

	// A region of a molecule that binds a TF complex [GO:0005667].
	sofa["TF_BINDING_SITE"] = sofa["SO:0000235"] = TF_BINDING_SITE;

	// The inframe interval between the stop codons of a reading frame which when read as sequential triplets, has the potential of encoding a sequential string of amino acids. TER(NNN)nTER.
	sofa["ORF"] = sofa["SO:0000236"] = ORF;

	// The DNA sequences extending on either side of a specific locus.
	sofa["FLANKING_REGION"] = sofa["SO:0000239"] = FLANKING_REGION;

	// RNA that comprises part of a ribosome, and that can provide both structural scaffolding and catalytic activity.
	sofa["RRNA"] = sofa["SO:0000252"] = RRNA;

	// Transfer RNA (tRNA) molecules are approximately 80 nucleotides in length. Their secondary structure includes four short double-helical elements and three loops (D, anti-codon, and T loops). Further hydrogen bonds mediate the characteristic L-shaped molecular structure. Transfer RNAs have two regions of fundamental functional importance: the anti-codon, which is responsible for specific mRNA codon recognition, and the 3' end, to which the tRNA's corresponding amino acid is attached (by aminoacyl-tRNA synthetases). Transfer RNAs cope with the degeneracy of the genetic code in two manners: having more than one tRNA (with a specific anti-codon) for a particular amino acid; and 'wobble' base-pairing, i.e. permitting non-standard base-pairing at the 3rd anti-codon position.
	sofa["TRNA"] = sofa["SO:0000253"] = TRNA;

	// A small nuclear RNA molecule involved in pre-mRNA splicing and processing.
	sofa["SNRNA"] = sofa["SO:0000274"] = SNRNA;

	// A snoRNA (small nucleolar RNA) is any one of a class of small RNAs that are associated with the eukaryotic nucleus as components of small nucleolar ribonucleoproteins. They participate in the processing or modifications of many RNAs, mostly ribosomal RNAs (rRNAs) though snoRNAs are also known to target other classes of RNA, including spliceosomal RNAs, tRNAs, and mRNAs via a stretch of sequence that is complementary to a sequence in the targeted RNA.
	sofa["SNORNA"] = sofa["SO:0000275"] = SNORNA;

	// Small, ~22-nt, RNA molecule that is the endogenous transcript of a miRNA gene. Micro RNAs are produced from precursor molecules (SO:0000647) that can form local hairpin structures, which ordinarily are processed (via the Dicer pathway) such that a single miRNA molecule accumulates from one arm of a hairpin precursor molecule. Micro RNAs may trigger the cleavage of their target molecules or act as translational repressors.
	sofa["MIRNA"] = sofa["SO:0000276"] = MIRNA;

	// A very short unit sequence of DNA (2 to 4 bp) that is repeated multiple times in tandem.
	sofa["MICROSATELLITE"] = sofa["SO:0000289"] = MICROSATELLITE;

	// The sequence is complementarily repeated on the opposite strand. It is a palindrome, and it may, or may not be hyphenated. Examples: GCTGATCAGC, or GCTGA-----TCAGC.
	sofa["INVERTED_REPEAT"] = sofa["SO:0000294"] = INVERTED_REPEAT;

	// The origin of replication; starting site for duplication of a nucleic acid molecule to give two identical copies.
	sofa["ORIGIN_OF_REPLICATION"] = sofa["SO:0000296"] = ORIGIN_OF_REPLICATION;

	// Part of the primary transcript that is clipped off during processing.
	sofa["CLIP"] = sofa["SO:0000303"] = CLIP;

	// A modified nucleotide, i.e. a nucleotide other than A, T, C. G or (in RNA) U.
	sofa["MODIFIED_BASE_SITE"] = sofa["SO:0000305"] = MODIFIED_BASE_SITE;

	// A nucleotide modified by methylation.
	sofa["METHYLATED_BASE_FEATURE"] = sofa["SO:0000306"] = METHYLATED_BASE_FEATURE;

	// Regions of a few hundred to a few thousand bases in vertebrate genomes that are relatively GC and CpG rich; they are typically unmethylated and often found near the 5' ends of genes.
	sofa["CPG_ISLAND"] = sofa["SO:0000307"] = CPG_ISLAND;

	// A repeat where the same sequence is repeated in the same direction. Example: GCTGA-----GCTGA.
	sofa["DIRECT_REPEAT"] = sofa["SO:0000314"] = DIRECT_REPEAT;

	// The first base where RNA polymerase begins to synthesize the RNA transcript.
	sofa["TRANSCRIPTION_START_SITE"] = sofa["SO:0000315"] = TRANSCRIPTION_START_SITE;

	// A contiguous sequence which begins with, and includes, a start codon and ends with, and includes, a stop codon.
	sofa["CDS"] = sofa["SO:0000316"] = CDS;

	// First codon to be translated by a ribosome.
	sofa["START_CODON"] = sofa["SO:0000318"] = START_CODON;

	// In mRNA, a set of three nucleotides that indicates the end of information for protein synthesis.
	sofa["STOP_CODON"] = sofa["SO:0000319"] = STOP_CODON;

	// A nucleotide sequence that may be used to identify a larger sequence.
	sofa["TAG"] = sofa["SO:0000324"] = TAG;

	// A short diagnostic sequence tag, serial analysis of gene expression (SAGE), that allows the quantitative and simultaneous analysis of a large number of transcripts.
	sofa["SAGE_TAG"] = sofa["SO:0000326"] = SAGE_TAG;

	// Region of sequence similarity by descent from a common ancestor.
	sofa["CONSERVED_REGION"] = sofa["SO:0000330"] = CONSERVED_REGION;

	// Short (typically a few hundred base pairs) DNA sequence that has a single occurrence in a genome and whose location and base sequence are known.
	sofa["STS"] = sofa["SO:0000331"] = STS;

	// Coding region of sequence similarity by descent from a common ancestor.
	sofa["CODING_CONSERVED_REGION"] = sofa["SO:0000332"] = CODING_CONSERVED_REGION;

	// The boundary between two exons in a processed transcript.
	sofa["EXON_JUNCTION"] = sofa["SO:0000333"] = EXON_JUNCTION;

	// Non-coding region of sequence similarity by descent from a common ancestor.
	sofa["NC_CONSERVED_REGION"] = sofa["SO:0000334"] = NC_CONSERVED_REGION;

	// A sequence that closely resembles a known functional gene, at another locus within a genome, that is non-functional as a consequence of (usually several) mutations that prevent either its transcription or translation (or both). In general, pseudogenes result from either reverse transcription of a transcript of their \"normal\" paralog (SO:0000043) (in which case the pseudogene typically lacks introns and includes a poly(A) tail) or from recombination (SO:0000044) (in which case the pseudogene is typically a tandem duplication of its \"normal\" paralog).
	sofa["PSEUDOGENE"] = sofa["SO:0000336"] = PSEUDOGENE;

	// A double stranded RNA duplex, at least 20bp long, used experimentally to inhibit gene function by RNA interference.
	sofa["RNAI_REAGENT"] = sofa["SO:0000337"] = RNAI_REAGENT;

	// Structural unit composed of a nucleic acid molecule which controls its own replication through the interaction of specific proteins at one or more origins of replication.
	sofa["CHROMOSOME"] = sofa["SO:0000340"] = CHROMOSOME;

	// A cytologically distinguishable feature of a chromosome, often made visible by staining, and usually alternating light and dark.
	sofa["CHROMOSOME_BAND"] = sofa["SO:0000341"] = CHROMOSOME_BAND;

	// A region of sequence, aligned to another sequence with some statistical significance, using an algorithm such as BLAST or SIM4.
	sofa["MATCH"] = sofa["SO:0000343"] = MATCH;

	// Region of a transcript that regulates splicing.
	sofa["SPLICE_ENHANCER"] = sofa["SO:0000344"] = SPLICE_ENHANCER;

	// Expressed Sequence Tag: The sequence of a single sequencing read from a cDNA clone or PCR product; typically a few hundred base pairs long.
	sofa["EST"] = sofa["SO:0000345"] = EST;

	// A match against a nucleotide sequence.
	sofa["NUCLEOTIDE_MATCH"] = sofa["SO:0000347"] = NUCLEOTIDE_MATCH;

	// A match against a protein sequence.
	sofa["PROTEIN_MATCH"] = sofa["SO:0000349"] = PROTEIN_MATCH;

	// A sequence of nucleotides that has been algorithmically derived from an alignment of two or more different sequences.
	sofa["ASSEMBLY"] = sofa["SO:0000353"] = ASSEMBLY;

	// A set of (usually) three nucleotide bases in a DNA or RNA sequence, which together signify a unique amino acid or the termination of translation.
	sofa["CODON"] = sofa["SO:0000360"] = CODON;

	// The junction where an insertion occurred.
	sofa["INSERTION_SITE"] = sofa["SO:0000366"] = INSERTION_SITE;

	// The junction in a genome where a transposable_element has inserted.
	sofa["TRANSPOSABLE_ELEMENT_INSERTION_SITE"] = sofa["SO:0000368"] = TRANSPOSABLE_ELEMENT_INSERTION_SITE;

	// A non-coding RNA, usually with a specific secondary structure, that acts to regulate gene expression.
	sofa["SMALL_REGULATORY_NCRNA"] = sofa["SO:0000370"] = SMALL_REGULATORY_NCRNA;

	// An RNA sequence that has catalytic activity with or without an associated ribonucleoprotein.
	sofa["ENZYMATIC_RNA"] = sofa["SO:0000372"] = ENZYMATIC_RNA;

	// An RNA with catalytic activity.
	sofa["RIBOZYME"] = sofa["SO:0000374"] = RIBOZYME;

	// 5.8S ribosomal RNA (5.8S rRNA) is a component of the large subunit of the eukaryotic ribosome. It is transcribed by RNA polymerase I as part of the 45S precursor that also contains 18S and 28S rRNA. Functionally, it is thought that 5.8S rRNA may be involved in ribosome translocation. It is also known to form covalent linkage to the p53 tumour suppressor protein. 5.8S rRNA is also found in archaea.
	sofa["RRNA_5.8S"] = sofa["SO:0000375"] = RRNA_5_8S;

	// A small catalytic RNA motif that catalyzes self-cleavage reaction. Its name comes from its secondary structure which resembles a carpenter's hammer. The hammerhead ribozyme is involved in the replication of some viroid and some satellite RNAs.
	sofa["HAMMERHEAD_RIBOZYME"] = sofa["SO:0000380"] = HAMMERHEAD_RIBOZYME;

	// The RNA molecule essential for the catalytic activity of RNase MRP, an enzymatically active ribonucleoprotein with two distinct roles in eukaryotes. In mitochondria it plays a direct role in the initiation of mitochondrial DNA replication. In the nucleus it is involved in precursor rRNA processing, where it cleaves the internal transcribed spacer 1 between 18S and 5.8S rRNAs.
	sofa["RNASE_MRP_RNA"] = sofa["SO:0000385"] = RNASE_MRP_RNA;

	// The RNA component of Ribonuclease P (RNase P), a ubiquitous endoribonuclease, found in archaea, bacteria and eukarya as well as chloroplasts and mitochondria. Its best characterised activity is the generation of mature 5 prime ends of tRNAs by cleaving the 5 prime leader elements of precursor-tRNAs. Cellular RNase Ps are ribonucleoproteins. RNA from bacterial RNase Ps retains its catalytic activity in the absence of the protein subunit, i.e. it is a ribozyme. Isolated eukaryotic and archaeal RNase P RNA has not been shown to retain its catalytic function, but is still essential for the catalytic activity of the holoenzyme. Although the archaeal and eukaryotic holoenzymes have a much greater protein content than the bacterial ones, the RNA cores from all the three lineages are homologous. Helices corresponding to P1, P2, P3, P4, and P10/11 are common to all cellular RNase P RNAs. Yet, there is considerable sequence variation, particularly among the eukaryotic RNAs.
	sofa["RNASE_P_RNA"] = sofa["SO:0000386"] = RNASE_P_RNA;

	// The RNA component of telomerase, a reverse transcriptase that synthesises telomeric DNA.
	sofa["TELOMERASE_RNA"] = sofa["SO:0000390"] = TELOMERASE_RNA;

	// U1 is a small nuclear RNA (snRNA) component of the spliceosome (involved in pre-mRNA splicing). Its 5' end forms complementary base pairs with the 5' splice junction, thus defining the 5' donor site of an intron. There are significant differences in sequence and secondary structure between metazoan and yeast U1 snRNAs, the latter being much longer (568 nucleotides as compared to 164 nucleotides in human). Nevertheless, secondary structure predictions suggest that all U1 snRNAs share a 'common core' consisting of helices I, II, the proximal region of III, and IV.
	sofa["U1_SNRNA"] = sofa["SO:0000391"] = U1_SNRNA;

	// U2 is a small nuclear RNA (snRNA) component of the spliceosome (involved in pre-mRNA splicing). Complementary binding between U2 snRNA (in an area lying towards the 5' end but 3' to hairpin I) and the branchpoint sequence (BPS) of the intron results in the bulging out of an unpaired adenine, on the BPS, which initiates a nucleophilic attack at the intronic 5' splice site, thus starting the first of two transesterification reactions that mediate splicing.
	sofa["U2_SNRNA"] = sofa["SO:0000392"] = U2_SNRNA;

	// U4 small nuclear RNA (U4 snRNA) is a component of the major U2-dependent spliceosome. It forms a duplex with U6, and with each splicing round, it is displaced from U6 (and the spliceosome) in an ATP-dependent manner, allowing U6 to refold and create the active site for splicing catalysis. A recycling process involving protein Prp24 re-anneals U4 and U6.
	sofa["U4_SNRNA"] = sofa["SO:0000393"] = U4_SNRNA;

	// An snRNA required for the splicing of the minor U12-dependent class of eukaryotic nuclear introns. It forms a base paired complex with U6atac_snRNA (SO:0000397).
	sofa["U4ATAC_SNRNA"] = sofa["SO:0000394"] = U4ATAC_SNRNA;

	// U5 RNA is a component of both types of known spliceosome. The precise function of this molecule is unknown, though it is known that the 5' loop is required for splice site selection and p220 binding, and that both the 3' stem-loop and the Sm site are important for Sm protein binding and cap methylation.
	sofa["U5_SNRNA"] = sofa["SO:0000395"] = U5_SNRNA;

	// U6 snRNA is a component of the spliceosome which is involved in splicing pre-mRNA. The putative secondary structure consensus base pairing is confined to a short 5' stem loop, but U6 snRNA is thought to form extensive base-pair interactions with U4 snRNA.
	sofa["U6_SNRNA"] = sofa["SO:0000396"] = U6_SNRNA;

	// U6atac_snRNA is an snRNA required for the splicing of the minor U12-dependent class of eukaryotic nuclear introns. It forms a base paired complex with U4atac_snRNA (SO:0000394).
	sofa["U6ATAC_SNRNA"] = sofa["SO:0000397"] = U6ATAC_SNRNA;

	// U11 snRNA plays a role in splicing of the minor U12-dependent class of eukaryotic nuclear introns, similar to U1 snRNA in the major class spliceosome it base pairs to the conserved 5' splice site sequence.
	sofa["U11_SNRNA"] = sofa["SO:0000398"] = U11_SNRNA;

	// The U12 small nuclear (snRNA), together with U4atac/U6atac, U5, and U11 snRNAs and associated proteins, forms a spliceosome that cleaves a divergent class of low-abundance pre-mRNA introns.
	sofa["U12_SNRNA"] = sofa["SO:0000399"] = U12_SNRNA;

	// U14 small nucleolar RNA (U14 snoRNA) is required for early cleavages of eukaryotic precursor rRNAs. In yeasts, this molecule possess a stem-loop region (known as the Y-domain) which is essential for function. A similar structure, but with a different consensus sequence, is found in plants, but is absent in vertebrates.
	sofa["U14_SNORNA"] = sofa["SO:0000403"] = U14_SNORNA;

	// A family of RNAs are found as part of the enigmatic vault ribonucleoprotein complex. The complex consists of a major vault protein (MVP), two minor vault proteins (VPARP and TEP1), and several small untranslated RNA molecules. It has been suggested that the vault complex is involved in drug resistance.
	sofa["VAULT_RNA"] = sofa["SO:0000404"] = VAULT_RNA;

	// Y RNAs are components of the Ro ribonucleoprotein particle (Ro RNP), in association with Ro60 and La proteins. The Y RNAs and Ro60 and La proteins are well conserved, but the function of the Ro RNP is not known. In humans the RNA component can be one of four small RNAs: hY1, hY3, hY4 and hY5. These small RNAs are predicted to fold into a conserved secondary structure containing three stem structures. The largest of the four, hY1, contains an additional hairpin.
	sofa["Y_RNA"] = sofa["SO:0000405"] = Y_RNA;

	// A large polynucleotide in eukaryotes, which functions as the small subunit of the ribosome.
	sofa["RRNA_18S"] = sofa["SO:0000407"] = RRNA_18S;

	// A region on the surface of a molecule that may interact with another molecule. When applied to polypeptides: Amino acids involved in binding or interactions. It can also apply to an amino acid bond which is represented by the positions of the two flanking amino acids.
	sofa["BINDING_SITE"] = sofa["SO:0000409"] = BINDING_SITE;

	// Any of the individual polynucleotide sequences produced by digestion of DNA with a restriction endonuclease.
	sofa["RESTRICTION_FRAGMENT"] = sofa["SO:0000412"] = RESTRICTION_FRAGMENT;

	// A region where the sequence differs from that of a specified sequence.
	sofa["SEQUENCE_DIFFERENCE"] = sofa["SO:0000413"] = SEQUENCE_DIFFERENCE;

	// The signal_peptide is a short region of the peptide located at the N-terminus that directs the protein to be secreted or part of membrane components.
	sofa["SIGNAL_PEPTIDE"] = sofa["SO:0000418"] = SIGNAL_PEPTIDE;

	// The extent of a polypeptide chain in the mature protein.
	sofa["MATURE_PROTEIN_REGION"] = sofa["SO:0000419"] = MATURE_PROTEIN_REGION;

	// A sequence that can autonomously replicate, as a plasmid, when transformed into a bacterial host.
	sofa["ARS"] = sofa["SO:0000436"] = ARS;

	// A small, 17-28-nt, small interfering RNA derived from transcripts of repetitive elements.
	sofa["RASIRNA"] = sofa["SO:0000454"] = RASIRNA;

	// A non-functional descendent of a functional entity.
	sofa["PSEUDOGENIC_REGION"] = sofa["SO:0000462"] = PSEUDOGENIC_REGION;

	// A non-functional descendant of an exon.
	sofa["DECAYED_EXON"] = sofa["SO:0000464"] = DECAYED_EXON;

	// One of the pieces of sequence that make up a golden path.
	sofa["GOLDEN_PATH_FRAGMENT"] = sofa["SO:0000468"] = GOLDEN_PATH_FRAGMENT;

	// A set of regions which overlap with minimal polymorphism to form a linear sequence.
	sofa["TILING_PATH"] = sofa["SO:0000472"] = TILING_PATH;

	// A piece of sequence that makes up a tiling_path (SO:0000472).
	sofa["TILING_PATH_FRAGMENT"] = sofa["SO:0000474"] = TILING_PATH_FRAGMENT;

	// A primary transcript that is never translated into a protein.
	sofa["NC_PRIMARY_TRANSCRIPT"] = sofa["SO:0000483"] = NC_PRIMARY_TRANSCRIPT;

	// The sequence of the 3' exon that is not coding.
	sofa["THREE_PRIME_CODING_EXON_NONCODING_REGION"] = sofa["SO:0000484"] = THREE_PRIME_CODING_EXON_NONCODING_REGION;

	// The sequence of the 5' exon preceding the start codon.
	sofa["FIVE_PRIME_CODING_EXON_NONCODING_REGION"] = sofa["SO:0000486"] = FIVE_PRIME_CODING_EXON_NONCODING_REGION;

	// A continuous piece of sequence similar to the 'virtual contig' concept of the Ensembl database.
	sofa["VIRTUAL_SEQUENCE"] = sofa["SO:0000499"] = VIRTUAL_SEQUENCE;

	// A region of sequence that is transcribed. This region may cover the transcript of a gene, it may emcompas the sequence covered by all of the transcripts of a alternately spliced gene, or it may cover the region transcribed by a polycistronic transcript. A gene may have 1 or more transcribed regions and a transcribed_region may belong to one or more genes.
	sofa["TRANSCRIBED_REGION"] = sofa["SO:0000502"] = TRANSCRIBED_REGION;

	// The recognition sequence necessary for endonuclease cleavage of an RNA transcript that is followed by polyadenylation; consensus=AATAAA.
	sofa["POLYA_SIGNAL_SEQUENCE"] = sofa["SO:0000551"] = POLYA_SIGNAL_SEQUENCE;

	// The site on an RNA transcript to which will be added adenine residues by post-transcriptional polyadenylation.
	sofa["POLYA_SITE"] = sofa["SO:0000553"] = POLYA_SITE;

	// A region of chromosome where the spindle fibers attach during mitosis and meiosis.
	sofa["CENTROMERE"] = sofa["SO:0000577"] = CENTROMERE;

	// A structure consisting of a 7-methylguanosine in 5'-5' triphosphate linkage with the first nucleotide of an mRNA. It is added post-transcriptionally, and is not encoded in the DNA.
	sofa["CAP"] = sofa["SO:0000581"] = CAP;

	// Group I catalytic introns are large self-splicing ribozymes. They catalyse their own excision from mRNA, tRNA and rRNA precursors in a wide range of organisms. The core secondary structure consists of 9 paired regions (P1-P9). These fold to essentially two domains, the P4-P6 domain (formed from the stacking of P5, P4, P6 and P6a helices) and the P3-P9 domain (formed from the P8, P3, P7 and P9 helices). Group I catalytic introns often have long ORFs inserted in loop regions.
	sofa["GROUP_I_INTRON"] = sofa["SO:0000587"] = GROUP_I_INTRON;

	// A self spliced intron.
	sofa["AUTOCATALYTICALLY_SPLICED_INTRON"] = sofa["SO:0000588"] = AUTOCATALYTICALLY_SPLICED_INTRON;

	// The signal recognition particle (SRP) is a universally conserved ribonucleoprotein. It is involved in the co-translational targeting of proteins to membranes. The eukaryotic SRP consists of a 300-nucleotide 7S RNA and six proteins: SRPs 72, 68, 54, 19, 14, and 9. Archaeal SRP consists of a 7S RNA and homologues of the eukaryotic SRP19 and SRP54 proteins. In most eubacteria, the SRP consists of a 4.5S RNA and the Ffh protein (a homologue of the eukaryotic SRP54 protein). Eukaryotic and archaeal 7S RNAs have very similar secondary structures, with eight helical elements. These fold into the Alu and S domains, separated by a long linker region. Eubacterial SRP is generally a simpler structure, with the M domain of Ffh bound to a region of the 4.5S RNA that corresponds to helix 8 of the eukaryotic and archaeal SRP S domain. Some Gram-positive bacteria (e.g. Bacillus subtilis), however, have a larger SRP RNA that also has an Alu domain. The Alu domain is thought to mediate the peptide chain elongation retardation function of the SRP. The universally conserved helix which interacts with the SRP54/Ffh M domain mediates signal sequence recognition. In eukaryotes and archaea, the SRP19-helix 6 complex is thought to be involved in SRP assembly and stabilizes helix 8 for SRP54 binding.
	sofa["SRP_RNA"] = sofa["SO:0000590"] = SRP_RNA;

	// A short 3'-uridylated RNA that can form a duplex (except for its post-transcriptionally added oligo_U tail (SO:0000609)) with a stretch of mature edited mRNA.
	sofa["GUIDE_RNA"] = sofa["SO:0000602"] = GUIDE_RNA;

	// Group II introns are found in rRNA, tRNA and mRNA of organelles in fungi, plants and protists, and also in mRNA in bacteria. They are large self-splicing ribozymes and have 6 structural domains (usually designated dI to dVI). A subset of group II introns also encode essential splicing proteins in intronic ORFs. The length of these introns can therefore be up to 3kb. Splicing occurs in almost identical fashion to nuclear pre-mRNA splicing with two transesterification steps. The 2' hydroxyl of a bulged adenosine in domain VI attacks the 5' splice site, followed by nucleophilic attack on the 3' splice site by the 3' OH of the upstream exon. Protein machinery is required for splicing in vivo, and long range intron-intron and intron-exon interactions are important for splice site positioning. Group II introns are further sub-classified into groups IIA and IIB which differ in splice site consensus, distance of bulged A from 3' splice site, some tertiary interactions, and intronic ORF phylogeny.
	sofa["GROUP_II_INTRON"] = sofa["SO:0000603"] = GROUP_II_INTRON;

	// The region between two known genes.
	sofa["INTERGENIC_REGION"] = sofa["SO:0000605"] = INTERGENIC_REGION;

	// Sequence of about 100 nucleotides of A added to the 3' end of most eukaryotic mRNAs.
	sofa["POLYA_SEQUENCE"] = sofa["SO:0000610"] = POLYA_SEQUENCE;

	// A pyrimidine rich sequence near the 3' end of an intron to which the 5'end becomes covalently bound during nuclear splicing. The resulting structure resembles a lariat.
	sofa["BRANCH_SITE"] = sofa["SO:0000611"] = BRANCH_SITE;

	// The polypyrimidine tract is one of the cis-acting sequence elements directing intron removal in pre-mRNA splicing.
	sofa["POLYPYRIMIDINE_TRACT"] = sofa["SO:0000612"] = POLYPYRIMIDINE_TRACT;

	// The base where transcription ends.
	sofa["TRANSCRIPTION_END_SITE"] = sofa["SO:0000616"] = TRANSCRIPTION_END_SITE;

	// A specific structure at the end of a linear chromosome, required for the integrity and maintenance of the end.
	sofa["TELOMERE"] = sofa["SO:0000624"] = TELOMERE;

	// Combination of short DNA sequence elements which suppress the transcription of an adjacent gene or genes.
	sofa["SILENCER"] = sofa["SO:0000625"] = SILENCER;

	// A trancriptional cis regulatory region that when located between a CM and a gene's promoter prevents the CRM from modulating that genes expression.
	sofa["INSULATOR"] = sofa["SO:0000627"] = INSULATOR;

	// 
	sofa["CHROMOSOMAL_STRUCTURAL_ELEMENT"] = sofa["SO:0000628"] = CHROMOSOMAL_STRUCTURAL_ELEMENT;

	// A repetitive sequence spanning 500 to 20,000 base pairs (a repeat unit is 5 - 30 base pairs).
	sofa["MINISATELLITE"] = sofa["SO:0000643"] = MINISATELLITE;

	// Antisense RNA is RNA that is transcribed from the coding, rather than the template, strand of DNA. It is therefore complementary to mRNA.
	sofa["ANTISENSE_RNA"] = sofa["SO:0000644"] = ANTISENSE_RNA;

	// The reverse complement of the primary transcript.
	sofa["ANTISENSE_PRIMARY_TRANSCRIPT"] = sofa["SO:0000645"] = ANTISENSE_PRIMARY_TRANSCRIPT;

	// A small RNA molecule that is the product of a longer exogenous or endogenous dsRNA, which is either a bimolecular duplex or very long hairpin, processed (via the Dicer pathway) such that numerous siRNAs accumulate from both strands of the dsRNA. SRNAs trigger the cleavage of their target molecules.
	sofa["SIRNA"] = sofa["SO:0000646"] = SIRNA;

	// Non-coding RNAs of about 21 nucleotides in length that regulate temporal development; first discovered in C. elegans.
	sofa["STRNA"] = sofa["SO:0000649"] = STRNA;

	// Ribosomal RNA transcript that structures the small subunit of the ribosome.
	sofa["SMALL_SUBUNIT_RRNA"] = sofa["SO:0000650"] = SMALL_SUBUNIT_RRNA;

	// Ribosomal RNA transcript that structures the large subunit of the ribosome.
	sofa["LARGE_SUBUNIT_RRNA"] = sofa["SO:0000651"] = LARGE_SUBUNIT_RRNA;

	// 5S ribosomal RNA (5S rRNA) is a component of the large ribosomal subunit in both prokaryotes and eukaryotes. In eukaryotes, it is synthesised by RNA polymerase III (the other eukaryotic rRNAs are cleaved from a 45S precursor synthesised by RNA polymerase I). In Xenopus oocytes, it has been shown that fingers 4-7 of the nine-zinc finger transcription factor TFIIIA can bind to the central region of 5S RNA. Thus, in addition to positively regulating 5S rRNA transcription, TFIIIA also stabilises 5S rRNA until it is required for transcription.
	sofa["RRNA_5S"] = sofa["SO:0000652"] = RRNA_5S;

	// A component of the large ribosomal subunit.
	sofa["RRNA_28S"] = sofa["SO:0000653"] = RRNA_28S;

	// An RNA transcript that does not encode for a protein rather the RNA molecule is the gene product.
	sofa["NCRNA"] = sofa["SO:0000655"] = NCRNA;

	// A region of sequence containing one or more repeat units.
	sofa["REPEAT_REGION"] = sofa["SO:0000657"] = REPEAT_REGION;

	// A repeat that is located at dispersed sites in the genome.
	sofa["DISPERSED_REPEAT"] = sofa["SO:0000658"] = DISPERSED_REPEAT;

	// An intron which is spliced by the spliceosome.
	sofa["SPLICEOSOMAL_INTRON"] = sofa["SO:0000662"] = SPLICEOSOMAL_INTRON;

	// A region of sequence that has been inserted.
	sofa["INSERTION"] = sofa["SO:0000667"] = INSERTION;

	// A match against an EST sequence.
	sofa["EST_MATCH"] = sofa["SO:0000668"] = EST_MATCH;

	// An RNA synthesized on a DNA or RNA template by an RNA polymerase.
	sofa["TRANSCRIPT"] = sofa["SO:0000673"] = TRANSCRIPT;

	// A region of nucleotide sequence targeting by a nuclease enzyme.
	sofa["NUCLEASE_SENSITIVE_SITE"] = sofa["SO:0000684"] = NUCLEASE_SENSITIVE_SITE;

	// The space between two bases in a sequence which marks the position where a deletion has occurred.
	sofa["DELETION_JUNCTION"] = sofa["SO:0000687"] = DELETION_JUNCTION;

	// A set of subregions selected from sequence contigs which when concatenated form a nonredundant linear sequence.
	sofa["GOLDEN_PATH"] = sofa["SO:0000688"] = GOLDEN_PATH;

	// A match against cDNA sequence.
	sofa["CDNA_MATCH"] = sofa["SO:0000689"] = CDNA_MATCH;

	// SNPs are single base pair positions in genomic DNA at which different sequence alternatives (alleles) exist in normal individuals in some population(s), wherein the least frequent allele has an abundance of 1% or greater.
	sofa["SNP"] = sofa["SO:0000694"] = SNP;

	// A sequence used in experiment.
	sofa["REAGENT"] = sofa["SO:0000695"] = REAGENT;

	// A short oligonucleotide sequence, of length on the order of 10's of bases; either single or double stranded.
	sofa["OLIGO"] = sofa["SO:0000696"] = OLIGO;

	// A sequence_feature with an extent of zero.
	sofa["JUNCTION"] = sofa["SO:0000699"] = JUNCTION;

	// A comment about the sequence.
	sofa["REMARK"] = sofa["SO:0000700"] = REMARK;

	// A region of sequence where the validity of the base calling is questionable.
	sofa["POSSIBLE_BASE_CALL_ERROR"] = sofa["SO:0000701"] = POSSIBLE_BASE_CALL_ERROR;

	// A region of sequence where there may have been an error in the assembly.
	sofa["POSSIBLE_ASSEMBLY_ERROR"] = sofa["SO:0000702"] = POSSIBLE_ASSEMBLY_ERROR;

	// A region of sequence implicated in an experimental result.
	sofa["EXPERIMENTAL_RESULT_REGION"] = sofa["SO:0000703"] = EXPERIMENTAL_RESULT_REGION;

	// A region (or regions) that includes all of the sequence elements necessary to encode a functional transcript. A gene may include regulatory regions, transcribed regions and/or other functional sequence regions.
	sofa["GENE"] = sofa["SO:0000704"] = GENE;

	// Two or more adjacent copies of a DNA sequence.
	sofa["TANDEM_REPEAT"] = sofa["SO:0000705"] = TANDEM_REPEAT;

	// The process that produces mature transcripts by combining exons of independent pre-mRNA molecules. The acceptor site lies on the 3' of these molecules.
	sofa["TRANS_SPLICE_ACCEPTOR_SITE"] = sofa["SO:0000706"] = TRANS_SPLICE_ACCEPTOR_SITE;

	// A region of nucleotide sequence corresponding to a known motif.
	sofa["NUCLEOTIDE_MOTIF"] = sofa["SO:0000714"] = NUCLEOTIDE_MOTIF;

	// A nucleic acid sequence that when read as sequential triplets, has the potential of encoding a sequential string of amino acids. It need not contain the start or stop codon.
	sofa["READING_FRAME"] = sofa["SO:0000717"] = READING_FRAME;

	// An ordered and oriented set of scaffolds based on somewhat weaker sets of inferential evidence such as one set of mate pair reads together with supporting evidence from ESTs or location of markers from SNP or microsatellite maps, or cytogenetic localization of contained markers.
	sofa["ULTRACONTIG"] = sofa["SO:0000719"] = ULTRACONTIG;

	// A region of a DNA molecule where transfer is initiated during the process of conjugation or mobilization.
	sofa["ORIT"] = sofa["SO:0000724"] = ORIT;

	// The transit_peptide is a short region at the N-terminus of the peptide that directs the protein to an organelle (chloroplast, mitochondrion, microbody or cyanelle).
	sofa["TRANSIT_PEPTIDE"] = sofa["SO:0000725"] = TRANSIT_PEPTIDE;

	// A gap in the sequence of known length. The unknown bases are filled in with N's.
	sofa["GAP"] = sofa["SO:0000730"] = GAP;

	// 
	sofa["GENE_GROUP_REGULATORY_REGION"] = sofa["SO:0000752"] = GENE_GROUP_REGULATORY_REGION;

	// A non functional descendent of an rRNA.
	sofa["PSEUDOGENIC_RRNA"] = sofa["SO:0000777"] = PSEUDOGENIC_RRNA;

	// A non functional descendent of a tRNA.
	sofa["PSEUDOGENIC_TRNA"] = sofa["SO:0000778"] = PSEUDOGENIC_TRNA;

	// A region of a chromosome.
	sofa["CHROMOSOME_PART"] = sofa["SO:0000830"] = CHROMOSOME_PART;

	// A region of a gene.
	sofa["GENE_MEMBER_REGION"] = sofa["SO:0000831"] = GENE_MEMBER_REGION;

	// A region of a transcript.
	sofa["TRANSCRIPT_REGION"] = sofa["SO:0000833"] = TRANSCRIPT_REGION;

	// A region of a processed transcript.
	sofa["PROCESSED_TRANSCRIPT_REGION"] = sofa["SO:0000834"] = PROCESSED_TRANSCRIPT_REGION;

	// A region of a primary transcript.
	sofa["PRIMARY_TRANSCRIPT_REGION"] = sofa["SO:0000835"] = PRIMARY_TRANSCRIPT_REGION;

	// 
	sofa["MRNA_REGION"] = sofa["SO:0000836"] = MRNA_REGION;

	// A region of UTR.
	sofa["UTR_REGION"] = sofa["SO:0000837"] = UTR_REGION;

	// Biological sequence region that can be assigned to a specific subsequence of a protein.
	sofa["POLYPEPTIDE_REGION"] = sofa["SO:0000839"] = POLYPEPTIDE_REGION;

	// A region within an intron.
	sofa["SPLICEOSOMAL_INTRON_REGION"] = sofa["SO:0000841"] = SPLICEOSOMAL_INTRON_REGION;

	// 
	sofa["GENE_COMPONENT_REGION"] = sofa["SO:0000842"] = GENE_COMPONENT_REGION;

	// 
	sofa["CDS_REGION"] = sofa["SO:0000851"] = CDS_REGION;

	// A large polynucleotide in Bacteria and Archaea, which functions as the small subunit of the ribosome.
	sofa["RRNA_16S"] = sofa["SO:0001000"] = RRNA_16S;

	// A large polynucleotide in Bacteria and Archaea, which functions as the large subunit of the ribosome.
	sofa["RRNA_23S"] = sofa["SO:0001001"] = RRNA_23S;

	// A large polynucleotide which functions as\npart of the large subunit of the ribosome in some eukaryotes.
	sofa["RRNA_25S"] = sofa["SO:0001002"] = RRNA_25S;

	// A DNA sequence that controls the expression of a gene.
	sofa["REGULATORY_REGION"] = sofa["SO:0005836"] = REGULATORY_REGION;

	// A collection of related genes.
	sofa["GENE_GROUP"] = sofa["SO:0005855"] = GENE_GROUP;

	// Any change in genomic DNA caused by a single event.
	sofa["SUBSTITUTION"] = sofa["SO:1000002"] = SUBSTITUTION;

	// When no simple or well defined DNA mutation event describes the observed DNA change, the keyword \"complex\" should be used. Usually there are multiple equally plausible explanations for the change.
	sofa["COMPLEX_SUBSTITUTION"] = sofa["SO:1000005"] = COMPLEX_SUBSTITUTION;

	// A single nucleotide change which has occurred at the same position of a corresponding nucleotide in a reference sequence.
	sofa["POINT_MUTATION"] = sofa["SO:1000008"] = POINT_MUTATION;

	// A continuous nucleotide sequence is inverted in the same position.
	sofa["INVERSION"] = sofa["SO:1000036"] = INVERSION;

	// A group of genes, whether linked as a cluster or not, that respond to a common regulatory signal.
	sofa["REGULON"] = sofa["SO:1001284"] = REGULON;

	// The sequence referred to by an entry in a databank such as Genbank or SwissProt.
	sofa["DATABANK_ENTRY"] = sofa["SO:2000061"] = DATABANK_ENTRY;

};
