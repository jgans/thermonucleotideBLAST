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

// Local functions
int read_gbk_key(ifstream &fin);
bool read_locus_GBK(ifstream &fin);
void read_accession_GBK(ifstream &fin, string &m_accession);
string read_source_GBK(ifstream &fin);
SEQPTR read_sequence_GBK(ifstream &m_fin, SEQPTR m_seq, unsigned int &m_seq_len);

SEQPTR read_base_count_GBK(ifstream &m_fin, unsigned int &m_seq_len);

unsigned int count_bases_GBK(ifstream &m_fin);

void load_custom_color_gbk(GeneAnnotation &m_annot, const string &m_data);

int next_key_GBK(ifstream &m_fin, const bool &m_clear_line = true);

int parse_gene_GBK(ifstream &m_fin, GeneAnnotation &m_gene);
int parse_rna_GBK(ifstream &m_fin, GeneAnnotation &m_rna);
int parse_rna_GBK(ifstream &m_fin, GeneAnnotation &m_rna, GeneAnnotation &m_gene, bool &m_add_rna);
int parse_trna_GBK(ifstream &m_fin, GeneAnnotation &m_trna);
int parse_trna_GBK(ifstream &m_fin, GeneAnnotation &m_trna, GeneAnnotation &m_gene, bool &m_add_rna);
int parse_imp_GBK(ifstream &m_fin, GeneAnnotation &m_imp);
int parse_user_GBK(ifstream &m_fin, GeneAnnotation &m_sig);
int parse_cds_GBK(ifstream &m_fin, GeneAnnotation &m_cds);
int parse_cds_GBK(ifstream &m_fin, GeneAnnotation &m_cds, GeneAnnotation &m_gene, bool &m_add_gene);

int parse_field_GBK(ifstream &m_fin, pair<string, string> &m_field);

// Enumerate all possible GBK keys
enum{	GBK_EOF = 0,
		GBK_NO_KEY,
		GBK_UNKNOWN_KEY,
		GBK_LOCUS,
		GBK_ACCESSION,
		GBK_VERSION,
		GBK_SOURCE,
		GBK_FEATURES,
		GBK_ORIGIN,
		GBK_CONTIG,
		GBK_BASE_COUNT
};

// Enumerate all annotation keys
enum {
	GBK_ANNOT_END = 0,
	GBK_ANNOT_SOURCE,
	GBK_ANNOT_GENE,
	GBK_ANNOT_CDS,
	GBK_ANNOT_RNA,
	GBK_ANNOT_tRNA,
	GBK_ANNOT_IMP,
	GBK_ANNOT_USER,
	GBK_ANNOT_UNKNOWN,
	GBK_ANNOT_NONE
};

// Track the line number for error reporting!
extern unsigned long int line_number;

bool DNAMol::loadGBK(const std::string &m_filename, streampos &m_pos)
{
	// Read a Genebank flat file

	// First, open the file [in binary mode to allow use of
	// read()].
	ifstream fin(m_filename.c_str(), ios::binary);

	if(!fin){
		throw error_msg("Unable to open Genbank Flat File (gbk)");
	}

	// Are we reading from the body of this file?
	if(m_pos > 0){
		fin.seekg(m_pos);
	}
	else{
		// If we're reading from the head of the file, reset the line counter
		line_number = 1;
	}

	// Free any exisiting sequence data	
	if(seq){
				
		delete [] seq;

		seq = NULL;
		seq_len = 0;
	}

	int key;
	string line;

	// Some defaults for GBK files (or until I find out how to 
	// parse these entries!)
	info_map[SOURCE] = "Unknown";

	while( (key = read_gbk_key(fin)) != GBK_EOF){
		switch(key){
			case GBK_NO_KEY:
				// Read and throw away the line
				getline(fin, line);
				line_number ++;
				break;
			case GBK_UNKNOWN_KEY:
				// Read and throw away the line
				getline(fin, line);
				line_number ++;
				break;
			case GBK_LOCUS:
				read_locus_GBK(fin);
				break;
			case GBK_ACCESSION:
				// Load the NCBI accesion as a SeqIdPtr
				read_accession_GBK(fin, accession);
				break;
			case GBK_VERSION:
				// The version is not currently stored
				break;
			case GBK_SOURCE:
				info_map[TAXA_NAME] = read_source_GBK(fin);
				break;
			case GBK_FEATURES:
				loadGBKFeatures(fin);
				break;
			case GBK_ORIGIN:
				
				seq = read_sequence_GBK(fin, seq, seq_len);

				processGeneList(true /* Loading this data for the first time */);
				
				m_pos = fin.tellg();

				// All done. Is there more data to read?
				return !fin.eof();
			case GBK_CONTIG:
				
				// Read and throw away the CONTIG record, which may be a multi-line
				// record. 
				while( getline(fin, line) ){
					
					// Is the last, non-white space character a comma? If so,
					// then this is a multi-line CONTIG record, and we need to
					// keep reading.
					bool is_multi_line = false;
					
					for(string::const_reverse_iterator i = line.rbegin();i != line.rend();++i){
						
						if(*i == ','){
						
							is_multi_line = true;
							break;
						}
						
						if( isalnum(*i) ){
							break;
						}
					}
					
					if(!is_multi_line){
						break;
					}
				}
				
				line_number ++;
				break;
			case GBK_BASE_COUNT:
			
				seq = read_base_count_GBK(fin, seq_len);
				break;

			default:
				throw error_msg("loadGBK: Unknown key encountered");
		};
	}

	return false;
}

void DNAMol::loadGBKFeatures(ifstream &m_fin)
{
	// Read and process all of the feature elements in a GBK file
	// Skip to the next line
	string buffer;

	getline(m_fin, buffer);
	line_number ++;
	
	int annot_key = next_key_GBK(m_fin);
	int last_annot_key = GBK_ANNOT_NONE;

	GeneAnnotation tmp_gene;

	while(annot_key != GBK_ANNOT_END){
		
		int cur_annot_key = annot_key;

		switch(annot_key){
			case GBK_ANNOT_END:
				return;
			case GBK_ANNOT_NONE:
				// Didn't find a key -- keep reading
				annot_key = next_key_GBK(m_fin);
				break;
			case GBK_ANNOT_SOURCE:
				// Skip the source feature for now. We'll need to parse
				// this feature to extract the taxon id).
				getline(m_fin, buffer);
				line_number ++;
				
				annot_key = next_key_GBK(m_fin);
				break;
			case GBK_ANNOT_GENE:
				annot_key = parse_gene_GBK(m_fin, tmp_gene);

				// Save this gene
				gene_list.push_back(tmp_gene);
				break;
			case GBK_ANNOT_CDS:
				if(last_annot_key == GBK_ANNOT_GENE){
					bool add_gene = true;

					annot_key = parse_cds_GBK(m_fin, tmp_gene, gene_list.back(), add_gene);

					if(add_gene){
						// Save this CDS
						gene_list.push_back(tmp_gene);
					}
					// else {we copied the protein records into the last gene read}
				}
				else{
					annot_key = parse_cds_GBK(m_fin, tmp_gene);

					// Save this CDS
					gene_list.push_back(tmp_gene);
				}

				break;
			case GBK_ANNOT_RNA:
				if(last_annot_key == GBK_ANNOT_GENE){
					bool add_rna = true;

					annot_key = parse_rna_GBK(m_fin, tmp_gene, gene_list.back(), add_rna);

					if(add_rna){
						// Save this RNA
						gene_list.push_back(tmp_gene);
					}
					// else {we copied the protein records into the last gene read}
				}
				else{
					annot_key = parse_rna_GBK(m_fin, tmp_gene);

					// Save this RNA
					gene_list.push_back(tmp_gene);
				}

				break;
			case GBK_ANNOT_tRNA:
				if(last_annot_key == GBK_ANNOT_GENE){
					bool add_trna = true;

					annot_key = parse_trna_GBK(m_fin, tmp_gene, gene_list.back(), add_trna);

					if(add_trna){
						// Save this tRNA
						gene_list.push_back(tmp_gene);
					}
					// else {we copied the protein records into the last gene read}
				}
				else{
					annot_key = parse_trna_GBK(m_fin, tmp_gene);

					// Save this tRNA
					gene_list.push_back(tmp_gene);
				}

				break;
			case GBK_ANNOT_IMP:
				annot_key = parse_imp_GBK(m_fin, tmp_gene);

				// Save this gene
				gene_list.push_back(tmp_gene);
				break;
			case GBK_ANNOT_USER:
				annot_key = parse_user_GBK(m_fin, tmp_gene);

				// Save this gene
				gene_list.push_back(tmp_gene);
				break;
			case GBK_ANNOT_UNKNOWN:
				// Do nothing for now
				break;
		};

		last_annot_key = cur_annot_key;
	}
}

int parse_cds_GBK(ifstream &m_fin, GeneAnnotation &m_cds)
{
	// Clear any existing info
	m_cds.clear();

	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	m_cds.type(GeneAnnotation::CDS);

	m_cds.is_complement(read_range(m_fin, range, seg_list));

	// Define multiple annotation segments if needed
	if( !seg_list.empty() ){
		m_cds.segments(seg_list);
	}

	m_cds.start(range.first);
	m_cds.stop(range.second);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 1024;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "gene"){
			m_cds.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			m_cds.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			m_cds.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			// Promote this CDS to a gene
			m_cds.type(GeneAnnotation::GENE);

			m_cds.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "EC_number"){
			m_cds.info(GeneAnnotation::EC, field.second);
		} else if(field.first == "protein_id"){
			// Set the accession
			m_cds.add_seqid(field.second);
		} else if(field.first == "db_xref"){
			m_cds.add_seqid(field.second);
		} else if(field.first == "pseudo"){
			// Turn this record into a pseduo gene
			m_cds.type(GeneAnnotation::PSEUDO_GENE);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(m_cds, field.second);
		}
	}

	return annot_key;
}

int parse_cds_GBK(ifstream &m_fin, GeneAnnotation &m_cds, GeneAnnotation &m_gene, bool &m_add_gene)
{
	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	bool is_comp = read_range(m_fin, range, seg_list);

	// Does this CDS record match the last gene read?
	if((range.first == m_gene.start()) && (range.second == m_gene.stop())){
		// We have a match!
		m_add_gene = false;
	}
	else{
		// No match -- this CDS record does not correspond to the last gene read
		m_add_gene = true;

		// Clear any existing info
		m_cds.clear();
		
		// Define multiple annotation segments if needed
		if( !seg_list.empty() ){
			m_cds.segments(seg_list);
		}

		m_cds.start(range.first);
		m_cds.stop(range.second);
		m_cds.is_complement(is_comp);
		m_cds.type(GeneAnnotation::CDS);
	}

	GeneAnnotation &gene_ref = ((m_add_gene == true) ? m_cds : m_gene);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 1024;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "gene"){
			gene_ref.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			gene_ref.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			gene_ref.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			gene_ref.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "EC_number"){
			gene_ref.info(GeneAnnotation::EC, field.second);
		} else if(field.first == "protein_id"){
			// Set the accession
			gene_ref.add_seqid(field.second);
		} else if(field.first == "db_xref"){
			gene_ref.add_seqid(field.second);
		} else if(field.first == "pseudo"){
			// Turn this record into a pseduo gene
			gene_ref.type(GeneAnnotation::PSEUDO_GENE);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(m_cds, field.second);
		}
	}

	return annot_key;
}


int parse_gene_GBK(ifstream &m_fin, GeneAnnotation &m_gene)
{
	// Clear any existing info
	m_gene.clear();

	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	m_gene.type(GeneAnnotation::GENE);

	m_gene.is_complement(read_range(m_fin, range, seg_list));

	// Define multiple annotation segments if needed
	if( !seg_list.empty() ){
		m_gene.segments(seg_list);
	}

	m_gene.start(range.first);
	m_gene.stop(range.second);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 32;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "gene"){
			m_gene.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			m_gene.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			m_gene.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			m_gene.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "pseudo"){
			// Turn this record into a pseudo gene
			m_gene.type(GeneAnnotation::PSEUDO_GENE);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(m_gene, field.second);
		}
	}

	return annot_key;
}

int parse_rna_GBK(ifstream &m_fin, GeneAnnotation &m_rna)
{
	// Clear any existing info
	m_rna.clear();

	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	m_rna.type(GeneAnnotation::RNA);

	m_rna.is_complement(read_range(m_fin, range, seg_list));

	// Define multiple annotation segments if needed
	if( !seg_list.empty() ){
		m_rna.segments(seg_list);
	}

	m_rna.start(range.first);
	m_rna.stop(range.second);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 32;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "gene"){
			m_rna.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			m_rna.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			m_rna.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			m_rna.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(m_rna, field.second);
		}
	}

	return annot_key;
}

int parse_rna_GBK(ifstream &m_fin, GeneAnnotation &m_rna, GeneAnnotation &m_gene, bool &m_add_rna)
{
	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	bool is_comp = read_range(m_fin, range, seg_list);

	// Does this RNA record match the last gene read?
	if((range.first == m_gene.start()) && (range.second == m_gene.stop())){
		// We have a match!
		m_add_rna = false;
	}
	else{
		// No match -- this RNA record does not correspond to the last gene read
		m_add_rna = true;

		// Clear any existing info
		m_rna.clear();
		
		// Define multiple annotation segments if needed
		if( !seg_list.empty() ){
			m_rna.segments(seg_list);
		}

		m_rna.start(range.first);
		m_rna.stop(range.second);
		m_rna.is_complement(is_comp);
	}	

	GeneAnnotation &rna_ref = ((m_add_rna == true) ? m_rna : m_gene);

	rna_ref.type(GeneAnnotation::RNA);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 32;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "gene"){
			rna_ref.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			rna_ref.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			rna_ref.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			rna_ref.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(rna_ref, field.second);
		}
	}

	return annot_key;
}

int parse_trna_GBK(ifstream &m_fin, GeneAnnotation &m_trna)
{
	// Clear any existing info
	m_trna.clear();

	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	m_trna.type(GeneAnnotation::tRNA);

	m_trna.is_complement(read_range(m_fin, range, seg_list));

	// Define multiple annotation segments if needed
	if( !seg_list.empty() ){
		m_trna.segments(seg_list);
	}

	m_trna.start(range.first);
	m_trna.stop(range.second);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 32;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "gene"){
			m_trna.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			m_trna.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			m_trna.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			m_trna.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(m_trna, field.second);
		}
	}

	return annot_key;
}

int parse_trna_GBK(ifstream &m_fin, GeneAnnotation &m_trna, GeneAnnotation &m_gene, bool &m_add_trna)
{
	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	bool is_comp = read_range(m_fin, range, seg_list);

	// Does this RNA record match the last gene read?
	if((range.first == m_gene.start()) && (range.second == m_gene.stop())){
		// We have a match!
		m_add_trna = false;
	}
	else{
		// No match -- this tRNA record does not correspond to the last gene read
		m_add_trna = true;

		// Clear any existing info
		m_trna.clear();
		
		// Define multiple annotation segments if needed
		if( !seg_list.empty() ){
			m_trna.segments(seg_list);
		}

		m_trna.start(range.first);
		m_trna.stop(range.second);
		m_trna.is_complement(is_comp);
	}	

	GeneAnnotation &trna_ref = ((m_add_trna == true) ? m_trna : m_gene);

	trna_ref.type(GeneAnnotation::tRNA);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 32;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "gene"){
			trna_ref.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			trna_ref.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			trna_ref.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			trna_ref.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(trna_ref, field.second);
		}
	}

	return annot_key;
}


int parse_imp_GBK(ifstream &m_fin, GeneAnnotation &m_imp)
{
	// Clear any existing info
	m_imp.clear();

	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	m_imp.type(GeneAnnotation::IMP);

	m_imp.is_complement(read_range(m_fin, range, seg_list));

	// Define multiple annotation segments if needed
	if( !seg_list.empty() ){
		m_imp.segments(seg_list);
	}

	m_imp.start(range.first);
	m_imp.stop(range.second);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 128;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "note"){
			m_imp.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			m_imp.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "standard_name"){
			m_imp.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "db_xref"){
			m_imp.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(m_imp, field.second);
		}
	}

	return annot_key;
}

int parse_user_GBK(ifstream &m_fin, GeneAnnotation &m_sig)
{
	// Clear any existing info
	m_sig.clear();

	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	m_sig.type(GeneAnnotation::USER);

	m_sig.is_complement(read_range(m_fin, range, seg_list));

	// Define multiple annotation segments if needed
	if( !seg_list.empty() ){
		m_sig.segments(seg_list);
	}

	m_sig.start(range.first);
	m_sig.stop(range.second);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 128;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "note"){
			m_sig.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "gene"){
			m_sig.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			m_sig.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "product"){
			m_sig.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "standard_name"){
			m_sig.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "db_xref"){
			m_sig.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(m_sig, field.second);
		}
	}

	return annot_key;
}

int parse_field_GBK(ifstream &m_fin, pair<string, string> &m_field)
{
	// Check for a possible annotation key
	const int annot_key = next_key_GBK(m_fin, false /* Don't clear the line */);

	if(annot_key != GBK_ANNOT_NONE){
		return annot_key;
	}

	// Read a field pair in one of two forms:
	//
	// Multiline:
	// /key="string stuff here
	//		over multiple lines"
	// Single Line
	// /key=number
	// Boolean
	// /key

	const unsigned int buffer_size = 96;
	char buffer[buffer_size];
	char *start_ptr;
	char *stop_ptr;
	char *ptr;
	
	// Count the number of matching '(' and ')'
	int paren_count = 0;

	m_fin.getline(buffer, buffer_size);
	line_number ++;
	
	// How many characters did we read?
	unsigned int len = m_fin.gcount();

	// Read the key from buffer
	start_ptr = strchr(buffer, '/');

	if(start_ptr == NULL){
		throw error_msg(":parse_field: Unable to find key start");
	}
	
	// Skip the '/' character
	start_ptr++;

	stop_ptr = strchr(start_ptr, '=');

	if(stop_ptr == NULL){
		// Note that keys do NOT have to have associated fields!
		stop_ptr = start_ptr;
		
		while((*stop_ptr != '\0') && !isspace(*stop_ptr)){
			stop_ptr++;
		}
		
		*stop_ptr = '\0';

		// Just return ...
		m_field.first = start_ptr;

		m_field.second = "";

		return annot_key;
	}

	// Save the key
	*stop_ptr = '\0';
	
	ptr = stop_ptr - 1;

	while(ptr != start_ptr){
		if(isspace(*ptr)){
			*ptr = '\0';

			ptr --;

			continue;
		}
		
		break;
	}
	
	m_field.first = start_ptr;

	// Read the field value
	start_ptr = stop_ptr + 1;

	while(isspace(*start_ptr)){
		start_ptr++;
	}

	if(*start_ptr == '('){
		paren_count ++;
	}

	// Do we have a single line field?
	if((paren_count == 0) && *start_ptr != '"'){
		// Yes -- this is a single line field
		stop_ptr = buffer + len - 2;

		while(isspace(*stop_ptr)){
			*stop_ptr = '\0';

			stop_ptr --;
		}
		
		m_field.second = start_ptr;

		return annot_key;
	}


	// We have a (possibly) multiline field.
	// skip the '"' character [but not the '(' character]
	if(paren_count == 0){
		start_ptr++;
	}

	// Clear the field variable
	m_field.second = "";

	stop_ptr = buffer + len - 2;

	while(true){
		
		while( (*stop_ptr == '\r') || isspace(*stop_ptr) ){
		
			*stop_ptr = '\0';
			stop_ptr --;
		}
		
		// The "(stop_ptr >= start_ptr)" test is to make sure that
		// a field that starts with a '"' on a single line does
		// not derail the parser: i.e.
		// /note="
		//   Highly similar to CBO0900 (85.5 38d), CBO2273 (85.5 38d)
		//   and CBO2267 (93.7 0d)"
		// (this example if from NC_009495.gbk). 
		if( (*stop_ptr == '"') && (stop_ptr >= start_ptr) ){
		
			*stop_ptr = '\0';
			
			// all done
			m_field.second += start_ptr;

			return annot_key;
		}

		if( (paren_count != 0) && (*stop_ptr == ')') ){
		
			// Count the number of parens so far
			ptr = start_ptr + 1; // Skip the first '('

			while(*ptr != '\0'){
				if(*ptr == '('){
					paren_count++;
				}

				if(*ptr == ')'){
					paren_count--;
				}

				ptr++;
			}

			if(paren_count == 0){
				// all done
				m_field.second += start_ptr;

				return annot_key;
			}
		}

		// Save the current line (see the comment above about
		// testing for the pointer ordering)
		if(stop_ptr >= start_ptr){
		
			m_field.second += start_ptr;
	
			// Add a space between multiple lines for easy reading
			m_field.second += ' ';
		}

		// Read another line
		m_fin.getline(buffer, buffer_size);
		line_number ++;
		
		len = m_fin.gcount();
		
		// Empty lines are not allowed!
		if(len == 0){
			throw error_msg("Unexpected blank line or end of file encountered");
		}

		start_ptr = buffer;

		while(isspace(*start_ptr)){
			start_ptr++;
		}

		stop_ptr = buffer + len - 2;
	}

	return annot_key;
}

int next_key_GBK(ifstream &m_fin, const bool &m_clear_line /* = true */)
{
	const int buffer_size = 21;
	
	char buffer[buffer_size + 1];
	char *start_ptr = buffer;

	// Terminate the array
	buffer[buffer_size] = '\0';

	if( !m_fin.read(buffer, buffer_size) ){
		throw error_msg("Unable to read next annotation key");
	}
	
	// Find the start of the string
	while(isspace(*start_ptr)){
		start_ptr++;
	}

	// Is this an empty string?
	if(*start_ptr == '\0'){
		if(m_clear_line){
			// Throw away the rest of the line
			m_fin.ignore(80, '\n');
			line_number ++;
		}

		return GBK_ANNOT_NONE;
	}

	// Make the sting upper case
	for(char *ptr = start_ptr;*ptr != '\0';ptr++){
		*ptr = toupper(*ptr);
	}

	// Have we read into the base count section?
	if(strncmp(start_ptr, "BASE", 4 /*strlen("BASE")*/) == 0){
	
		// Rewind the stream by buffer_size characters
		m_fin.seekg(-buffer_size, ios::cur);

		return GBK_ANNOT_END;
	}

	if(strncmp(start_ptr, "CONTIG", 6 /*strlen("CONTIG")*/) == 0){
	
		// Rewind the stream by buffer_size characters
		m_fin.seekg(-buffer_size, ios::cur);

		return GBK_ANNOT_END;
	}
	
	// Added the test for "//" to prevent a missing nucleic acid sequence
	// record from derailing the parsing of the entire GBK file (2/17/2017).
	if( (strncmp(start_ptr, "ORIGIN", 6 /*strlen("ORIGIN")*/) == 0) ||
	    (strncmp(start_ptr, "//", 2 /*strlen("//")*/) == 0) ){
		// Rewind the stream by buffer_size characters
		m_fin.seekg(-buffer_size, ios::cur);

		return GBK_ANNOT_END;
	}

	if(strncmp(start_ptr, "CDS", 3 /*strlen("CDS")*/) == 0){
		return GBK_ANNOT_CDS;
	}

	if(strncmp(start_ptr, "SOURCE", 6 /*strlen("SOURCE")*/) == 0){
		return GBK_ANNOT_SOURCE;
	}

	if(strncmp(start_ptr, "GENE", 4 /*strlen("GENE")*/) == 0){
		return GBK_ANNOT_GENE;
	}

	if(strncmp(start_ptr, "TRNA", 4 /*strlen("TRNA")*/) == 0){
		return GBK_ANNOT_tRNA;
	}

	// Does the buffer contain the string "RNA"?
	if(strstr(start_ptr, "RNA") != NULL){
		return GBK_ANNOT_RNA;
	}
	
	if(strncmp(start_ptr, "USER", 4 /*strlen("USER")*/) == 0){
		return GBK_ANNOT_USER;
	}
	
	// Did not match this key
	return GBK_ANNOT_IMP;
}

SEQPTR read_base_count_GBK(ifstream &m_fin, unsigned int &m_seq_len)
{
	// Read the number of bases to expect in the sequence
	const char buffer_size = 126;
	char buffer[buffer_size];
	char *start_ptr, *stop_ptr;

	m_fin.getline(buffer, buffer_size);
	line_number ++;
	
	///////////// Read the number of a's /////////////
	// 1) skip leading spaces
	start_ptr = buffer;

	while( !isdigit(*start_ptr) ){
		start_ptr++;
	}

	// 2) Find the end of the "a" counts
	for(stop_ptr = start_ptr;isdigit(*stop_ptr);stop_ptr++){
		// Do nothing!
	}
	
	*stop_ptr = '\0';

	const int a_count = atoi(start_ptr);

	start_ptr = stop_ptr + 1;

	///////////// Read the number of c's /////////////
	// 1) skip leading spaces
	while( !isdigit(*start_ptr) ){
		start_ptr++;
	}

	// 2) Find the end of the "c" counts
	for(stop_ptr = start_ptr;isdigit(*stop_ptr);stop_ptr++){
		// Do nothing!
	}
	
	*stop_ptr = '\0';

	const int c_count = atoi(start_ptr);

	start_ptr = stop_ptr + 1;

	//////////// Read the number of g's /////////////
	// 1) skip leading spaces
	while( !isdigit(*start_ptr) ){
		start_ptr++;
	}

	// 2) Find the end of the "g" counts
	for(stop_ptr = start_ptr;isdigit(*stop_ptr);stop_ptr++){
		// Do nothing!
	}
	
	*stop_ptr = '\0';

	const int g_count = atoi(start_ptr);

	start_ptr = stop_ptr + 1;
	
	//////////// Read the number of t's /////////////
	// 1) skip leading spaces
	while( !isdigit(*start_ptr) ){
		start_ptr++;
	}

	// 2) Find the end of the "t" counts
	for(stop_ptr = start_ptr;isdigit(*stop_ptr);stop_ptr++){
		// Do nothing!
	}
	
	*stop_ptr = '\0';

	const int t_count = atoi(start_ptr);

	start_ptr = stop_ptr + 1;

	//////////// Read the number of Others's /////////////
	// 1) skip leading spaces but allow for optional "others" count
	while((*start_ptr != '\0') && !isdigit(*start_ptr)){
		start_ptr++;
	}
	
	int other_count = 0;

	if(*start_ptr != '\0'){
		// 2) Find the end of the "others" counts
		for(stop_ptr = start_ptr;isdigit(*stop_ptr);stop_ptr++){
			// Do nothing!
		}
		
		*stop_ptr = '\0';

		other_count = atoi(start_ptr);
	}

	// Set the sequence size
	m_seq_len = a_count + t_count + g_count + c_count + other_count;
	
	SEQPTR seq = new SEQBASE [m_seq_len + SEQ_HEADER_SIZE];
	
	if(!seq){
		throw __FILE__ ": Unable to allocate memory for sequence data";
	}
	
	return seq;
}

// Read nucleotide sequence and return the gc content
SEQPTR read_sequence_GBK(ifstream &m_fin, SEQPTR m_seq, unsigned int &m_seq_len)
{
	streampos pos = m_fin.tellg();
	SEQPTR seq = m_seq;
		
	// Have we read the number of bases in the sequence?
	if(m_seq_len == 0){

		// We don't know the size of the DNA molecule
		// Count the number of bases and ALSO the GC content!
		m_seq_len = count_bases_GBK(m_fin);
		
		seq = new SEQBASE [m_seq_len + SEQ_HEADER_SIZE];
		
		if(!seq){
			throw __FILE__ ": Unable to allocate memory for sequence data";
		}

		// Restore the state of the input file
		m_fin.clear();

		m_fin.seekg(pos);
	}

	// The buffer size is a parameter that must be tuned for the file IO of 
	// a given machine.
	const unsigned int buffer_size = 2046;
	char buffer[buffer_size];
	char *ptr;
	
	unsigned int base_count = 0;
	
	SEQPTR iter = seq;
	
	// Initialize the sequence length component of the header
	memcpy( iter, &m_seq_len, sizeof(unsigned int) );
	iter += sizeof(unsigned int);
	
	// First, throw away the line that contains "ORIGIN"
	m_fin.getline(buffer, buffer_size);
	line_number ++;
	
	// Track the current position
	pos += strlen(buffer);

	// Read until we hit a "/" symbol or reach the end of the file
	while( !m_fin.eof() ){
	
		m_fin.read(buffer, buffer_size);
		ptr = buffer;
		
		const unsigned int len = m_fin.gcount();

		for(unsigned int i = 0;i < len;i++, ptr++){
			
			if(*ptr == '/'){
			
				// all done!
				if( base_count != m_seq_len ){					
					throw error_msg("Did not read enough bases");
				}
				
				pos += i;

				m_fin.clear();
				m_fin.seekg(pos);
				
				return seq;
			}
			
			const char local_base = toupper(*ptr);
			
			// Skip numbers and white space
			if( (local_base >= 'A') && (local_base <= 'Z') ){
				
				switch( toupper(local_base) ){
					case 'A':
						*(iter++) = DB_A;
						break;
					case 'T':
						*(iter++) = DB_T;
						break;
					case 'G':
						*(iter++) = DB_G;
						break;
					case 'C':
						*(iter++) = DB_C;
						break;
					case 'I':
						*(iter++) = DB_I;
						break;
					case 'M':
						*(iter++) = DB_M;
						break;
					case 'R':
						*(iter++) = DB_R;
						break;
					case 'S':
						*(iter++) = DB_S;
						break;
					case 'V':
						*(iter++) = DB_V;
						break;
					case 'W':
						*(iter++) = DB_W;
						break;
					case 'Y':
						*(iter++) = DB_Y;
						break;
					case 'H':
						*(iter++) = DB_H;
						break;
					case 'K':
						*(iter++) = DB_K;
						break;
					case 'D':
						*(iter++) = DB_D;
						break;
					case 'B':
						*(iter++) = DB_B;
						break;
					case 'N':
						*(iter++) = DB_N;
						break;
					default:
						*(iter++) = DB_UNKNOWN;
						break;
				};
				
				// Keep track of the number of bases read
				base_count++;
			}
		}

		// Update the current position
		pos += len;
	}

	throw error_msg(":read_sequence: Could not find end-of-sequence terminator");
	
	return NULL;
}

unsigned int count_bases_GBK(ifstream &m_fin)
{
	// The buffer size is a parameter that must be tuned for the file IO of 
	// a given machine.
	const unsigned int buffer_size = 2048;
	char buffer[buffer_size];
	char *ptr;
	
	unsigned int base_count = 0;
	
	// First, throw away the line that contains "ORIGIN"
	m_fin.getline(buffer, buffer_size);
	
	// We're going to rewind, so don't increment the line number counter
	
	// Read until we hit a "/" symbol or reach the end of the file
	while( !m_fin.eof() ){
	
		m_fin.read(buffer, buffer_size);
		ptr = buffer;
		
		const unsigned int len = m_fin.gcount();

		for(unsigned int i = 0;i < len;i++, ptr++){

			// We expect upper case letters, so test for these first
			if( (*ptr >= 'A') && (*ptr <= 'Z') ){
			
				base_count++;
				continue;
			}
			
			if(*ptr == '/'){

				// all done!
				return base_count;
			}
			
			// Test for lower case letters
			if( (*ptr >= 'a') && (*ptr <= 'z') ){
			
				base_count++;
				continue;
			}
					
			#ifdef OPTIMIZE
			 
			switch(buffer[i]){
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
				case ' ':
				case '\t':
				case '\n':
				case '\r':
					// Skip digits and white space
					break;
				case '/':
					// all done!
					return base_count;
				default:
					base_count++;
					break;
			};
			
			#endif // OPTIMIZE
		}
	}

	throw error_msg(":count_bases: Could not find end-of-sequence terminator");

	return 0;
}

string read_source_GBK(ifstream &fin)
{
	string taxa;
	string buffer;

	getline(fin, buffer);
	line_number ++;

	stringstream ss(buffer);

	while(ss >> buffer){
		if(taxa.empty()){
			taxa = buffer;
		}
		else{
			taxa += " " + buffer;
		}
	}

	return taxa;
}

void read_accession_GBK(ifstream &fin, string &m_accession)
{
	fin >> m_accession;
}

bool read_locus_GBK(ifstream &fin)
{
	string line;

	getline(fin, line);
	line_number ++;

	string::iterator iter;

	for(iter = line.begin();iter != line.end(); iter++){
		*iter = tolower(*iter);
	}

	return (string::npos != line.find("circular"));
}

int read_gbk_key(ifstream &fin)
{
	const int buffer_size = 12;
	
	string key;
	
	// Allow for malformed GBK files that do *not* have the correct space-padding after
	// key words
	for(int i = 0;i < buffer_size;++i){
		
		const int c = fin.get();
		
		if(c == '\r'){
			continue;
		}
		
		if(c == '\n'){
			break;
		}
		
		if(c == EOF){
			return GBK_EOF;
		}
		
		if( !isspace(c) ){
			key.push_back(c);
		}
	}

	if(key.empty()){
		return GBK_NO_KEY;
	}

	string::iterator iter;
	

	for(iter = key.begin();iter != key.end();iter++){
		*iter = toupper(*iter);
	}

	if(key == "LOCUS"){
		return GBK_LOCUS;
	}

	if(key == "ACCESSION"){
		return GBK_ACCESSION;
	}

	if(key == "VERSION"){
		return GBK_VERSION;
	}

	if(key == "SOURCE"){
		return GBK_SOURCE;
	}

	if(key == "FEATURES"){
		return GBK_FEATURES;
	}

	if(key == "CONTIG"){
		return GBK_CONTIG;
	}
	
	if(key == "ORIGIN"){
		return GBK_ORIGIN;
	}
	
	if(key == "BASE"){
		return GBK_BASE_COUNT;
	}

	// Did not match this key
	return GBK_UNKNOWN_KEY;
}

void load_custom_color_gbk(GeneAnnotation &m_annot, const string &m_data)
{
	// The load_custom_color_gbk function is only included for code-base
	// compatibility. No changes are made to m_annot and this function
	// will remain an empty stub.
}
