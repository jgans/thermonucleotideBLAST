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

#include "annotation.h"
#include <string.h>

using namespace std;

SeqIdPtr write_accession_GFF3(SeqIdPtr &m_sip, const string &m_accession);
SeqIdPtr write_gi_GFF3(SeqIdPtr &m_sip, const unsigned int &m_gi);

void DNAMol::load(const GFF3File &m_fin, const string &m_source)
{
	info_map[SOURCE] = m_source;
	
	// Is there sequence associated with this source?
	if(m_fin.has_sequence(m_source) == true){
		
		// Load the sequence
		const string &seq_ref = m_fin.sequence(m_source);
		
		seq_len = seq_ref.size();
		
		seq = new SEQBASE [seq_len + SEQ_HEADER_SIZE];
		
		if(!seq){
			throw __FILE__ ": Unable to allocate memory for sequence data";
		}
		
		SEQPTR iter = seq;
	
		// Initialize the sequence length component of the header
		memcpy( iter, &seq_len, sizeof(unsigned int) );
		iter += sizeof(unsigned int);

		for(unsigned int i = 0;i < seq_len;i++){
		
			switch(seq_ref[i]){
				case 'A':
				case 'a':
					*(iter++) = DB_A;
					break;
				case 'T':
				case 't':
					*(iter++) = DB_T;
					break;
				case 'G':
				case 'g':
					*(iter++) = DB_G;
					break;
				case 'C':
				case 'c':
					*(iter++) = DB_C;
					break;
				case 'I':
				case 'i':
					*(iter++) = DB_I;
					break;
				case 'M':
				case 'm':
					*(iter++) = DB_M;
					break;
				case 'R':
				case 'r':
					*(iter++) = DB_R;
					break;
				case 'S':
				case 's':
					*(iter++) = DB_S;
					break;
				case 'V':
				case 'v':
					*(iter++) = DB_V;
					break;
				case 'W':
				case 'w':
					*(iter++) = DB_W;
					break;
				case 'Y':
				case 'y':
					*(iter++) = DB_Y;
					break;
				case 'H':
				case 'h':
					*(iter++) = DB_H;
					break;
				case 'K':
				case 'k':
					*(iter++) = DB_K;
					break;
				case 'D':
				case 'd':
					*(iter++) = DB_D;
					break;
				case 'B':
				case 'b':
					*(iter++) = DB_B;
					break;
				case 'N':
				case 'n':
					*(iter++) = DB_N;
					break;
				case '-':
					*(iter++) = DB_GAP;
					break;
				default:
					*(iter++) = DB_UNKNOWN;
					break;
			};
		}
	}
	
	// Are there features associated with this source?
	if(m_fin.has_features(m_source) == false){
		return;
	}
	
	// There *are* associated features. Load `em up!
	// feature_id is a list of all feature ids for this source and
	// annot are the corresponding features.
	// Currently, parent-child relationships are ignored (because I'm not sure how to
	// represent a child with multiple parents -- which is allowed according to the GFF3
	// file format.
	const vector<string> feature_id = m_fin.feature_id(m_source);
	const map<string, GFF3Record> &annot = m_fin.feature_map(m_source);
	
	for(vector<string>::const_iterator id = feature_id.begin();id != feature_id.end();id++){
		
		map<string, GFF3Record>::const_iterator f = annot.find(*id);
		
		if( f == annot.end() ){
			throw "DNAMol::load: Unable to find features associated with GFF3 id";
		}
		
		const GFF3Record& annot = f->second;
		
		GeneAnnotation tmp;
		pair<GFF3Record::const_iterator, GFF3Record::const_iterator> attrib;
		bool valid_annot = true;

		// Mapping from the GFF3 ontology (SOFA) to the terms used in the GeneAnnotation class
		switch( annot.feature_type() ){
			case GFF3Record::SEQUENCE_ONTOLOGY:
				valid_annot = false;
				break;
			case GFF3Record::CONTIG:
				valid_annot = false;

				// Has organism_name been specified?
				attrib = annot("organism_name");

				if(attrib.first != attrib.second){
					info_map[TAXA_NAME] = attrib.first->second;
				}

				// Has an accession or gi been specified?
				attrib = annot("Dbxref");
				
				for(GFF3Record::const_iterator i = attrib.first;i != attrib.second;i++){

					string::size_type loc = i->second.find("GenBank:");
					
					if(loc != string::npos){

						loc += 8;
						sip = write_accession_GFF3( sip, i->second.substr(loc, i->second.size() - loc) );
					}

					loc = i->second.find("gi:");
					
					if(loc != string::npos){

						loc += 3;
						sip = write_gi_GFF3( sip, 
							atoi( i->second.substr(loc, i->second.size() - loc).c_str() ) );
					}
				}

				break;
			case GFF3Record::CDS:
				tmp.type(GeneAnnotation::CDS);
				break;
			case GFF3Record::GENE:
				tmp.type(GeneAnnotation::GENE);
				break;
			case GFF3Record::PSEUDOGENE:
			case GFF3Record::PSEUDOGENIC_REGION:
			case GFF3Record::PSEUDOGENIC_RRNA:
			case GFF3Record::PSEUDOGENIC_TRNA:
				tmp.type(GeneAnnotation::PSEUDO_GENE);
				break;
			case GFF3Record::TRNA:
				tmp.type(GeneAnnotation::tRNA);
				break;
			case GFF3Record::SCRNA:
			case GFF3Record::MRNA:
			case GFF3Record::RRNA:
			case GFF3Record::SNRNA:
			case GFF3Record::SNORNA:
			case GFF3Record::MIRNA:
			case GFF3Record::SMALL_REGULATORY_NCRNA:
			case GFF3Record::ENZYMATIC_RNA:
			case GFF3Record::RIBOZYME:
			case GFF3Record::RRNA_5_8S:
			case GFF3Record::HAMMERHEAD_RIBOZYME:
			case GFF3Record::RNASE_MRP_RNA:
			case GFF3Record::RNASE_P_RNA:
			case GFF3Record::TELOMERASE_RNA:
			case GFF3Record::U1_SNRNA:
			case GFF3Record::U2_SNRNA:
			case GFF3Record::U4_SNRNA:
			case GFF3Record::U4ATAC_SNRNA:
			case GFF3Record::U5_SNRNA:
			case GFF3Record::U6_SNRNA:
			case GFF3Record::U6ATAC_SNRNA:
			case GFF3Record::U11_SNRNA:
			case GFF3Record::U12_SNRNA:
			case GFF3Record::U14_SNORNA:
			case GFF3Record::VAULT_RNA:
			case GFF3Record::Y_RNA:
			case GFF3Record::RRNA_18S:
			case GFF3Record::RASIRNA:
			case GFF3Record::SRP_RNA:
			case GFF3Record::GUIDE_RNA:
			case GFF3Record::ANTISENSE_RNA:
			case GFF3Record::SIRNA:
			case GFF3Record::STRNA:
			case GFF3Record::SMALL_SUBUNIT_RRNA:
			case GFF3Record::LARGE_SUBUNIT_RRNA:
			case GFF3Record::RRNA_5S:
			case GFF3Record::RRNA_28S:
			case GFF3Record::NCRNA:
			case GFF3Record::MRNA_REGION:
			case GFF3Record::RRNA_16S:
			case GFF3Record::RRNA_23S:
			case GFF3Record::RRNA_25S:
				tmp.type(GeneAnnotation::RNA);
				break;
			case GFF3Record::TF_BINDING_SITE:
			case GFF3Record::BINDING_SITE:
				tmp.type(GeneAnnotation::TFBS);
				break;
			case GFF3Record::PRIMER:
				tmp.type(GeneAnnotation::PRIMER);
				break;
			default:
				tmp.type(GeneAnnotation::IMP);
				break;
		};
		
		if(!valid_annot){
			continue;
		}

		switch( annot.feature_strand() ){
			case GFF3Record::GFF3_PLUS_STRAND:
				tmp.strand(Seq_strand_plus);
				break;
			case GFF3Record::GFF3_MINUS_STRAND:
				tmp.strand(Seq_strand_minus);
				break;
			default:
				tmp.strand(Seq_strand_both);
				break;
		};
		
		tmp.start( annot.feature_start() );
		tmp.stop( annot.feature_stop() );

		// Import attributes that are used by the GeneAnnotation class
		attrib = annot("gene_symbol");
		if(attrib.first != attrib.second){
			tmp.info(GeneAnnotation::LOCUS, attrib.first->second);
		}
		else{
			
			// Try to load name if gene_symbol is not found
			attrib = annot("Name");
			if(attrib.first != attrib.second){
				tmp.info(GeneAnnotation::LOCUS, attrib.first->second);
			}
		}

		attrib = annot("locus_tag");
		if(attrib.first != attrib.second){
			tmp.info(GeneAnnotation::LOCUS_TAG, attrib.first->second);
		}
		else{

			attrib = annot("Name");
			if(attrib.first != attrib.second){
				tmp.info(GeneAnnotation::LOCUS_TAG, attrib.first->second);
			}
		}

		attrib = annot("description");
		if(attrib.first != attrib.second){
			tmp.info(GeneAnnotation::PRODUCT, attrib.first->second);
		}

		SeqIdPtr local_sip = NULL;

		attrib = annot("Dbxref");
		for(GFF3Record::const_iterator i = attrib.first;i != attrib.second;i++){

			string::size_type loc = i->second.find("EC:");

			if(loc != string::npos){

				loc += 3;
				tmp.info( GeneAnnotation::EC, i->second.substr(loc, i->second.size() - loc) );
			}

			loc = i->second.find("protein_id:");

			if(loc != string::npos){

				loc += 11;
				local_sip = write_accession_GFF3( local_sip, i->second.substr(loc, i->second.size() - loc) );
			}
		}

		// Set the SeqId
		if(local_sip){
			tmp.seqid(local_sip);

			local_sip = SeqIdSetFree(local_sip);
		}

		// Is this annotation segmented?
		const list<GFFSegment> &annot_seg = annot.feature_seg();
		
		if(annot_seg.size() > 1){
			
			// Convert the list of GFFSegment to a list of pair<unsigned int, unsigned int>
			list< pair<unsigned int, unsigned int> > seg;
			
			for(list<GFFSegment>::const_iterator i = annot_seg.begin();i != annot_seg.end();i++){
				seg.push_back(i->range);
			}
			
			tmp.segments(seg);
		}

		//tmp.codon_start(annot_seg.front().phase);

		// Save this annotation
		gene_list.push_back(tmp);
	}
	
	processGeneList(true /* Loading this data for the first time */);
}

// Write an accession to the given SeqIdPtr. If an accession entry
// already exists, this code will overwrite it! The updated SedIdPtr is
// returned and the input pointer is invalid.
SeqIdPtr write_accession_GFF3(SeqIdPtr &m_sip, const string &m_accession)
{
	SeqIdPtr sip = NULL;
	SeqIdPtr tmp_new = NULL;
	SeqIdPtr tmp_old = NULL;

	// Is this a properly formatted accession?
	if(m_accession.find('|') != string::npos){
		// Yes
		sip = SeqIdParse ((char*)m_accession.c_str());
	}
	else{
		// No
		sip = SeqIdParse( (char*)string("gb|" + m_accession).c_str() );
	}

	if(m_sip == NULL){
		return sip;
	}

	if(sip == NULL){
		throw error_msg(":write_accession: Unable to parse SeqId");
	}

	tmp_new = sip;
	tmp_old = m_sip;

	while(tmp_old != NULL){
		if(tmp_old->choice != SEQID_GENBANK){
			tmp_new->next = SeqIdDup(tmp_old);
			tmp_new = tmp_new->next;
		}

		tmp_old = tmp_old->next;
	}

	// Free the old ptr
	m_sip = SeqIdSetFree(m_sip);

	return sip;
}

// Write a gi to the given SeqIdPtr. If a gi entry
// already exists, this code will overwrite it! The updated SedIdPtr is
// returned and the input pointer is invalid.
SeqIdPtr write_gi_GFF3(SeqIdPtr &m_sip, const unsigned int &m_gi)
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
