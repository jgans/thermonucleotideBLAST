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

#ifdef USE_NCBI

extern "C" {
	#include "subutil.h" // For definition of TOPOLOGY_LINEAR, TOPOLOGY_CIRCULAR, etc.
	#include "accentr.h"
	#include "accutils.h"
	#include "sqnutils.h"
	#include "pmfapi.h"
	#include "ent2api.h"
	#include "explore.h"
}

#endif // USE_NCBI

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <string.h>

using namespace std;

#ifdef USE_NCBI

// Local functions:
void BuildBioseqList(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent);
void FindBioseq(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent);
const char* GetAsnFileType(const string &name);

NLM_EXTERN void GetFeatListCallback(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent);
NLM_EXTERN void ApplyOffsetCallback(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent);

static void touch_bioseq(BioseqPtr m_bsp, Pointer m_userdata);
static Boolean LIBCALLBACK load_bioseq_data(
   SeqLocPtr m_slp, SeqMgrSegmentContextPtr m_context);

void add_loc_offset(SeqLocPtr anp, const Int4 &m_offset);

////////////////////////////////////////////////////////////////////////////
// Global variables for use with NCBI Callback functions
SeqIdPtr global_search_id = NULL;
list<SeqEntryPtr> global_far_sep;
bool is_fetch_enabled = false;

#endif // USE_NCBI

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
		case ASN_1:
			#ifdef USE_NCBI
			ret = loadASN(m_filename, m_pos);
			#endif // USE_NCBI
			break;
		case PTT:
			// Multiple genomes are not currently supported for ASN
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

// Load an entry from the NCBI database
bool DNAMol::load(const int &m_gi)
{
	#ifdef USE_NCBI
	return loadNetASN(m_gi);
	#else
	return false;
	#endif // USE_NCBI
}

#ifdef USE_NCBI
// This function is based on the code contained in ncbi/demo/seqtest.c
// Note that the calling function is responsible for freeing the memory
// associated with "SeqEntryPtr sep".
bool DNAMol::loadASN(SeqEntryPtr sep, SeqIdPtr sid)
{
	BioseqPtr genome_bioseq = NULL;
	SeqPortPtr     spp;
	Uint1          residue;
	
	// Free any exisiting sequence data	
	if(seq){
				
		delete [] seq;

		seq = NULL;
		seq_len = 0;
	}
	
	// Handle ASN.1 data with external links
	const bool fetch_state = is_fetch_enabled;
		
	global_far_sep.clear();		
	Uint2 entityID = ObjMgrGetEntityIDForChoice(sep);
	SeqMgrIndexFeatures(entityID, NULL);
	VisitBioseqsInSep(sep, NULL, touch_bioseq);
	
	if(!fetch_state && is_fetch_enabled){
		
		// If PubSeqFetchEnable started off disabled but finished enabled,
		// then we've used the web to download additional data to supplement an
		// ASN file. Turn off Seq-fetching ...
		
		PubSeqFetchDisable();
		is_fetch_enabled = false;
	}

	// Assume that SeqEntryLoad has already been called!
	//if (! SeqEntryLoad()) {
	//	throw __FILE__ ":DNAMol::loadASN: SeqEntryLoad failed";
	//}

	if(NULL == sep){
		throw "loadASN got a NULL SeqEntryPtr";
	}
			
	// Get taxanomic info about this DNAMol
	taxaInfo( SeqEntryGetSeqDescr(sep, Seq_descr_source, NULL) );
	
	global_search_id = sid;
	BioseqExplore(sep, (Pointer) &genome_bioseq, FindBioseq);
	global_search_id = NULL;
	
	// The DNA molecule of interest is the first bioseq in the file
	if(NULL == genome_bioseq){
		throw "loadASN: Unable to find a single DNA molecule!";
	}

	// Save the SeqId
	sip = SeqIdSetDup(genome_bioseq->id);
	
	// Save the sequence data
	spp = SeqPortNew(genome_bioseq, 0, -1, Seq_strand_plus, Seq_code_iupacna);
	
	if(spp == NULL){
		throw "Didn't read any sequence data. "
			"If downloading, obtain the file directly from NCBI";
	}
	
	//seq = DNA3Seq(BioseqGetLen(genome_bioseq));
	seq_len = BioseqGetLen(genome_bioseq);
	seq = new SEQBASE [seq_len + SEQ_HEADER_SIZE];
	
	if(!seq){
		throw __FILE__ ": Unable to allocate memory for sequence data";
	}
	
	SEQPTR iter = seq;
	
	// Initialize the sequence header
	memcpy( iter, &seq_len, sizeof(unsigned int) );
	iter += sizeof(unsigned int);
	
	while ((residue = SeqPortGetResidue(spp)) != SEQPORT_EOF) {

		switch(residue){
			case 'A': case 'a': 
				*(iter++) = DB_A;
				break;
			case 'T': case 't': 
				*(iter++) = DB_T;
				break;
			case 'G': case 'g': 
				*(iter++) = DB_G;
				break;
			case 'C': case 'c':
				*(iter++) = DB_C;
				break;

			// G or T or C
			case 'B':	case 'b':
				*(iter++) = DB_B;
				break;
			// G or A or T
			case 'D':	case 'd':
				*(iter++) = DB_D;
				break;
			// A or C or T
			case 'H':	case 'h':
				*(iter++) = DB_H;
				break;
			// G or T
			case 'K':	case 'k':
				*(iter++) = DB_K;
				break;
			// A or C
			case 'M':	case 'm':
				*(iter++) = DB_M;
				break;
			// A or C or G or T
			case 'N':	case 'n':
				*(iter++) = DB_N;
				break;
			// G or A
			case 'R':	case 'r':
				*(iter++) = DB_R;
				break;
			// G or C
			case 'S':	case 's':
				*(iter++) = DB_S;
				break;
			// G or C or A
			case 'V':	case 'v':
				*(iter++) = DB_V;
				break;
			// A or T
			case 'W':	case 'w':
				*(iter++) = DB_W;
				break;
			// T or C
			case 'Y':	case 'y':
				*(iter++) = DB_Y;
				break;
			case SEQPORT_VIRT:
				*(iter++) = DB_GAP;
				break;
			default:
				*(iter++) = DB_UNKNOWN;
				break;
		};
	}
				
	spp = SeqPortFree(spp);
		
	list<SeqFeatPtr> feat_list[SEQFEAT_MAX];
	list<SeqFeatPtr>::iterator feat_iter;
	
	// Extract all of the relevant feature pointers in a single pass
	global_search_id = sid;
	SeqEntryExplore(sep, (Pointer)(feat_list), GetFeatListCallback);
	global_search_id = NULL;
	
	if(global_far_sep.empty() == false){
		
		list<SeqEntryPtr>::const_iterator far_iter;
		
		for(far_iter = global_far_sep.begin();far_iter != global_far_sep.end();far_iter++){
			SeqEntryExplore(*far_iter, (Pointer)(feat_list), GetFeatListCallback);
		}
		
		// We're done saving the feature pointers in the array of feature lists
		global_far_sep.clear();
	}
	
	for(feat_iter = feat_list[SEQFEAT_CDREGION].begin();
		feat_iter != feat_list[SEQFEAT_CDREGION].end();feat_iter++){

		GeneAnnotation tmp_gene;
		gene_list.push_back(tmp_gene);

		gene_list.back() = *feat_iter;

		// Did we parse the entry correctly?
		if(gene_list.back().type() == GeneAnnotation::NONE){
			gene_list.pop_back();
		}
	}
	
	for(feat_iter = feat_list[SEQFEAT_GENE].begin();
		feat_iter != feat_list[SEQFEAT_GENE].end();feat_iter++){

		GeneAnnotation tmp_gene;
		
		tmp_gene = *feat_iter;
		
		// Check to see if this gene record is already present in the
		// form of a CDS region
		list<GeneAnnotation>::iterator iter = 
			find(gene_list.begin(), gene_list.end(), tmp_gene);
			
		if(iter == gene_list.end()){

			// Couldn't find the gene
			gene_list.push_back(tmp_gene);

			// Did we parse the entry correctly?
			if(gene_list.back().type() == GeneAnnotation::NONE){
				gene_list.pop_back();
			}
		}
		else{

			// Change the annotation type from CDS to GENE.
			// Don't change PSEUDO_GENES to GENES -- they stay pseudo!
			if(iter->type() == GeneAnnotation::CDS){
				iter->type(tmp_gene.type());
			}
			
			// Add the gene locus and/or locus tag
			iter->info(GeneAnnotation::LOCUS, tmp_gene.info(GeneAnnotation::LOCUS));
			iter->info(GeneAnnotation::LOCUS_TAG, tmp_gene.info(GeneAnnotation::LOCUS_TAG));
		}
	}
		
	for(feat_iter = feat_list[SEQFEAT_RNA].begin();
		feat_iter != feat_list[SEQFEAT_RNA].end();feat_iter++){

		GeneAnnotation tmp_gene;

		tmp_gene = *feat_iter;

		// Check to see if this gene record is already present in a 
		// previously recorded region
		list<GeneAnnotation>::iterator iter = 
			find(gene_list.begin(), gene_list.end(), tmp_gene);
			
		if(iter == gene_list.end()){

			// Couldn't find the gene
			gene_list.push_back(tmp_gene);

			// Did we parse the entry correctly?
			if(gene_list.back().type() == GeneAnnotation::NONE){
				gene_list.pop_back();
			}
		}
		else{

			// Change the annotation type to (t)RNA
			iter->type(tmp_gene.type());
			
			// DON'T add the RNA locus and/or locus tag (since it has none!)
		}
	}
	
	for(feat_iter = feat_list[SEQFEAT_IMP].begin();
		feat_iter != feat_list[SEQFEAT_IMP].end();feat_iter++){

		GeneAnnotation tmp_gene;
		gene_list.push_back(tmp_gene);

		gene_list.back() = *feat_iter;

		// Did we parse the entry correctly?
		if(gene_list.back().type() == GeneAnnotation::NONE){
			gene_list.pop_back();
		}
	}	
	
	for(feat_iter = feat_list[SEQFEAT_PROT].begin();
		feat_iter != feat_list[SEQFEAT_PROT].end();feat_iter++){
	
		if(!(*feat_iter)->location){
			continue;
		}
		
		list<GeneAnnotation>::iterator g_iter = find(gene_list.begin(),
				gene_list.end(), SeqLocId((*feat_iter)->location));
				
		if(g_iter == gene_list.end()){
			continue;
		}
		
		// Assume that sfp->data.choice == SEQFEAT_PROT
		ProtRefPtr prp = (ProtRefPtr)(*feat_iter)->data.value.ptrvalue;
		
		// Set the product annotation
		if(prp->name && prp->name->data.ptrvalue){
			g_iter->info(GeneAnnotation::PRODUCT, (char*)(prp->name->data.ptrvalue));
		}
		else{
			if(prp->desc){
				g_iter->info(GeneAnnotation::PRODUCT, prp->desc);
			}
		}
		
		// Is there an ec value?
		if(prp->ec && prp->ec->data.ptrvalue){
			g_iter->info(GeneAnnotation::EC, (char*)(prp->ec->data.ptrvalue));
		}
	}
	
	processGeneList(true /* Loading this data for the first time */);
	
	// Clean up any objects stored by the NCBI toolkit
	ObjMgrFreeCache(0);
	
	return true;
}

#endif // USE_NCBI

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

#ifdef USE_NCBI

bool DNAMol::loadASN(const std::string &m_filename, streampos &m_pos)
{
	AsnIoPtr		aip;
	AsnTypePtr	atp_se;
	AsnTypePtr	atp;
	AsnModulePtr	amp;
	SeqEntryPtr	sep;

	// Load SeqEntry object loader and sequence alphabets
	if( !SeqEntryLoad() ){
		throw "loadASN: SeqEntryLoad failed";
	}
	
	// Determine the file type (binary or ASCII) by probing the file
	// with GetAsnFileType.
	if ( (aip = AsnIoOpen( (char*)m_filename.c_str(), (char*)GetAsnFileType(m_filename) ) ) == NULL){
		throw "loadASN: Unable to open input file";
	}
	
	// Jump to the specified location in the file
	AsnIoSeek(aip, m_pos);
	
	// We restrict ourselves to reading ASN.1 files that contain annotation and
	// sequence data!
	amp = AsnAllModPtr();
	atp_se = AsnFind("Seq-entry");
	atp = AsnReadId(aip, amp, NULL);

	// The following test only works with ASCII text files
	if( (aip->type & ASNIO_TEXT) && (atp != atp_se) ){
	
		aip = AsnIoClose(aip);

		// If we were in the process of reading mulitple genomes from a single
		// ASN file, then we have read the last genome. Return false to stop
		// checking this file for more data.
		if( (long)m_pos != 0 ){
			return false;
		}
	
		throw "This ASN.1 file does not contain a Seq-entry";
	}

	// Rewind the file
	AsnIoSeek(aip, m_pos);

	// Read in the whole entry into the Sequence Entry Pointer, sep.
	// Close the ASN stream, which in turn closes the input file.
	sep = SeqEntryAsnRead(aip, NULL);
	
	if(sep == NULL){
		// There is no genome to load
		return false;
	}
	
	// Get the current location
	m_pos = AsnIoTell(aip);
	
	aip = AsnIoClose(aip);
			
	loadASN(sep);
	
	sep = SeqEntryFree(sep);
	
	return true;
}

// This function is based on ncbi\demo\getfasta.c
bool DNAMol::loadNetASN(const unsigned int &m_gi)
{
	SeqEntryPtr sep = NULL;
	
	// Load SeqEntry object loader and sequence alphabets
	if (! SeqEntryLoad()) {
		throw "loadNetASN: SeqEntryLoad failed";
	}

	is_fetch_enabled = PubSeqFetchEnable();

	try{
		
		sep = PubSeqSynchronousQuery(m_gi, 0, 0);

		if (sep == NULL){
			throw "Error downloading: sep was NULL";
		}
		
		// This SeqEntryPtr may contain multiple entires. Extract
		// the entry that corresponds to m_gi
		// Search for the desired entry
		SeqId sid;
		
		sid.choice = SEQID_GI;
		sid.data.intvalue = m_gi;
		sid.next = NULL;
		
		loadASN(sep, &sid);		
	}
	catch(const char *error){
	
		if(sep){
			sep = SeqEntryFree(sep);
		}

		PubSeqFetchDisable();
		is_fetch_enabled = false;
		
		throw error;
	}
	catch(...){
	
		if(sep){
			sep = SeqEntryFree(sep);
		}

		PubSeqFetchDisable();
		is_fetch_enabled = false;
		
		throw ":loadNetASN: Unknown error occured parsing ASN data";
	}

	if(sep){
		sep = SeqEntryFree(sep);
	}

	PubSeqFetchDisable();
	is_fetch_enabled = false;
	
	return true;
}


bool DNAMol::loadNetASN(const unsigned int &m_gi, const string &m_filename)
{
	SeqEntryPtr sep = NULL;
	
	// Load SeqEntry object loader and sequence alphabets
	if (! SeqEntryLoad()) {
		throw "loadNetASN: SeqEntryLoad failed";
	}

	is_fetch_enabled = PubSeqFetchEnable();

	try{
		
		sep = PubSeqSynchronousQuery(m_gi, 0, 0);
		
		if (sep == NULL){
			throw "Error downloading: sep was NULL";
		}
		
		// This SeqEntryPtr may contain multiple entires. Extract
		// the entry that corresponds to m_gi
		// Search for the desired entry
		SeqId sid;
		
		sid.choice = SEQID_GI;
		sid.data.intvalue = m_gi;
		sid.next = NULL;
		
		loadASN(sep, &sid);
	}
	catch(const char *error){
		if(sep){
			sep = SeqEntryFree(sep);
		}

		PubSeqFetchDisable();
		is_fetch_enabled = false;
		
		throw error;
	}
	catch(...){
		if(sep){
			sep = SeqEntryFree(sep);
		}

		PubSeqFetchDisable();
		is_fetch_enabled = false;
		
		throw "Unknown error occured parsing ASN data";
	}
			
	// Save this sep to disk
	if(sep == NULL){
		throw "Unable to save ASN data to disk";
	}

	// Save as a text file to allow user inspection (for now)
	AsnIoPtr aip = AsnIoOpen((char*)m_filename.c_str(), "w");
		
	if(aip == NULL){
	
		sep = SeqEntryFree(sep);
		throw "Unable to save file to disk. Please check file name.";
	}
		
	SeqEntryAsnWrite(sep, aip, NULL);
	
	if(aip){
		aip = AsnIoClose(aip);
	}
	
	if(sep){
		sep = SeqEntryFree(sep);
	}
	
	PubSeqFetchDisable();
	is_fetch_enabled = false;
	
	return true;
}

#endif // USE_NCBI

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

#ifdef USE_NCBI
void DNAMol::taxaInfo(ValNodePtr m_vnp)
{	
	if(!m_vnp){
		return;
	}
	
	BioSourcePtr bsp;
	OrgRefPtr orp;
	OrgNamePtr onp;

	bsp = (BioSourcePtr)m_vnp->data.ptrvalue;
	
	if(!bsp){
		return;
	}

	switch(bsp->genome){
		case 0:
			info_map[SOURCE] = "Unknown";
			break;
		case 1:
			info_map[SOURCE] = "genomic";
			break;
		case 2:
			info_map[SOURCE] = "chloroplast";
			break;
		case 3:
			info_map[SOURCE] = "chromoplast";
			break;
		case 4:
			info_map[SOURCE] = "kinetoplast";
			break;
		case 5:
			info_map[SOURCE] = "mitochondrion";
			break;
		case 6:
			info_map[SOURCE] = "plastid";
			break;
		case 7:
			info_map[SOURCE] = "macronuclear";
			break;
		case 8:
			info_map[SOURCE] = "extrachrom";
			break;
		case 9:
			info_map[SOURCE] = "plasmid";
			break;
		case 10:
			info_map[SOURCE] = "transposon";
			break;
		case 11:
			info_map[SOURCE] = "insertion-seq";
			break;
		case 12:
			info_map[SOURCE] = "cyanelle";
			break;
		case 13:
			info_map[SOURCE] = "proviral";
			break;
		case 14:
			info_map[SOURCE] = "viral";
			break;
		case 15:
			info_map[SOURCE] = "nucleomorph";
			break;
		case 16:
			info_map[SOURCE] = "apicoplast";
			break;
		case 17:
			info_map[SOURCE] = "leucoplast";
			break;
		case 18:
			info_map[SOURCE] = "proplastid";
			break;
		case 19:
			info_map[SOURCE] = "endogenous_virus";
			break;
		default:
			throw "taxaInfo: Unknown biosource!";
	};
	
	orp = bsp->org;

	if(orp->taxname){
		info_map[TAXA_NAME] = orp->taxname;
	}
	else{
		if(orp->common){
			info_map[TAXA_NAME] = orp->taxname;
		}
	}

	onp = orp->orgname;

	if(!onp){
		return;
	}

	if(onp){
		if(onp->lineage){
			info_map[LINEAGE] = onp->lineage;
		}
		
		BinomialOrgNamePtr bonp;
		
		switch(onp->choice){
			case 1: // binomial
				bonp = (BinomialOrgNamePtr)onp->data;
				
				if(!bonp){
					return;
				}
				
				if(bonp->genus){
					info_map[GENUS] = bonp->genus;
				}
				
				if(bonp->species){
					info_map[SPECIES] = bonp->species;
				}
				
				if(bonp->subspecies){
					info_map[SUBSPECIES] = bonp->subspecies;
				}
				
				break;
		}
	}
}

#endif // USE_NCBI

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
		
		if( (buffer[i] == '{') || (buffer[i] == '}') ){
		
			// ASCII Asn.1
			ret = DNAMol::ASN_1;
			break;
		}

		if(buffer[i] < 0){
		
			// Binary Asn.1
			ret = DNAMol::ASN_1;
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

#ifdef USE_NCBI
GeneAnnotation& GeneAnnotation::operator=(const SeqFeatPtr &sfp)
{

	if(!sfp){
		throw "GeneAnnotation::=: NULL sfp";
	}

	if(!sfp->location){
		throw "GeneAnnotation::=: NULL sfp->location";
	}
	
	int tmp_start = SeqLocStart(sfp->location);
	int tmp_stop = SeqLocStop(sfp->location);

	if( (tmp_start < 0) || (tmp_stop < 0)){
		gene_start = 0;
		gene_stop = 0;
	}
	else{
		gene_start = (unsigned int)tmp_start;
		gene_stop = (unsigned int)tmp_stop;
	}

	// Is this one contiguous SeqLoc entry ?
	if(IS_one_loc(sfp->location, true) == false){

		// Look through the list of SeqLocs to extract the info
		// we need!
		SeqLocPtr slp = NULL;

		while( (slp = SeqLocFindNext(sfp->location, slp)) != NULL){
			// Save all of the segments
			tmp_start = SeqLocStart(slp);
			tmp_stop = SeqLocStop(slp);
			
			if( (tmp_start >= 0) && (tmp_stop >= 0)){
				// There is a bug in SeqLocStart/SeqLocStop that can cause
				// incorrect gene_start/gene_stop values (see NC_006151.asn for an example;
				// Note the NULL value in the seqloc).
				gene_start = min((int)gene_start, tmp_start);
				gene_stop = max((int)gene_stop, tmp_stop);

				seg_list.push_back( make_pair((unsigned int)tmp_start, (unsigned int)tmp_stop) );
			}
		}

		// Sort the seg list (important for annotations that span the origin)
		seg_list.sort();
	}

	complement = (SeqLocStrand(sfp->location) == Seq_strand_minus);
	
	if(sip){
		sip = SeqIdSetFree(sip);
	}
	
	if(sfp->product){
		sip = SeqIdSetDup(SeqLocId(sfp->product));
	}
	
	if(sfp->comment){
		info_map[NOTE] = sfp->comment;
	}
	
	GeneRefPtr grp;
	RnaRefPtr rrp;
	ImpFeatPtr ifp;
	CdRegionPtr crp;

	switch(sfp->data.choice){
		case SEQFEAT_GENE: // Gene
			grp = (GeneRefPtr)sfp->data.value.ptrvalue;
			
			// Is this a gene or pseudo gene?
			gene_type = (sfp->pseudo) ? PSEUDO_GENE : GENE;

			if(grp->locus){
				info_map[LOCUS] = grp->locus;
			}

			if(grp->locus_tag){
				info_map[LOCUS_TAG] = grp->locus_tag;
			}

			break;
		case SEQFEAT_CDREGION: // CD-region
			crp = (CdRegionPtr)sfp->data.value.ptrvalue;

			// Is this a CDS or pseudo gene?
			gene_type = (sfp->pseudo) ? PSEUDO_GENE : CDS;
			
			{
				// Create a fake gene locus that provides a helpful name that
				// will be displayed on the CDS
				info_map[LOCUS] = seq_id_str();
				
				string::size_type loc = info_map[LOCUS].find('|');

				if(loc != string::npos){
					loc += 1;

					info_map[LOCUS] = info_map[LOCUS].substr(loc, info_map[LOCUS].size() - loc);
				}
			}

			break;
		case SEQFEAT_PROT: // Prot-ref
			throw "GeneAnnotation::=: Protein record assigned to gene!";
			
			break;
		case SEQFEAT_RNA: // RNA
			rrp = (RnaRefPtr)sfp->data.value.ptrvalue;
			
			if(rrp == NULL){
				throw "GeneAnnotation::=: NULL rrp";
			}
			
			// Note that there is also an mRNA type (that has been 
			// associated with packed-int SeqLoc structures). The non-local
			// nature of the packed-int SeqLoc throws off the parsing scheme, so
			// just skip these entries (they will have gene types == NONE and hence
			// be removed).
			switch(rrp->ext.choice){
				case 0: // Fall through
				case 1:
					gene_type = RNA;

					if(rrp->ext.value.ptrvalue != NULL){
						info_map[PRODUCT] = (char*)rrp->ext.value.ptrvalue;
					}
					break;
				case 2:
					// Codon is stored in:
					// tRNAPtr trp = (tRNAPtr)rrp->ext.value.ptrvalue;
					// trp->aa
					gene_type = tRNA;
					info_map[PRODUCT] = "Codon";
					break;
			};
			
			break;
		case SEQFEAT_IMP: // Imp
			ifp = (ImpFeatPtr)sfp->data.value.ptrvalue;
			
			gene_type = IMP;
			
			if(ifp == NULL){
				throw "GeneAnnotation::gene: NULL ifp";
			}
			
			info_map[PRODUCT] = ifp->key;
			
			break;
	}
	
	return *this;
}	

// This SeqEntry exploration function copies the current pointer position in the
// the Bioseq entry to a list of Bioseq pointers
void BuildBioseqList(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	((BioseqPtr PNTR) data)[index] = (BioseqPtr)sep->data.ptrvalue;
	return;
}

// For use with BioseqExplore: Only copy the location of the first bioseq
void FindBioseq(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	// If global_search_id is NULL, return the first bioseq in the record
	// (note that global_search_id is a global variable!)
	if(global_search_id == NULL){
	
		if(0 == index){
			*((BioseqPtr PNTR) data) = (BioseqPtr)sep->data.ptrvalue;
		}
	}
	else{
	
		// If global_search_id is *not* NULL, return the bioseq record that
		// matches global_search_id.
		if( SeqIdIn(global_search_id, BioseqPtr(sep->data.ptrvalue)->id) ){
		
			// We found a match!
			*((BioseqPtr PNTR) data) = (BioseqPtr)sep->data.ptrvalue;
		}
	}
}

NLM_EXTERN void GetFeatListCallback(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	
	if(sep == NULL){
		return;
	}
	
	BioseqPtr bsp = NULL;
	BioseqSetPtr bssp = NULL;
	SeqAnnotPtr sap = NULL;
	SeqFeatPtr sfp = NULL;
	SeqIdPtr sid = NULL;
	list<SeqFeatPtr> *feat_list = (list<SeqFeatPtr> *)data;

	if (IS_Bioseq(sep))
	{
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		sid = bsp->id;
		sap = bsp->annot;
	}
	else
	{
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		sap = bssp->annot;
	}

	while (sap != NULL)
	{
		if (sap->type == 1)  /* feature table */
		{
			// If the global variable global_search_id is NULL, load all
			// entires. If global_search_id == the current SeqId, load
			// the corresponding features
			if( (global_search_id == NULL) || SeqIdIn(global_search_id, sid) ){
			
				for (sfp = (SeqFeatPtr)(sap->data); sfp != NULL; sfp = sfp->next){
					feat_list[sfp->data.choice].push_back(sfp);
				}
			}
			else{
				
				if(global_search_id != NULL){
				
					for (sfp = (SeqFeatPtr)(sap->data); sfp != NULL; sfp = sfp->next){
						
						if( SeqIdIn( global_search_id, SeqLocId(sfp->location) ) ){
							feat_list[sfp->data.choice].push_back(sfp);
						}
					}
				}
			}
		}

		sap = sap->next;
	}
}

NLM_EXTERN void ApplyOffsetCallback(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	
	if(sep == NULL){
		return;
	}
	
	BioseqPtr bsp = NULL;
	SeqAnnotPtr sap = NULL;
	SeqFeatPtr sfp = NULL;
	Int4 from = *((Int4 *)data);
	
	if (IS_Bioseq(sep) )
	{
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		sap = bsp->annot;		
	}
	else
	{
		return;
	}

	while (sap != NULL)
	{
		if (sap->type == 1)  /* feature table */
		{
			for (sfp = (SeqFeatPtr)(sap->data); sfp != NULL; sfp = sfp->next){
				
				// Add an offset to the feature
				add_loc_offset(sfp->location, from);
			}
		}

		sap = sap->next;
	}
}

void add_loc_offset(SeqLocPtr anp, const Int4 &m_offset)
{
	SeqIdPtr sip = NULL;
	
	if(anp == NULL){
		return;
	}

	switch(anp->choice){
		case SEQLOC_INT:
			
			if((SeqIntPtr)anp->data.ptrvalue == NULL){
				throw "Unable to apply annotation offset";
			}
			
			( (SeqIntPtr)anp->data.ptrvalue )->from += m_offset;
			( (SeqIntPtr)anp->data.ptrvalue )->to += m_offset;
			break;
		case SEQLOC_PNT:
			
			if((SeqIntPtr)anp->data.ptrvalue == NULL){
				throw "Unable to apply annotation offset";
			}
			
			( (SeqPntPtr)anp->data.ptrvalue )->point += m_offset;
			break;
		case SEQLOC_MIX:    /* mix -- more than one seq */
		case SEQLOC_EQUIV:    /* equiv -- ditto */
        	case SEQLOC_PACKED_INT:    /* packed int */

			sip = SeqLocId(anp);

			if (sip != NULL){      /* all on one Bioseq */

				SeqLocPtr slp = (SeqLocPtr)anp->data.ptrvalue;

				while (slp != NULL){
				
					add_loc_offset(slp, m_offset);	
					slp = slp->next;
				}
			}

			break;
	};
}

#endif // USE_NCBI

// Test the input file to determine if it is in either binary or ASCII format.
const char* GetAsnFileType(const string &name)
{
        ifstream fin(name.c_str());
 
        if(!fin.is_open()){
                throw "GetAsnFileType:Unable to open file for testing.";
        }
 
        const int sample_size = 15;
        char sample[sample_size];
 
        fin.read(sample, sample_size);
 
        if(fin.gcount() != sample_size){
	   
                throw "GetAsnFileType: Unable to read the required numnber of samples";
        }
 
        for(int i = 0;i < sample_size;i++){
	   
                // The method of detecting a binary file is not very sophisticated: look for
                // any negative 8 bit values (i.e. those values that do not map to ASCII
                // characters).
                if(sample[i] < 0){
			 
                        // We found a non-ASCII character
                        fin.close();
 
                        return "rb";
                }
        }
 
        fin.close();
 
        return "r";
}

string wrap_string(const string &m_str, const unsigned int &m_len)
{
	if(m_len == 0){
		throw __FILE__ ":wrap_string: m_len == 0";
	}
	
	stringstream sout;
	const unsigned int len = m_str.size();
	unsigned int count = 1;
	
	for(unsigned int i = 0;i < len;i++, count++){

		if( (count >= m_len) && isspace(m_str[i]) ){
		
			sout << ENDL;
			count = 1;
		}
		else{
			sout << m_str[i];
		}
	}
	
	return sout.str();
}

bool integer(const string &m_int)
{
	for(string::const_iterator iter = m_int.begin();iter != m_int.end();iter++){
		if( !isdigit(*iter) ){
			return false;
		}
	}
	
	return true;
}

string trim_space(const string &m_str)
{
	int len = m_str.size();
	int start = 0;
	
	if(len == 0){
		return "";
	}
	
	int stop = len - 1;
	
	while( (start < len) && isspace(m_str[start]) ){
		start++;
	}
	
	while( (start <= stop) && isspace(m_str[stop]) ){
		stop--;
	}
	
	return m_str.substr(start, stop - start + 1);
}

// Remove any end of line characters
string pack_string(const string &m_str)
{
	// Push back onto a list of char (instead of string) to retain
	// compatability with the old windows version of std::string
	// [which does not allow string::push_back()].
	list<char> buffer;
	char last_char = -1;
	
	for(string::const_iterator iter = m_str.begin();iter != m_str.end();iter++){
	
		if(*iter == '\n'){
		
			if( (last_char != '\n') && (last_char != '\r') ){
				buffer.push_back(' ');
			}
			
			last_char = *iter;
			continue;
		}
		
		if(*iter == '\r'){
		
			if( (last_char != '\n') && (last_char != '\r') ){
				buffer.push_back(' ');
			}
			
			last_char = *iter;
			continue;
		}

		buffer.push_back(*iter);
		
		last_char = *iter;
	}

	const unsigned int len = buffer.size();
	
	string tmp(len, 'A');
	unsigned int j = 0;
	
	for(list<char>::const_iterator i = buffer.begin();i != buffer.end();i++, j++){
		tmp[j] = *i;
	}
	
	return tmp;
}

#ifdef USE_NCBI
static void touch_bioseq(BioseqPtr m_bsp, Pointer m_userdata)
{
	if(m_bsp == NULL){
		return;
	}
	
	if( (m_bsp->repr == Seq_repr_delta) || (m_bsp->repr == Seq_repr_seg) ){
		
		// If we're loading an ASN file from disk, we may need to go to NCBI over the 
		// net to get additional data
		if(is_fetch_enabled == false){
		
			is_fetch_enabled = PubSeqFetchEnable();
		}
		
		SeqMgrExploreSegments(m_bsp, NULL, load_bioseq_data);
   	} 
}

static Boolean LIBCALLBACK load_bioseq_data(
   SeqLocPtr m_slp, SeqMgrSegmentContextPtr m_context)
{
	BioseqPtr  bsp = NULL;
   	Uint2      entityID;
   	SeqLocPtr  loc;
   	SeqIdPtr   sip;
	Int4 	 from;
	Int4 	 to;
	SeqEntryPtr sep = NULL;
  
   	if( (m_slp == NULL) || (m_context == NULL) ){
		return FALSE;
	}

   	from = m_context->cumOffset;
   	to = from + m_context->to - m_context->from;

   	sip = SeqLocId(m_slp);
	
   	if(sip == NULL){
	
		loc = SeqLocFindNext(m_slp, NULL);
     	
		if(loc != NULL){
       		sip = SeqLocId(loc);
     	}
   	}
	
	if(sip == NULL){
   		return TRUE;
	}
	
	if(sip->choice != SEQID_GI){
		return TRUE;
	}
	
	// Remote fetch genome component if not already in memory
	sep = PubSeqSynchronousQuery(sip->data.intvalue, 0, 0);

	if( (sep == NULL) || (IS_Bioseq(sep) == false) ){
		return TRUE;
	}
	
	// Save this gi in a global list so we know to extract features from it
	global_far_sep.push_back(sep);
	
	bsp = (BioseqPtr)sep->data.ptrvalue;
	
	//bsp = BioseqLockById(sip);

   	if(bsp == NULL){
		return TRUE;
	}

   	entityID = ObjMgrGetEntityIDForPointer(bsp);

   	if(entityID != m_context->entityID){
		
		//if segment not packaged in record, may need to feature index it
		if(SeqMgrFeaturesAreIndexed(entityID) == 0){
			
			SeqMgrIndexFeatures(entityID, NULL);
		}
   	}
	
	// Correct all sequence features by adding the cumulative offset calculated above
	SeqEntryExplore(sep, (Pointer)(&from), ApplyOffsetCallback);	
		
   	return TRUE;
}
#endif // USE_NCBI
