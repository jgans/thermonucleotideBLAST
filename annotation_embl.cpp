#include "annotation.h"

#include <sstream>
#include <string.h>
#include <zlib.h>

using namespace std;

// Local functions
int read_EMBL_key(gzFile fin);
bool read_locus_EMBL(gzFile fin);
void read_accession_EMBL(gzFile fin, string &m_str);
string read_source_EMBL(gzFile fin);
SEQPTR read_sequence_EMBL(gzFile m_fin, unsigned int &m_seq_len);

int next_key_EMBL(gzFile m_fin, const bool &m_clear_line = true);

int parse_gene_EMBL(gzFile m_fin, GeneAnnotation &m_gene);
int parse_rna_EMBL(gzFile m_fin, GeneAnnotation &m_rna);
int parse_rna_EMBL(gzFile m_fin, GeneAnnotation &m_rna, GeneAnnotation &m_gene, bool &m_add_rna);
int parse_trna_EMBL(gzFile m_fin, GeneAnnotation &m_trna);
int parse_trna_EMBL(gzFile m_fin, GeneAnnotation &m_trna, GeneAnnotation &m_gene, bool &m_add_rna);
int parse_imp_EMBL(gzFile m_fin, GeneAnnotation &m_imp);
int parse_cds_EMBL(gzFile m_fin, GeneAnnotation &m_cds);
int parse_cds_EMBL(gzFile m_fin, GeneAnnotation &m_cds, GeneAnnotation &m_gene, bool &m_add_gene);

int parse_field_EMBL(gzFile m_fin, pair<string, string> &m_field);

// Enumerate all possible EMBL keys
enum{	EMBL_EOF = 0,
		EMBL_NO_KEY,
		EMBL_UNKNOWN_KEY,
		EMBL_LOCUS,
		EMBL_ACCESSION,
		EMBL_VERSION,
		EMBL_SOURCE,
		EMBL_FEATURES,
		EMBL_ORIGIN
};

// Enumerate all annotation keys
enum {
	EMBL_ANNOT_END = 0,
	EMBL_ANNOT_SOURCE,
	EMBL_ANNOT_GENE,
	EMBL_ANNOT_CDS,
	EMBL_ANNOT_RNA,
	EMBL_ANNOT_tRNA,
	EMBL_ANNOT_IMP,
	EMBL_ANNOT_UNKNOWN,
	EMBL_ANNOT_NONE
};

#define	MAX_LINE_LEN	1024

bool DNAMol::loadEMBL(gzFile m_fin, size_t &m_pos)
{
	// Read a EMBL file

	// Are we reading from the body of this file?
	if(m_pos > 0){
		gzseek(m_fin, m_pos, SEEK_SET);
	}

	int key;
	char line[MAX_LINE_LEN];

	// Some defaults for EMBL files (or until I find out how to 
	// parse these entries!)
	info_map[SOURCE] = "Unknown";

	while( (key = read_EMBL_key(m_fin)) != EMBL_EOF){
		switch(key){
			case EMBL_NO_KEY:
				// Read and throw away the line
				if( (gzgets(m_fin, line, MAX_LINE_LEN) == NULL) || !strip_eol(line, MAX_LINE_LEN) ){
					throw __FILE__ ":DNAMol::loadEMBL: Unable to read EMBL_NO_KEY";
				}
				break;
			case EMBL_UNKNOWN_KEY:
				// Read and throw away the line
				if( (gzgets(m_fin, line, MAX_LINE_LEN) == NULL) || !strip_eol(line, MAX_LINE_LEN) ){
					throw __FILE__ ":DNAMol::loadEMBL: Unable to read EMBL_UNKNOWN_KEY";
				}
				break;
			case EMBL_LOCUS:
				read_locus_EMBL(m_fin);
				break;
			case EMBL_ACCESSION:
				// Load the NCBI accesion as a SeqIdPtr
				read_accession_EMBL(m_fin, accession);
				break;
			case EMBL_VERSION:
				// The version is not currently stored
				break;
			case EMBL_SOURCE:
				info_map[TAXA_NAME] = read_source_EMBL(m_fin);
				break;
			case EMBL_FEATURES:
				loadEMBLFeatures(m_fin);
				break;
			case EMBL_ORIGIN:
			
				// The read_sequence function has the option to 
				// set the gc_content
				if(seq){
				
					delete [] seq;
					
					seq = NULL;
					seq_len = 0;
				}
				
				seq = read_sequence_EMBL(m_fin, seq_len);

				processGeneList(true /* Loading this data for the first time */);
				
				m_pos = gztell(m_fin);

				// All done. Is there more data to read?
				return (gzeof(m_fin) == 0);
			default:
				throw "loadEMBL: Unknown key encountered!";
		};
	}

	return false;
}

void DNAMol::loadEMBLFeatures(gzFile m_fin)
{
	// Read and process all of the feature elements in an EMBL file
	// Skip to the next line
	char buffer[MAX_LINE_LEN];

	size_t pos = gztell(m_fin);
	
	if( (gzgets(m_fin, buffer, MAX_LINE_LEN) == NULL) || !strip_eol(buffer, MAX_LINE_LEN) ){
		throw __FILE__ ":DNAMol::loadEMBLFeatures: Error reading the first line";
	}

	
	if(strstr(buffer, "FT") != NULL){

		// This is an Artemis file -- don't skip this line!
		gzseek(m_fin, pos, SEEK_SET);
	}
	
	int annot_key = next_key_EMBL(m_fin);
	int last_annot_key = EMBL_ANNOT_NONE;

	GeneAnnotation tmp_gene;

	while(annot_key != EMBL_ANNOT_END){
		
		int cur_annot_key = annot_key;

		switch(annot_key){
			case EMBL_ANNOT_END:
				return;
			case EMBL_ANNOT_NONE:
				// Didn't find a key -- keep reading
				annot_key = next_key_EMBL(m_fin);
				break;
			case EMBL_ANNOT_SOURCE:
				// Skip the source feature for now. We'll need to parse
				// this feature to extract the taxon id).
				if( (gzgets(m_fin, buffer, MAX_LINE_LEN) == NULL) || !strip_eol(buffer, MAX_LINE_LEN) ){
					throw __FILE__ ":DNAMol::loadEMBLFeatures: Error reading EMBL_ANNOT_SOURCE";
				}

				annot_key = next_key_EMBL(m_fin);
				break;
			case EMBL_ANNOT_GENE:
				annot_key = parse_gene_EMBL(m_fin, tmp_gene);

				// Save this gene
				gene_list.push_back(tmp_gene);
				break;
			case EMBL_ANNOT_CDS:
				if(last_annot_key == EMBL_ANNOT_GENE){
					bool add_gene = true;

					annot_key = parse_cds_EMBL(m_fin, tmp_gene, gene_list.back(), add_gene);

					if(add_gene){
						// Save this CDS
						gene_list.push_back(tmp_gene);
					}
					// else {we copied the protein records into the last gene read}
				}
				else{
					annot_key = parse_cds_EMBL(m_fin, tmp_gene);

					// Save this CDS
					gene_list.push_back(tmp_gene);
				}

				break;
			case EMBL_ANNOT_RNA:
				if(last_annot_key == EMBL_ANNOT_GENE){
					bool add_rna = true;

					annot_key = parse_rna_EMBL(m_fin, tmp_gene, gene_list.back(), add_rna);

					if(add_rna){
						// Save this RNA
						gene_list.push_back(tmp_gene);
					}
					// else {we copied the protein records into the last gene read}
				}
				else{
					annot_key = parse_rna_EMBL(m_fin, tmp_gene);

					// Save this RNA
					gene_list.push_back(tmp_gene);
				}

				break;
			case EMBL_ANNOT_tRNA:
				if(last_annot_key == EMBL_ANNOT_GENE){
					bool add_trna = true;

					annot_key = parse_trna_EMBL(m_fin, tmp_gene, gene_list.back(), add_trna);

					if(add_trna){
						// Save this tRNA
						gene_list.push_back(tmp_gene);
					}
					// else {we copied the protein records into the last gene read}
				}
				else{
					annot_key = parse_trna_EMBL(m_fin, tmp_gene);

					// Save this tRNA
					gene_list.push_back(tmp_gene);
				}

				break;
			case EMBL_ANNOT_IMP:
				annot_key = parse_imp_EMBL(m_fin, tmp_gene);

				// Save this gene
				gene_list.push_back(tmp_gene);
				break;
			case EMBL_ANNOT_UNKNOWN:
				// Do nothing for now
				break;
		};

		last_annot_key = cur_annot_key;
	}
}

int parse_cds_EMBL(gzFile m_fin, GeneAnnotation &m_cds)
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

	while( (annot_key = parse_field_EMBL(m_fin, field)) == EMBL_ANNOT_NONE){
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
			
			string::size_type pos = field.second.find(':');

			if(pos != string::npos){

				pos ++;
				m_cds.add_seqid( field.second.substr(pos, field.second.size() - pos) );
			}
		} else if(field.first == "pseudo"){
			// Change this CDS into a pseduo-gene
			m_cds.type(GeneAnnotation::PSEUDO_GENE);
		}
	}

	return annot_key;
}

int parse_cds_EMBL(gzFile m_fin, GeneAnnotation &m_cds, GeneAnnotation &m_gene, bool &m_add_gene)
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

	while( (annot_key = parse_field_EMBL(m_fin, field)) == EMBL_ANNOT_NONE){
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
			
			string::size_type pos = field.second.find(':');

			if(pos != string::npos){

				pos ++;
				gene_ref.add_seqid( field.second.substr(pos, field.second.size() - pos) );
			}
		} else if(field.first == "pseudo"){
			// Change this CDS into a pseduo-gene
			gene_ref.type(GeneAnnotation::PSEUDO_GENE);
		}
	}

	return annot_key;
}

int parse_gene_EMBL(gzFile m_fin, GeneAnnotation &m_gene)
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

	while( (annot_key = parse_field_EMBL(m_fin, field)) == EMBL_ANNOT_NONE){
		if(field.first == "gene"){
			m_gene.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			m_gene.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			m_gene.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			m_gene.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "pseudo"){
			// Turn this record into a pseduo gene
			m_gene.type(GeneAnnotation::PSEUDO_GENE);
		}
	}

	return annot_key;
}

int parse_rna_EMBL(gzFile m_fin, GeneAnnotation &m_rna)
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

	while( (annot_key = parse_field_EMBL(m_fin, field)) == EMBL_ANNOT_NONE){
		if(field.first == "gene"){
			m_rna.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			m_rna.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			m_rna.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			m_rna.info(GeneAnnotation::PRODUCT, field.second);
		}
	}

	return annot_key;
}

int parse_rna_EMBL(gzFile m_fin, GeneAnnotation &m_rna, GeneAnnotation &m_gene, bool &m_add_rna)
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

	while( (annot_key = parse_field_EMBL(m_fin, field)) == EMBL_ANNOT_NONE){
		if(field.first == "gene"){
			rna_ref.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			rna_ref.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			rna_ref.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			rna_ref.info(GeneAnnotation::PRODUCT, field.second);
		}
	}

	return annot_key;
}

int parse_trna_EMBL(gzFile m_fin, GeneAnnotation &m_trna)
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

	while( (annot_key = parse_field_EMBL(m_fin, field)) == EMBL_ANNOT_NONE){
		if(field.first == "gene"){
			m_trna.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			m_trna.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			m_trna.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			m_trna.info(GeneAnnotation::PRODUCT, field.second);
		}
	}

	return annot_key;
}

int parse_trna_EMBL(gzFile m_fin, GeneAnnotation &m_trna, GeneAnnotation &m_gene, bool &m_add_trna)
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

	while( (annot_key = parse_field_EMBL(m_fin, field)) == EMBL_ANNOT_NONE){
		if(field.first == "gene"){
			trna_ref.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			trna_ref.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			trna_ref.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			trna_ref.info(GeneAnnotation::PRODUCT, field.second);
		}
	}

	return annot_key;
}


int parse_imp_EMBL(gzFile m_fin, GeneAnnotation &m_imp)
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

	while( (annot_key = parse_field_EMBL(m_fin, field)) == EMBL_ANNOT_NONE){
		if(field.first == "note"){
			m_imp.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			m_imp.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "standard_name"){
			m_imp.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "db_xref"){
			m_imp.info(GeneAnnotation::LOCUS_TAG, field.second);
		}
	}

	return annot_key;
}

int parse_field_EMBL(gzFile m_fin, pair<string, string> &m_field)
{
	// Check for a possible annotation key
	const int annot_key = next_key_EMBL(m_fin, false /* Don't clear the line */);

	if(annot_key != EMBL_ANNOT_NONE){
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

	//const unsigned int buffer_size = 96;
	const unsigned int buffer_size = 102;
	char buffer[buffer_size];
	char *start_ptr;
	char *stop_ptr;
	char *ptr;
	
	// Count the number of matching '(' and ')'
	int paren_count = 0;

	if( (gzgets(m_fin, buffer, buffer_size) == NULL) || !strip_eol(buffer, buffer_size) ){
		throw __FILE__ ":parse_field_EMBL: Error reading line (1)";
	}
	
	// How many characters did we read?
	unsigned int len = strlen(buffer);

	// Read the key from buffer
	start_ptr = strchr(buffer, '/');

	if(start_ptr == NULL){
		throw ":parse_field_EMBL: Unable to find key start";
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
		stop_ptr = buffer + (len - 1);

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
		
		while((*stop_ptr == '\r') || isspace(*stop_ptr)){
			*stop_ptr = '\0';
			stop_ptr --;
		}
		
		// The "(stop_ptr > start_ptr)" test is to make sure that
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

		// Read another line, skipping the first two characters!
		gzgetc(m_fin);
		gzgetc(m_fin);

		if( (gzgets(m_fin, buffer, buffer_size) == NULL) || !strip_eol(buffer, buffer_size) ){
			throw __FILE__ ":parse_field_EMBL: Unable to read line (2)";
		}
		
		len = strlen(buffer);

		// Empty lines are not allowed!
		if(len == 0){
			throw "Unexpected blank line or end of file encountered!";
		}

		start_ptr = buffer;

		while(isspace(*start_ptr)){
			start_ptr++;
		}

		stop_ptr = buffer + (len - 1);
	}

	return annot_key;
}


int next_key_EMBL(gzFile m_fin, const bool &m_clear_line /* = true */)
{
	char buffer[MAX_LINE_LEN];
	char *start_ptr = buffer;

	// First read the two letter code
	if( gzread(m_fin, buffer, 2) ){
		throw ":next_key_EMBL: Unable to read next annotation key";
	}

	// Check the two letter record code
	if(buffer[0] == 'X'){

		// Throw away the rest of the line
		if( (gzgets(m_fin, buffer, MAX_LINE_LEN) == NULL) || !strip_eol(buffer, MAX_LINE_LEN) ){
			throw ":next_key_EMBL: Unable to read remaining characters in line";
		}

		return EMBL_ANNOT_END;
	}
	
	if(buffer[0] == 'S'){

		// This is an Artemis file, rewind by two bases
		gzungetc(buffer[1], m_fin);
		gzungetc(buffer[0], m_fin);
		
		return EMBL_ANNOT_END;
	}

	if((buffer[0] == 'F') && (buffer[1] == 'H')){

		// Throw away the rest of the line
		if( (gzgets(m_fin, buffer, MAX_LINE_LEN) == NULL) || !strip_eol(buffer, MAX_LINE_LEN) ){
			throw __FILE__ ":next_key_EMBL: Unable to skip remaining characters in line";
		}

		return EMBL_ANNOT_NONE;
	}
	
	// We had better have an "FT" !
	if((buffer[0] != 'F') && (buffer[1] != 'T')){
		throw __FILE__ ":next_key_EMBL: Premature end of file or blank line encountered";
	}

	if( (gzgets(m_fin, buffer, MAX_LINE_LEN) == NULL) || !strip_eol(buffer, MAX_LINE_LEN) ){
		throw __FILE__ ":next_key_EMBL: Unable to read next annotation key";
	}

	// Find the start of the string
	while(isspace(*start_ptr)){
		start_ptr++;
	}

	// Is this an empty string?
	if(*start_ptr == '\0'){

		if(m_clear_line){
			// Throw away the rest of the line
			if( (gzgets(m_fin, buffer, MAX_LINE_LEN) == NULL) || !strip_eol(buffer, MAX_LINE_LEN) ){
				throw __FILE__ ":next_key_EMBL: Unable to discard remaining characters in line";
			}
		}

		return EMBL_ANNOT_NONE;
	}

	// Make the sting upper case
	for(char *ptr = start_ptr;*ptr != '\0';ptr++){
		*ptr = toupper(*ptr);
	}

	if(strncmp(start_ptr, "CDS", 3 /*strlen("CDS")*/) == 0){
		return EMBL_ANNOT_CDS;
	}

	if(strncmp(start_ptr, "SOURCE", 6 /*strlen("SOURCE")*/) == 0){
		return EMBL_ANNOT_SOURCE;
	}

	if(strncmp(start_ptr, "GENE", 4 /*strlen("GENE")*/) == 0){
		return EMBL_ANNOT_GENE;
	}

	if(strncmp(start_ptr, "TRNA", 4 /*strlen("TRNA")*/) == 0){
		return EMBL_ANNOT_tRNA;
	}

	// Does the buffer contain the string "RNA"?
	if(strstr(start_ptr, "RNA") != NULL){
		return EMBL_ANNOT_RNA;
	}

	// Did not match this key
	return EMBL_ANNOT_IMP;
}

// EMBL files use a two letter key code
int read_EMBL_key(gzFile m_fin)
{
	const int buffer_size = 2;
	
	char buffer[buffer_size + 1];

	// Terminate the string and zero the array
	memset(buffer, '\0', buffer_size);

	if(gzread(m_fin, buffer, buffer_size) != buffer_size){
		return EMBL_EOF;
	}

	buffer[0] = toupper(buffer[0]);
	buffer[1] = toupper(buffer[1]);

	if((buffer[0] == 'X') && (buffer[1] == 'X')){
		return EMBL_NO_KEY;
	}

	if((buffer[0] == 'I') && (buffer[1] == 'D')){
		return EMBL_LOCUS;
	}

	if((buffer[0] == 'A') && (buffer[1] == 'C')){
		return EMBL_ACCESSION;
	}

	if((buffer[0] == 'S') && (buffer[1] == 'V')){
		return EMBL_VERSION;
	}

	if((buffer[0] == 'O') && (buffer[1] == 'S')){
		return EMBL_SOURCE;
	}

	if((buffer[0] == 'F') && (buffer[1] == 'H')){
		return EMBL_FEATURES;
	}

	if((buffer[0] == 'S') && (buffer[1] == 'Q')){
		return EMBL_ORIGIN;
	}
	
	if((buffer[0] == 'F') && (buffer[1] == 'T')){
		
		// Rewind by two bases
		gzungetc(buffer[1], m_fin);
		gzungetc(buffer[0], m_fin);

		return EMBL_FEATURES;
	}
	
	// Did not match this key
	return EMBL_UNKNOWN_KEY;
}

bool read_locus_EMBL(gzFile m_fin)
{
	char line[MAX_LINE_LEN];

	if( (gzgets(m_fin, line, MAX_LINE_LEN) == NULL) || !strip_eol(line, MAX_LINE_LEN) ){
		throw __FILE__ ":read_locus_EMBL: Unable to read line";
	}
	
	const int len = strlen(line);

	for(int i = 0;i < len;++i){
		line[i] = tolower(line[i]);
	}

	return (strstr(line, "circular") != NULL);
}

void read_accession_EMBL(gzFile m_fin, string &m_accession)
{
	m_accession.clear();

	char c;

	while( ( c = gzgetc(m_fin) ) != -1 ){

		if( !m_accession.empty() && isspace(c) ){

			gzungetc(c, m_fin);
			return;
		}

		m_accession.push_back(c);
	}
}

string read_source_EMBL(gzFile m_fin)
{
	string taxa;
	string buffer;

	char line[MAX_LINE_LEN];

	if( (gzgets(m_fin, line, MAX_LINE_LEN) == NULL) || !strip_eol(line, MAX_LINE_LEN) ){
		throw __FILE__ ":read_source_EMBL: Unable to read line";
	}

	buffer = line;

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

// Read nucleotide sequence and return the gc content
SEQPTR read_sequence_EMBL(gzFile m_fin, unsigned int &m_seq_len)
{
	// The buffer size is a parameter that must be tuned for the file IO of 
	// a given machine.
	const unsigned int buffer_size = 2046;
	char buffer[buffer_size];

	size_t pos = gztell(m_fin);

	// First, read the base counts from the Sequence line, i.e.:
	// Sequence 70159 BP; 19675 A; 15631 C; 15833 G; 19020 T; 0 other;
	if( (gzgets(m_fin, buffer, buffer_size) == NULL) || !strip_eol(buffer, buffer_size) ){
		throw __FILE__ ":read_sequence_EMBL: Unable to read line (1)";
	}
	
	// Track the current position
	pos += strlen(buffer);

	int a_count = -1;
	int t_count = -1;
	int g_count = -1;
	int c_count = -1;
	int other_count = -1;
	
	char *ptr = strchr(buffer, ';') + 1;
	list<char> tmp_num;

	while(*ptr != '\0'){
	
		if( isspace(*ptr) ){
			ptr++;
			continue;
		}

		if( isdigit(*ptr) ){
		
			tmp_num.push_back(*ptr - '0');

			ptr++;
			continue;
		}

		switch(*ptr){
			case 'A':
			case 'a':
				a_count = list_to_int(tmp_num);
				break;
			case 'T':
			case 't':
				t_count = list_to_int(tmp_num);
				break;
			case 'G':
			case 'g':
				g_count = list_to_int(tmp_num);
				break;
			case 'C':
			case 'c':
				c_count = list_to_int(tmp_num);
				break;
			case 'O':
			case 'o':
				other_count = list_to_int(tmp_num);

				ptr += 3; //sizeof("the")
				break;
		};

		ptr ++;
	}

	if(a_count < 0){
		throw __FILE__ ":read_sequence_EMBL: Unable to count number of A's";
	}

	if(t_count < 0){
		throw __FILE__ ":read_sequence_EMBL: Unable to count number of T's";
	}

	if(g_count < 0){
		throw __FILE__ ":read_sequence_EMBL: Unable to count number of G's";
	}

	if(c_count < 0){
		throw __FILE__ ":read_sequence_EMBL: Unable to count number of C's";
	}

	if(other_count < 0){
		throw __FILE__ ":read_sequence_EMBL: Unable to count number of \"others\"'s";
	}
	
	// Set the sequence size
	m_seq_len = a_count + t_count + g_count + c_count + other_count;

	// Allocate memory for the sequence data
	SEQPTR seq = new SEQBASE [m_seq_len + SEQ_HEADER_SIZE];
	
	if(!seq){
		throw __FILE__ ":read_sequence_EMBL: Unable to allocate memory for sequence data";
	}

	unsigned int base_count = 0;
	SEQPTR iter = seq;
	
	// Initialize the sequence header
	memcpy( iter, &m_seq_len, sizeof(unsigned int) );
	iter += sizeof(unsigned int);
	
	// Read until we hit a "/" symbol or reach the end of the file
	while( gzeof(m_fin) == 0 ){
	
		if( (gzgets(m_fin, buffer, buffer_size) == NULL) || !strip_eol(buffer, buffer_size) ){
			throw __FILE__ ":read_sequence_EMBL: Unable to read line (2)";
		}

		ptr = buffer;
		
		const unsigned int len = strlen(buffer);

		for(unsigned int i = 0;i < len;i++, ptr++){
			
			if(*ptr == '/'){
			
				// all done!
				if(base_count != m_seq_len){
					throw "Did not read enough bases!";
				}
				
				pos += i;

				gzclearerr(m_fin);
				gzseek(m_fin, pos, SEEK_SET);
				
				return seq;
			}
			
			const char local_base = toupper(*ptr);
			
			// Skip numbers and white space
			if( (local_base >= 'A') && (local_base <= 'Z') ){
				
				*(iter++) = ascii_to_hash_base(local_base);
				
				// Keep track of the number of bases read
				base_count++;
			}
		}

		// Update the current position
		pos += len;
	}

	throw __FILE__ ":read_sequence_EMBL: Could not find end-of-sequence terminator!";
	
	return NULL;
}

