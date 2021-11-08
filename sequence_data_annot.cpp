#include "sequence_data.h"

// Needed for sort()
#include <algorithm>

using namespace std;

void sequence_data::load_gbk(const std::string &m_filename)
{
	format = GBK;
	
	// Clear any existing annotation data
	mol.clear();
	size_t pos = 0;
	
	gzFile fin = gzopen(m_filename.c_str(), "r");

	// Increase the size of the internal zlib buffer used for decompression
	gzbuffer(fin, 32768);

	if(fin == NULL){
		throw __FILE__ ":sequence_data::load_gbk: Unable to open input file";
	}

	while(true){

		mol.push_back( DNAMol() );
		
		if(mol.back().load(fin, DNAMol::GBK, pos) == false){
			break;
		}
	}
	
	gzclose(fin);

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
}

void sequence_data::load_embl(const std::string &m_filename)
{
	format = EMBL;
	
	// Clear any existing annotation data
	mol.clear();
	size_t pos = 0;
	
	gzFile fin = gzopen(m_filename.c_str(), "r");

	// Increase the size of the internal zlib buffer used for decompression
	gzbuffer(fin, 32768);

	if(fin == NULL){
		throw __FILE__ ":sequence_data::load_embl: Unable to open input file";
	}

	while(true){
		
		mol.push_back( DNAMol() );
		
		if(mol.back().load(fin, DNAMol::EMBL, pos) == false){
			break;
		}
	}
	
	gzclose(fin);

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
}
