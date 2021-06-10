#include "sequence_data.h"

// Needed for sort()
#include <algorithm>

using namespace std;

void sequence_data::load_gbk(const std::string &m_filename)
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
}

void sequence_data::load_embl(const std::string &m_filename)
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
}

void sequence_data::load_gff3(const std::string &m_filename)
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
}

void sequence_data::load_ptt(const std::string &m_filename)
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
}
