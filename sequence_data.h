#ifndef __SEQUENCE_DATA
#define __SEQUENCE_DATA

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include <deque>
#include <fstream>
#include <string>

#include <zlib.h>

#include "seq_hash.h"

#include "annotation.h"

#ifdef USE_BLAST_DB
#include <objtools/blast/seqdb_reader/seqdb.hpp>
#endif // USE_BLAST_DB

// Use a 512 kb buffer for reading fasta files from disk
// (only applies when we're not memory mapping files). Note that
// this buffer size should never be less than the maximum allowed
// defline size (in bytes)
#define	FASTA_BUFFER_SIZE	524288

typedef unsigned long long int file_index;

class sequence_data
{
	public:
	
	enum {FASTA_SLOW, FASTQ_SLOW, 
		NCBI, REMOTE,

		// The file formats with annotation start here
		GBK, EMBL, NONE};
	
	private:
	
	#ifdef USE_BLAST_DB
	NCBI_NS_NCBI::CSeqDB *blast_db_ptr;
	#endif // USE_BLAST_DB
	
	unsigned char format;
		
	unsigned int read_bio_seq_fasta_slow(std::pair<std::string, SEQPTR> &m_seq, 
		const unsigned int &m_index) const;

	unsigned int read_bio_seq_fasta_slow(std::pair<std::string, SEQPTR> &m_seq, 
		const unsigned int &m_index, const int &m_start, 
		const int &m_stop) const;
		
	unsigned int read_bio_seq_fastq_slow(std::pair<std::string, SEQPTR> &m_seq, 
		const unsigned int &m_index) const;

	unsigned int read_bio_seq_fastq_slow(std::pair<std::string, SEQPTR> &m_seq, 
		const unsigned int &m_index, const int &m_start, 
		const int &m_stop) const;
		
	unsigned int read_bio_seq_annot(std::pair<std::string, SEQPTR> &m_seq, 
		const unsigned int &m_index) const;
	
	unsigned int read_bio_seq_annot(std::pair<std::string, SEQPTR> &m_seq, 
		const unsigned int &m_index, const int &m_start, 
		const int &m_stop) const;
		
	#ifdef USE_BLAST_DB
	unsigned int read_bio_seq_ncbi(std::pair<std::string, SEQPTR> &m_seq,
		const unsigned int &m_index) const;
		
	unsigned int read_bio_seq_ncbi(std::pair<std::string, SEQPTR> &m_seq,
		unsigned int m_index, const int &m_start, 
		const int &m_stop) const;
	#endif // USE_BLAST_DB
	
	#ifdef USE_MPI
	unsigned int read_bio_seq_remote(std::pair<std::string, SEQPTR> &m_seq, 
		const unsigned int &m_index) const;
		
	unsigned int read_bio_seq_remote(std::pair<std::string, SEQPTR> &m_seq, 
		const unsigned int &m_index, const int &m_start, 
		const int &m_stop) const;
	#endif // USE_MPI
	
	// If annotations have been provided, this vector of DNAMol 
	// stores all sequence and associated annotation information.
	// Note that a DNAMol does *not* have a valid copy constructor
	// (so we can not use a vector).
	std::list<DNAMol> mol;
	
	// Map a sequence number to the location in a file
	std::deque<file_index> seq_index;
	
	// Store (sequence length, index) pairs, sorted by length
	std::vector< std::pair<unsigned int, unsigned int> > seq_length;
	
	gzFile fasta_in;
	
	bool _verbose;
		
	void load_fasta(const std::string &m_filename);
	void load_fastq(const std::string &m_filename);
	void load_gbk(const std::string &m_filename);
	void load_embl(const std::string &m_filename);
	
	public:

	sequence_data()
	{
		format = REMOTE;
		
		#ifdef USE_BLAST_DB
		blast_db_ptr = NULL;
		#endif // USE_BLAST_DB
		
		fasta_in = NULL;
		
		_verbose = true;
	};
	
	~sequence_data()
	{
		close();
	};
	
	void close()
	{
		
		// Purge any existing sequence indicies
		seq_index.clear();
		
		#ifdef USE_BLAST_DB
		if(blast_db_ptr != NULL){

			delete blast_db_ptr;
			blast_db_ptr = NULL;
		}
		#endif // USE_BLAST_DB
		
		if(fasta_in != NULL){
			
			gzclose(fasta_in);
			fasta_in = NULL;
		}
	};
	
	// The number of target sequences
	size_t size() const;
	
	// The number of target sequences taking into account "virtual"
	// targets created by target sequence fragmentation
	size_t effective_size(const unsigned int &m_max_len) const;
	
	void open(const std::string &m_filename, const std::vector<std::string> &m_blast_include,
		const std::vector<std::string> &m_blast_exclude);
	
	unsigned int read_bio_seq(std::pair<std::string, SEQPTR> &m_seq, const unsigned int &m_index) const;
		
	unsigned int read_bio_seq(std::pair<std::string, SEQPTR> &m_seq, const unsigned int &m_index,
		const unsigned int &m_start, const unsigned int &m_stop) const;
	
	inline int file_format() const
	{
		return format;
	};
	
	// Does the file format contain annotation information?
	inline bool is_annot_format() const
	{
		return ( (format == GBK) || (format == EMBL) );
	};
		
	// A helper function. Does the specified file format use indicies?
	inline bool wants_indicies(const int &m_format) const
	{
		return (m_format == FASTA_SLOW);
	};
	
	// Use previously computed indicies for loading sequence records from fasta files
	void indicies(const std::deque<file_index> &m_indicies)
	{
		seq_index = m_indicies;
	};
	
	// Return the current set of sequence record indicies
	inline const std::deque<file_index>& indicies() const
	{
		return seq_index;
	};
	
	inline void verbose(const bool &m_verbose)
	{
		_verbose = m_verbose;
	};
	
	const DNAMol& annot(const size_t &m_index) const
	{
		if( m_index > mol.size() ){
			throw "annot: m_index > mol.size()";
		}
		
		std::list<DNAMol>::const_iterator iter = mol.begin();
		
		for(size_t i = 0;i < m_index;i++){
			iter++;
		}
		
		return *iter;
	};
	
	// Return the *approximate* size (in bases) of the m_i^th sequence or zero
	// if there are no sequences in the database
	inline size_t approx_seq_len(const size_t &m_i) const
	{
		return ( (seq_length.size() <= m_i) ? 0 : seq_length[m_i].first);
	};
};

struct sequence_order
{
	inline bool operator()(const std::pair<unsigned int, size_t> &m_a, const std::pair<unsigned int, size_t> &m_b) const
	{
		// Descending order
		//return m_a.first > m_b.first;
		
		// Ascending order
		return m_a.first < m_b.first;
	};
};

std::pair<unsigned int, unsigned int> seq_len_increment(const unsigned int &m_len, 
	const unsigned int &m_max_len);

#endif // __SEQUENCE_DATA
