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

#ifndef __SEQUENCE_DATA
#define __SEQUENCE_DATA

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include <deque>
#include <fstream>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>

#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <sys/mman.h>
#endif // WIN32

#include <fcntl.h>

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
	
	// Note that PTT must be the last annotation type in the enumeration
	enum {FASTA_MMAP = 0, FASTA_SLOW, FASTQ_MMAP, FASTQ_SLOW, 
		NCBI, REMOTE,

		// The file formats with annotation start here
		GBK, EMBL, GFF3, PTT, NONE};
	
	private:
	
	#ifdef USE_BLAST_DB
	NCBI_NS_NCBI::CSeqDB *blast_db_ptr;
	#endif // USE_BLAST_DB
	
	unsigned char format;
	
	unsigned int read_bio_seq_fasta_mmap(std::pair<std::string, SEQPTR> &m_seq, 
		const unsigned int &m_index) const;
	
	unsigned int read_bio_seq_fasta_mmap(std::pair<std::string, SEQPTR> &m_seq, 
		const unsigned int &m_index, const int &m_start, 
		const int &m_stop) const;
		
	unsigned int read_bio_seq_fasta_slow(std::pair<std::string, SEQPTR> &m_seq, 
		const unsigned int &m_index) const;

	unsigned int read_bio_seq_fasta_slow(std::pair<std::string, SEQPTR> &m_seq, 
		const unsigned int &m_index, const int &m_start, 
		const int &m_stop) const;
	
	unsigned int read_bio_seq_fastq_mmap(std::pair<std::string, SEQPTR> &m_seq, 
		const unsigned int &m_index) const;
	
	unsigned int read_bio_seq_fastq_mmap(std::pair<std::string, SEQPTR> &m_seq, 
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
		const unsigned int &m_index, const int &m_start, 
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
	std::deque< std::pair<unsigned int, unsigned int> > seq_length;
	
	#ifdef WIN32
	HANDLE fasta_in;
	HANDLE h_map_file;
	HANDLE buffer;
	#else
	int fasta_in;
	unsigned char *buffer;
	#endif // WIN32
	
	file_index buffer_size;
	
	bool _verbose;
		
	void load_fasta(const std::string &m_filename, const bool &m_sort_by_len,
		const bool &m_allow_fasta_mmap);
	void load_fastq(const std::string &m_filename, const bool &m_sort_by_len,
		const bool &m_allow_fastq_mmap);
	void load_gbk(const std::string &m_filename, const bool &m_sort_by_len);
	void load_embl(const std::string &m_filename, const bool &m_sort_by_len);
	void load_ptt(const std::string &m_filename, const bool &m_sort_by_len);
	void load_gff3(const std::string &m_filename, const bool &m_sort_by_len);
	
	public:

	sequence_data()
	{
		format = REMOTE;
		
		#ifdef USE_BLAST_DB
		blast_db_ptr = NULL;
		#endif // USE_BLAST_DB
		
		#ifdef WIN32
		fasta_in = NULL;
		h_map_file = NULL;
		#else
		fasta_in = -1;
		#endif // WIN32

		buffer = NULL;
		buffer_size = 0;
		
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
		
		#ifdef WIN32
		if(fasta_in != NULL){

			if(h_map_file != NULL){
				UnmapViewOfFile(h_map_file);
				CloseHandle(h_map_file);
			}
			
			buffer = NULL;
			buffer_size = 0;
			CloseHandle(fasta_in);
			fasta_in = NULL;
		}
		#else
		if(fasta_in != -1){

			if(buffer != NULL){
				munmap(buffer, buffer_size);
			}
			
			buffer = NULL;
			buffer_size = 0;
			
			::close(fasta_in);
			fasta_in = -1;
		}
		#endif // WIN32
	};
	
	// The number of target sequences
	size_t size() const;
	
	// The number of target sequences taking into account "virtual"
	// targets created by target sequence fragmentation
	size_t effective_size(const unsigned int &m_max_len) const;
	
	void open(const std::string &m_filename, const bool &m_allow_fasta_mmap,
		const bool &m_sort_by_len);
	
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
		return ( (format >= GBK) && (format <= PTT) );
	};
		
	// A helper function. Does the specified file format use indicies?
	inline bool wants_indicies(const int &m_format) const
	{
		return ( (m_format == FASTA_MMAP) || (m_format == FASTA_SLOW) );
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
		// Remap the index in case we're using length sorted
		// indicies. If we're not using length sorted indicies,
		// then the value of the index will not be changed.
		const size_t local_index = length_sorted_index(m_index);
		
		if( local_index > mol.size() ){
			throw "annot: m_index > mol.size()";
		}
		
		std::list<DNAMol>::const_iterator iter = mol.begin();
		
		for(size_t i = 0;i < local_index;i++){
			iter++;
		}
		
		return *iter;
	};
	
	inline size_t length_sorted_index(const size_t &m_i) const
	{
		return ( (seq_length.size() <= m_i) ? m_i : seq_length[m_i].second);
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
