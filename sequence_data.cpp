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

// For compatibility with windows, we must include the
// tntblast.h header before we include the sequence_data.h header.
// This insures that the proper version of std::min and std::max is
// used in the hybrid_sig.h file (which will otherwise be confused by
// the min and max macros in windows.h that is included by sequence_data.h).
#include "tntblast.h"
#include "sequence_data.h"

#include <errno.h>

// WIN32 only?
#include <time.h>

using namespace std;

#ifdef USE_NCBI
unsigned int bioseq_to_SEQPTR(SEQPTR &m_seq, const BioseqPtr &m_bsp, const int &m_start, const int &m_stop);
#endif // USE_NCBI

void sequence_data::open(const string &m_filename, const bool &m_allow_fasta_mmap,
	const bool &m_sort_by_len)
{
	#ifdef USE_NCBI

	// Does this filename refer to a blast database?
	rdbfp = readdb_new( (char*)m_filename.c_str(), false /* Not a protein database*/ );
			
	// If this is a valid pointer, return right away
	if(rdbfp){
	
		format = NCBI;		
			
		const size_t num_seq = size();

		seq_length.resize(num_seq);

		for(size_t i = 0;i < num_seq;i++){

			// Don't use readdb_get_sequence_length -- it's too slow on large databases
			const unsigned int seq_len = readdb_get_sequence_length_approx(rdbfp, i);

			seq_length[i] = make_pair(seq_len, i);
		}
		
		if(m_sort_by_len == true){
			sort( seq_length.begin(), seq_length.end(), sequence_order() );
		}
		
		return;
	}
	#endif // USE_NCBI
	
	// Is this file a recognized annotation file (i.e. GBK, ASN.1, etc.)?
	switch( file_type(m_filename) ){
	
		case DNAMol::FASTA:
			load_fasta(m_filename, m_sort_by_len, m_allow_fasta_mmap);
			return;
		case DNAMol::FASTQ:
			load_fastq(m_filename, m_sort_by_len, m_allow_fasta_mmap);
			return;
		case DNAMol::ASN_1:
			#ifdef USE_NCBI
			load_asn(m_filename, m_sort_by_len);
			#else
			throw "NCBI ASN.1 input file detected. Please compile thermonucleotideBLAST with the NCBI toolkit";
			#endif // USE_NCBI
			return;
		case DNAMol::GBK:
			load_gbk(m_filename, m_sort_by_len);
			return;
		case DNAMol::EMBL:
			load_embl(m_filename, m_sort_by_len);
			return;
		case DNAMol::GFF3:
			load_gff3(m_filename, m_sort_by_len);
			return;
		case DNAMol::PTT:
			throw "The PTT file format is not currently supported";
			//load_ptt(m_filename, m_sort_by_len);
			return;
		default:
			throw "File not found or unrecognized file type (not fasta, gbk, asn.1, ptt or embl)";
	};
}

unsigned int sequence_data::read_bio_seq(std::pair<std::string, SEQPTR> &m_seq, 	
	const unsigned int &m_index) const
{

	unsigned int seq_len = 0;
	
	switch(format){
		case FASTA_MMAP:
			seq_len = read_bio_seq_fasta_mmap(m_seq, m_index);
			break;
		case FASTA_SLOW:
			seq_len = read_bio_seq_fasta_slow(m_seq, m_index);
			break;
		case FASTQ_MMAP:
			seq_len = read_bio_seq_fastq_mmap(m_seq, m_index);
			break;
		case FASTQ_SLOW:
			seq_len = read_bio_seq_fastq_slow(m_seq, m_index);
			break;
		#ifdef USE_NCBI			
		case NCBI:
			seq_len = read_bio_seq_ncbi(m_seq, m_index);
		#endif // USE_NCBI
			break;
		#ifdef USE_MPI
		case REMOTE:
			seq_len = read_bio_seq_remote(m_seq, m_index);
			break;
		#endif // USE_MPI
		case ASN_1:
		case GBK:
		case EMBL:
		case GFF3:
		case PTT:
			seq_len = read_bio_seq_annot(m_seq, m_index);
			break;
		default:
			throw __FILE__ ":sequence_data::read_bio_seq: Unknown database format";
	};
	
	return seq_len;
}

unsigned int sequence_data::read_bio_seq(std::pair<std::string, SEQPTR> &m_seq, 	
	const unsigned int &m_index, const unsigned int &m_start, 
	const unsigned int &m_stop) const
{

	unsigned int seq_len = 0;
	
	switch(format){
		case FASTA_MMAP:
			seq_len = read_bio_seq_fasta_mmap(m_seq, m_index, m_start, m_stop);
			break;
		case FASTA_SLOW:
			seq_len = read_bio_seq_fasta_slow(m_seq, m_index, m_start, m_stop);
			break;
		case FASTQ_MMAP:
			seq_len = read_bio_seq_fastq_mmap(m_seq, m_index, m_start, m_stop);
			break;
		case FASTQ_SLOW:
			seq_len = read_bio_seq_fastq_slow(m_seq, m_index, m_start, m_stop);
			break;
		#ifdef USE_NCBI			
		case NCBI:
			seq_len = read_bio_seq_ncbi(m_seq, m_index, m_start, m_stop);
		#endif // USE_NCBI
			break;
		#ifdef USE_MPI
		case REMOTE:
			seq_len = read_bio_seq_remote(m_seq, m_index, m_start, m_stop);
			break;
		#endif // USE_MPI
		case ASN_1:
		case GBK:
		case EMBL:
		case GFF3:
		case PTT:
			seq_len = read_bio_seq_annot(m_seq, m_index, m_start, m_stop);
			break;
		default:
			throw __FILE__ ":sequence_data::read_bio_seq: Unknown database format";
	};
	
	return seq_len;
}

#ifdef USE_MPI
unsigned int sequence_data::read_bio_seq_remote(pair<string, SEQPTR> &m_seq, 
	const unsigned int &m_index) const
{
	return read_bio_seq_remote(m_seq, m_index, 0, -1);
}

unsigned int sequence_data::read_bio_seq_remote(pair<string, SEQPTR> &m_seq, 
	const unsigned int &m_index, const int &m_start, const int &m_stop) const
{
	MPI_Status status;
		
	// Request this sequence from the master
	if(MPI_Send( (unsigned int*)&m_index, 1, MPI_UNSIGNED, 0 /* master */, 
	   SEQ_REQUEST, MPI_COMM_WORLD ) != MPI_SUCCESS){

	   throw __FILE__ ":read_bio_seq_remote: Error sending SEQ_REQUEST";
	}

	int range[2];
	
	range[0] = m_start;
	range[1] = m_stop;
	
	if(MPI_Send( range, 2, MPI_INT, 0 /* master */, 
	   SEQ_REQUEST, MPI_COMM_WORLD ) != MPI_SUCCESS){

	   throw __FILE__ ":read_bio_seq_remote: Error sending SEQ_REQUEST";
	}
	
	MPI_Probe(0 /* master */, SEQ_REQUEST, MPI_COMM_WORLD, &status);

	int buffer_size;

	MPI_Get_count(&status, MPI_BYTE, &buffer_size);

	if(buffer_size == 0){
		throw __FILE__ ":read_bio_seq_remote: Error buffer size is 0";
	}

	unsigned char* buffer = new unsigned char [buffer_size];

	if(buffer == NULL){
		throw __FILE__ ":read_bio_seq_remote: Error allocating receive buffer";
	}

	if(MPI_Recv(buffer, buffer_size, MPI_BYTE, 0 /* master */, 
		SEQ_REQUEST, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS){

		throw __FILE__ ":read_bio_seq_remote: Error receiving buffer";
	}

	unsigned char *ptr = buffer;

	// Unpack the defline
	m_seq.first = (char*)ptr;
	
	const unsigned int defline_size = m_seq.first.size() + 1;
	
	ptr += defline_size;
	
	const unsigned int seq_size = buffer_size - defline_size;
	
	m_seq.second = new SEQBASE [seq_size];

	if(m_seq.second == NULL){
		throw __FILE__ ":read_bio_seq_remote: Error allocating SEQBASE buffer";
	}

	memcpy(m_seq.second, ptr, seq_size);

	// Clean up the send buffer. The sequence buffer must be deallocated by
	// the calling function
	delete [] buffer;
	
	// Return the number of bases (excluding separators and terminator)  
	return SEQ_SIZE(m_seq.second);

}
#endif // USE_MPI

#ifdef USE_NCBI
unsigned int sequence_data::read_bio_seq_ncbi(pair<string, SEQPTR> &m_seq, 
	const unsigned int &m_index) const
{
	return read_bio_seq_ncbi(m_seq, m_index, 0, -1);
}

unsigned int sequence_data::read_bio_seq_ncbi(pair<string, SEQPTR> &m_seq, 
	const unsigned int &m_index, const int &m_start, const int &m_stop) const
{
	BioseqPtr bsp = NULL;
	
	// Read this sequence from the database file. The "omp critical"
	// is needed to prevent multiple openMP threads from encountering
	// a race condition in the NCBI toolbox
	#pragma omp critical (NcbiToolbox)
	{
		
		// Note that neither readdb_get_bioseq nor readdb_get_bioseq_ex is threadsafe
		// (they both call functions like SeqLocAsnLoad() that modify global library variables
		// which can lead to race conditions). However, by using readdb_get_bioseq_ex, we can
		// explicitly avoid the overhead associated with the NCBI ObjMgr code and also avoid
		// the need for an additional critical section when deleting the bsp-associated memory.
		
		//bsp = readdb_get_bioseq(rdbfp, m_index);
		bsp = readdb_get_bioseq_ex(rdbfp, m_index, FALSE /* don't use the ObjMgr */, FALSE);
		
		if(bsp == NULL){
			throw __FILE__ ":read_bio_seq_ncbi: NULL bsp";
		}
	
		// Save the description string as a defline so we can tell the 
		// user something about this sequence
		char def_buffer[DEFLINE_BUFFER_SIZE];

		// Make sure that we always have a valid string	
		def_buffer[DEFLINE_BUFFER_SIZE - 1] = '\0';
	
		// Pass the buffer length as def_len - 1 to insure that our '\0'
		// character does not get clobbered!
		BioseqLabel(bsp, def_buffer, DEFLINE_BUFFER_SIZE - 1, OM_LABEL_CONTENT);

		m_seq.first = def_buffer;

		// Add a separator
		def_buffer[0] = ' ';

		// Pass the buffer length as def_len - 1 - 1 to insure that our '\0'
		// character does not get clobbered!
		SeqDescLabel(bsp->descr, def_buffer + 1, DEFLINE_BUFFER_SIZE - 2, OM_LABEL_CONTENT);
		
		m_seq.first += def_buffer;
	}
	
	// Convert this BioseqPtr into a SEQPTR buffer
	const unsigned int seq_len = bioseq_to_SEQPTR(m_seq.second, bsp, m_start, m_stop);
		
	// By using readdb_get_bioseq_ex (see above), we can explicitly avoid the NCBI ObjMgr
	// (which is not thread safe). This allows us to free the bsp associated memory
	// locally in each thread (and avoid an additional omp critical section).
	//#pragma omp critical (NcbiToolbox)
	//{
		//bsp = BioseqFree(bsp);
		SeqIdSetFree(bsp->id);
		BioseqFreeComponents(bsp);
		MemFree(bsp);
	//}
	
	return seq_len;
}
#endif // USE_NCBI

unsigned int sequence_data::read_bio_seq_annot(pair<string, SEQPTR> &m_seq, 
		const unsigned int &m_index) const
{
	return read_bio_seq_annot(m_seq, m_index, 0, -1);
}
		
unsigned int sequence_data::read_bio_seq_annot(pair<string, SEQPTR> &m_seq, 
		const unsigned int &m_index, const int &m_start, const int &m_stop) const
{
	if( m_index >= mol.size() ){
		throw __FILE__ ":read_bio_seq_annot: index out of bounds";
	}
	
	list<DNAMol>::const_iterator iter = mol.begin();

	for(unsigned int i = 0;i < m_index;i++){
		iter++;
	}
	
	m_seq.first = iter->seq_id_str() + " " + iter->info(DNAMol::TAXA_NAME);
		
	// Old version
	//m_seq.second = iter->sequence();
	unsigned int seq_len = iter->num_bases();
	
	const unsigned int start = m_start;
	
	unsigned int stop;
	
	if( (m_stop < 0) || (m_stop >= int(seq_len) ) ){
		stop = seq_len - 1;
	}
	else{
	
		stop = m_stop;
	}
	
	// Recompute the sequence size, allowing for the possibility that start > stop
	seq_len = (start > stop) ? 0 : stop - start + 1;
		
	// Allocate space to hold the sub-sequence
	m_seq.second = new SEQBASE [seq_len + SEQ_HEADER_SIZE];
	
	if(m_seq.second == NULL){
		throw __FILE__ ":sequence_data::read_bio_seq_annot: Error allocating SEQBASE buffer";
	}
	
	
	SEQPTR ptr = m_seq.second;
	
	// Initialize the sequence header
	memcpy( ptr, &seq_len, sizeof(unsigned int) );
	ptr += sizeof(unsigned int);
	
	if(seq_len > 0){
	
		memcpy( ptr, SEQ_START( iter->sequence() ) + start, seq_len*sizeof(unsigned char) );
		ptr += seq_len*sizeof(unsigned char);
	}
	
	return seq_len;
}

size_t sequence_data::size() const
{
	// Note that we can't just return seq_length.size() here, since this
	// function is called before seq_length is initialized!
	
	switch(format){
		case FASTA_MMAP:
		case FASTA_SLOW:
		case FASTQ_MMAP:
		case FASTQ_SLOW:
			return seq_index.size();
		case NCBI:
			#ifdef USE_NCBI
			{
				if(!rdbfp){
					throw __FILE__ ":sequence_data::size: NCBI database is not open";
				}
				
				Int8 total_len = 0L;
				Int4 num_seq = 0;

				// Get the total length of all sequences and the number of sequences
				readdb_get_totals(rdbfp, &total_len, &num_seq);
				
				return size_t(num_seq);
			}
			#endif // USE_NCBI
			
			return 0;
		case REMOTE:
			throw __FILE__ ":sequence_data::size: size is not defined for remote databases";
		case ASN_1:
		case GBK:
		case EMBL:
		case PTT:
			return mol.size();
		default:
			throw __FILE__ ":sequence_data::size: size is not defined for uninitialized databases";
	};
}

size_t sequence_data::effective_size(const unsigned int &m_max_len) const
{
	
	const size_t len = size();
	size_t effective_len = 0;
	
	for(size_t i = 0;i < len;i++){
		
		if(m_max_len >= seq_length[i].first){
		
			effective_len ++;
			continue;
		}
		
		// If we get here, then m_max_len < seq_length[i].first
		effective_len += seq_len_increment(seq_length[i].first, m_max_len).second;
	}
	
	return effective_len;
}

pair<unsigned int, unsigned int> seq_len_increment(const unsigned int &m_len, const unsigned int &m_max_len)
{
	if(m_len <= m_max_len){
		return make_pair(m_len - 1, 1);
	}
	
	// Find the smallest number of paritions, n, such that m_len/n <= m_max_len
	unsigned int n = 2;
	
	while(m_len > n*m_max_len){
		n++;
	}
	
	// If there's a remainder, add one to the length
	return make_pair(m_len/n + ( (m_len%n != 0) ? 1 : 0), n);
}

#ifdef USE_NCBI
unsigned int bioseq_to_SEQPTR(SEQPTR &m_seq, const BioseqPtr &m_bsp, const int &m_start, const int &m_stop)
{
	unsigned int seq_len = BioseqGetLen(m_bsp);
	
	const unsigned int start = m_start;
	
	unsigned int stop;
	
	if( (m_stop < 0) || (m_stop >= int(seq_len) ) ){
		stop = seq_len - 1;
	}
	else{
	
		stop = m_stop;
	}
	
	// Recompute the sequence size, allowing for the possibility that start > stop
	seq_len = (start > stop) ? 0 : stop - start + 1;
	
	SeqPortPtr spp = NULL;
	Uint1 residue;
	
	// Allocate space for the new sequence
	m_seq = new SEQBASE [seq_len + SEQ_HEADER_SIZE];
	
	if(m_seq == NULL){
		throw __FILE__ ":bioseq_to_SEQPTR: Error allocating SEQBASE buffer";
	}

	SEQPTR iter = m_seq;
	
	// Initialize the sequence header
	memcpy( iter, &seq_len, sizeof(unsigned int) );
	iter += sizeof(unsigned int);
	
	if(seq_len > 0){
	
		// Save the sequence data
		spp = SeqPortNew(m_bsp, start, stop, Seq_strand_plus, Seq_code_iupacna);

		while( (residue = SeqPortGetResidue(spp) ) != SEQPORT_EOF){

			switch(residue)
			{
				case 'A': case 'a':
					*iter = DB_A;
					break;
				case 'T': case 't':
				case 'U': case 'u': // Map RNA -> DNA
					*iter = DB_T;
					break;
				case 'G': case 'g':
					*iter = DB_G;
					break;
				case 'C': case 'c':
					*iter = DB_C;
					break;
				// G or T or C
				case 'B':	case 'b':
					*iter = DB_B;
					break;
				// G or A or T
				case 'D':	case 'd':
					*iter = DB_D;
					break;
				// A or C or T
				case 'H':	case 'h':
					*iter = DB_H;
					break;
				// G or T
				case 'K':	case 'k':
					*iter = DB_K;
					break;
				// A or C
				case 'M':	case 'm':
					*iter = DB_M;
					break;
				// A or C or G or T
				case 'N':	case 'n':
					*iter = DB_N;
					break;
				// G or A
				case 'R':	case 'r':
					*iter = DB_R;
					break;
				// G or C
				case 'S':	case 's':
					*iter = DB_S;
					break;
				// G or C or A
				case 'V':	case 'v':
					*iter = DB_V;
					break;
				// A or T
				case 'W':	case 'w':
					*iter = DB_W;
					break;
				// T or C
				case 'Y':	case 'y':
					*iter = DB_Y;
					break;
				case SEQPORT_VIRT:
					*iter = DB_GAP;
					break;
				default:
					*iter = DB_UNKNOWN;
					break;
			};

			iter++;
		}
		
		spp = SeqPortFree(spp);
	}
	
	return seq_len;
}

#endif // USE_NCBI
