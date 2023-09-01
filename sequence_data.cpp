#include "tntblast.h"
#include "sequence_data.h"

#include <errno.h>

// WIN32 only?
#include <time.h>

#ifdef USE_BLAST_DB
#include <objects/seq/Seq_descr.hpp>
#endif // USE_BLAST_DB

using namespace std;

// Does the input string only contains digits?
inline bool is_number(const string &m_str)
{
	for(string::const_iterator i = m_str.begin();i != m_str.end();++i){
		if( !isdigit(*i) ){
			return false;
		}
	}

	return true;
}

void sequence_data::open(const string &m_filename, const vector<string> &m_blast_include,
	const vector<string> &m_blast_exclude)
{
	#ifdef USE_BLAST_DB

	// Does this filename refer to a blast database?
	try{

		blast_db_ptr = new NCBI_NS_NCBI::CSeqDB(m_filename, NCBI_NS_NCBI::CSeqDB::eNucleotide);

		if(blast_db_ptr == NULL){
        	throw __FILE__ ":sequence_data::open: Unable to allocate CSeqDB";
        }

		format = NCBI;

		deque<string> accession_include;
		deque<string> accession_exclude;

		set<NCBI_NS_NCBI::TTaxId> taxid_include;
		set<NCBI_NS_NCBI::TTaxId> taxid_exclude;

		// The m_blast_include and m_blast_exclude vectors can contain either accessions or
		// taxid values (stored as strings). If only numbers appear in the string, then we
		// assume that the value refers to a TaxId.
		for(vector<string>::const_iterator i = m_blast_include.begin();i != m_blast_include.end();++i){

			if( is_number(*i) ){
				taxid_include.insert( atoi( i->c_str() ) );
			}
			else{
				accession_include.push_back(*i);
			}
		}

		for(vector<string>::const_iterator i = m_blast_exclude.begin();i != m_blast_exclude.end();++i){

			if( is_number(*i) ){
				taxid_exclude.insert( atoi( i->c_str() ) );
			}
			else{
				accession_exclude.push_back(*i);
			}
		}

		// How can we transparantly restrict the database search by TaxIDs that are *above* the
		// level of species? See the error below:
		// 		BLAST Database error: Taxonomy ID(s) not found. This could be because the ID(s) provided are 
		// 		not at or below the species level. Please use get_species_taxids.sh to get taxids for nodes 
		// 		higher than species (see https://www.ncbi.nlm.nih.gov/books/NBK546209/).

		deque<unsigned int> oid_include;
		deque<unsigned int> oid_exclude;

		// For the specified *included* accessions, compute the corresponding OID values
		for(deque<string>::const_iterator i = accession_include.begin();i != accession_include.end();++i){

			vector<int> oid;

			try{

				blast_db_ptr->AccessionToOids(*i, oid);

				if( oid.empty() ){
					throw;
				}
			}
			catch(...){

				cerr << "Unable to find accession " << *i << " in BLAST database" << endl;
				throw "Unable to find an included accession in BLAST database";
			}

			for(vector<int>::const_iterator j = oid.begin();j != oid.end();++j){

				// Make sure we can convert from int -> unsigned int
				if(*j < 0){
					throw __FILE__ ":sequence_data::open: Invalid (negative) TaxID";
				}

				oid_include.push_back(*j);
			}
		}

		// For the specified *excluded* accessions, compute the corresponding OID values
		for(deque<string>::const_iterator i = accession_exclude.begin();i != accession_exclude.end();++i){

			vector<int> oid;

			try{

				blast_db_ptr->AccessionToOids(*i, oid);

				if( oid.empty() ){
					throw;
				}
			}
			catch(...){

				cerr << "Unable to find accession " << *i << " in BLAST database" << endl;
				throw "Unable to find an excluded accession in BLAST database";
			}

			for(vector<int>::const_iterator j = oid.begin();j != oid.end();++j){

				// Make sure we can convert from int -> unsigned int
				if(*j < 0){
					throw __FILE__ ":sequence_data::open: Invalid (negative) TaxID";
				}

				oid_exclude.push_back(*j);
			}
		}

		// For the specified *included* TaxIDs, compute the corresponding OID values
		if( !taxid_include.empty() ){

			vector< NCBI_NS_NCBI::blastdb::TOid > tmp;

			try{
				blast_db_ptr->TaxIdsToOids(taxid_include, tmp);
			}
			catch(exception &error){

				cerr << error.what() << endl;
				throw "Unable to find an included TaxID in BLAST database";
			}

			for(vector< NCBI_NS_NCBI::blastdb::TOid >::const_iterator i = tmp.begin();i != tmp.end();++i){

				// Make sure we can convert from int -> unsigned int
				if(*i < 0){
					throw __FILE__ ":sequence_data::open: Invalid (negative) TaxID";
				}

				oid_include.push_back(*i);
			}
		}

		// For the specified *excluded* TaxIDs, compute the corresponding OID values
		if( !taxid_exclude.empty() ){

			vector< NCBI_NS_NCBI::blastdb::TOid > tmp;

			try{
				blast_db_ptr->TaxIdsToOids(taxid_exclude, tmp);
			}
			catch(exception &error){

				cerr << error.what() << endl;
				throw "Unable to find an excluded TaxID in BLAST database";
			}

			for(vector< NCBI_NS_NCBI::blastdb::TOid >::const_iterator i = tmp.begin();i != tmp.end();++i){

				// Make sure we can convert from int -> unsigned int
				if(*i < 0){
					throw __FILE__ ":sequence_data::open: Invalid (negative) TaxID";
				}

				oid_exclude.push_back(*i);
			}
		}

		// Make the oid values to include and exclude unique. We could have used std::set, but
		// wanted to be prepared for the case when we have huge numbers of unique OID values
		// to process.
		sort( oid_include.begin(), oid_include.end() );
		oid_include.erase( unique( oid_include.begin(), oid_include.end() ), oid_include.end() );

		sort( oid_exclude.begin(), oid_exclude.end() );
		oid_exclude.erase( unique( oid_exclude.begin(), oid_exclude.end() ), oid_exclude.end() );

		const unsigned int num_seq = blast_db_ptr->GetNumSeqs();

		seq_length.reserve(num_seq);

		for(unsigned int i = 0;i < num_seq;i++){

			if( !oid_include.empty() ){

				deque<unsigned int>::const_iterator iter = 
					lower_bound(oid_include.begin(), oid_include.end(), i);

				if( (iter == oid_include.end() ) || (*iter != i) ){

					// Skip OIDs that are *not* in the include list
					continue;
				}
			}

			if( !oid_exclude.empty() ){

				deque<unsigned int>::const_iterator iter = 
					lower_bound(oid_exclude.begin(), oid_exclude.end(), i);

				if( (iter != oid_exclude.end() ) && (*iter == i) ){

					// Skip OIDs that *are* the exclude list
					continue;
				}
			}

			// Don't use readdb_get_sequence_length -- it's too slow on large databases
			const unsigned int seq_len = blast_db_ptr->GetSeqLengthApprox(i);

			seq_length.push_back( make_pair(seq_len, i) );
		}

		return;
	}
	catch(...){
		// This is either *not* a BLAST database, or there is an error in the BLAST database
	}
	#endif // USE_BLAST_DB

	// Is this file a recognized annotation file (i.e. GBK, EMBL, etc.)?
	switch( file_type(m_filename) ){
	
		case DNAMol::FASTA:
			load_fasta(m_filename);
			return;
		case DNAMol::FASTQ:
			load_fastq(m_filename);
			return;
		case DNAMol::GBK:
			load_gbk(m_filename);
			return;
		case DNAMol::EMBL:
			load_embl(m_filename);
			return;
		default:
			throw "File not found, unrecognized file type, or error reading BLAST database";
	};
}

unsigned int sequence_data::read_bio_seq(std::pair<std::string, SEQPTR> &m_seq, 	
	const unsigned int &m_index) const
{
	unsigned int seq_len = 0;
	
	switch(format){
		case FASTA_SLOW:
			seq_len = read_bio_seq_fasta_slow(m_seq, m_index);
			break;
		case FASTQ_SLOW:
			seq_len = read_bio_seq_fastq_slow(m_seq, m_index);
			break;
		#ifdef USE_BLAST_DB			
		case NCBI:
			seq_len = read_bio_seq_ncbi(m_seq, m_index);
		#endif // USE_BLAST_DB
			break;
		#ifdef USE_MPI
		case REMOTE:
			seq_len = read_bio_seq_remote(m_seq, m_index);
			break;
		#endif // USE_MPI
		case GBK:
		case EMBL:
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
		case FASTA_SLOW:
			seq_len = read_bio_seq_fasta_slow(m_seq, m_index, m_start, m_stop);
			break;
		case FASTQ_SLOW:
			seq_len = read_bio_seq_fastq_slow(m_seq, m_index, m_start, m_stop);
			break;
		#ifdef USE_BLAST_DB			
		case NCBI:
			seq_len = read_bio_seq_ncbi(m_seq, m_index, m_start, m_stop);
		#endif // USE_BLAST_DB
			break;
		#ifdef USE_MPI
		case REMOTE:
			seq_len = read_bio_seq_remote(m_seq, m_index, m_start, m_stop);
			break;
		#endif // USE_MPI
		case GBK:
		case EMBL:
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

	int local_buffer_size;

	MPI_Get_count(&status, MPI_BYTE, &local_buffer_size);

	if(local_buffer_size == 0){
		throw __FILE__ ":read_bio_seq_remote: Error buffer size is 0";
	}

	unsigned char* buffer = new unsigned char [local_buffer_size];

	if(buffer == NULL){
		throw __FILE__ ":read_bio_seq_remote: Error allocating receive buffer";
	}

	if(MPI_Recv(buffer, local_buffer_size, MPI_BYTE, 0 /* master */, 
		SEQ_REQUEST, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS){

		throw __FILE__ ":read_bio_seq_remote: Error receiving buffer";
	}

	unsigned char *ptr = buffer;

	// Unpack the defline
	m_seq.first = (char*)ptr;
	
	const unsigned int defline_size = m_seq.first.size() + 1;
	
	ptr += defline_size;
	
	const unsigned int seq_size = local_buffer_size - defline_size;
	
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

#ifdef USE_BLAST_DB
unsigned int sequence_data::read_bio_seq_ncbi(pair<string, SEQPTR> &m_seq, 
	const unsigned int &m_index) const
{
	return read_bio_seq_ncbi(m_seq, m_index, 0, -1);
}

unsigned int sequence_data::read_bio_seq_ncbi(pair<string, SEQPTR> &m_seq, 
	unsigned int m_index, const int &m_start, const int &m_stop) const
{

	// Map the sequence index to the OID of the BLAST database
	if( m_index > seq_length.size() ){
		throw __FILE__ ":sequence_data::read_bio_seq_ncbi: Error looking up BLAST OID";
	}

	m_index = seq_length[m_index].second;

	unsigned int seq_len = 0;
	
	if(blast_db_ptr == NULL){
		throw __FILE__ ":sequence_data::read_bio_seq_ncbi: Invalid blast_db_ptr pointer";
	}

	// Read this sequence from the database file. The "omp critical" is needed to prevent multiple openMP threads 
	// from encountering a race condition in the NCBI toolbox
	// Update August 31, 2023: This race condition may be related to the NCBI toolbox race condition mentioned in 
	// https://www.repeatmasker.org/ (see RMBlast 2.6.0 BUGFIX). Rather than requiring users to patch the NCBI toolbox
	// code, it is easier to protect access with an OpenMP critical section (as long as most of the compute time is spent
	// aligning sequences).
	#pragma omp critical (NcbiToolbox)
	{
		// There is at least one example of a blast 'nt' database (from May 10, 2021) for which at least
		// one sequence index throws an error when we attempt to load the CBioseq. If this happens, return
		// a zero length sequence.
		try{

			ncbi::CRef<ncbi::objects::CBioseq> bs = blast_db_ptr->GetBioseqNoData(m_index);

			string title;
			string accession;

			const list< ncbi::CRef< ncbi::CSeq_id > >& seqIds = bs->GetId();
			
			for (list< ncbi::CRef< ncbi::CSeq_id > >::const_iterator cit = seqIds.begin();cit != seqIds.end(); cit++)
			{
				const ncbi::CTextseq_id* textId = (*cit)->GetTextseq_Id();

				if(textId && textId->CanGetAccession() ){

					if( textId->CanGetAccession() ){
						accession = textId->GetAccession();
					}

					if( !accession.empty() ){
						break;
					}
				}
			}

			ITERATE( ncbi::CSeq_descr::Tdata, desc, bs->GetDescr().Get() ) {
				if( (*desc)->Which() == ncbi::CSeqdesc::e_Title ) {
					title = (*desc)->GetTitle();
					break;
				}
			}

			if( accession.empty() ){
				m_seq.first = title;
			}
			else{

				if( title.empty() ){
					m_seq.first = accession;
				}
				else{
					m_seq.first = accession + " " + title;
				}
			}

			// Extract the sequence data
			const char* seq_buffer = NULL;
			const int start = m_start;

			// The stop value that we pass to GetAmbigSeq is not
			// included in the extracted bases, i.e. [start, stop).
			int stop = min( m_stop + 1, blast_db_ptr->GetSeqLength(m_index) );
			
			if(start > stop){

				stop = blast_db_ptr->GetAmbigSeq(m_index,
					&seq_buffer,
					ncbi::kSeqDBNuclNcbiNA8);
			}
			else{
				blast_db_ptr->GetAmbigSeq(m_index,
					&seq_buffer,
					ncbi::kSeqDBNuclNcbiNA8,
					start, stop);
			}

			seq_len = stop - start;

			// Allocate space to hold the sub-sequence
			m_seq.second = new SEQBASE [seq_len + SEQ_HEADER_SIZE];
			
			if(m_seq.second == NULL){
				throw __FILE__ ":sequence_data::read_bio_seq_annot: Error allocating SEQBASE buffer";
			}
			
			SEQPTR out_ptr = m_seq.second;
			
			// Initialize the sequence header
			memcpy( out_ptr, &seq_len, sizeof(unsigned int) );
			out_ptr += sizeof(unsigned int);
			
			const char *in_ptr = seq_buffer;

			for(int i = start;i < stop;++i,++out_ptr,++in_ptr){

				// kSeqDBNuclNcbiNA8
				#define BLAST_DB_A	1
				#define BLAST_DB_C 	2
				#define BLAST_DB_G	4
				#define BLAST_DB_T	8

				switch(*in_ptr){
					case BLAST_DB_A:
						*out_ptr = DB_A;
						break;
					case BLAST_DB_T:
						*out_ptr = DB_T;
						break;
					case BLAST_DB_G:
						*out_ptr = DB_G;
						break;
					case BLAST_DB_C:
						*out_ptr = DB_C;
						break;
					// G or T or C
					case (BLAST_DB_G | BLAST_DB_T | BLAST_DB_C):
						*out_ptr = DB_B;
						break;
					// G or A or T
					case (BLAST_DB_G | BLAST_DB_A | BLAST_DB_T):
						*out_ptr = DB_D;
						break;
					// A or C or T
					case (BLAST_DB_A | BLAST_DB_C | BLAST_DB_T):
						*out_ptr = DB_H;
						break;
					// G or T
					case (BLAST_DB_G | BLAST_DB_T):
						*out_ptr = DB_K;
						break;
					// A or C
					case (BLAST_DB_A | BLAST_DB_C):
						*out_ptr = DB_M;
						break;
					// A or C or G or T
					case (BLAST_DB_A | BLAST_DB_C | BLAST_DB_G | BLAST_DB_T):
						*out_ptr = DB_N;
						break;
					// G or A
					case (BLAST_DB_G | BLAST_DB_A):
						*out_ptr = DB_R;
						break;
					// G or C
					case (BLAST_DB_G | BLAST_DB_C):
						*out_ptr = DB_S;
						break;
					// G or C or A
					case (BLAST_DB_G | BLAST_DB_C | BLAST_DB_A):
						*out_ptr = DB_V;
						break;
					// A or T
					case (BLAST_DB_A | BLAST_DB_T):
						*out_ptr = DB_W;
						break;
					// T or C
					case (BLAST_DB_T | BLAST_DB_C):
						*out_ptr = DB_Y;
						break;
					default:
						*out_ptr = DB_UNKNOWN;
						break;
				};
			}

			blast_db_ptr->RetAmbigSeq(&seq_buffer);
		}
		catch(...){

			m_seq.first = "Invalid";
			m_seq.second = NULL;
		}
	}
	
	return seq_len;
}
#endif // USE_BLAST_DB

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
		case FASTA_SLOW:
		case FASTQ_SLOW:
			// If we have read any data, |seq_index| = num_seq + 1.
			// If have have yet to read data, then |seq_index| == 0.
			return seq_index.empty() ? 0 : seq_index.size() - 1;
		case NCBI:
			return seq_length.size();
		case REMOTE:
			throw __FILE__ ":sequence_data::size: size is not defined for remote databases";
		case GBK:
		case EMBL:
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
