#include "sequence_data.h"

#include "errno.h"

// Needed for sort()
#include <algorithm>

#include <limits.h>

using namespace std;

void sequence_data::load_fasta(const std::string &m_filename)
{
	// Try to open the file as a (potentially compressed) fasta file
	fasta_in = gzopen( m_filename.c_str(), "r");

	if(fasta_in == NULL){
		throw "Unable to open fasta sequence file for reading";
	}
	
	format = FASTA_SLOW;
	
	if(seq_index.empty() == true){
	
		// Parse through the fasta file
		// 1) count the number of sequences
		// 2) the file location of each defline (the '>' symbol)

		file_index buffer_size = 0; // <-- global location in file

		char file_buffer[FASTA_BUFFER_SIZE];
		int bytes_read = 0;

		bool read_fasta = false;

		while( (bytes_read = gzread(fasta_in, file_buffer, FASTA_BUFFER_SIZE) ) > 0 ){

			for(int i = 0;i < bytes_read;++i, ++buffer_size){

				if( !read_fasta && (file_buffer[i] == '>') ){
					
					read_fasta = true;
					seq_index.push_back(buffer_size); // Save the location of this fasta header
				}
				else if(file_buffer[i] == '\n'){
					
					// Reset the read_fasta flag. This is needed to properly handle fasta headers that
					// contain multiple '>' symbols.
					read_fasta = false;
				}
			}
		}

		seq_index.push_back(buffer_size); // Save the location of the end of the fasta file

		gzrewind(fasta_in);
		
		const size_t num_seq = size();

		if(num_seq != 0){

			seq_length.resize(num_seq);

			for(size_t i = 0;i < num_seq;i++){

				// This is an overestimate of the sequence length since
				// the defline sizes are also included (as well as any white space)
				const unsigned int seq_len = seq_index[i + 1] - seq_index[i];

				seq_length[i] = make_pair(seq_len, i);
			}

			// The last sequences size is determined by its start position and the number of bytes
			// in the file
			seq_length[num_seq - 1] = make_pair(seq_index[num_seq] - seq_index[num_seq - 1], num_seq - 1);
		}
	}
}

void sequence_data::load_fastq(const std::string &m_filename)
{
	// Try to open the file as a fastq file
	fasta_in = gzopen( m_filename.c_str(), "r");

	if(fasta_in == NULL){
		throw "Unable to open fastq sequence file for reading";
	}
	
	format = FASTQ_SLOW;
	
	// If we have not been given fastq indicies (by the load_fasta_indicies() function),
	// then we need to load them now
	if(seq_index.empty() == true){
	
		// Parse through the fastq file
		// 1) count the number of sequences
		// 2) the file location of each defline (the '@' symbol)

		file_index buffer_size = 0; // <-- global location in file

		char file_buffer[FASTA_BUFFER_SIZE];
		int bytes_read = 0;
		char last_header = '-';
		bool read_eol = true;

		while( (bytes_read = gzread(fasta_in, file_buffer, FASTA_BUFFER_SIZE) ) > 0 ){

			// Since the '@' symbol can appear in quality scores, only test for the
			// '@' as the first character
			for(int i = 0;i < bytes_read;++i, ++buffer_size){

				switch(file_buffer[i]){
					case '@':
						if(read_eol){
					
							if(last_header != '+'){
								seq_index.push_back(buffer_size);
							}

							last_header = '@';
						}

						read_eol = false;
						
						break;
					case '+':
						if(read_eol){
					
							// Handle the special case of the quality
							// line header followed by a '+' score as the first
							// value on the next line
							if(last_header == '+'){
								last_header = '-';
							}
							else{
								last_header = '+';
							}
						}

						read_eol = false;
						break;
					case ' ':
					case '\t':
						// Skip spaces and tabs
						break;
					case '\n':
					case '\r':
						read_eol = true;
						break;
					default:
						if(read_eol){
							last_header = '-';
						}
						read_eol = false;
						break;
				};
			}
		}

		seq_index.push_back(buffer_size); // Save the location of the end of the fastq file
	
		gzrewind(fasta_in);
			
		const size_t num_seq = size();

		if(num_seq != 0){

			seq_length.resize(num_seq);

			for(size_t i = 0;i < num_seq;i++){

				// This is an overestimate of the sequence length since
				// the defline sizes are also included (as well as any white space)
				const unsigned int seq_len = seq_index[i + 1] - seq_index[i];

				seq_length[i] = make_pair(seq_len, i);
			}

			// The last sequences size is determined by its start position and the number of bytes
			// in the file
			seq_length[num_seq - 1] = make_pair(seq_index[num_seq] - seq_index[num_seq - 1], num_seq - 1);
		}
	}
}

unsigned int sequence_data::read_bio_seq_fasta_slow(pair<string, SEQPTR> &m_seq, 
	const unsigned int &m_index) const
{
	return read_bio_seq_fasta_slow(m_seq, m_index, 0, -1);
}
	
unsigned int sequence_data::read_bio_seq_fasta_slow(pair<string, SEQPTR> &m_seq, 
	const unsigned int &m_index, const int &m_start, const int &m_stop) const
{
	if(format != FASTA_SLOW){
		throw ":sequence_data::read_bio_seq_fasta_slow: Database is not in fasta format!";
	}
	
	if( m_index >= seq_index.size() ){
		throw ":sequence_data::read_bio_seq_fasta_slow: Index out of bounds";
	}

	// Move to the start of the fasta defline (first base will be the '>' symbol)
	//lseek(fasta_in, seq_index[m_index], SEEK_SET);
	// Note that lseek is *not* thread safe (i.e. not openMP safe). Track the
	// position in the file buffer "by hand".
	file_index file_offset = seq_index[m_index];

	// Buffer the reading of fasta records in chunks of FASTA_BUFFER_SIZE or size of the 
	// record -- whichever is smaller.
	file_index record_size = seq_index[m_index + 1] - seq_index[m_index];
	
	const file_index file_buffer_size = min(record_size, file_index(FASTA_BUFFER_SIZE) );
	
	char *file_buffer = new char[file_buffer_size];
	
	if(file_buffer == NULL){
		throw "Unable to allocate fasta file buffer";
	}
	
	file_index total_bytes_read = 0;
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	// Read the defline
	///////////////////////////////////////////////////////////////////////////////////////////////
	ssize_t bytes_read = 0;

	#pragma omp critical
	{
		// The zlib functions are not thread safe!
		gzseek(fasta_in, file_offset, SEEK_SET);

		bytes_read = gzread(fasta_in, file_buffer, file_buffer_size);
	}
	
	if(bytes_read == -1){
	
		delete [] file_buffer;
		throw "Error while reading fasta file";
	}
	
	file_offset += bytes_read;

	total_bytes_read += file_index(bytes_read);

	// The first illegal buffer location
	char* last = file_buffer + bytes_read;

	// Add an offset of 1 to skip the '>' symbol
	char *ptr = file_buffer + 1;

	// Skip any leading white spaces
	while( (last > ptr) && isspace(*ptr) ){
		ptr++;
	}

	if(ptr == last){
	
		delete [] file_buffer;
		throw "Truncated fasta file detected!";
	}	
	
	char* defline_start = ptr;

	while( (last > ptr) && (*ptr != '\n') && (*ptr != '\r') ){
		++ptr;
	}

	if(ptr == last){

		delete [] file_buffer;
		throw "Truncated fasta file detected!";
	}

	m_seq.first = string(defline_start, ptr - defline_start);

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Read the sequence
	///////////////////////////////////////////////////////////////////////////////////////////////
	// (over)-estimate the number of bases in the sequence data.
	// Add the number of bytes in the record minus the number of bytes scanned so far (as part of
	// the defline).
	file_index seq_size = SEQ_HEADER_SIZE + record_size - (ptr - file_buffer);
	
	const unsigned int start = m_start;
	
	unsigned int stop;
	
	if( (m_stop < 0) || (m_stop >= int(seq_size) ) ){
		stop = seq_size - 1;
	}
	else{
	
		stop = m_stop;
	}
	
	// Recompute the sequence size, allowing for the possibility that start > stop
	seq_size = (start > stop) ? 0 : stop - start + 1;
	seq_size += SEQ_HEADER_SIZE;
	
	// Allocate memory to hold the selected sequence. This memory must be deallocated by the calling function
	// when the sequence is no longer needed.
	m_seq.second = new SEQBASE[seq_size];
	
	if(m_seq.second == NULL){
	
		delete [] file_buffer;
		throw __FILE__ ":sequence_data::read_bio_seq_fasta_slow: Unable to allocate memory for sequence";
	}
	
	// A pointer to the num_base header that preceeds the sequence data
	unsigned int *num_base = (unsigned int*)(m_seq.second);
	
	SEQPTR seq_ptr = m_seq.second;
	
	// Initialize the sequence header
	memset( seq_ptr, 0, sizeof(unsigned int) );
	seq_ptr += sizeof(unsigned int);
	
	unsigned int index = 0;
	
	while(true){
		
		if(index > stop){
			break;
		}
		
		if(ptr == last){
			
			// We've exhausted the buffer -- can we read more?
			if(total_bytes_read >= record_size){
			
				// we're done!
				break;
			}

			#pragma omp critical
			{
				// The zlib functions are not thread safe!
				gzseek(fasta_in, file_offset, SEEK_SET);

				bytes_read = gzread( fasta_in, file_buffer, 
					min(file_buffer_size, record_size - total_bytes_read) );
			}
			
			if(bytes_read == -1){
			
				delete [] file_buffer;
				delete [] m_seq.second;
				
				throw "Error while reading fasta file";
			}
			
			file_offset += bytes_read;

			total_bytes_read += file_index(bytes_read);

			// The first illegal buffer location
			last = file_buffer + bytes_read;

			ptr = file_buffer;
		}
		
		// Ignore '*' and spaces when parsing sequences
		if( !isspace(*ptr) && (*ptr != '*') && (*ptr != '\r') && (index++ >= start) ){

			*(seq_ptr++) = ascii_to_hash_base(*ptr);
			++(*num_base);
		}
		
		++ptr;
	}
	
	// Clean up the file buffer
	delete [] file_buffer;
	
	// Return the number of bases (excluding separators and terminator). The actual size
	// of the sequence memory returned may be greater than seq_size (due to end-of-line and
	// other white space characters).
	return SEQ_SIZE(m_seq.second);
}

unsigned int sequence_data::read_bio_seq_fastq_slow(pair<string, SEQPTR> &m_seq, 
	const unsigned int &m_index) const
{
	return read_bio_seq_fastq_slow(m_seq, m_index, 0, -1);
}
	
unsigned int sequence_data::read_bio_seq_fastq_slow(pair<string, SEQPTR> &m_seq, 
	const unsigned int &m_index, const int &m_start, const int &m_stop) const
{
	if(format != FASTQ_SLOW){
		throw ":sequence_data::read_bio_seq_fastq_slow: Database is not in fastq format!";
	}
	
	if( m_index >= seq_index.size() ){
		throw ":sequence_data::read_bio_seq_fastq_slow: Index out of bounds";
	}

	// Move to the start of the fasta defline (first base will be the '>' symbol)
	//lseek(fasta_in, seq_index[m_index], SEEK_SET);
	// Note that lseek is *not* thread safe (i.e. not openMP safe). Track the
	// position in the file buffer "by hand".
	file_index file_offset = seq_index[m_index];

	// Buffer the reading of fasta records in chunks of FASTA_BUFFER_SIZE or size of the 
	// record -- whichever is smaller.
	file_index record_size = seq_index[m_index + 1] - seq_index[m_index];
	
	const file_index file_buffer_size = min(record_size, file_index(FASTA_BUFFER_SIZE) );
	
	char *file_buffer = new char[file_buffer_size];
	
	if(file_buffer == NULL){
		throw "Unable to allocate fastq file buffer";
	}
	
	file_index total_bytes_read = 0;
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	// Read the defline
	///////////////////////////////////////////////////////////////////////////////////////////////
	ssize_t bytes_read = 0;
	
	#pragma omp critical
	{
		// The zlib functions are not thread safe!
		gzseek(fasta_in, file_offset, SEEK_SET);

		bytes_read = gzread(fasta_in, file_buffer, file_buffer_size);
	}

	if(bytes_read == -1){
	
		delete [] file_buffer;
		throw "Error while reading fastq file";
	}
	
	file_offset += bytes_read;

	total_bytes_read += file_index(bytes_read);

	// The first illegal buffer location
	char* last = file_buffer + bytes_read;

	// Add an offset of 1 to skip the '@' symbol
	char *ptr = file_buffer + 1;

	// Skip any leading white spaces
	while( (last > ptr) && isspace(*ptr) ){
		ptr++;
	}

	if(ptr == last){
	
		delete [] file_buffer;
		throw "Truncated fastq file detected!";
	}	
	
	char* defline_start = ptr;

	while( (last > ptr) && (*ptr != '\n') && (*ptr != '\r') ){
		ptr++;
	}

	if(ptr == last){
	
		delete [] file_buffer;
		throw "Truncated fastq file detected!";
	}

	m_seq.first = string(defline_start, ptr - defline_start);
	
	// Skip the '\n' or '\r' character
	ptr++;
	
	if(ptr == last){
		throw "Truncated fastq file detected!";
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	// Read the sequence
	///////////////////////////////////////////////////////////////////////////////////////////////
	// (over)-estimate the number of bases in the sequence data.
	// Add the number of bytes in the record minus the number of bytes scanned so far (as part of
	// the defline).
	file_index seq_size = SEQ_HEADER_SIZE + record_size - (ptr - file_buffer);
	
	const unsigned int start = m_start;
	
	unsigned int stop;
	
	if( (m_stop < 0) || (m_stop >= int(seq_size) ) ){
		stop = seq_size - 1;
	}
	else{
	
		stop = m_stop;
	}
	
	// Recompute the sequence size, allowing for the possibility that start > stop
	seq_size = (start > stop) ? 0 : stop - start + 1;
	seq_size += SEQ_HEADER_SIZE;
	
	// Allocate memory to hold the selected sequence. This memory must be deallocated by the calling function
	// when the sequence is no longer needed.
	m_seq.second = new SEQBASE[seq_size];
	
	if(m_seq.second == NULL){
	
		delete [] file_buffer;
		throw __FILE__ ":sequence_data::read_bio_seq_fastq_slow: Unable to allocate memory for sequence";
	}
	
	unsigned int *num_base = (unsigned int*)(m_seq.second);
	
	SEQPTR seq_ptr = m_seq.second;
	
	// Initialize the sequence header
	memset( seq_ptr, 0, sizeof(unsigned int) );
	seq_ptr += sizeof(unsigned int);
	
	unsigned int index = 0;
	
	while(true){
		
		if(index > stop){
			break;
		}
		
		if(ptr == last){
			
			// We've exhausted the buffer -- can we read more?
			if(total_bytes_read >= record_size){
			
				// we're done!
				break;
			}
			
			#pragma omp critical
			{
				// The zlib functions are not thread safe!
				gzseek(fasta_in, file_offset, SEEK_SET);

				bytes_read = gzread( fasta_in, file_buffer, 
					min(file_buffer_size, record_size - total_bytes_read) );
			}
			
			if(bytes_read == -1){
			
				delete [] file_buffer;
				delete [] m_seq.second;
				
				throw "Error while reading fastq file";
			}
			
			file_offset += bytes_read;

			total_bytes_read += file_index(bytes_read);

			// The first illegal buffer location
			last = file_buffer + bytes_read;

			ptr = file_buffer;
		}
		
		// The FASTQ requires that all sequence be on a single line
		if( (*ptr == '\n') || (*ptr == '\r') ){
			break;
		}
		
		// Ignore '*' and spaces when parsing sequences
		if( !isspace(*ptr) && (*ptr != '*') && (index++ >= start) ){

			*(seq_ptr++) = ascii_to_hash_base(*ptr);
			++(*num_base);
		}
		
		++ptr;
	}
	
	// Clean up the file buffer
	delete [] file_buffer;
	
	// Return the number of bases (excluding separators and terminator). The actual size
	// of the sequence memory returned may be greater than seq_size (due to end-of-line and
	// other white space characters).
	return SEQ_SIZE(m_seq.second);
}
