#include "sequence_data.h"

#include "errno.h"

// Needed for sort()
#include <algorithm>

#include <limits.h>

using namespace std;

void sequence_data::load_fasta(const std::string &m_filename, const bool &m_sort_by_len,
	const bool &m_allow_fasta_mmap)
{
	bool allow_fasta_mmap = m_allow_fasta_mmap;
	
	// Try to open the file as a fasta file
	#ifdef WIN32

	// Note when compiling: for recent versions of Visual Studio, you will need
	// to set the option "Configuration Properties->General->Project Defaults->Character set"
	// to "Use Multi-Byte charater set" instead of the default "Use Unicode character set".
	// Otherwise, CreateFile() will expect Unicode strings (instead of char*)
	fasta_in = CreateFile(
		m_filename.c_str(),
		GENERIC_READ,
		FILE_SHARE_READ,
		NULL, // security descriptor
		OPEN_EXISTING,
		FILE_ATTRIBUTE_READONLY,
		NULL);

	if(fasta_in == INVALID_HANDLE_VALUE){
		throw "Unable to open fasta sequence file for reading";
	}

	#else
	fasta_in = ::open( m_filename.c_str(), O_RDONLY);

	if(fasta_in == -1){
		throw "Unable to open fasta sequence file for reading";
	}

	#endif// WIN32
	
	struct stat file_info;
	
	if(stat(m_filename.c_str(), &file_info) != 0){
		throw "Unable to read fasta file size";
	}
	
	buffer_size = file_info.st_size;
	
	// Don't try to memory map files larger than 2^32 - 1
	if(buffer_size > UINT_MAX){
		allow_fasta_mmap = false;
	}
	
	if(allow_fasta_mmap == true){
	
		format = FASTA_MMAP;
		
		#ifdef WIN32
		h_map_file = CreateFileMapping(fasta_in,
				NULL, // If the security descriptor is null, this handle can not be inherited
				PAGE_READONLY,
				0, 0,
				NULL);

		if(h_map_file == NULL){
			
			if(_verbose){
				cerr << "Unable to memory map the fasta file -- insufficient memory" << endl;
			}
			
			format = FASTA_SLOW;
		}
		else{

			buffer = MapViewOfFile(h_map_file,
				FILE_MAP_READ,
				0, // high order bits of offset
        			0, // low order bits of offset
        			0); // Bytes to map

			if(buffer == NULL){
				
				if(_verbose){
					cerr << "Unable to memory map the fasta file -- insufficient memory" << endl;
				}
				
				format = FASTA_SLOW;
			}
		}

		#else

		buffer = (unsigned char*)mmap(0, buffer_size, PROT_READ, MAP_SHARED, fasta_in, 0);

		if(buffer == MAP_FAILED){

			switch(errno){
				case EACCES:
					throw "Specified fasta file is a non-regular file -- please check for filename typos";
				case EAGAIN:
					throw "Specified fasta file has been locked";
				case ENOMEM:
					
					if(_verbose){
						cerr << "Unable to memory map the fasta file -- insufficient memory" << endl;
					}
					
					format = FASTA_SLOW;
					buffer = NULL;
					break;
				case ENODEV:
				
					if(_verbose){
						cerr << "File system does not support memory mapping" << endl;
					}
					
					format = FASTA_SLOW;
					buffer = NULL;
					break;
				default:
				
					if(_verbose){
						cerr << "Unable to memory map the fasta file" << endl;
					}
					
					format = FASTA_SLOW;
					buffer = NULL;
					break;
			};
			
			if(_verbose){
				cerr << "Fasta file will be read using standard disk access" << endl;
			}
		}
		#endif // WIN32
	}
	else{
		format = FASTA_SLOW;
	}
	
	// If we have not been given fasta indicies (by the load_fasta_indicies() function),
	// then we need to load them now
	if(seq_index.empty() == true){
	
		// Parse through the fasta file
		// 1) count the number of sequences
		// 2) the file location of each defline (the '>' symbol)

		file_index i;

		if(format == FASTA_MMAP){

			// Make sure that multiple '>' symbols on the same line don't derail
			// the parser (this has been observed for sequences from the RDP).
			bool read_eol = true;
			
			for(i = 0;i < buffer_size;i++){

				switch( ( (char*)(buffer) )[i] ){
					case '>':
						if(read_eol){
							seq_index.push_back(i);
							read_eol = false;
						}
						break;
					case '\n':
						read_eol = true;
						break;
				};				
			}
		}
		else{ // format == FASTA_SLOW

			char *file_buffer = new char [FASTA_BUFFER_SIZE];

			if(file_buffer == NULL){
				throw "Unable to allocate memory for indexing fasta file";
			}

			i = 0;

			// Keep the user informed
			if(_verbose){
				cerr << "Indexing fasta file, please wait:";
			}
			
			time_t profile = time(NULL);
			
			const char waiting[] = "|/-\\";
			const unsigned int waiting_size = strlen(waiting);
			unsigned int waiting_index = 0;
			
			const unsigned int waiting_every = 100;
			unsigned int waiting_count = 0;
			
			if(_verbose){
				cerr << waiting[waiting_index++];
			}
			
			// Make sure that multiple '>' symbols on the same line don't derail
			// the parser (this has been observed for sequences from the RDP).
			bool read_eol = true;
			
			while(i < buffer_size){
				
				if(_verbose){
				
					// Only update the waiting indicator evey waiting_every
					// passes through the white loop.
					if(waiting_count%waiting_every == 0){
						
						cerr << '\b';
						cerr << waiting[waiting_index++];
						waiting_index = waiting_index%waiting_size;
					}
					
					waiting_count++;
				}
				
				#ifdef WIN32

					DWORD bytes_read = 0;

					if( !ReadFile(fasta_in, file_buffer, FASTA_BUFFER_SIZE, &bytes_read, NULL) ){

						delete [] file_buffer;
						throw "Error while reading fasta file";
					}

				#else

					ssize_t bytes_read = read(fasta_in, file_buffer, FASTA_BUFFER_SIZE);

					if(bytes_read == -1){

						delete [] file_buffer;
						throw "Error while reading fasta file";
					}
				#endif // WIN32

				for(file_index j = 0;j < file_index(bytes_read);j++, i++){
					
					switch(file_buffer[j]){
						case '>':
							if(read_eol){
								seq_index.push_back(i);
								read_eol = false;
							}
							break;
						case '\n':
							read_eol = true;
							break;
					};	
				}
			}

			delete [] file_buffer;
			
			if(_verbose){
				
				profile = time(NULL) - profile;
				cerr << "\bindexing complete (" << profile << " sec)" << endl;
			}
		}
			
		const size_t num_seq = size();

		if(num_seq != 0){

			const size_t num_seq_m_1 = num_seq - 1;

			seq_length.resize(num_seq);

			for(size_t i = 0;i < num_seq_m_1;i++){

				// This is an overestimate of the sequence length since
				// the defline sizes are also included (as well as any white space)
				const unsigned int seq_len = seq_index[i + 1] - seq_index[i];

				seq_length[i] = make_pair(seq_len, i);
			}

			// The last sequences size is determined by its start position and the number of bytes
			// in the file
			seq_length[num_seq_m_1] = make_pair(buffer_size - seq_index[num_seq_m_1], num_seq_m_1);
		}
		
		if(m_sort_by_len == true){
			sort( seq_length.begin(), seq_length.end(), sequence_order() );
		}
	}
}

void sequence_data::load_fastq(const std::string &m_filename, const bool &m_sort_by_len,
	const bool &m_allow_fastq_mmap)
{
	bool allow_fastq_mmap = m_allow_fastq_mmap;
	
	// Try to open the file as a fastq file
	#ifdef WIN32

	// Note when compiling: for recent versions of Visual Studio, you will need
	// to set the option "Configuration Properties->General->Project Defaults->Character set"
	// to "Use Multi-Byte charater set" instead of the default "Use Unicode character set".
	// Otherwise, CreateFile() will expect Unicode strings (instead of char*)
	fasta_in = CreateFile(
		m_filename.c_str(),
		GENERIC_READ,
		FILE_SHARE_READ,
		NULL, // security descriptor
		OPEN_EXISTING,
		FILE_ATTRIBUTE_READONLY,
		NULL);

	if(fasta_in == INVALID_HANDLE_VALUE){
		throw "Unable to open fastq sequence file for reading";
	}

	#else
	fasta_in = ::open( m_filename.c_str(), O_RDONLY);

	if(fasta_in == -1){
		throw "Unable to open fastq sequence file for reading";
	}

	#endif// WIN32
	
	struct stat file_info;
	
	if(stat(m_filename.c_str(), &file_info) != 0){
		throw "Unable to read fastq file size";
	}
	
	buffer_size = file_info.st_size;
	
	// Don't try to memory map files larger than 2^32 - 1
	if(buffer_size > UINT_MAX){
		allow_fastq_mmap = false;
	}
	
	if(allow_fastq_mmap == true){
	
		format = FASTQ_MMAP;
		
		#ifdef WIN32
		h_map_file = CreateFileMapping(fasta_in,
				NULL, // If the security descriptor is null, this handle can not be inherited
				PAGE_READONLY,
				0, 0,
				NULL);

		if(h_map_file == NULL){
			
			if(_verbose){
				cerr << "Unable to memory map the fastq file -- insufficient memory" << endl;
			}
			
			format = FASTQ_SLOW;
		}
		else{

			buffer = MapViewOfFile(h_map_file,
				FILE_MAP_READ,
				0, // high order bits of offset
        			0, // low order bits of offset
        			0); // Bytes to map

			if(buffer == NULL){
				
				if(_verbose){
					cerr << "Unable to memory map the fastq file -- insufficient memory" << endl;
				}
				
				format = FASTQ_SLOW;
			}
		}

		#else

		buffer = (unsigned char*)mmap(0, buffer_size, PROT_READ, MAP_SHARED, fasta_in, 0);

		if(buffer == MAP_FAILED){

			switch(errno){
				case EACCES:
					throw "Specified fastq file is a non-regular file -- please check for filename typos";
				case EAGAIN:
					throw "Specified fastq file has been locked";
				case ENOMEM:
					
					if(_verbose){
						cerr << "Unable to memory map the fastq file -- insufficient memory" << endl;
					}
					
					format = FASTQ_SLOW;
					buffer = NULL;
					break;
				case ENODEV:
				
					if(_verbose){
						cerr << "File system does not support memory mapping" << endl;
					}
					
					format = FASTQ_SLOW;
					buffer = NULL;
					break;
				default:
				
					if(_verbose){
						cerr << "Unable to memory map the fastq file" << endl;
					}
					
					format = FASTQ_SLOW;
					buffer = NULL;
					break;
			};
			
			if(_verbose){
				cerr << "Fastq file will be read using standard disk access" << endl;
			}
		}
		#endif // WIN32
	}
	else{
		format = FASTQ_SLOW;
	}
	
	// If we have not been given fastq indicies (by the load_fasta_indicies() function),
	// then we need to load them now
	if(seq_index.empty() == true){
	
		// Parse through the fastq file
		// 1) count the number of sequences
		// 2) the file location of each defline (the '@' symbol)

		file_index i;

		if(format == FASTQ_MMAP){

			bool read_eol = true;			
			char last_header = '-';
			
			for(i = 0;i < buffer_size;i++){
				
				switch( ( (char*)(buffer) )[i] ){
					case '@':
						if(read_eol){
						
							if(last_header != '+'){
								seq_index.push_back(i);
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
		else{ // format == FASTQ_SLOW

			char *file_buffer = new char [FASTA_BUFFER_SIZE];

			if(file_buffer == NULL){
				throw "Unable to allocate memory for indexing fastq file";
			}

			i = 0;

			// Keep the user informed
			if(_verbose){
				cerr << "Indexing fastq file, please wait:";
			}
			
			time_t profile = time(NULL);
			
			const char waiting[] = "|/-\\";
			const unsigned int waiting_size = strlen(waiting);
			unsigned int waiting_index = 0;
			
			const unsigned int waiting_every = 1000;
			unsigned int waiting_count = 0;
			
			if(_verbose){
				cerr << waiting[waiting_index++];
			}
			
			bool read_eol = true;
			char last_header = '-';
			
			while(i < buffer_size){
				
				if(_verbose){
				
					// Only update the waiting indicator evey waiting_every
					// passes through the white loop.
					if(waiting_count%waiting_every == 0){
						
						cerr << '\b';
						cerr << waiting[waiting_index++];
						waiting_index = waiting_index%waiting_size;
					}
					
					waiting_count++;
				}
				
				#ifdef WIN32

					DWORD bytes_read = 0;

					if( !ReadFile(fasta_in, file_buffer, FASTA_BUFFER_SIZE, &bytes_read, NULL) ){

						delete [] file_buffer;
						throw "Error while reading fastq file";
					}

				#else

					ssize_t bytes_read = read(fasta_in, file_buffer, FASTA_BUFFER_SIZE);

					if(bytes_read == -1){

						delete [] file_buffer;
						throw "Error while reading fastq file";
					}
				#endif // WIN32

				for(file_index j = 0;j < file_index(bytes_read);j++, i++){
					
					switch(file_buffer[j]){
						case '@':
							if(read_eol){
						
								if(last_header != '+'){
									seq_index.push_back(i);
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

			delete [] file_buffer;
			
			if(_verbose){
				
				profile = time(NULL) - profile;
				cerr << "\bindexing complete (" << profile << " sec)" << endl;
			}
		}
			
		const size_t num_seq = size();

		if(num_seq != 0){

			const size_t num_seq_m_1 = num_seq - 1;

			seq_length.resize(num_seq);

			for(size_t i = 0;i < num_seq_m_1;i++){

				// This is an overestimate of the sequence length since
				// the defline sizes are also included (as well as any white space)
				const unsigned int seq_len = seq_index[i + 1] - seq_index[i];

				seq_length[i] = make_pair(seq_len, i);
			}

			// The last sequences size is determined by its start position and the number of bytes
			// in the file
			seq_length[num_seq_m_1] = make_pair(buffer_size - seq_index[num_seq_m_1], num_seq_m_1);
		}
		
		if(m_sort_by_len == true){
			sort( seq_length.begin(), seq_length.end(), sequence_order() );
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
	
	#ifdef WIN32

		// Low 32 bits of the new file location
		//LONG low_order_loc = LONG(seq_index[m_index] & 0xffffffff);

		// High 32 bits of the new file location. Right shifting a unigned int by 32 bits
		// does not seem to work on 32-bit windows (VS 6.0)
		//LONG high_order_loc = LONG( (seq_index[m_index] >> 32) & 0xffffffff );

		LARGE_INTEGER li;
		
		memset( &li, 0, sizeof(LARGE_INTEGER) );

		memcpy( &li, &(seq_index[m_index]), sizeof(seq_index[m_index]) );

		SetFilePointer( 
		  fasta_in, 
		  li.LowPart, 
		  &li.HighPart, 
		  FILE_BEGIN
		); 
	#else

		// Move to the start of the fasta defline (first base will be the '>' symbol)
		//lseek(fasta_in, seq_index[m_index], SEEK_SET);
		// Note that lseek is *not* thread safe (i.e. not openMP safe). Track the
		// position in the file buffer "by hand".
		file_index file_offset = seq_index[m_index];
	#endif // WIN32

	// Buffer the reading of fasta records in chunks of FASTA_BUFFER_SIZE or size of the 
	// record -- whichever is smaller.
	file_index record_size = (m_index == seq_index.size() - 1) ? buffer_size - seq_index[m_index] : 
		seq_index[m_index + 1] - seq_index[m_index];
	
	const file_index file_buffer_size = min(record_size, file_index(FASTA_BUFFER_SIZE) );
	
	char *file_buffer = new char[file_buffer_size];
	
	if(file_buffer == NULL){
		throw "Unable to allocate fasta file buffer";
	}
	
	file_index total_bytes_read = 0;
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	// Read the defline
	///////////////////////////////////////////////////////////////////////////////////////////////
	#ifdef WIN32

		DWORD bytes_read = 0;

		if( !ReadFile(fasta_in, file_buffer, file_buffer_size, &bytes_read, NULL) ){
			throw "Error while reading fasta file";
		}

	#else

		// We can't use read(), since lseek is not thread safe!
		//ssize_t bytes_read = read(fasta_in, file_buffer, file_buffer_size);

		ssize_t bytes_read = pread(fasta_in, file_buffer, file_buffer_size, file_offset);
		
		if(bytes_read == -1){
		
			delete [] file_buffer;
			throw "Error while reading fasta file";
		}
		
		file_offset += bytes_read;

	#endif

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
		ptr++;
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
			
			#ifdef WIN32

				DWORD bytes_read = 0;

				if( !ReadFile(fasta_in, file_buffer, 
						min(file_buffer_size, record_size - total_bytes_read), 
						&bytes_read, NULL) ){

					delete [] file_buffer;
					delete [] m_seq.second;
					throw "Error while reading fasta file";
				}

			#else
				// Read more data from the file, but don't use read(), since
				// lseek (which would be needed to position the file pointer) is
				// not thread safe
				//bytes_read = read(fasta_in, file_buffer, 
				//	min(file_buffer_size, record_size - total_bytes_read) );
				
				bytes_read = pread(fasta_in, file_buffer, 
					min(file_buffer_size, record_size - total_bytes_read), file_offset);
				
				if(bytes_read == -1){
				
					delete [] file_buffer;
					delete [] m_seq.second;
					
					throw "Error while reading fasta file";
				}
				
				file_offset += bytes_read;
			#endif

			total_bytes_read += file_index(bytes_read);

			// The first illegal buffer location
			last = file_buffer + bytes_read;

			ptr = file_buffer;
		}
		
		// Make sure that all bases are valid
		switch(*ptr){		
			case 'A': case 'a':
				if(index++ >= start){
					*(seq_ptr++) = DB_A;
					(*num_base)++;
				}
				break;
			case 'T': case 't': 
			case 'U': case 'u': // Map RNA -> DNA
				if(index++ >= start){
					*(seq_ptr++) = DB_T;
					(*num_base)++;
				}
				break;
			case 'G': case 'g': 
				if(index++ >= start){
					*(seq_ptr++) = DB_G;
					(*num_base)++;
				}
				break;
			case 'C': case 'c':
				if(index++ >= start){
					*(seq_ptr++) = DB_C;
					(*num_base)++;
				}
				break;
				
			// G or T or C
			case 'B':	case 'b':
				if(index++ >= start){
					*(seq_ptr++) = DB_B;
					(*num_base)++;
				}
				break;
			// G or A or T
			case 'D':	case 'd':
				if(index++ >= start){
					*(seq_ptr++) = DB_D;
					(*num_base)++;
				}
				break;
			// A or C or T
			case 'H':	case 'h':
				if(index++ >= start){
					*(seq_ptr++) = DB_H;
					(*num_base)++;
				}
				break;
			// G or T
			case 'K':	case 'k':
				if(index++ >= start){
					*(seq_ptr++) = DB_K;
					(*num_base)++;
				}
				break;
			// A or C
			case 'M':	case 'm':
				if(index++ >= start){
					*(seq_ptr++) = DB_M;
					(*num_base)++;
				}
				break;
			// A or C or G or T
			case 'N':	case 'n':
				if(index++ >= start){
					*(seq_ptr++) = DB_N;
					(*num_base)++;
				}
				break;
			// G or A
			case 'R':	case 'r':
				if(index++ >= start){
					*(seq_ptr++) = DB_R;
					(*num_base)++;
				}
				break;
			// G or C
			case 'S':	case 's':
				if(index++ >= start){
					*(seq_ptr++) = DB_S;
					(*num_base)++;
				}
				break;
			// G or C or A
			case 'V':	case 'v':
				if(index++ >= start){
					*(seq_ptr++) = DB_V;
					(*num_base)++;
				}
				break;
			// A or T
			case 'W':	case 'w':
				if(index++ >= start){
					*(seq_ptr++) = DB_W;
					(*num_base)++;
				}
				break;
			// T or C
			case 'Y':	case 'y':
				if(index++ >= start){
					*(seq_ptr++) = DB_Y;
					(*num_base)++;
				}
				break;
			
			case 'X':	case 'x':
				if(index++ >= start){
					*(seq_ptr++) = DB_UNKNOWN;
					(*num_base)++;
				}
				break;
				
			// Gap
			case '-':
				if(index++ >= start){
					*(seq_ptr++) = DB_GAP;
					(*num_base)++;
				}
				break;
				
			// White space and '*' symbols
			case ' ': case '\t': case '\n': case '\r': case '*':
				break;
			default:
			
				cerr << "\nFound the illegal base \"" << char(*ptr) << "\" in fasta record:" << endl;
				cerr << m_seq.first << endl;
				
				// Don't panic -- treat the base as an unknown and keep going
				//delete [] file_buffer;
				//delete [] m_seq.second;
				
				//throw "Illegal base in fasta file";
				if(index++ >= start){
					*(seq_ptr++) = DB_UNKNOWN;
					(*num_base)++;
				}
				break;
		};
		
		ptr++;
	}
	
	// Clean up the file buffer
	delete [] file_buffer;
	
	// Return the number of bases (excluding separators and terminator). The actual size
	// of the sequence memory returned may be greater than seq_size (due to end-of-line and
	// other white space characters).
	return SEQ_SIZE(m_seq.second);
}

unsigned int sequence_data::read_bio_seq_fasta_mmap(pair<string, SEQPTR> &m_seq, 
	const unsigned int &m_index) const
{
	return read_bio_seq_fasta_mmap(m_seq, m_index, 0, -1);
}

unsigned int sequence_data::read_bio_seq_fasta_mmap(pair<string, SEQPTR> &m_seq, 
	const unsigned int &m_index, const int &m_start, const int &m_stop) const
{
	if(format != FASTA_MMAP){
		throw ":sequence_data::read_bio_seq_fasta_mmap: Database is not in fasta format!";
	}
	
	if( m_index >= seq_index.size() ){
		throw ":sequence_data::read_bio_seq_fasta_mmap: Index out of bounds";
	}
	
	if(buffer == NULL){
		throw ":sequence_data::read_bio_seq_fasta_mmap: Invalid buffer";
	}
	
	// The first illegal buffer location
	const char* last = (char*)buffer + buffer_size;
	
	// Add an offset of 1 to skip the '>' symbol
	char *ptr = (char*)buffer + seq_index[m_index] + 1;
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	// Read the defline
	///////////////////////////////////////////////////////////////////////////////////////////////
	
	// Skip any leading white spaces
	while( (last > ptr) && isspace(*ptr) ){
		ptr++;
	}
	
	if(ptr == last){
		throw "Truncated fasta file detected!";
	}
	
	char* defline_start = ptr;
	
	while( (last > ptr) && (*ptr != '\n') && (*ptr != '\r') ){
		ptr++;
	}
	
	if(ptr == last){
		throw "Truncated fasta file detected!";
	}
	
	m_seq.first = string(defline_start, ptr - defline_start);
	
	// Count the number of bases in the sequence data
	char* seq_start = ptr;
		
	unsigned int num_base = 0;
	
	while( (last > ptr) && (*ptr != '>') ){
		
		// Make sure that all bases are valid
		switch(*ptr){
			case 'A': case 'a': case 'T': case 't': case 'G': case 'g': case 'C': case 'c':
			
			// Map RNA -> DNA
			case 'U': case 'u':
			
			// G or T or C
			case 'B':	case 'b':
			
			// G or A or T
			case 'D':	case 'd':
			
			// A or C or T
			case 'H':	case 'h':
			
			// G or T
			case 'K':	case 'k':
			
			// A or C
			case 'M':	case 'm':
			
			// A or C or G or T
			case 'N':	case 'n':
			
			// G or A
			case 'R':	case 'r':
			
			// G or C
			case 'S':	case 's':
			
			// G or C or A
			case 'V':	case 'v':
			
			// A or T
			case 'W':	case 'w':
			
			// T or C
			case 'Y':	case 'y':
			
			// Unknown base
			case 'X':	case 'x':
			
			// Gap
			case '-':
			
				num_base ++;
				break;
			
			// White space and '*' symbols
			case ' ': case '\t': case '\n': case '\r': case '*':
				break;
			default:
				cerr << "\nFound the illegal base \"" << char(*ptr) << "\" in fasta record:" << endl;
				cerr << m_seq.first << endl;
				
				// Don't panic, treat the base as an unknown and keep going
				//throw "Illegal base in fasta file";
				num_base ++;
				break;
		};
		
		ptr++;
	}	
	
	if(m_start < 0){
		throw __FILE__ ":sequence_data::read_bio_seq_fasta_mmap: m_start is out of bounds";
	}
	
	const unsigned int start = m_start;
	
	unsigned int stop;
	
	if( (m_stop < 0) || (m_stop >= int(num_base) ) ){
		stop = num_base - 1;
	}
	else{
	
		stop = m_stop;		
	}
		
	// Recompute the number of bases, allowing for the possibility that start > stop  
	num_base = (start > stop) ? 0 : stop - start + 1;
	
	// Allocate memory to hold the selected sequence. This memory must be deallocated by the calling function
	// when the sequence is no longer needed.
	m_seq.second = new SEQBASE[num_base + SEQ_HEADER_SIZE];
	
	if(m_seq.second == NULL){
		throw __FILE__ ":sequence_data::read_bio_seq_fasta_mmap: Unable to allocate memory for sequence";
	}
	
	ptr = seq_start;
	SEQPTR seq_ptr = m_seq.second;
	
	// Initialize the sequence header
	memcpy( seq_ptr, &num_base, sizeof(unsigned int) );
	seq_ptr += sizeof(unsigned int);
		
	unsigned int index = 0;

	if(num_base > 0){
	
		while( (last > ptr) && (*ptr != '>') ){

			if(index > stop){
				break;
			}

			// Make sure that all bases are valid
			switch(*ptr){
				case 'A': case 'a': 
					if(index++ >= start){
						*(seq_ptr++) = DB_A;
					}
					break;
				case 'T': case 't': 
				case 'U': case 'u': 
					if(index++ >= start){
						*(seq_ptr++) = DB_T;
					}
					break;
				case 'G': case 'g': 
					if(index++ >= start){
						*(seq_ptr++) = DB_G;
					}
					break;
				case 'C': case 'c':
					if(index++ >= start){
						*(seq_ptr++) = DB_C;
					}
					break;

				// G or T or C
				case 'B':	case 'b':
					if(index++ >= start){
						*(seq_ptr++) = DB_B;
					}
					break;
				// G or A or T
				case 'D':	case 'd':
					if(index++ >= start){
						*(seq_ptr++) = DB_D;
					}
					break;
				// A or C or T
				case 'H':	case 'h':
					if(index++ >= start){
						*(seq_ptr++) = DB_H;
					}
					break;
				// G or T
				case 'K':	case 'k':
					if(index++ >= start){
						*(seq_ptr++) = DB_K;
					}
					break;
				// A or C
				case 'M':	case 'm':
					if(index++ >= start){
						*(seq_ptr++) = DB_M;
					}
					break;
				// A or C or G or T
				case 'N':	case 'n':
					if(index++ >= start){
						*(seq_ptr++) = DB_N;
					}
					break;
				// G or A
				case 'R':	case 'r':
					if(index++ >= start){
						*(seq_ptr++) = DB_R;
					}
					break;
				// G or C
				case 'S':	case 's':
					if(index++ >= start){
						*(seq_ptr++) = DB_S;
					}
					break;
				// G or C or A
				case 'V':	case 'v':
					if(index++ >= start){
						*(seq_ptr++) = DB_V;
					}
					break;
				// A or T
				case 'W':	case 'w':
					if(index++ >= start){
						*(seq_ptr++) = DB_W;
					}
					break;
				// T or C
				case 'Y':	case 'y':
					if(index++ >= start){
						*(seq_ptr++) = DB_Y;
					}
					break;

				case 'X':	case 'x':
					if(index++ >= start){
						*(seq_ptr++) = DB_UNKNOWN;
					}
					break;

				// Gap
				case '-':
					if(index++ >= start){
						*(seq_ptr++) = DB_GAP;
					}
					break;

				// White space and '*' symbols
				case ' ': case '\t': case '\n': case '\r': case '*':
					break;
				default:
					// Don't panic, treat the base as an unknown and keep going
					//cerr << "Found the illegal base \"" << char(*ptr) << "\" in fasta record:" << endl;
					//cerr << m_seq.first << endl;
					//throw "Illegal base in fasta file";
					if(index++ >= start){
						*(seq_ptr++) = DB_UNKNOWN;
					}
					break;
			};

			ptr++;
		}
	}
	
	// Return the number of bases (excluding separators and terminator)  
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
	
	#ifdef WIN32

		// Low 32 bits of the new file location
		//LONG low_order_loc = LONG(seq_index[m_index] & 0xffffffff);

		// High 32 bits of the new file location. Right shifting a unigned int by 32 bits
		// does not seem to work on 32-bit windows (VS 6.0)
		//LONG high_order_loc = LONG( (seq_index[m_index] >> 32) & 0xffffffff );

		LARGE_INTEGER li;
		
		memset( &li, 0, sizeof(LARGE_INTEGER) );

		memcpy( &li, &(seq_index[m_index]), sizeof(seq_index[m_index]) );

		SetFilePointer( 
		  fasta_in, 
		  li.LowPart, 
		  &li.HighPart, 
		  FILE_BEGIN
		); 
	#else

		// Move to the start of the fasta defline (first base will be the '>' symbol)
		//lseek(fasta_in, seq_index[m_index], SEEK_SET);
		// Note that lseek is *not* thread safe (i.e. not openMP safe). Track the
		// position in the file buffer "by hand".
		file_index file_offset = seq_index[m_index];
	#endif // WIN32

	// Buffer the reading of fasta records in chunks of FASTA_BUFFER_SIZE or size of the 
	// record -- whichever is smaller.
	file_index record_size = (m_index == seq_index.size() - 1) ? buffer_size - seq_index[m_index] : 
		seq_index[m_index + 1] - seq_index[m_index];
	
	const file_index file_buffer_size = min(record_size, file_index(FASTA_BUFFER_SIZE) );
	
	char *file_buffer = new char[file_buffer_size];
	
	if(file_buffer == NULL){
		throw "Unable to allocate fastq file buffer";
	}
	
	file_index total_bytes_read = 0;
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	// Read the defline
	///////////////////////////////////////////////////////////////////////////////////////////////
	#ifdef WIN32

		DWORD bytes_read = 0;

		if( !ReadFile(fasta_in, file_buffer, file_buffer_size, &bytes_read, NULL) ){
			throw "Error while reading fastq file";
		}

	#else

		// We can't use read(), since lseek is not thread safe!
		//ssize_t bytes_read = read(fasta_in, file_buffer, file_buffer_size);

		ssize_t bytes_read = pread(fasta_in, file_buffer, file_buffer_size, file_offset);
		
		if(bytes_read == -1){
		
			delete [] file_buffer;
			throw "Error while reading fastq file";
		}
		
		file_offset += bytes_read;

	#endif

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
			
			#ifdef WIN32

				DWORD bytes_read = 0;

				if( !ReadFile(fasta_in, file_buffer, 
						min(file_buffer_size, record_size - total_bytes_read), 
						&bytes_read, NULL) ){

					delete [] file_buffer;
					delete [] m_seq.second;
					throw "Error while reading fastq file";
				}

			#else
				// Read more data from the file, but don't use read(), since
				// lseek (which would be needed to position the file pointer) is
				// not thread safe
				//bytes_read = read(fasta_in, file_buffer, 
				//	min(file_buffer_size, record_size - total_bytes_read) );
				
				bytes_read = pread(fasta_in, file_buffer, 
					min(file_buffer_size, record_size - total_bytes_read), file_offset);
				
				if(bytes_read == -1){
				
					delete [] file_buffer;
					delete [] m_seq.second;
					
					throw "Error while reading fastq file";
				}
				
				file_offset += bytes_read;
			#endif

			total_bytes_read += file_index(bytes_read);

			// The first illegal buffer location
			last = file_buffer + bytes_read;

			ptr = file_buffer;
		}
		
		// The FASTQ requires that all sequence be on a single line
		if( (*ptr == '\n') || (*ptr == '\r') ){
			break;
		}
		
		// Make sure that all bases are valid
		switch(*ptr){		
			case 'A': case 'a':
				if(index++ >= start){
					*(seq_ptr++) = DB_A;
					(*num_base)++;
				}
				break;
			case 'T': case 't': 
			case 'U': case 'u': // Map RNA -> DNA
				if(index++ >= start){
					*(seq_ptr++) = DB_T;
					(*num_base)++;
				}
				break;
			case 'G': case 'g': 
				if(index++ >= start){
					*(seq_ptr++) = DB_G;
					(*num_base)++;
				}
				break;
			case 'C': case 'c':
				if(index++ >= start){
					*(seq_ptr++) = DB_C;
					(*num_base)++;
				}
				break;
				
			// G or T or C
			case 'B':	case 'b':
				if(index++ >= start){
					*(seq_ptr++) = DB_B;
					(*num_base)++;
				}
				break;
			// G or A or T
			case 'D':	case 'd':
				if(index++ >= start){
					*(seq_ptr++) = DB_D;
					(*num_base)++;
				}
				break;
			// A or C or T
			case 'H':	case 'h':
				if(index++ >= start){
					*(seq_ptr++) = DB_H;
					(*num_base)++;
				}
				break;
			// G or T
			case 'K':	case 'k':
				if(index++ >= start){
					*(seq_ptr++) = DB_K;
					(*num_base)++;
				}
				break;
			// A or C
			case 'M':	case 'm':
				if(index++ >= start){
					*(seq_ptr++) = DB_M;
					(*num_base)++;
				}
				break;
			// A or C or G or T
			case 'N':	case 'n':
				if(index++ >= start){
					*(seq_ptr++) = DB_N;
					(*num_base)++;
				}
				break;
			// G or A
			case 'R':	case 'r':
				if(index++ >= start){
					*(seq_ptr++) = DB_R;
					(*num_base)++;
				}
				break;
			// G or C
			case 'S':	case 's':
				if(index++ >= start){
					*(seq_ptr++) = DB_S;
					(*num_base)++;
				}
				break;
			// G or C or A
			case 'V':	case 'v':
				if(index++ >= start){
					*(seq_ptr++) = DB_V;
					(*num_base)++;
				}
				break;
			// A or T
			case 'W':	case 'w':
				if(index++ >= start){
					*(seq_ptr++) = DB_W;
					(*num_base)++;
				}
				break;
			// T or C
			case 'Y':	case 'y':
				if(index++ >= start){
					*(seq_ptr++) = DB_Y;
					(*num_base)++;
				}
				break;
			
			case 'X':	case 'x':
				if(index++ >= start){
					*(seq_ptr++) = DB_UNKNOWN;
					(*num_base)++;
				}
				break;
				
			// Gap
			case '-':
				if(index++ >= start){
					*(seq_ptr++) = DB_GAP;
					(*num_base)++;
				}
				break;
				
			// Skip white space and '*' symbols (don't include '\n' or '\r' as white space as fastq files
			// have all sequence on a single line).
			case ' ': case '\t': case '*':
				break;
			default:
			
				cerr << "\nFound the illegal base \"" << char(*ptr) << "\" in fastq record:" << endl;
				cerr << m_seq.first << endl;
				
				// Don't panic -- treat the base as an unknown and keep going
				//delete [] file_buffer;
				//delete [] m_seq.second;
				
				//throw "Illegal base in fasta file";
				if(index++ >= start){
					*(seq_ptr++) = DB_UNKNOWN;
					(*num_base)++;
				}
				break;
		};
		
		ptr++;
	}
	
	// Clean up the file buffer
	delete [] file_buffer;
	
	// Return the number of bases (excluding separators and terminator). The actual size
	// of the sequence memory returned may be greater than seq_size (due to end-of-line and
	// other white space characters).
	return SEQ_SIZE(m_seq.second);
}

unsigned int sequence_data::read_bio_seq_fastq_mmap(pair<string, SEQPTR> &m_seq, 
	const unsigned int &m_index) const
{
	return read_bio_seq_fastq_mmap(m_seq, m_index, 0, -1);
}

unsigned int sequence_data::read_bio_seq_fastq_mmap(pair<string, SEQPTR> &m_seq, 
	const unsigned int &m_index, const int &m_start, const int &m_stop) const
{
	if(format != FASTQ_MMAP){
		throw ":sequence_data::read_bio_seq_fastq_mmap: Database is not in fastq format!";
	}
	
	if( m_index >= seq_index.size() ){
		throw ":sequence_data::read_bio_seq_fastq_mmap: Index out of bounds";
	}
	
	if(buffer == NULL){
		throw ":sequence_data::read_bio_seq_fastq_mmap: Invalid buffer";
	}
	
	// The first illegal buffer location
	const char* last = (char*)buffer + buffer_size;
	
	// Add an offset of 1 to skip the '@' symbol
	char *ptr = (char*)buffer + seq_index[m_index] + 1;
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	// Read the defline
	///////////////////////////////////////////////////////////////////////////////////////////////
	
	// Skip any leading white spaces
	while( (last > ptr) && isspace(*ptr) ){
		ptr++;
	}
	
	if(ptr == last){
		throw "Truncated fastq file detected!";
	}
	
	char* defline_start = ptr;
	
	while( (last > ptr) && (*ptr != '\n') && (*ptr != '\r') ){
		ptr++;
	}
	
	if(ptr == last){
		throw "Truncated fastq file detected!";
	}
	
	m_seq.first = string(defline_start, ptr - defline_start);
	
	// Skip the '\n' or '\r' character
	ptr++;
	
	if(ptr == last){
		throw "Truncated fastq file detected!";
	}
	
	// Count the number of bases in the sequence data
	char* seq_start = ptr;
		
	unsigned int num_base = 0;
	
	while( (last > ptr) && (*ptr != '\n') && (*ptr != '\r') ){

		// Make sure that all bases are valid
		switch(*ptr){
			case 'A': case 'a': case 'T': case 't': case 'G': case 'g': case 'C': case 'c':
			
			// Map RNA -> DNA
			case 'U': case 'u':
			
			// G or T or C
			case 'B':	case 'b':
			
			// G or A or T
			case 'D':	case 'd':
			
			// A or C or T
			case 'H':	case 'h':
			
			// G or T
			case 'K':	case 'k':
			
			// A or C
			case 'M':	case 'm':
			
			// A or C or G or T
			case 'N':	case 'n':
			
			// G or A
			case 'R':	case 'r':
			
			// G or C
			case 'S':	case 's':
			
			// G or C or A
			case 'V':	case 'v':
			
			// A or T
			case 'W':	case 'w':
			
			// T or C
			case 'Y':	case 'y':
			
			// Unknown base
			case 'X':	case 'x':
			
			// Gap
			case '-':
			
				num_base ++;
				break;
			
			// Ignore white space and '*' symbols (don't include '\n' or '\r' as white space as fastq files
			// have all sequence on a single line).
			case ' ': case '\t': case '*':
				break;
			default:
				cerr << "\nFound the illegal base \"" << char(*ptr) << "\" in fastq record:" << endl;
				cerr << m_seq.first << endl;
				
				// Don't panic, treat the base as an unknown and keep going
				//throw "Illegal base in fasta file";
				num_base ++;
				break;
		};
		
		ptr++;
	}
	
	if(m_start < 0){
		throw __FILE__ ":sequence_data::read_bio_seq_fastq_mmap: m_start is out of bounds";
	}
	
	const unsigned int start = m_start;
	
	unsigned int stop;
	
	if( (m_stop < 0) || (m_stop >= int(num_base) ) ){
		stop = num_base - 1;
	}
	else{
	
		stop = m_stop;		
	}
		
	// Recompute the number of bases, allowing for the possibility that start > stop  
	num_base = (start > stop) ? 0 : stop - start + 1;
	
	// Allocate memory to hold the selected sequence. This memory must be deallocated by the calling function
	// when the sequence is no longer needed.
	m_seq.second = new SEQBASE[num_base + SEQ_HEADER_SIZE];
	
	if(m_seq.second == NULL){
		throw __FILE__ ":sequence_data::read_bio_seq_fastq_mmap: Unable to allocate memory for sequence";
	}
	
	ptr = seq_start;
	SEQPTR seq_ptr = m_seq.second;
	
	// Initialize the sequence header
	memcpy( seq_ptr, &num_base, sizeof(unsigned int) );
	seq_ptr += sizeof(unsigned int);
		
	unsigned int index = 0;

	if(num_base > 0){
	
		while( (last > ptr) && (*ptr != '\n') && (*ptr != '\r') ){

			if(index > stop){
				break;
			}

			// Make sure that all bases are valid
			switch(*ptr){
				case 'A': case 'a': 
					if(index++ >= start){
						*(seq_ptr++) = DB_A;
					}
					break;
				case 'T': case 't': 
				case 'U': case 'u': 
					if(index++ >= start){
						*(seq_ptr++) = DB_T;
					}
					break;
				case 'G': case 'g': 
					if(index++ >= start){
						*(seq_ptr++) = DB_G;
					}
					break;
				case 'C': case 'c':
					if(index++ >= start){
						*(seq_ptr++) = DB_C;
					}
					break;

				// G or T or C
				case 'B':	case 'b':
					if(index++ >= start){
						*(seq_ptr++) = DB_B;
					}
					break;
				// G or A or T
				case 'D':	case 'd':
					if(index++ >= start){
						*(seq_ptr++) = DB_D;
					}
					break;
				// A or C or T
				case 'H':	case 'h':
					if(index++ >= start){
						*(seq_ptr++) = DB_H;
					}
					break;
				// G or T
				case 'K':	case 'k':
					if(index++ >= start){
						*(seq_ptr++) = DB_K;
					}
					break;
				// A or C
				case 'M':	case 'm':
					if(index++ >= start){
						*(seq_ptr++) = DB_M;
					}
					break;
				// A or C or G or T
				case 'N':	case 'n':
					if(index++ >= start){
						*(seq_ptr++) = DB_N;
					}
					break;
				// G or A
				case 'R':	case 'r':
					if(index++ >= start){
						*(seq_ptr++) = DB_R;
					}
					break;
				// G or C
				case 'S':	case 's':
					if(index++ >= start){
						*(seq_ptr++) = DB_S;
					}
					break;
				// G or C or A
				case 'V':	case 'v':
					if(index++ >= start){
						*(seq_ptr++) = DB_V;
					}
					break;
				// A or T
				case 'W':	case 'w':
					if(index++ >= start){
						*(seq_ptr++) = DB_W;
					}
					break;
				// T or C
				case 'Y':	case 'y':
					if(index++ >= start){
						*(seq_ptr++) = DB_Y;
					}
					break;

				case 'X':	case 'x':
					if(index++ >= start){
						*(seq_ptr++) = DB_UNKNOWN;
					}
					break;

				// Gap
				case '-':
					if(index++ >= start){
						*(seq_ptr++) = DB_GAP;
					}
					break;

				// Skip white space and '*' symbols (don't include '\n' or '\r' as white space as fastq files
				// have all sequence on a single line).
				case ' ': case '\t': case '*':
					break;
				default:
					// Don't panic, treat the base as an unknown and keep going
					//cerr << "Found the illegal base \"" << char(*ptr) << "\" in fasta record:" << endl;
					//cerr << m_seq.first << endl;
					//throw "Illegal base in fasta file";
					if(index++ >= start){
						*(seq_ptr++) = DB_UNKNOWN;
					}
					break;
			};

			ptr++;
		}
	}
	
	// Return the number of bases (excluding separators and terminator)  
	return SEQ_SIZE(m_seq.second);
}
