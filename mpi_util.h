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

#ifndef __MPI_UTIL
#define __MPI_UTIL

#include <mpi.h>
#include <vector>
#include <deque>
#include <string>

//template<typename _tmp>
//void send(const std::vector<_tmp> &m_array);

//template<typename _tmp>
//void receive(std::vector<_tmp> &m_array);

// Send from rank 0 to all other ranks
template<typename _tmp>
inline void send(const std::vector<_tmp> &m_array)
{
	unsigned int num_elem = m_array.size();
	unsigned int buffer_len = num_elem*sizeof(_tmp);
	unsigned char *buffer = new unsigned char [buffer_len];
	
	if(buffer == NULL){
		throw __FILE__ ":send: Unable to allocate memory for send buffer";
	}
	
	unsigned char *ptr = buffer;
	typename std::vector<_tmp>::const_iterator iter;
	
	for(iter = m_array.begin();iter != m_array.end();iter++){
		memcpy(ptr, &(*iter), sizeof(_tmp) );
		
		ptr += sizeof(_tmp);
	}
	
	MPI_Bcast(&num_elem, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast(buffer, buffer_len, MPI_BYTE, 0, MPI_COMM_WORLD);
	
	delete [] buffer;
}

// Receive from rank 0
template<typename _tmp>
inline void receive(std::vector<_tmp> &m_array)
{
	unsigned int num_elem;
	
	MPI_Bcast(&num_elem, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	
	unsigned char *buffer = new unsigned char [num_elem*sizeof(_tmp)];
	
	if(buffer == NULL){
		throw __FILE__ ":receive: Unable to allocate memory for buffer";
	}
	
	MPI_Bcast(buffer, num_elem*sizeof(_tmp), MPI_BYTE, 0, MPI_COMM_WORLD);
	
	m_array.resize(num_elem);
	
	unsigned char *ptr = buffer;
	typename std::vector<_tmp>::iterator iter;
	
	for(iter = m_array.begin();iter != m_array.end();iter++){

		*iter = *( (_tmp *)ptr );
		
		ptr += sizeof(_tmp);
	}
	
	delete [] buffer;
}

inline void send(const std::string &m_str)
{
	// Include the '\0' string terminator
	unsigned int len = m_str.size() + 1;
	
	MPI_Bcast(&len, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast((char*)m_str.c_str(), len, MPI_CHAR, 0, MPI_COMM_WORLD);
}

// Receive from rank 0
inline void receive(std::string &m_str)
{
	unsigned int len;
	
	MPI_Bcast(&len, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	
	char *buffer = new char [len];
	
	if(buffer == NULL){
		throw __FILE__ ":receive: Unable to allocate memory for buffer";
	}
	
	MPI_Bcast(buffer,len, MPI_CHAR, 0, MPI_COMM_WORLD);
	
	m_str = buffer;
	
	delete [] buffer;
}

template<typename _tmp>
inline void send(const std::deque<_tmp> &m_array)
{
	unsigned int num_elem = m_array.size();
	unsigned int buffer_len = num_elem*sizeof(_tmp);
	unsigned char *buffer = new unsigned char [buffer_len];
	
	if(buffer == NULL){
		throw __FILE__ ":send: Unable to allocate memory for send buffer";
	}
	
	unsigned char *ptr = buffer;
	typename std::deque<_tmp>::const_iterator iter;
	
	for(iter = m_array.begin();iter != m_array.end();iter++){
		memcpy(ptr, &(*iter), sizeof(_tmp) );
		
		ptr += sizeof(_tmp);
	}
	
	MPI_Bcast(&num_elem, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast(buffer, buffer_len, MPI_BYTE, 0, MPI_COMM_WORLD);
	
	delete [] buffer;
}

// Receive from rank 0
template<typename _tmp>
inline void receive(std::deque<_tmp> &m_array)
{
	unsigned int num_elem;
	
	MPI_Bcast(&num_elem, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	
	unsigned char *buffer = new unsigned char [num_elem*sizeof(_tmp)];
	
	if(buffer == NULL){
		throw __FILE__ ":receive: Unable to allocate memory for buffer";
	}
	
	MPI_Bcast(buffer, num_elem*sizeof(_tmp), MPI_BYTE, 0, MPI_COMM_WORLD);
	
	m_array.resize(num_elem);
	
	unsigned char *ptr = buffer;
	typename std::deque<_tmp>::iterator iter;
	
	for(iter = m_array.begin();iter != m_array.end();iter++){

		*iter = *( (_tmp *)ptr );
		
		ptr += sizeof(_tmp);
	}
	
	delete [] buffer;
}
#endif // __MPI_UTIL
