#ifndef __MPI_UTIL
#define __MPI_UTIL

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include <limits.h>
#include <deque>
#include <vector>
#include <string>
#include <string.h> // memcpy
#include <unordered_map>

// Forward function definitions for containers (needed to be able to transport nested C++ structures):
template<class T> size_t mpi_size(const std::deque<T> &m_obj);
template<class T> unsigned char* mpi_pack(unsigned char* m_ptr, const std::deque<T> &m_obj);
template<class T> unsigned char* mpi_unpack(unsigned char* m_ptr, std::deque<T> &m_obj);

template<class T> size_t mpi_size(const std::vector<T> &m_obj);
template<class T> unsigned char* mpi_pack(unsigned char* m_ptr, const std::vector<T> &m_obj);
template<class T> unsigned char* mpi_unpack(unsigned char* m_ptr, std::vector<T> &m_obj);

template<class A, class B> size_t mpi_size(const std::pair<A, B> &m_obj);
template<class A, class B> unsigned char* mpi_pack(unsigned char* m_ptr, const std::pair<A, B> &m_obj);
template<class A, class B> unsigned char* mpi_unpack(unsigned char* m_ptr, std::pair<A, B> &m_obj);

// Use a template for simple objects. Specialize as needed for more complex types
template <class T>
size_t mpi_size(const T &m_obj)
{
	// Force a *compile* time test of whether this is a native or derived type
	static_assert(std::is_fundamental<T>::value || std::is_enum<T>::value,
		":mpi_size: Non-fundamental or non-enum type passed as template");
		
	return sizeof(m_obj);
}

template <class T>
unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj)
{
	// Force a *compile* time test of whether this is a native or derived type
	static_assert(std::is_fundamental<T>::value || std::is_enum<T>::value,
		":mpi_pack: Non-fundamental or non-enum type passed as template");
		
	memcpy( m_ptr, &m_obj, sizeof(m_obj) );
	m_ptr += sizeof(m_obj);
	
	return m_ptr;
}

template <class T>
unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj)
{
	// Force a *compile* time test of whether this is a native or derived type
	static_assert(std::is_fundamental<T>::value || std::is_enum<T>::value,
		":mpi_unpack: Non-fundamental or non-enum type passed as template");
		
	memcpy( &m_obj, m_ptr, sizeof(m_obj) );
	m_ptr += sizeof(m_obj);
	
	return m_ptr;
}

// Specialization for string
template<>
	size_t mpi_size(const std::string &m_str);
template<>
	unsigned char* mpi_pack(unsigned char* m_ptr, const std::string &m_str);
template<>
	unsigned char* mpi_unpack(unsigned char* m_ptr, std::string &m_str);

/////////////////////////////////////////////////////////////////////////////////////////
// Overload for std::deque
/////////////////////////////////////////////////////////////////////////////////////////
template<class T>
size_t mpi_size(const std::deque<T> &m_obj)
{
	size_t len = sizeof(size_t);
	
	for(typename std::deque<T>::const_iterator i = m_obj.begin();i != m_obj.end();++i){
		len += mpi_size(*i);
	}
	
	return len;
}

template<class T>
unsigned char* mpi_pack(unsigned char* m_ptr, const std::deque<T> &m_obj)
{
	size_t len = m_obj.size();
	
	memcpy( m_ptr, &len, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	for(typename std::deque<T>::const_iterator i = m_obj.begin();i != m_obj.end();++i){
		m_ptr = mpi_pack(m_ptr, *i);
	}
	
	return m_ptr;
}

template<class T>
unsigned char* mpi_unpack(unsigned char* m_ptr, std::deque<T> &m_obj)
{
	size_t len;
	
	memcpy( &len, m_ptr, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	m_obj.resize(len);
	
	for(size_t i = 0;i < len;++i){
		m_ptr = mpi_unpack(m_ptr, m_obj[i]);
	}
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Overload for std::vector
/////////////////////////////////////////////////////////////////////////////////////////
template<class T>
size_t mpi_size(const std::vector<T> &m_obj)
{
	size_t len = sizeof(size_t);
	
	for(typename std::vector<T>::const_iterator i = m_obj.begin();i != m_obj.end();++i){
		len += mpi_size(*i);
	}
	
	return len;
}

template<class T>
unsigned char* mpi_pack(unsigned char* m_ptr, const std::vector<T> &m_obj)
{
	size_t len = m_obj.size();
	
	memcpy( m_ptr, &len, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	for(typename std::vector<T>::const_iterator i = m_obj.begin();i != m_obj.end();++i){
		m_ptr = mpi_pack(m_ptr, *i);
	}
	
	return m_ptr;
}

template<class T>
unsigned char* mpi_unpack(unsigned char* m_ptr, std::vector<T> &m_obj)
{
	size_t len;
	
	memcpy( &len, m_ptr, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	m_obj.resize(len);
	
	for(size_t i = 0;i < len;++i){
		m_ptr = mpi_unpack(m_ptr, m_obj[i]);
	}
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Overload for std::pair
/////////////////////////////////////////////////////////////////////////////////////////
template<class A, class B>
size_t mpi_size(const std::pair<A, B> &m_obj)
{
	return mpi_size(m_obj.first) + mpi_size(m_obj.second);
}

template<class A, class B>
unsigned char* mpi_pack(unsigned char* m_ptr, const std::pair<A, B> &m_obj)
{
	m_ptr = mpi_pack(m_ptr, m_obj.first);
	m_ptr = mpi_pack(m_ptr, m_obj.second);
        
        return m_ptr;
}

template<class A, class B>
unsigned char* mpi_unpack(unsigned char* m_ptr, std::pair<A, B> &m_obj)
{       
        m_ptr = mpi_unpack(m_ptr, m_obj.first);
        m_ptr = mpi_unpack(m_ptr, m_obj.second);
        
        return m_ptr;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Generic broadcast (from rank m_src_rank to all other ranks)
//////////////////////////////////////////////////////////////////////////////////////////////
template <class T>
void broadcast(T &m_obj, const int &m_my_rank, const int &m_src_rank)
{
	size_t len = (m_my_rank == m_src_rank) ? mpi_size(m_obj) : 0;
	
	#ifdef USE_MPI
	MPI_Bcast(&len, sizeof(len), MPI_BYTE, m_src_rank, MPI_COMM_WORLD);
	#endif // USE_MPI

	unsigned char *buffer = new unsigned char[len];
	
	if(buffer == NULL){
		throw __FILE__ ":broadcast: Unable to allocate buffer";
	}
	
	if(m_my_rank == m_src_rank){
		mpi_pack(buffer, m_obj);		
	}

	// MPI_Bcast can send memory buffers with at most INT_MAX elements.
	// If we need to send more than INT_MAX elements, we will need to
	// break the buffer into chunks. To reduce the total amount of memory
	// used, we can make the max chunk size even smaller than INT_MAX.
	
	unsigned char *ptr = buffer;
	const size_t max_buffer_size = 1073741824UL; //Conserve memory with a max buffer size of 1 GB
	
	while(len > 0){
		
		const size_t curr_len = (len > max_buffer_size) ? max_buffer_size : len;
		
		#ifdef USE_MPI
		MPI_Bcast(ptr, curr_len, MPI_BYTE, m_src_rank, MPI_COMM_WORLD);
		#endif // USE_MPI

		len -= curr_len;
		ptr += curr_len;
	}
	
	if(m_my_rank != m_src_rank){
		mpi_unpack(buffer, m_obj);
	}
		
	delete [] buffer;
}

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
	
	#ifdef USE_MPI
	MPI_Bcast(&num_elem, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast(buffer, buffer_len, MPI_BYTE, 0, MPI_COMM_WORLD);
	#endif // USE_MPI

	delete [] buffer;
}

// Receive from rank 0
template<typename _tmp>
inline void receive(std::vector<_tmp> &m_array)
{
	unsigned int num_elem;
	
	#ifdef USE_MPI
	MPI_Bcast(&num_elem, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	#endif // USE_MPI

	unsigned char *buffer = new unsigned char [num_elem*sizeof(_tmp)];
	
	if(buffer == NULL){
		throw __FILE__ ":receive: Unable to allocate memory for buffer";
	}
	
	#ifdef USE_MPI
	MPI_Bcast(buffer, num_elem*sizeof(_tmp), MPI_BYTE, 0, MPI_COMM_WORLD);
	#endif // USE_MPI

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
	#ifdef USE_MPI
	// Include the '\0' string terminator
	unsigned int len = m_str.size() + 1;
	
	MPI_Bcast(&len, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast((char*)m_str.c_str(), len, MPI_CHAR, 0, MPI_COMM_WORLD);
	#endif // USE_MPI
}

// Receive from rank 0
inline void receive(std::string &m_str)
{
	unsigned int len = 0;
	
	#ifdef USE_MPI
	MPI_Bcast(&len, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	#endif // USE_MPI

	char *buffer = new char [len];
	
	if(buffer == NULL){
		throw __FILE__ ":receive: Unable to allocate memory for buffer";
	}
	
	#ifdef USE_MPI
	MPI_Bcast(buffer,len, MPI_CHAR, 0, MPI_COMM_WORLD);
	#endif // USE_MPI

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
	
	#ifdef USE_MPI
	MPI_Bcast(&num_elem, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast(buffer, buffer_len, MPI_BYTE, 0, MPI_COMM_WORLD);
	#endif // USE_MPI

	delete [] buffer;
}

// Receive from rank 0
template<typename _tmp>
inline void receive(std::deque<_tmp> &m_array)
{
	unsigned int num_elem = 0;
	
	#ifdef USE_MPI
	MPI_Bcast(&num_elem, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	#endif // USE_MPI

	unsigned char *buffer = new unsigned char [num_elem*sizeof(_tmp)];
	
	if(buffer == NULL){
		throw __FILE__ ":receive: Unable to allocate memory for buffer";
	}
	
	#ifdef USE_MPI
	MPI_Bcast(buffer, num_elem*sizeof(_tmp), MPI_BYTE, 0, MPI_COMM_WORLD);
	#endif // USE_MPI

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
