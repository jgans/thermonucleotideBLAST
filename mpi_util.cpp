#include "mpi_util.h"
#include "options.h"
#include "hybrid_sig.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for std::string
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<string>(const string &m_str)
{
	return sizeof(size_t) + m_str.size();
}

template<>
unsigned char* mpi_pack<string>(unsigned char* m_ptr, const string &m_str)
{
	size_t len = m_str.size();
	
	memcpy( m_ptr, &len, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	memcpy(m_ptr, m_str.c_str(), len);
	m_ptr += len;
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack<string>(unsigned char* m_ptr, string &m_str)
{
	size_t len;
	
	memcpy( &len, m_ptr, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	m_str.assign( (char*)m_ptr, len );
	m_ptr += len;
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for Options
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const Options &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		OPTIONS_MEMBERS
	#undef VARIABLE

	return ret;	
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, Options &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
		OPTIONS_MEMBERS
    #undef VARIABLE
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const Options &m_obj)
{

	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
		OPTIONS_MEMBERS
	#undef VARIABLE

	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for hybrid_sig
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const hybrid_sig &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		HYBRID_SIG_MEMBERS
	#undef VARIABLE

	return ret;	
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, hybrid_sig &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
		HYBRID_SIG_MEMBERS
    #undef VARIABLE
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const hybrid_sig &m_obj)
{

	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
		HYBRID_SIG_MEMBERS
	#undef VARIABLE

	return m_ptr;
}