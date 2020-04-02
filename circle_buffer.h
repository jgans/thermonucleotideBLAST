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

// A circular buffer
// v 1.0
// J. D. Gans
// Los Alamos National Lab
// 4/28/2005
//
//	6/12/06	Added indexing via operator[].

#ifndef __CIRCLE_BUFFER
#define __CIRCLE_BUFFER

#include <stdlib.h>

typedef unsigned int CircleBuffer_index;

// Forward declaration for the iterators
template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE>
class CircleBuffer_iterator;

template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE>
class CircleBuffer_const_iterator;

template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE>
class CircleBuffer_reverse_iterator;

template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE>
class CircleBuffer_const_reverse_iterator;
		
// A circular buffer for more efficient memory management
template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE>
class CircleBuffer {

	private:
		T data[__INTERNAL_BUFFER_SIZE];
		CircleBuffer_index head, tail;

	public:
		
		friend class CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE>;
		friend class CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE>;

		friend class CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE>;
		friend class CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE>;
	
		typedef CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE> iterator;
		typedef CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE> const_iterator;
		
		typedef CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> reverse_iterator;
		typedef CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> const_reverse_iterator;

		CircleBuffer()
		{
			head = tail = 0;
		};
		
		iterator begin();
		
		const_iterator begin() const;
		
		iterator end();

		const_iterator end() const;

		reverse_iterator rbegin();

		const_reverse_iterator rbegin() const;

		reverse_iterator rend();

		const_reverse_iterator rend() const;

		inline void push_back(const T &m_elem)
		{
			data[tail] = m_elem;

			tail = (tail + 1) % __INTERNAL_BUFFER_SIZE;
		};
		
		inline void pop_back()
		{
			tail = (tail == 0) ? __INTERNAL_BUFFER_SIZE - 1 : tail - 1;
		};

		inline void push_front(const T &m_elem)
		{
			head = (head == 0) ? __INTERNAL_BUFFER_SIZE - 1 : head - 1;

			data[head] = m_elem;
		};

		inline void pop_front()
		{
			head = (head + 1) % __INTERNAL_BUFFER_SIZE;
		};

		inline unsigned int size() const
		{
			return (tail >= head) ? tail - head : (__INTERNAL_BUFFER_SIZE - head) + tail;
		};

		inline bool empty() const
		{
			return (head == tail);
		};

		inline void clear()
		{
			head = tail = 0;
		};

		inline T & front()
		{
			return data[head];
		};

		inline T & back()
		{
			return (tail == 0) ? data[__INTERNAL_BUFFER_SIZE - 1] : data[tail - 1];
		};
		
		inline T front() const
		{
			return data[head];
		};

		inline T back() const
		{
			return (tail == 0) ? data[__INTERNAL_BUFFER_SIZE - 1] : data[tail - 1];
		};
		
		// No bounds checking!
		inline T & operator[](const unsigned int &m_index)
		{
			return data[(head + m_index)%__INTERNAL_BUFFER_SIZE];
		};
		
		inline const T & operator[](const unsigned int &m_index) const
		{
			return data[(head + m_index)%__INTERNAL_BUFFER_SIZE];
		};
};

//////////////////////////////////////////////////////////////////////////////////////////////
// Forward iterator

template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE>
class CircleBuffer_iterator {
	friend class CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE>;
	
	private:
		CircleBuffer<T, __INTERNAL_BUFFER_SIZE> *ptr;
		CircleBuffer_index cur;
				
	public:
		inline CircleBuffer_iterator() : ptr(NULL), cur(0)
		{
		};
		
		inline CircleBuffer_iterator(const CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy) : 
			ptr(m_copy.ptr), cur(m_copy.cur)
		{
		};

		inline CircleBuffer_iterator(CircleBuffer<T, __INTERNAL_BUFFER_SIZE> &m_parent, 
			const CircleBuffer_index &m_cur) : ptr(&m_parent), cur(m_cur)
		{
		};
				
		inline bool operator != (const CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy) const
		{
			return (ptr != m_copy.ptr) || (cur != m_copy.cur);
		};
		
		inline bool operator == (const CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy) const
		{
			return (ptr == m_copy.ptr) && (cur == m_copy.cur);
		};

		inline CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE> & operator = (const CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy)
		{
			ptr = m_copy.ptr;
			cur = m_copy.cur;

			return (*this);
		};
		
		// Prefix
		inline CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE> & operator ++ ()
		{
			cur = (cur != ptr->tail) ? (cur + 1) % __INTERNAL_BUFFER_SIZE : cur;

			return (*this);
		};

		// Postfix
		inline CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE> operator ++ (int)
		{
			CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE> ret(*this);

			cur = (cur != ptr->tail) ? (cur + 1) % __INTERNAL_BUFFER_SIZE : cur;

			return ret;
		};
		
		// Prefix
		inline CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE> & operator -- ()
		{
			if(cur == 0){
				cur = (cur != ptr->head) ? __INTERNAL_BUFFER_SIZE - 1: cur;
			}
			else{
				cur = (cur != ptr->head) ? (cur - 1) % __INTERNAL_BUFFER_SIZE : cur;
			}

			return (*this);
		};

		// Postfix
		inline CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE> operator -- (int)
		{
			CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE> ret(*this);

			if(cur == 0){
				cur = (cur != ptr->head) ? __INTERNAL_BUFFER_SIZE - 1: cur;
			}
			else{
				cur = (cur != ptr->head) ? (cur - 1) % __INTERNAL_BUFFER_SIZE : cur;
			}

			return ret;
		};

		inline T & operator * ()
		{
			return ptr->data[cur];
		};
		
		inline T * operator -> ()
		{
			return &(ptr->data[cur]);
		};
};

template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE>
class CircleBuffer_const_iterator {
	private:
		const CircleBuffer<T, __INTERNAL_BUFFER_SIZE> *ptr;
		CircleBuffer_index cur;
				
	public:
		inline CircleBuffer_const_iterator() : ptr(NULL), cur(0)
		{
		};
		
		inline CircleBuffer_const_iterator(const CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy) : 
			ptr(m_copy.ptr), cur(m_copy.cur)
		{
		};
		
		inline CircleBuffer_const_iterator(const CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy) : 
			ptr(m_copy.ptr), cur(m_copy.cur)
		{
		};
		
		inline CircleBuffer_const_iterator(const CircleBuffer<T, __INTERNAL_BUFFER_SIZE> &m_parent, 
			const CircleBuffer_index &m_cur) : ptr(&m_parent), cur(m_cur)
		{
		};
		
		inline bool operator != (const CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy) const
		{
			return (ptr != m_copy.ptr) || (cur != m_copy.cur);
		};
		
		inline bool operator == (const CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy) const
		{
			return (ptr == m_copy.ptr) && (cur == m_copy.cur);
		};

		inline CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE> & operator = (const CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy)
		{
			ptr = m_copy.ptr;
			cur = m_copy.cur;

			return (*this);
		};

		inline CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE> & operator = (const CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy)
		{
			ptr = m_copy.ptr;
			cur = m_copy.cur;

			return (*this);
		};
		
		// Prefix
		inline CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE> & operator ++ ()
		{
			cur = (cur != ptr->tail) ? (cur + 1) % __INTERNAL_BUFFER_SIZE : cur;

			return (*this);
		};

		// Postfix
		inline CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE> operator ++ (int)
		{
			CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE> ret(*this);

			cur = (cur != ptr->tail) ? (cur + 1) % __INTERNAL_BUFFER_SIZE : cur;

			return ret;
		};

		// Prefix
		inline CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE> & operator -- ()
		{
			if(cur == 0){
				cur = (cur != ptr->head) ? __INTERNAL_BUFFER_SIZE - 1: cur;
			}
			else{
				cur = (cur != ptr->head) ? (cur - 1) % __INTERNAL_BUFFER_SIZE : cur;
			}

			return (*this);
		};

		// Postfix
		inline CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE> operator -- (int)
		{
			CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE> ret(*this);

			if(cur == 0){
				cur = (cur != ptr->head) ? __INTERNAL_BUFFER_SIZE - 1: cur;
			}
			else{
				cur = (cur != ptr->head) ? (cur - 1) % __INTERNAL_BUFFER_SIZE : cur;
			}

			return ret;
		};

		inline T operator * ()
		{
			return ptr->data[cur];
		};
		
		inline T * operator -> ()
		{
			return &(ptr->data[cur]);
		};
};

////////////////////////////////////////////////////////////////////////////////////////
// Reverse iterator

template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE>
class CircleBuffer_reverse_iterator {
	friend class CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE>;
	
	private:
		CircleBuffer<T, __INTERNAL_BUFFER_SIZE> *ptr;
		CircleBuffer_index cur;
				
	public:
		inline CircleBuffer_reverse_iterator() : ptr(NULL), cur(0)
		{
		};
		
		inline CircleBuffer_reverse_iterator(const CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy) : 
			ptr(m_copy.ptr), cur(m_copy.cur)
		{	
		};

		inline CircleBuffer_reverse_iterator(CircleBuffer<T, __INTERNAL_BUFFER_SIZE> &m_parent, 
			const CircleBuffer_index &m_cur) : ptr(&m_parent), cur(m_cur)
		{	
		};
				
		inline bool operator != (const CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy) const
		{
			return (ptr != m_copy.ptr) || (cur != m_copy.cur);
		};
		
		inline bool operator == (const CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy) const
		{
			return (ptr == m_copy.ptr) && (cur == m_copy.cur);
		};

		inline CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> & operator = (const CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy)
		{
			ptr = m_copy.ptr;
			cur = m_copy.cur;

			return (*this);
		};

		// Prefix
		inline CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> & operator -- ()
		{
			cur = (cur != ptr->tail) ? (cur + 1) % __INTERNAL_BUFFER_SIZE : cur;

			return (*this);
		};

		// Postfix
		inline CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> operator -- (int)
		{
			CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> ret(*this);

			cur = (cur != ptr->tail) ? (cur + 1) % __INTERNAL_BUFFER_SIZE : cur;

			return ret;
		};
		
		// Prefix
		inline CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> & operator ++ ()
		{
			if(cur == 0){
				cur = ( (cur != ptr->head) && (cur != __INTERNAL_BUFFER_SIZE) ) ? 
					__INTERNAL_BUFFER_SIZE - 1: __INTERNAL_BUFFER_SIZE;
			}
			else{
				cur = ( (cur != ptr->head) && (cur != __INTERNAL_BUFFER_SIZE) ) ? 
					(cur - 1) % __INTERNAL_BUFFER_SIZE : __INTERNAL_BUFFER_SIZE;
			}

			return (*this);
		};

		// Postfix
		inline CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> operator ++ (int)
		{
			CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> ret(*this);

			if(cur == 0){
				cur = ( (cur != ptr->head) && (cur != __INTERNAL_BUFFER_SIZE) ) ? 
					__INTERNAL_BUFFER_SIZE - 1: __INTERNAL_BUFFER_SIZE;
			}
			else{
				cur = ( (cur != ptr->head) && (cur != __INTERNAL_BUFFER_SIZE) ) ? 
					(cur - 1) % __INTERNAL_BUFFER_SIZE : __INTERNAL_BUFFER_SIZE;
			}

			return ret;
		};

		inline T & operator * ()
		{
			return ptr->data[cur];
		};
		
		inline T * operator -> ()
		{
			return &(ptr->data[cur]);
		};
};

template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE>
class CircleBuffer_const_reverse_iterator {
	
	private:
		const CircleBuffer<T, __INTERNAL_BUFFER_SIZE> *ptr;
		CircleBuffer_index cur;
				
	public:
		inline CircleBuffer_const_reverse_iterator() : ptr(NULL), cur(0)
		{
		};
		
		inline CircleBuffer_const_reverse_iterator(const CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy) : 
			ptr(m_copy.ptr), cur(m_copy.cur)
		{
		};
		
		inline CircleBuffer_const_reverse_iterator(const CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy) : 
			ptr(m_copy.ptr), cur(m_copy.cur)
		{
		};
		
		inline CircleBuffer_const_reverse_iterator(const CircleBuffer<T, __INTERNAL_BUFFER_SIZE> &m_parent, 
			const CircleBuffer_index &m_cur) : ptr(&m_parent), cur(m_cur)
		{
		};
		
		inline bool operator != (const CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy) const
		{
			return (ptr != m_copy.ptr) || (cur != m_copy.cur);
		};
		
		inline bool operator == (const CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy) const
		{
			return (ptr == m_copy.ptr) && (cur == m_copy.cur);
		};

		inline CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> & operator = (const CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy)
		{
			ptr = m_copy.ptr;
			cur = m_copy.cur;

			return (*this);
		};

		inline CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> & operator = (const CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> &m_copy)
		{
			ptr = m_copy.ptr;
			cur = m_copy.cur;

			return (*this);
		};
		
		// Prefix
		inline CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> & operator -- ()
		{
			cur = (cur != ptr->tail) ? (cur + 1) % __INTERNAL_BUFFER_SIZE : cur;

			return (*this);
		};

		// Postfix
		inline CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> operator -- (int)
		{
			CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> ret(*this);

			cur = (cur != ptr->tail) ? (cur + 1) % __INTERNAL_BUFFER_SIZE : cur;

			return ret;
		};

		// Prefix
		inline CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> & operator ++ ()
		{
			if(cur == 0){
				cur = ( (cur != ptr->head) && (cur != __INTERNAL_BUFFER_SIZE) ) ? 
					__INTERNAL_BUFFER_SIZE - 1: __INTERNAL_BUFFER_SIZE;
			}
			else{
				cur = ( (cur != ptr->head) && (cur != __INTERNAL_BUFFER_SIZE) ) ? (cur - 1) % 
					__INTERNAL_BUFFER_SIZE : __INTERNAL_BUFFER_SIZE;
			}

			return (*this);
		};

		// Postfix
		inline CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> operator ++ (int)
		{
			CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE> ret(*this);

			if(cur == 0){
				cur = ( (cur != ptr->head) && (cur != __INTERNAL_BUFFER_SIZE) ) ? 
					__INTERNAL_BUFFER_SIZE - 1: __INTERNAL_BUFFER_SIZE;
			}
			else{
				cur = ( (cur != ptr->head) && (cur != __INTERNAL_BUFFER_SIZE) ) ? 
					(cur - 1) % __INTERNAL_BUFFER_SIZE : __INTERNAL_BUFFER_SIZE;
			}

			return ret;
		};

		inline T operator * ()
		{
			return ptr->data[cur];
		};
		
		inline T * operator -> ()
		{
			return &(ptr->data[cur]);
		};
};

template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE> 
CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE> 
CircleBuffer<T, __INTERNAL_BUFFER_SIZE>::begin()
{
	return iterator(*this, head);
};

template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE>
CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE>
CircleBuffer<T, __INTERNAL_BUFFER_SIZE>::begin() const
{
	return const_iterator(*this, head);
};

template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE>
CircleBuffer_iterator<T, __INTERNAL_BUFFER_SIZE>
CircleBuffer<T, __INTERNAL_BUFFER_SIZE>::end()
{
	return iterator(*this, tail);
};

template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE>
CircleBuffer_const_iterator<T, __INTERNAL_BUFFER_SIZE>
CircleBuffer<T, __INTERNAL_BUFFER_SIZE>::end() const
{
	return const_iterator(*this, tail);
};

template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE>
CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE>
CircleBuffer<T, __INTERNAL_BUFFER_SIZE>::rbegin()
{	
	return (tail == 0) ? 
		reverse_iterator(*this, __INTERNAL_BUFFER_SIZE - 1) :
		reverse_iterator(*this, tail - 1);
};

template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE>
CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE>
CircleBuffer<T, __INTERNAL_BUFFER_SIZE>::rbegin() const
{
	return (tail == 0) ? 
		const_reverse_iterator(*this, __INTERNAL_BUFFER_SIZE - 1) :
		const_reverse_iterator(*this, tail - 1);
};

template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE>
CircleBuffer_reverse_iterator<T, __INTERNAL_BUFFER_SIZE>
CircleBuffer<T, __INTERNAL_BUFFER_SIZE>::rend()
{
	return reverse_iterator(*this, __INTERNAL_BUFFER_SIZE);
};

template <class T, CircleBuffer_index __INTERNAL_BUFFER_SIZE>
CircleBuffer_const_reverse_iterator<T, __INTERNAL_BUFFER_SIZE>
CircleBuffer<T, __INTERNAL_BUFFER_SIZE>::rend() const
{
	return const_reverse_iterator(*this, __INTERNAL_BUFFER_SIZE);
};

#endif // __CIRCLE_BUFFER
