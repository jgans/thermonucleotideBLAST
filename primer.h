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

#ifndef __PCR_PRIMER_HEURISTIC
#define	__PCR_PRIMER_HEURISTIC


#include "circle_buffer.h"
#include <vector>

// This class encapsulates heuristic rules for PCR primer design
// There are currently 5 heuristic criteria that are considered:
//
//	a) Avoid runs of 3 or more G's or C's at the 3' prime end.			POLY_3_GC
//	b) The 5 bases at the 5' end should contain no more than			MULTI_5_GC
//			3 G's or C's if no two pyrimidines (T|C) are adjacent
//			2 G's or C's otherwise
//	c) Polyprimidine (T|C) and polypurine (A|G) runs should be			NO_POLY_RUNS
//		 avoided.
//	d) Avoid 3' terminal T's							NO_3_T
//	e) GC_min <= %GC content <= GC_max						GC_CONTENT
// ------------------------------------------------------------------------------
// 	TaqMan (tm) Probe Design rules
//	f) No 5' G 									NO_5_G
//	g) Avoid runs of identical nuleotides						NO_IDENTICAL_RUNS
//	h) No 5' G AND 3' C (i.e. we will be unable to avoid				NO_5_G_3_C
//		 a 5'G by switching strands)
//	i) No penultimate 5' G (i.e. G in second 5' position)				NO_5_PENULTIMATE_G
//		From the ABI PrimerExpress manual
// ------------------------------------------------------------------------------
// 	TaqMan (tm) MGB Probe Design rules
//	j) No 3' G residues 9 (from the ABI PrimerExpress manual)			NO_3_POLY_G
//		Avoid 5'-...GGG-MGB-3' and 5'-...GGAG-MGB-3'
//	k) No CC dinucleotides in the middle of the probe				NO_MIDDLE_CC
//		(from the ABI PrimerExpress manual). Consider the middle 50%
//		of the probe.

#define	BAD_BASE				-1
#define	PCR_VALID				0

// Please note that the bit variables that will be passed to the 
// PCRPrimer constructor need to be unsigned int's (to prevent the
// compiler from getting confused about ambigous typecasting)
#define	POLY_3_GC				1u
#define	MULTI_5_GC				(1u << 1)
#define	NO_POLY_RUNS				(1u << 2)
#define	NO_3_T					(1u << 3)
#define	GC_CONTENT				(1u << 4)
#define	NO_5_G					(1u << 5)
#define	NO_5_G_3_C				(1u << 6)
#define	NO_IDENTICAL_RUNS			(1u << 7)
#define	NO_5_PENULTIMATE_G			(1u << 8)
#define	NO_3_POLY_G				(1u << 9)
#define	NO_MIDDLE_CC				(1u << 10)

#define	HEURISTIC_PCR_GC_MIN	0.40f
#define	HEURISTIC_PCR_GC_MAX	0.60f
#define	HEURISTIC_PCR_RUN		4

#define	MAX_PCR_LENGTH			128

// For NO_MIDDLE_CC
#define	MIDDLE_LOWER_BOUND	0.25f
#define	MIDDLE_UPPER_BOUND	0.75f

// JDG 2/3/10
// The NCBI C toolkit has started using a "PCRPrimer" name. In order to remain compatible with
// the NCBI toolkit, I am hiding the PCRPrimer class name (in this file) behind a unique namespace

namespace ASSY_HEURISTIC{

class PCRPrimer {
private:
		// What options are we testing?
		unsigned int mask;

		// Allowed ranges for GC conent
		float gc_min;
		float gc_max;

		// How many bases constitute a "run" for NO_POLY_RUNS ?
		unsigned int run_len;

		CircleBuffer<unsigned char, MAX_PCR_LENGTH> primer_buffer;

		int operator() (const CircleBuffer<unsigned char, MAX_PCR_LENGTH> &m_primer, 
			const unsigned int &m_mask) const;
		
		bool report_verbose;
public:
	enum {A = 0, T, G, C};

	PCRPrimer(const unsigned int &m_mask, 
		const float &m_min_gc = HEURISTIC_PCR_GC_MIN, 
		const float &m_max_gc = HEURISTIC_PCR_GC_MAX) :
		mask(m_mask)
	{
		run_len = HEURISTIC_PCR_RUN;
		gc_range(m_min_gc, m_max_gc);
		report_verbose = true;
	};
	
	PCRPrimer(const float &m_min_gc = HEURISTIC_PCR_GC_MIN, 
		const float &m_max_gc = HEURISTIC_PCR_GC_MAX)
	{
		mask = 0;
		run_len = HEURISTIC_PCR_RUN;
		gc_range(m_min_gc, m_max_gc);
		report_verbose = true;
	};
	
	inline void verbose(const bool &m_verbose)
	{
		report_verbose = m_verbose;
	};
	
	inline unsigned int options() const
	{
		return mask;
	};

	inline void options(const unsigned int &m_mask)
	{
		mask = m_mask;
	};

	inline void run(const unsigned int &m_run_len)
	{
		run_len = m_run_len;
	};

	inline void gc_range(const float &m_min_gc, const float &m_max_gc)
	{
		if(m_min_gc > 1.0f){
			throw "GC min out of bounds (> 1.0)";
		}

		if(m_min_gc < 0.0f){
			throw "GC min out of bounds (< 0.0)";
		}
		
		if(m_max_gc > 1.0f){
			throw "GC max out of bounds (> 1.0)";
		}

		if(m_max_gc < 0.0f){
			throw "GC max out of bounds (< 0.0)";
		}
		
		if(m_max_gc < m_min_gc){
			throw "GC range out of bounds (max < min)";
		}

		gc_min = m_min_gc;
		gc_max = m_max_gc;
	};
	
	inline int operator() (const std::string &m_primer) const
	{
		return (*this)(m_primer, mask);
	};
	
	int operator() (const std::string &m_primer, const unsigned int &m_mask) const;
	
	// Use the internal buffer
	int operator() () const;
	
	// Use the internal buffer and the specified mask
	int operator() (const unsigned int &m_mask) const;
	
	inline void push_back(const char &m_base)
	{
		switch(m_base){
			case 'A':
			case 'a':
				primer_buffer.push_back(A);
				break;
			case 'T':
			case 't':
				primer_buffer.push_back(T);
				break;
			case 'G':
			case 'g':
				primer_buffer.push_back(G);
				break;
			case 'C':
			case 'c':
				primer_buffer.push_back(C);
				break;
			default:
				// This is an illegal base!
				primer_buffer.clear();
		};
	};
	
	inline void push_back_native(const unsigned char &m_base)
	{
		primer_buffer.push_back(m_base);
	};

	inline void push_front(const char &m_base)
	{
		switch(m_base){
			case 'A':
			case 'a':
				primer_buffer.push_front(A);
				break;
			case 'T':
			case 't':
				primer_buffer.push_front(T);
				break;
			case 'G':
			case 'g':
				primer_buffer.push_front(G);
				break;
			case 'C':
			case 'c':
				primer_buffer.push_front(C);
				break;
			default:
				// This is an illegal base!
				primer_buffer.clear();
		};
	};
	
	inline void push_front_native(const unsigned char &m_base)
	{
		primer_buffer.push_front(m_base);
	};

	inline void pop_front()
	{
		if(primer_buffer.empty() == false){
			primer_buffer.pop_front();
		}
	};

	inline void pop_back()
	{
		if(primer_buffer.empty() == false){
			primer_buffer.pop_back();
		}
	};

	inline void clear()
	{
		primer_buffer.clear();
	};

	inline unsigned int size() const
	{
		return primer_buffer.size();
	};

	std::string str() const;
	
	std::string error(const int &m_code) const;
	
	// A helper function
	std::pair<unsigned int, unsigned int> num_g_num_c(const std::string &m_seq) const;
};

} // namespace ASSY_HEURISTIC

unsigned int parse_pcr_filter_options(const std::vector<std::string> &m_opt);

#endif // __PCR_PRIMER_HEURISTIC
