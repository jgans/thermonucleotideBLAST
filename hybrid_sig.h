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

#ifndef __SIGNATURE_CLASS
#define __SIGNATURE_CLASS

#include <string.h> // for memcpy

#include <string>
#include <algorithm> // for std::pair
#include <vector>
#include <sstream>

// The symbol defined by COMMENT_SYMBOL is used to comment out
// lines in the assya input file
#define	COMMENT_SYMBOL						'#'

// What assay are we simulating?
enum {ASSAY_PCR, ASSAY_PROBE, ASSAY_PADLOCK, ASSAY_AFFYMETRIX, ASSAY_NONE};

// A class for storing hybridization signatures
class hybrid_sig {

	private:
		int id; // An identifier that is unique to an assay
		int degen_id; // An identifier that is unique to the assays generated 
					  // by the expansion of degenerate assays. If an assay does not
					  // contain any degeneracy, then degen_id == id.
		int seq_index; // The sequence index that was matched by this record
		
		inline int def_to_gi(const std::string &m_defline) const
		{
			
			size_t start = m_defline.find("gi|");
			
			if(start == std::string::npos){
				return -1;
			}
		
			start += 3; // strlen("gi|")
			
			const size_t len = m_defline.size();
			
			size_t stop = start;
			
			while( (stop < len) && isdigit(m_defline[stop]) ){
				stop ++;
			}
			
			return atoi( m_defline.substr(start, stop - start).c_str() );
		};

	public:
		// Enumerate the DNA strand orientation
		// A value of PLUS for probe_strand means the probe binds to the 
		// forward strand, while a value of MINUS means the probe binds to
		// the reverse strand
		enum {PLUS = 0, MINUS};
		
		hybrid_sig() // default constructor
		{
			id = degen_id = -1;
			seq_index = -1;
			
			amplicon_range = std::make_pair(0, 0);
			probe_range = std::make_pair(0, 0);
			
			forward_tm = -1.0f;
			reverse_tm = -1.0f;
			probe_tm = -1.0f;
			
			forward_hairpin_tm = -1.0f;
			reverse_hairpin_tm = -1.0f;
			forward_dimer_tm = -1.0f;
			reverse_dimer_tm = -1.0f;
			primer_dimer_tm = -1.0f;
		
			probe_hairpin_tm = -1.0f;
			probe_dimer_tm = -1.0f;
		
			forward_dH = 100.0f;
			forward_dS = 0.0f;
		
			reverse_dH = 100.0f;
			reverse_dS = 0.0f;
		
			probe_dH = 100.0f;
			probe_dS = 0.0f;
		
			primer_strand = PLUS;
			probe_strand = PLUS;

			forward_degen = 1;
			reverse_degen = 1;
			probe_degen = 1;
			
			forward_mm = -1;
			forward_gap = -1;
			reverse_mm = -1;
			reverse_gap = -1;
			probe_mm = -1;
			probe_gap = -1;
			
			ct = -1.0f;
		};
		
		// forward + reverse + probe
		hybrid_sig(const std::string &m_name, const std::string &m_forward, 
			const std::string &m_reverse, const std::string &m_probe, 
			const unsigned int &m_id)
		{
			id = degen_id = m_id;
			seq_index = -1;
			
			name = m_name;
			forward_oligo = m_forward;
			reverse_oligo = m_reverse;
			probe_oligo = m_probe;
			
			amplicon_range = std::make_pair(0, 0);
			probe_range = std::make_pair(0, 0);
			
			forward_tm = -1.0f;
			reverse_tm = -1.0f;
			probe_tm = -1.0f;
			
			forward_hairpin_tm = -1.0f;
			reverse_hairpin_tm = -1.0f;
			forward_dimer_tm = -1.0f;
			reverse_dimer_tm = -1.0f;
			primer_dimer_tm = -1.0f;
		
			probe_hairpin_tm = -1.0f;
			probe_dimer_tm = -1.0f;
			
			forward_dH = 100.0f;
			forward_dS = 0.0f;
		
			reverse_dH = 100.0f;
			reverse_dS = 0.0f;
		
			probe_dH = 100.0f;
			probe_dS = 0.0f;
			
			primer_strand = PLUS;
			probe_strand = PLUS;
			
			forward_primer_clamp = -1;
			reverse_primer_clamp = -1;

			forward_degen = 1;
			reverse_degen = 1;
			probe_degen = 1;
			
			forward_mm = -1;
			forward_gap = -1;
			reverse_mm = -1;
			reverse_gap = -1;
			probe_mm = -1;
			probe_gap = -1;
			
			ct = -1.0f;
		}
		
		// forward + reverse
		hybrid_sig(const std::string &m_name, const std::string &m_forward, 
			const std::string &m_reverse, const unsigned int &m_id)
		{
			id = degen_id = m_id;
			seq_index = -1;
			
			name = m_name;
			forward_oligo = m_forward;
			reverse_oligo = m_reverse;
			
			amplicon_range = std::make_pair(0, 0);
			probe_range = std::make_pair(0, 0);
			
			forward_tm = -1.0f;
			reverse_tm = -1.0f;
			probe_tm = -1.0f;
			
			forward_hairpin_tm = -1.0f;
			reverse_hairpin_tm = -1.0f;
			forward_dimer_tm = -1.0f;
			reverse_dimer_tm = -1.0f;
			primer_dimer_tm = -1.0f;
		
			probe_hairpin_tm = -1.0f;
			probe_dimer_tm = -1.0f;
			
			forward_dH = 100.0f;
			forward_dS = 0.0f;
		
			reverse_dH = 100.0f;
			reverse_dS = 0.0f;
		
			probe_dH = 100.0f;
			probe_dS = 0.0f;
			
			primer_strand = PLUS;
			probe_strand = PLUS;
			
			forward_primer_clamp = -1;
			reverse_primer_clamp = -1;

			forward_degen = 1;
			reverse_degen = 1;
			probe_degen = 1;
			
			forward_mm = -1;
			forward_gap = -1;
			reverse_mm = -1;
			reverse_gap = -1;
			probe_mm = -1;
			probe_gap = -1;
			
			ct = -1.0f;
		}
		
		// probe
		hybrid_sig(const std::string &m_name, const std::string &m_probe, 
			const unsigned int &m_id)
		{
			id = degen_id = m_id;
			seq_index = -1;
			
			name = m_name;
			probe_oligo = m_probe;
			
			amplicon_range = std::make_pair(0, 0);
			probe_range = std::make_pair(0, 0);
			
			forward_tm = -1.0f;
			reverse_tm = -1.0f;
			probe_tm = -1.0f;
			
			forward_hairpin_tm = -1.0f;
			reverse_hairpin_tm = -1.0f;
			forward_dimer_tm = -1.0f;
			reverse_dimer_tm = -1.0f;
			primer_dimer_tm = -1.0f;
		
			probe_hairpin_tm = -1.0f;
			probe_dimer_tm = -1.0f;
			
			forward_dH = 100.0f;
			forward_dS = 0.0f;
		
			reverse_dH = 100.0f;
			reverse_dS = 0.0f;
		
			probe_dH = 100.0f;
			probe_dS = 0.0f;
			
			primer_strand = PLUS;
			probe_strand = PLUS;
			
			forward_primer_clamp = -1;
			reverse_primer_clamp = -1;

			forward_degen = 1;
			reverse_degen = 1;
			probe_degen = 1;
			
			forward_mm = -1;
			forward_gap = -1;
			reverse_mm = -1;
			reverse_gap = -1;
			probe_mm = -1;
			probe_gap = -1;
			
			ct = -1.0f;
		}
		
		~hybrid_sig()
		{
			// Do nothing!
		};
		
		std::string name;
		std::string forward_oligo;
		std::string reverse_oligo;
		std::string probe_oligo;
		
		std::string amplicon_def;
		std::string amplicon;
		
		std::string forward_align;
		std::string reverse_align;
		std::string probe_align;
		
		std::pair<int, int> amplicon_range;
		std::pair<int, int> probe_range;
		
		float forward_tm;
		float reverse_tm;
		float probe_tm;
		
		float forward_hairpin_tm;
		float reverse_hairpin_tm;
		float forward_dimer_tm;
		float reverse_dimer_tm;
		float primer_dimer_tm;
		float probe_hairpin_tm;
		float probe_dimer_tm;
		
		float forward_dH;
		float forward_dS;
		
		float reverse_dH;
		float reverse_dS;
		
		float probe_dH;
		float probe_dS;
		
		int primer_strand;
		int probe_strand;
		
		// The smallest number of exact 3' matches for the forward
		// and reverse primer
		int forward_primer_clamp;
		int reverse_primer_clamp;
		
		// The level of degeneracy for each assay oligo
		int forward_degen;
		int reverse_degen;
		int probe_degen;
		
		// Mismatch and gap values for each oligo
		int forward_mm;
		int forward_gap;
		int reverse_mm;
		int reverse_gap;
		int probe_mm;
		int probe_gap;
		
		// The cycle threshold value for real-time PCR
		float ct;
		
		inline size_t mpi_size() const
		{
			return 	sizeof(id) + sizeof(degen_id) + 
					sizeof(seq_index) +
					name.size() + 1 + // Include '\0'
					forward_oligo.size() + 1 + // Include '\0'
					reverse_oligo.size() + 1 + // Include '\0'
					probe_oligo.size() + 1 + // Include '\0'
					amplicon_def.size() + 1 + // Include '\0'
					amplicon.size() + 1 + // Include '\0'
					forward_align.size() + 1 + // Include '\0'
					reverse_align.size() + 1 + // Include '\0'
					probe_align.size() + 1	+ // Include '\0'
					2*sizeof(int) +	// amplicon_range
					2*sizeof(int) +	// probe_range
					3*sizeof(float) + // forward, reverse and probe Tms
					7*sizeof(float) + // oligo dimer and hairpin Tms
					6*sizeof(float) + // forward, reverse and probe dH and dS
					sizeof(primer_strand) +
					sizeof(probe_strand) +
					sizeof(forward_primer_clamp) +
					sizeof(reverse_primer_clamp) +
					3*sizeof(int) +	// The assay oligo degeneracy levels
					6*sizeof(int) + // The mismatch and gap parameters
					sizeof(float); // Ct
		
		};
		
		inline unsigned char* mpi_pack(unsigned char* m_ptr) const
		{
			memcpy( m_ptr, &id, sizeof(id) );
			m_ptr += sizeof(id);
			
			memcpy( m_ptr, &degen_id, sizeof(degen_id) );
			m_ptr += sizeof(degen_id);

			memcpy( m_ptr, &seq_index, sizeof(seq_index) );
			m_ptr += sizeof(seq_index);
			
			memcpy( m_ptr, name.c_str(), name.size() + 1 );
			m_ptr += name.size() + 1;
			
			memcpy( m_ptr, forward_oligo.c_str(), forward_oligo.size() + 1 );
			m_ptr += forward_oligo.size() + 1;
			
			memcpy( m_ptr, reverse_oligo.c_str(), reverse_oligo.size() + 1 );
			m_ptr += reverse_oligo.size() + 1;
			
			memcpy( m_ptr, probe_oligo.c_str(), probe_oligo.size() + 1 );
			m_ptr += probe_oligo.size() + 1;
			
			memcpy( m_ptr, amplicon_def.c_str(), amplicon_def.size() + 1 );
			m_ptr += amplicon_def.size() + 1;
			
			memcpy( m_ptr, amplicon.c_str(), amplicon.size() + 1 );
			m_ptr += amplicon.size() + 1;
			
			memcpy( m_ptr, forward_align.c_str(), forward_align.size() + 1 );
			m_ptr += forward_align.size() + 1;
			
			memcpy( m_ptr, reverse_align.c_str(), reverse_align.size() + 1 );
			m_ptr += reverse_align.size() + 1;
			
			memcpy( m_ptr, probe_align.c_str(), probe_align.size() + 1 );
			m_ptr += probe_align.size() + 1;
			
			memcpy( m_ptr, &(amplicon_range.first), sizeof(int) );
			m_ptr += sizeof(amplicon_range.first);
			
			memcpy( m_ptr, &(amplicon_range.second), sizeof(int) );
			m_ptr += sizeof(amplicon_range.second);
			
			memcpy( m_ptr, &(probe_range.first), sizeof(int) );
			m_ptr += sizeof(probe_range.first);
			
			memcpy( m_ptr, &(probe_range.second), sizeof(int) );
			m_ptr += sizeof(probe_range.second);
			
			///////////////////////////////////////////////////
			// Melting temperatures
			memcpy( m_ptr, &forward_tm, sizeof(forward_tm) );
			m_ptr += sizeof(forward_tm);
			
			memcpy( m_ptr, &reverse_tm, sizeof(reverse_tm) );
			m_ptr += sizeof(reverse_tm);
			
			memcpy( m_ptr, &probe_tm, sizeof(probe_tm) );
			m_ptr += sizeof(probe_tm);
			
			///////////////////////////////////////////////////
			// Dimer and hairpin melting temperatures
			memcpy( m_ptr, &forward_hairpin_tm, sizeof(forward_hairpin_tm) );
			m_ptr += sizeof(forward_hairpin_tm);
			
			memcpy( m_ptr, &reverse_hairpin_tm, sizeof(reverse_hairpin_tm) );
			m_ptr += sizeof(reverse_hairpin_tm);
			
			memcpy( m_ptr, &forward_dimer_tm, sizeof(forward_dimer_tm) );
			m_ptr += sizeof(forward_dimer_tm);
			
			memcpy( m_ptr, &reverse_dimer_tm, sizeof(reverse_dimer_tm) );
			m_ptr += sizeof(reverse_dimer_tm);
			
			memcpy( m_ptr, &primer_dimer_tm, sizeof(primer_dimer_tm) );
			m_ptr += sizeof(primer_dimer_tm);
			
			memcpy( m_ptr, &probe_hairpin_tm, sizeof(probe_hairpin_tm) );
			m_ptr += sizeof(probe_hairpin_tm);
			
			memcpy( m_ptr, &probe_dimer_tm, sizeof(probe_dimer_tm) );
			m_ptr += sizeof(probe_dimer_tm);
			
			///////////////////////////////////////////////////
			// Thermodynamic parameters
			memcpy( m_ptr, &forward_dH, sizeof(forward_dH) );
			m_ptr += sizeof(forward_dH);
			
			memcpy( m_ptr, &forward_dS, sizeof(forward_dS) );
			m_ptr += sizeof(forward_dS);
			
			memcpy( m_ptr, &reverse_dH, sizeof(reverse_dH) );
			m_ptr += sizeof(reverse_dH);
			
			memcpy( m_ptr, &reverse_dS, sizeof(reverse_dS) );
			m_ptr += sizeof(reverse_dS);
			
			memcpy( m_ptr, &probe_dH, sizeof(probe_dH) );
			m_ptr += sizeof(probe_dH);
			
			memcpy( m_ptr, &probe_dS, sizeof(probe_dS) );
			m_ptr += sizeof(probe_dS);
			
			///////////////////////////////////////////////////
			// Strand concentrations
			memcpy( m_ptr, &primer_strand, sizeof(primer_strand) );
			m_ptr += sizeof(primer_strand);
			
			memcpy( m_ptr, &probe_strand, sizeof(probe_strand) );
			m_ptr += sizeof(probe_strand);
			
			///////////////////////////////////////////////////
			// Clamp values
			memcpy( m_ptr, &forward_primer_clamp, sizeof(forward_primer_clamp) );
			m_ptr += sizeof(forward_primer_clamp);
			
			memcpy( m_ptr, &reverse_primer_clamp, sizeof(reverse_primer_clamp) );
			m_ptr += sizeof(reverse_primer_clamp);
			
			///////////////////////////////////////////////////
			// Degeneracy levels
			memcpy( m_ptr, &forward_degen, sizeof(forward_degen) );
			m_ptr += sizeof(forward_degen);

			memcpy( m_ptr, &reverse_degen, sizeof(reverse_degen) );
			m_ptr += sizeof(reverse_degen);

			memcpy( m_ptr, &probe_degen, sizeof(probe_degen) );
			m_ptr += sizeof(probe_degen);

			///////////////////////////////////////////////////
			// Mismatch and gap parameters
			memcpy( m_ptr, &forward_mm, sizeof(forward_mm) );
			m_ptr += sizeof(forward_mm);
			
			memcpy( m_ptr, &forward_gap, sizeof(forward_gap) );
			m_ptr += sizeof(forward_gap);
			
			memcpy( m_ptr, &reverse_mm, sizeof(reverse_mm) );
			m_ptr += sizeof(reverse_mm);
			
			memcpy( m_ptr, &reverse_gap, sizeof(reverse_gap) );
			m_ptr += sizeof(reverse_gap);
			
			memcpy( m_ptr, &probe_mm, sizeof(probe_mm) );
			m_ptr += sizeof(probe_mm);
			
			memcpy( m_ptr, &probe_gap, sizeof(probe_gap) );
			m_ptr += sizeof(probe_gap);
			
			///////////////////////////////////////////////////
			// Ct
			memcpy( m_ptr, &ct, sizeof(ct) );
			m_ptr += sizeof(ct);
			
			return m_ptr;
		};
		
		inline unsigned char* mpi_unpack(unsigned char* m_ptr)
		{
			memcpy( &id, m_ptr, sizeof(id) );
			m_ptr += sizeof(id);

			memcpy( &degen_id, m_ptr, sizeof(degen_id) );
			m_ptr += sizeof(degen_id);
			
			memcpy( &seq_index, m_ptr, sizeof(seq_index) );
			m_ptr += sizeof(seq_index);
			
			name = (char*)m_ptr;
			m_ptr += name.size() + 1;
			
			forward_oligo = (char*)m_ptr;
			m_ptr += forward_oligo.size() + 1;
			
			reverse_oligo = (char*)m_ptr;
			m_ptr += reverse_oligo.size() + 1;
			
			probe_oligo = (char*)m_ptr;
			m_ptr += probe_oligo.size() + 1;
			
			amplicon_def = (char*)m_ptr;
			m_ptr += amplicon_def.size() + 1;
			
			amplicon = (char*)m_ptr;
			m_ptr += amplicon.size() + 1;
			
			forward_align = (char*)m_ptr;
			m_ptr += forward_align.size() + 1;
			
			reverse_align = (char*)m_ptr;
			m_ptr += reverse_align.size() + 1;
			
			probe_align = (char*)m_ptr;
			m_ptr += probe_align.size() + 1;
			
			memcpy( &(amplicon_range.first), m_ptr, sizeof(int) );
			m_ptr += sizeof(amplicon_range.first);
			
			memcpy( &(amplicon_range.second), m_ptr, sizeof(int) );
			m_ptr += sizeof(amplicon_range.second);
			
			memcpy( &(probe_range.first), m_ptr, sizeof(int) );
			m_ptr += sizeof(probe_range.first);
			
			memcpy( &(probe_range.second), m_ptr, sizeof(int) );
			m_ptr += sizeof(probe_range.second);
			
			///////////////////////////////////////////////////
			// Melting temperatures
			memcpy( &forward_tm, m_ptr, sizeof(forward_tm) );
			m_ptr += sizeof(forward_tm);
			
			memcpy( &reverse_tm, m_ptr, sizeof(reverse_tm) );
			m_ptr += sizeof(reverse_tm);
			
			memcpy( &probe_tm, m_ptr, sizeof(probe_tm) );
			m_ptr += sizeof(probe_tm);
			
			///////////////////////////////////////////////////
			// Dimer and hairpin melting temperatures
			memcpy( &forward_hairpin_tm, m_ptr, sizeof(forward_hairpin_tm) );
			m_ptr += sizeof(forward_hairpin_tm);
			
			memcpy( &reverse_hairpin_tm, m_ptr, sizeof(reverse_hairpin_tm) );
			m_ptr += sizeof(reverse_hairpin_tm);
			
			memcpy( &forward_dimer_tm, m_ptr, sizeof(forward_dimer_tm) );
			m_ptr += sizeof(forward_dimer_tm);
			
			memcpy( &reverse_dimer_tm, m_ptr, sizeof(reverse_dimer_tm) );
			m_ptr += sizeof(reverse_dimer_tm);
			
			memcpy( &primer_dimer_tm, m_ptr, sizeof(primer_dimer_tm) );
			m_ptr += sizeof(primer_dimer_tm);
			
			memcpy( &probe_hairpin_tm, m_ptr, sizeof(probe_hairpin_tm) );
			m_ptr += sizeof(probe_hairpin_tm);
			
			memcpy( &probe_dimer_tm, m_ptr, sizeof(probe_dimer_tm) );
			m_ptr += sizeof(probe_dimer_tm);
			
			///////////////////////////////////////////////////
			// Thermodynamic parameters
			memcpy( &forward_dH, m_ptr, sizeof(forward_dH) );
			m_ptr += sizeof(forward_dH);
			
			memcpy( &forward_dS, m_ptr, sizeof(forward_dS) );
			m_ptr += sizeof(forward_dS);
			
			memcpy( &reverse_dH, m_ptr, sizeof(reverse_dH) );
			m_ptr += sizeof(reverse_dH);
			
			memcpy( &reverse_dS, m_ptr, sizeof(reverse_dS) );
			m_ptr += sizeof(reverse_dS);
			
			memcpy( &probe_dH, m_ptr, sizeof(probe_dH) );
			m_ptr += sizeof(probe_dH);
			
			memcpy( &probe_dS, m_ptr, sizeof(probe_dS) );
			m_ptr += sizeof(probe_dS);
			
			///////////////////////////////////////////////////
			// Strand concentrations
			memcpy( &primer_strand, m_ptr, sizeof(primer_strand) );
			m_ptr += sizeof(primer_strand);
			
			memcpy( &probe_strand, m_ptr, sizeof(probe_strand) );
			m_ptr += sizeof(probe_strand);
			
			///////////////////////////////////////////////////
			// Clamp values
			memcpy( &forward_primer_clamp, m_ptr, sizeof(forward_primer_clamp) );
			m_ptr += sizeof(forward_primer_clamp);
			
			memcpy( &reverse_primer_clamp, m_ptr, sizeof(reverse_primer_clamp) );
			m_ptr += sizeof(reverse_primer_clamp);
			
			///////////////////////////////////////////////////
			// Degeneracy levels
			memcpy( &forward_degen, m_ptr, sizeof(forward_degen) );
			m_ptr += sizeof(forward_degen);

			memcpy( &reverse_degen, m_ptr, sizeof(reverse_degen) );
			m_ptr += sizeof(reverse_degen);

			memcpy( &probe_degen, m_ptr, sizeof(probe_degen) );
			m_ptr += sizeof(probe_degen);

			///////////////////////////////////////////////////
			// Mismatch and gap parameters
			memcpy( &forward_mm, m_ptr, sizeof(forward_mm) );
			m_ptr += sizeof(forward_mm);
			
			memcpy( &forward_gap, m_ptr, sizeof(forward_gap) );
			m_ptr += sizeof(forward_gap);
			
			memcpy( &reverse_mm, m_ptr, sizeof(reverse_mm) );
			m_ptr += sizeof(reverse_mm);
			
			memcpy( &reverse_gap, m_ptr, sizeof(reverse_gap) );
			m_ptr += sizeof(reverse_gap);
			
			memcpy( &probe_mm, m_ptr, sizeof(probe_mm) );
			m_ptr += sizeof(probe_mm);
			
			memcpy( &probe_gap, m_ptr, sizeof(probe_gap) );
			m_ptr += sizeof(probe_gap);
			
			///////////////////////////////////////////////////
			// Ct
			memcpy( &ct, m_ptr, sizeof(ct) );
			m_ptr += sizeof(ct);
			
			return m_ptr;
		};
		
		inline int my_id() const
		{
			return id;
		};
		
		inline int my_degen_id() const
		{
			return degen_id;
		};

		inline void my_degen_id(const int &m_id)
		{
			// Id's must be >= 0
			if(m_id < 0){
				throw __FILE__ ":my_degen_id: m_id < 0";
			}
			
			degen_id = m_id;
		};
		
		inline void seq_id(const int &m_id)
		{
			seq_index = m_id;
		};
		
		inline int seq_id() const
		{
			return seq_index;
		};
		
		inline bool has_primers() const
		{
			return ( !forward_oligo.empty() && !reverse_oligo.empty() );
		};
		
		inline bool has_probe() const
		{
			return ( !probe_oligo.empty() );
		};
		
		inline bool probe_overlap(const hybrid_sig &m_sig) const
		{			
			// Overlapping probes must be on the same strand
			if(probe_strand != m_sig.probe_strand){
				return false;
			}
			
			// Is the first edge of this probe within the bounds of m_sig?
			if( (probe_range.first >= m_sig.probe_range.first) && 
			    (probe_range.first <= m_sig.probe_range.second) ){
			    
			    return true;
			}
			
			// Is the second edge of this probe within the bounds of m_sig?
			if( (probe_range.second >= m_sig.probe_range.first) && 
			    (probe_range.second <= m_sig.probe_range.second) ){
			    
			    return true;
			}
			
			// Does this probe contain m_sig?
			if( (probe_range.first <= m_sig.probe_range.first) && 
			    (probe_range.second >= m_sig.probe_range.second) ){
			    
			    return true;
			}
			
			// If we get here, the probes don't overlap
			return false;
		};

		inline bool operator>(const hybrid_sig &m_rhs) const
		{
			if(id == m_rhs.id){
			
				// Sort (in descending order) by:
				// 1) minimum primer melting temperature
				// 2) probe melting temperature
				// 3) maximum primer melting temperature
				// 4) target sequence index
				if( min_primer_tm() == m_rhs.min_primer_tm() ){
					
					if(probe_tm == m_rhs.probe_tm){
						
						if( max_primer_tm() == m_rhs.max_primer_tm() ){
						
							// Sort by target sequence
							return seq_index > m_rhs.seq_index;
						}
						
						return ( max_primer_tm() < m_rhs.max_primer_tm() );
					}
					
					return (probe_tm < m_rhs.probe_tm);
				}
				
				return ( min_primer_tm() < m_rhs.min_primer_tm() );
			}
			
			return (id > m_rhs.id);
		};

		inline bool operator<(const hybrid_sig &m_rhs) const
		{
			if(id == m_rhs.id){
			
				// Sort (in descending order) by:
				// 1) minimum primer melting temperature
				// 2) probe melting temperature
				// 3) maximum primer melting temperature
				// 4) target sequence index
				if( min_primer_tm() == m_rhs.min_primer_tm() ){
				
					if(probe_tm == m_rhs.probe_tm){
						
						if( max_primer_tm() == m_rhs.max_primer_tm() ){
						
							// Sort by target sequence
							return seq_index < m_rhs.seq_index;
						}
						
						return ( max_primer_tm() > m_rhs.max_primer_tm() );
					}
					
					return (probe_tm > m_rhs.probe_tm);
				}
				
				return ( min_primer_tm() > m_rhs.min_primer_tm() );
			}
			
			return (id < m_rhs.id);
		};
		
		inline int min_primer_clamp() const
		{
			return std::min(forward_primer_clamp, reverse_primer_clamp);
		};
		
		inline int max_primer_clamp() const
		{
			return std::max(forward_primer_clamp, reverse_primer_clamp);
		};
		
		inline float min_primer_tm() const
		{
			// Clamp the min Tm at zero
			return std::max(0.0f, std::min(forward_tm, reverse_tm) );
		};
		
		inline float max_primer_tm() const
		{
			return std::max(forward_tm, reverse_tm);
		};
		
		inline void offset_ranges(const int &m_off)
		{
		
			if(has_primers() == true){
			
				amplicon_range.first += m_off;
				amplicon_range.second += m_off;
			}
			
			if(has_probe() == true){
			
				probe_range.first += m_off;
				probe_range.second += m_off;
			}
		};
		
		inline bool start_overlap(const int &m_start) const
		{
			if(has_primers() == true){
			
				return (amplicon_range.first <= m_start);
			}
			else{ // has_probe == true
			
				return (probe_range.first <= m_start);
			}
		};
		
		inline bool stop_overlap(const int &m_stop) const
		{
			if(has_primers() == true){
			
				return (amplicon_range.second >= m_stop);
			}
			else{ // has_probe == true
			
				return (probe_range.second >= m_stop);
			}
		};
		
		inline std::string assay_string() const
		{
			return name;
			
			//std::stringstream ssout;
			//
			//ssout << name << '\t';
			//
			//if(forward_oligo.empty() == false){
			//	ssout << forward_oligo << '\t';
			//}
			//
			//if(reverse_oligo.empty() == false){
			//	ssout << reverse_oligo << '\t';
			//}
			//
			//if(probe_oligo.empty() == false){
			//	ssout << probe_oligo << '\t';
			//}
			//
			//return ssout.str();
		};
};

class sort_by_match{

	public:
		sort_by_match() {};
		~sort_by_match() {};

		inline bool operator()(const hybrid_sig &m_a, const hybrid_sig &m_b)
		{
			if( m_a.my_id() == m_b.my_id() ){
				return ( m_a.seq_id() < m_b.seq_id() );
			}
			
			return ( m_a.my_id() < m_b.my_id() );
		};
};

class sort_by_loc{

	public:
	
		sort_by_loc() {};
		~sort_by_loc() {};
		
		inline bool operator()(const hybrid_sig &m_a, const hybrid_sig &m_b)
		{
			// First sort by id
			if( m_a.my_id() == m_b.my_id() ){
			
				// Then by target sequence
				if( m_a.seq_id() == m_b.seq_id() ){
					
					// Then by range
					if( m_a.has_primers() == true){
						
						return m_a.amplicon_range < m_b.amplicon_range;
					}
					else{ // m_a.has_probe() == true
						
						return m_a.probe_range < m_b.probe_range;
					}
				}
				
				return ( m_a.seq_id() < m_b.seq_id() );
			}
			
			return ( m_a.my_id() < m_b.my_id() );
		};
};

class unique_by_loc{

	public:
	
		unique_by_loc() {};
		~unique_by_loc() {};
		
		inline bool operator()(const hybrid_sig &m_a, const hybrid_sig &m_b)
		{
			// First sort by id
			if( m_a.my_id() == m_b.my_id() ){
			
				// Then by target sequence
				if( m_a.seq_id() == m_b.seq_id() ){
					
					// Then by range
					if( m_a.has_primers() == true){
						
						return m_a.amplicon_range == m_b.amplicon_range;
					}
					else{ // m_a.has_probe() == true
						
						return m_a.probe_range == m_b.probe_range;
					}
				}
				
				return false;
			}
			
			return false;
		};
};

size_t read_input_file(const std::string &m_file, std::vector<hybrid_sig> &m_sig_list, 
	const bool &m_ignore_probe, const bool &m_force_probe);

#endif // __SIGNATURE_CLASS
