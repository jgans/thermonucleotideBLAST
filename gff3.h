// ThermonucleotideBLAST
// 
// Copyright (c) 2008, Los Alamos National Security, LLC
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
#ifndef	__GFF3_ANNOTATION
#define	__GFF3_ANNOTATION

#include <string>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <list>

#define	GFF3_DELIM				'\t'
#define	GFF3_ATTRIBUTE_DELIM	';'
#define	GFF3_TAG_VALUE_DELIM	'='
#define	GFF3_VALUE_DELIM		','

// Since all types are enumerated from 0, a value of -1 allows us to define an invalid type 
// that will not overlap future valid SOFA types.
#define INVALID_TYPE			-1

struct GFFSegment {

	std::pair<int, int> range;
	int phase;

	GFFSegment()
	{
		// Do nothing
	};

	GFFSegment(const GFFSegment &m_rhs)
	{
		range = m_rhs.range;
		phase = m_rhs.phase;
	};

	inline bool operator < (const GFFSegment &m_rhs) const
	{
		if(range.first < m_rhs.range.first){
			return true;
		}

		if(range.first > m_rhs.range.first){
			return false;
		}

		if(range.second < m_rhs.range.second){
			return true;
		}

		if(range.second > m_rhs.range.second){
			return false;
		}

		return (phase < m_rhs.phase);
	};

	inline bool operator == (const GFFSegment &m_rhs) const
	{
		return (range.first == m_rhs.range.first) && (range.second == m_rhs.range.second) &&
			(phase == m_rhs.phase);
	};
};

class GFF3Record {
	
	private:
		
		double rec_score;
		std::string rec_source;
		int rec_type;
		unsigned char rec_strand;
		std::list<GFFSegment> rec_seg; // Allow for annotations with multiple segments

		std::multimap<std::string, std::string> attrib;
		
		static std::string sofa_version;
		static std::map<std::string, int> sofa;
		
		std::string parse_attrib(const std::string &m_line_buffer);
		void parse_attrib_pair(const std::string &m_line_buffer, std::string &m_id);
		
		inline int parse_SOFA(std::string m_type) const
		{
			// Convert the search string to upper case
			for(std::string::iterator i = m_type.begin();i != m_type.end();i++){
				*i = toupper(*i);
			}
			
			std::map<std::string, int>::const_iterator iter = sofa.find(m_type);

			if( iter != sofa.end() ){
				return iter->second;
			}

			// Could not find a mapping, so return unknown
			return UNKNOWN;
		};
		
		void init();
	public:
		
		typedef std::multimap<std::string, std::string>::const_iterator const_iterator;
		typedef std::multimap<std::string, std::string>::iterator iterator;
		
		// Include the enumerated GFF3 annotation feature types
		#include "gff3_sofa.h"
			
		enum {GFF3_PLUS_STRAND, GFF3_MINUS_STRAND, GFF3_NO_STRAND, 
			GFF3_UKNOWN_STRAND};
			
		GFF3Record()
		{
			rec_score = 0.0;
			rec_type = INVALID_TYPE;
		};
		
		// Create a GFF3Record by parsing an input line
		GFF3Record(const std::string &m_line_buffer, 
			std::string &m_seqid, std::string &m_id)
		{
			rec_score = 0.0;
			rec_type = INVALID_TYPE;

			parse(m_line_buffer, m_seqid, m_id);
		};
		
		GFF3Record(const GFF3Record &m_copy)
		{
			*this = m_copy;
		};
		
		inline GFF3Record& operator=(const GFF3Record &m_rhs)
		{
			rec_score = m_rhs.rec_score;
			rec_source = m_rhs.rec_source;
			rec_type = m_rhs.rec_type;
			rec_seg = m_rhs.rec_seg;
			rec_strand = m_rhs.rec_strand;
		
			attrib = m_rhs.attrib;
			return *this;
		};
		
		// The "+=" operator has been overloaded to allow for concatination of multiline
		// records (i.e. records with the same id that have been split across multiple lines).
		GFF3Record& operator+=(const GFF3Record &m_rhs);

		void parse(const std::string &m_line_buffer, 
			std::string &m_seqid, std::string &m_id);
		
		inline std::pair<const_iterator, const_iterator> operator()(const std::string &m_key) const
		{
			return attrib.equal_range(m_key);
		};
		
		inline int feature_type() const
		{
			return rec_type;
		};
		
		inline double feature_score() const
		{
			return rec_score;
		};
		
		inline std::string feature_source() const
		{
			return rec_source;
		};
		
		inline unsigned char feature_strand() const
		{
			return rec_strand;
		};
		
		inline const std::list<GFFSegment>& feature_seg() const
		{
			return rec_seg;
		};
		
		inline std::pair<int, int> feature_range() const
		{
			if(rec_seg.empty() == true){
				throw "feature_range: empty seg list";
			}

			return std::make_pair( rec_seg.front().range.first, 
				rec_seg.back().range.second );
		};
		
		inline int feature_start() const
		{
			if(rec_seg.empty() == true){
				throw "feature_start: empty seg list";
			}

			return rec_seg.front().range.first;
		};
		
		inline int feature_stop() const
		{
			if(rec_seg.empty() == true){
				throw "feature_stop: empty seg list";
			}

			return rec_seg.back().range.second;
		};
		
		inline int feature_len() const
		{
			if(rec_seg.empty() == true){
				throw "feature_len: empty seg list";
			}

			return rec_seg.back().range.second - rec_seg.front().range.first + 1;
		};
		
		inline bool has_parent() const
		{
			const std::pair<const_iterator, const_iterator> range = 
				attrib.equal_range("Parent");
			
			return (range.first != range.second);
		};
		
		inline std::list<std::string> parent() const
		{
			const std::pair<const_iterator, const_iterator> range = 
				attrib.equal_range("Parent");
			
			std::list<std::string> ret;
			
			for(const_iterator i = range.first;i != range.second;i++){
				ret.push_back(i->second);
			}
			
			return ret;
		};
};

class GFF3File {
	private:
				
		std::string line_buffer;
		
		int pragma_version;
		std::string pragma_feature_ontology;
		std::string pragma_attribute_ontology;
		std::string pragma_source_ontology;
		std::string pragma_species;
		std::string pragma_genome_build_source;
		std::string pragma_genome_build_name;
		
		std::map< std::string, std::pair<int, int> > pragma_sequence_region;
		
		// The high level maps are indexed by seqid. The lower level multimap is
		// indexed by (unique) ID
		std::map<std::string, std::string> seq;
		std::map<std::string, std::map<std::string, GFF3Record> > features;	
		
		void parse(std::ifstream &m_fin);
		void parse_line(std::ifstream &m_fin);
		void parse_pragma(std::ifstream &m_fin);
		void parse_fasta(std::ifstream &m_fin);
		
	public:
	
	GFF3File()
	{
		pragma_version = -1;
	};
	
	GFF3File(const std::string &m_file)
	{
		pragma_version = -1;
		
		std::ifstream fin( m_file.c_str() );
		
		if(!fin){
			throw "Unable to open GFF3 file";
		}
		
		parse(fin);
		
		fin.close();
	};
	
	GFF3File(std::ifstream &m_fin)
	{
		pragma_version = -1;
		
		parse(m_fin);
	};
	
	inline void clear()
	{
		line_buffer.clear();
		seq.clear();
		features.clear();
	};
	
	int version() const
	{
		return pragma_version;
	};
	
	// Is this file valid?
	inline operator bool() const
	{
		// Only version 3 is currently defined
		return (pragma_version == 3);
	};
	
	std::vector<std::string> feature_source() const
	{
		std::set<std::string> tmp;
		
		typedef std::map<std::string, std::map<std::string, GFF3Record> >::const_iterator I;
		
		for(I iter = features.begin();iter != features.end();iter++){
			tmp.insert(iter->first);
		}
		
		std::vector<std::string> ret( tmp.size() );
		
		unsigned int index = 0;
		
		for(std::set<std::string>::const_iterator i = tmp.begin();i != tmp.end();i++){
			ret[index++] = *i;
		}
		
		return ret;
	};
	
	// Allow the user to test for the presence of associated sequence data
	inline bool has_sequence(const std::string &m_source) const
	{
		return ( seq.find(m_source) != seq.end() );
	};
	
	// Allow the user to test for the presence of associated feature (i.e. annotation) data
	inline bool has_features(const std::string &m_source) const
	{
		return ( features.find(m_source) != features.end() );
	};
	
	inline const std::string& sequence(const std::string &m_source) const
	{
		std::map<std::string, std::string>::const_iterator iter = seq.find(m_source);
		
		if( iter == seq.end() ){
			throw "No sequence found for specified source";
		}
		
		return iter->second;
	};
	
	inline const std::map<std::string, GFF3Record>& feature_map(const std::string &m_source) const
	{
		std::map<std::string, std::map<std::string, GFF3Record> >::const_iterator iter = features.find(m_source);
		
		if( iter == features.end() ){
			throw "No features found for specified source";
		}
		
		return iter->second;
	};
	
	inline std::vector<std::string> feature_id(const std::string &m_source) const
	{
		std::map<std::string, std::map<std::string, GFF3Record> >::const_iterator iter = features.find(m_source);
		
		if( iter == features.end() ){
			throw "No features found for specified source";
		}
		
		typedef std::map<std::string, GFF3Record>::const_iterator I;
		
		std::vector<std::string> ret( iter->second.size() );
		
		unsigned int index = 0;
		
		for(I i = iter->second.begin();i != iter->second.end();i++){
			ret[index++] = i->first;
		}
		
		return ret;
	};
};

#endif // __GFF3_ANNOTATION
