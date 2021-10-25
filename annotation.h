#ifndef __GENOME_ANNOTATION
#define	__GENOME_ANNOTATION

#include "seq.h"
#include <zlib.h>
#include <string>
#include <vector>
#include <deque>
#include <list>
#include <iostream>
#include <sstream>

// OS specific defines and includes
#ifdef WIN32

// Handle the annoying carrige return problem in windows
#define ENDL	"\r" << endl
#define EOL		"\r\n"

#else

#define ENDL	endl
#define EOL	"\n"

#endif // WIN32


class GeneAnnotation{
private:
	unsigned int gene_type;
	
	bool complement; // Is this gene on the opposite strand?
	unsigned int gene_start, gene_stop;
	
	// Is this gene divided into multiple regions?
	std::list< std::pair<unsigned int, unsigned int> > seg_list;
	
	// For many years tntblast used the NCBI C toolkit SeqIdPtr to track
	// sequence ids. Now, just store a deque of all of the accessions
	// associated with a given gene
	std::deque<std::string> id;
	
	std::vector<std::string> info_map;
	
public:

	enum {LOCUS = 0, LOCUS_TAG, PRODUCT, EC, NOTE, 
		PRIMER_TM, PRIMER_GC, PRIMER_HAIRPIN_TM, PRIMER_DIMER_TM, 
		COG_CODE, COG_ID,
		LAST_INFO}; // info
	enum {CDS = 0, GENE, PSEUDO_GENE, RNA, tRNA, IMP, PRIMER, TFBS, USER, NONE,
		LAST_ANNOT}; // gene_type
		
	GeneAnnotation()
	{
		complement = false;
		gene_type = NONE;
		gene_start = gene_stop = 0;
		id.clear();
		info_map = std::vector<std::string>(LAST_INFO);
	};
	
	GeneAnnotation(const GeneAnnotation &m_copy)
	{
		complement = false;
		gene_type = NONE;
		gene_start = gene_stop = 0;

		*this = m_copy;
	};
	
	~GeneAnnotation()
	{
		
	};
	
	inline void clear()
	{		
		complement = false;
		gene_type = NONE;
		gene_start = gene_stop = 0;
		id.clear();
		info_map = std::vector<std::string>(LAST_INFO);
		seg_list.clear();
	}

	GeneAnnotation& operator=(const GeneAnnotation& m_copy)
	{
		gene_type = m_copy.gene_type;
		complement = m_copy.complement;
		gene_start = m_copy.gene_start;
		gene_stop = m_copy.gene_stop;
		id = m_copy.id;
		info_map = m_copy.info_map;
		seg_list = m_copy.seg_list;

		return (*this);
	};

	inline bool operator<(const GeneAnnotation& m_copy) const
	{		
		const bool this_overlap = (gene_start > gene_stop);
		const bool copy_overlap = (m_copy.gene_start > m_copy.gene_stop);

		// Do either this annotation or the copy annotation (but not 
		// both) overlap the origin?
		if( (this_overlap || copy_overlap) && !(this_overlap && copy_overlap) ){

			// One of the two annotations overlaps the origin. By definition, this
			// annotation is "less than" the other.
			return this_overlap;
		}
		
		if(gene_start < m_copy.gene_start){
			return true;
		}
		
		if(gene_start > m_copy.gene_start){
			return false;
		}
		
		// start == m_copy.start
		return gene_stop < m_copy.gene_stop;
	};
	
	inline bool operator>(const GeneAnnotation& m_copy) const
	{		
		const bool this_overlap = (gene_start > gene_stop);
		const bool copy_overlap = (m_copy.gene_start > m_copy.gene_stop);

		// Do either this annotation or the copy annotation (but not 
		// both) overlap the origin?
		if( (this_overlap || copy_overlap) && !(this_overlap && copy_overlap) ){

			// One of the two annotations overlaps the origin. By definition, the
			// annotation that does NOT overlap the origin is "greater than" the other.
			return copy_overlap;
		}

		if(gene_start > m_copy.gene_start){
			return true;
		}
		
		if(gene_start < m_copy.gene_start){
			return false;
		}
		
		// start == m_copy.start
		return gene_stop > m_copy.gene_stop;
	};

	inline bool overlaps_origin() const
	{
		return (gene_start > gene_stop);
	};
	
	inline void add_seqid(const std::string &m_id)
	{
		id.push_back(m_id);
	};
	
	inline std::string seq_id_str() const
	{
		std::string ret;

		// Concatinate all of the accession strings into a single string
		for(std::deque<std::string>::const_iterator i = id.begin();i != id.end();++i){

			if( i == id.begin() ){
				ret = *i;
			}
			else{
				ret = ret + "|" + *i;
			}
		}

		return ret;
	};

	inline void info(unsigned int m_key, const char *m_str)
	{
		if(m_str != NULL){

			if(m_key >= LAST_INFO){
				throw "info:: key out of bounds";
			}

			info_map[m_key] = m_str;
		}
	};
	
	inline void info(unsigned int m_key, const std::string &m_str)
	{
		if(m_key >= LAST_INFO){
			throw "info:: key out of bounds";
		}

		info_map[m_key] = m_str;
	};
	
	inline std::string info(unsigned int m_key) const
	{	
		return info_map[m_key];
	};
	
	inline std::string name() const
	{
		return (info_map[LOCUS] != "") ? info_map[LOCUS] : info_map[LOCUS_TAG];
	};

	inline void primer_id(const unsigned int &m_id)
	{
		std::stringstream sout;

		sout << "P" << m_id;

		info_map[LOCUS] = sout.str();
	};

	inline void primer_tm(const float &m_tm)
	{
		std::stringstream sout;

		sout << "Tm = " << m_tm;

		info_map[PRIMER_TM] = sout.str();
	};

	inline void primer_dimer_tm(const float &m_tm)
	{
		std::stringstream sout;

		sout << "Dimer Tm = " << m_tm;

		info_map[PRIMER_DIMER_TM] = sout.str();
	};

	inline void primer_hairpin_tm(const float &m_tm)
	{
		std::stringstream sout;

		sout << "Hairpin Tm = " << m_tm;

		info_map[PRIMER_HAIRPIN_TM] = sout.str();
	};

	inline void primer_gc(const float &m_gc)
	{
		std::stringstream sout;

		sout << "%G+C = " << 100.0f*m_gc;

		info_map[PRIMER_GC] = sout.str();
	};

	// Does this annotation encode a protein?
	inline bool is_protein() const
	{
		return (gene_type < RNA);
	};

	inline bool is_gene() const
	{
		// This function assumes that all types greater than tRNA are
		// NOT genes!
		return (gene_type <= tRNA);
	};

	inline bool is_intergenic() const
	{
		return (gene_type == NONE);
	};
	
	inline bool has_id_str() const
	{
		return !((gene_type == IMP) || (gene_type == NONE));
	};

	inline bool is_complement() const
	{
		return complement;
	};

	inline void is_complement(const bool &m_comp)
	{
		complement = m_comp;
	};

	inline void start(const unsigned int &m_start)
	{
		gene_start = m_start;
	};
	
	inline  unsigned int start() const
	{
		return gene_start;
	};
	
	inline void stop(const unsigned int &m_stop)
	{
		gene_stop = m_stop;
	};
	
	inline unsigned int stop() const
	{
		return gene_stop;
	};

	inline unsigned int length(const unsigned int &m_seq_len) const
	{
		if(gene_start > gene_stop){
			
			if(m_seq_len < gene_start){
				throw __FILE__ ":length: gene_start is out of bounds!";
			}
			
			return (m_seq_len - gene_start) + gene_stop + 1;
		}

		return gene_stop - gene_start + 1;
	};
	
	inline void strand(const unsigned int &m_strand)
	{
		complement = (m_strand == Seq_strand_minus);
	};
	
	inline unsigned int strand() const
	{
		if(complement){
			return Seq_strand_minus;
		}
		
		return Seq_strand_plus;
	};
	
	inline void type(const unsigned int &m_type)
	{
		if(m_type > NONE){
			throw "GeneAnnotation: Unknown gene type";
		}
		
		gene_type = m_type;
	};
	
	inline unsigned int type() const
	{
		return gene_type;
	};
	
	inline void segments(const std::list< std::pair<unsigned int, unsigned int> > &m_seg_list)
	{
		seg_list = m_seg_list;
	};

	inline bool is_segmented() const
	{
		return (seg_list.empty() == false);
	};

	inline const std::list< std::pair<unsigned int, unsigned int> >& segments() const
	{
		return seg_list;
	};
	
	// Turn a range or segment string into a start/stop and segment list
	void range(const std::string &m_range);
	
	// Using the start/stop and segment list, return a valid range string
	std::string range() const;
	
	// Handle the special case of a gene that overlaps the genome start location
	inline bool handle_gene_start_overlap(const unsigned int &m_genome_len)
	{
		if( (gene_start == 0) && (seg_list.empty() == false) ){

			// Look for a segment stop that is equal to m_genome_len - 1
			const unsigned int tmp_len = m_genome_len - 1;
			std::list< std::pair<unsigned int, unsigned int> >::const_iterator iter;
			
			unsigned int tmp_start = 0;
			unsigned int tmp_stop = 0;

			for(iter = seg_list.begin();iter != seg_list.end();iter++){
				
				if(iter->first == 0){
					tmp_stop = iter->second;
				}

				if(iter->second == tmp_len){
					tmp_start = iter->first;
				}
			}

			if((tmp_start != 0) && (tmp_stop != 0)){
				// We found a match! Set the new start and stop files
				// and clear the seq_list.
				seg_list.clear();
				gene_start = tmp_start;
				gene_stop = tmp_stop;

				return true;
			}

		}
		
		if(gene_stop >= m_genome_len){
			gene_stop = gene_stop - m_genome_len;
			
			return true;
		}
		
		if(gene_start > gene_stop){
			return true;
		}
		
		return false;
	};
};

// Classes to manage an entire DNA molecule
class DNAMol{
private:
	std::string accession;

	// The sequence data
	SEQPTR seq;
	
	// The number of bases stored in "seq"
	// (note that this does *not* include any prefix or suffix
	// bases)
	unsigned int seq_len;
	
	// From the SeqLoc:	
	unsigned int start, stop;
	
	std::vector<std::string> info_map;
	
	std::list<GeneAnnotation> gene_list;
	
	// Load from a GBK file
	bool loadGBK(gzFile m_fin, size_t &m_pos);
	void loadGBKFeatures(gzFile m_fin);

	// Load from a EMBL file
	bool loadEMBL(gzFile m_fin, size_t &m_pos);
	void loadEMBLFeatures(gzFile m_fin);
	
public:
	enum {SOURCE, TAXA_NAME, LINEAGE, 
		GENUS, SPECIES, SUBSPECIES, LAST_INFO}; // info_map
	
	// Data input types:
	// GBK = ASCII Genbank annontation flat file (similar to EMBL)
	// EMBL = ASCII EMBL annotation file (similar to GBK)
	// FASTA = FASTA file
	// FASTQ = FASTQ file
	enum {GBK, EMBL, FASTA, FASTQ, LAST_ANNOT_FORMAT};
	
	DNAMol()
	{
		seq = NULL;
		seq_len = 0;
		start = stop = 0;
		info_map = std::vector<std::string>(LAST_INFO);
	};
	
	DNAMol(const DNAMol &m_copy)
	{
		seq = NULL;
		seq_len = 0;
		start = stop = 0;
		info_map = std::vector<std::string>(LAST_INFO);
		
		// Avoid shallow copy errors -- need to write a copy constructor
		if(m_copy.seq != NULL){
			throw ":DNAMol: Copy constructor has not been provided";
		}
	};
	
	~DNAMol()
	{
	};
	
	inline DNAMol& operator=(const DNAMol &m_copy)
	{
		throw ":DNAMol: Copy constructor has not been provided";
	};
	
	// A helper function for processing gene annotation lists
	// (i.e. sorting, adding intergenic space, counting the numner of 
	// annotations, etc.).
	void processGeneList(const bool &m_loading = false);

	inline std::string seq_id_str() const
	{
		return accession;
	};
	
	inline bool operator==(const std::string &m_accession)
	{		
		return (accession == m_accession);
	};

	inline std::list<GeneAnnotation>::iterator begin()
	{
		return gene_list.begin();
	};
	
	inline std::list<GeneAnnotation>::iterator end()
	{
		return gene_list.end();
	};

	inline std::list<GeneAnnotation>::const_iterator begin() const
	{
		return gene_list.begin();
	};
	
	inline std::list<GeneAnnotation>::const_iterator end() const
	{
		return gene_list.end();
	};

	inline bool empty() const
	{
		return (seq == NULL);
	};
	
	inline SEQPTR sequence() const
	{
		return seq;
	};

	inline unsigned int num_bases() const
	{
		return seq_len;
	};
	
	// Load and append the annotation file
	bool load(gzFile m_fin, const unsigned int &m_type, size_t &m_pos);
	
	bool import(const DNAMol &m_mol, const unsigned int &m_type);
	
	// Flush all data
	void clear()
	{
		accession.clear();

		start = stop = 0;
		info_map.clear();

		if(seq != NULL){

			delete [] seq;
			seq = NULL;
		}
		
		gene_list.clear();
	};
	
	inline std::string info(unsigned int m_key) const
	{
		if(m_key >= LAST_INFO){
			throw "info: key out of bounds";
		}

		return info_map[m_key];
	};
};

// Win32 STL implementation requires std::greater for std::list.sort
// I suspect that other implementations use std::less

namespace std{
	template <> 
	struct greater< list<GeneAnnotation>::iterator >
	{
		bool operator()(const list<GeneAnnotation>::iterator &m_a, 
			const list<GeneAnnotation>::iterator &m_b)
		{
			return greater<GeneAnnotation> () (*m_a, *m_b); 
		};
	};
}


// In annotation.cpp
int file_type(const std::string &m_filename);
std::string wrap_string(const std::string &m_str, const unsigned int &m_len);
bool integer(const std::string &m_int);
std::string trim_space(const std::string &m_str);
std::string pack_string(const std::string &m_str);

// In annotation_util.cpp
bool read_range(gzFile m_fin, std::pair<unsigned int, unsigned int> &m_range,
	std::list< std::pair<unsigned int, unsigned int> > &m_seg_list);
unsigned int list_to_int(std::list<char> &m_number);
char *error_msg(const char *m_error);
bool strip_eol(char *m_buffer, const size_t &m_max_len);

// In annotation_export.cpp
std::string itoa(const int &m_int);

#endif // __GENOME_ANNOTATION
