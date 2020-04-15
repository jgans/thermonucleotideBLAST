- Version 1.04
	- First public release
	- Added the NCBI toolkit to the Visual C++ project and exe (Windows)
- Version 1.05
	- Fixed bug in fasta output (alignments were being written)
	- Changed the sorting rule for formatting output. Results are now
	  sorted by: (1) minimum primer Tm, (2) probe TM, (3) index of the target
	  sequence in the target database.
	- Added the ability to specify both 5' (using --probe-clamp5) and 3' 
	  (using --probe-clamp3) clamp values for ligation based assays. This is
	  important since ligation based assays (like Padlock probes and MOL-PCR)
	  offer more SNP discrimination at the 3' end of the upstream probe that is
	  adjacent to the downstream probe.
	- Improved the clarity of the output text for Padlock probe queries
	- Added "-fopenmp" to the linker flags in Makefile.am
	- Rewrote thermodynamic alignment algorithm
	- Added the Dinkelbach fractional programming algorithm
- Version 1.06
	- Added GFF3 sequence format
	- Fixed bug in fasta parser
- Version 2.0
	- Modified the output order to be sorted by: (1) min primer Tm, (2) probe Tm,
	  (3) max primer Tm [this is the new criteria] and (4) target sequence index.
	- Map input 'U' residues (i.e. RNA) to 'T' residues (i.e. DNA)
	- Added the number of mismatches and number of gaps (for each oligo) to the standard 
	  output format.
	- Added an additional (long overdue) test to PCR assays with probes to insure that
	  the probe does not overlap the primer that binds to the same strand (as this
	  would prevent probe hydrolysis for TaqMan PCR assays).
	- Fixed bug in non-memory-mapped fasta file parsing.
	- Added the ability to read FASTQ files.
	- Modified the parser to skip '*' symbols in sequence data (treat as a space).
	- Fixed a bug when multiplexing PCR assays with probes (which caused *no* assays to be found!)
	- Fixed a bug when selecting the best matching assay (--best-match). When the minimum primer Tm
	  is the same for two assay variants, we now sort on the maximum primer Tm to select the variant
	  to keep.
	- Fixed a bug in DNAMol::annot (needed to use the length sorted index, if present).
	- Fixed the GBK and EMBL parsing code to read degenerate bases (in addition to just A, T, G and C).
	- Changed the behaviour when reporting the number of mismatches. Unaligned bases are no longer counted
	  as mismatches -- the actual number of mismatches (assuming a gap free alignment) is now reported.
- Version 2.01
	- Fixed compiler errors by adding "#include <string.h>" and "#include <limits.h>" to a number of source files
- Version 2.02
	- Stop sorting the database sequences by length. This stratgey was initially used to provide better load balancing
	  (by searching longer sequences first), but now only slows things down! It is much faster to read sequences in the
	  order in which they are stored in the database, rather than jumping around the database to read the larger sequences
	  first.
	- When tntblast is run without *any* command line arguments, it now prints a list of command line arguments and exists
	  (rather than complaining about a missing input file).
- Version 2.03
	- Fixed bug introduced when databases were no longer searched by length.
	- Updated annotation_gbk.cpp to enable parsing of GBK files that contain the "CONTIG" field.
- Version 2.1 (March 28, 2020)
	- Modified alignment style in output to show unaligned bases in the assay oligo and target sequences. This is
	 intended to assist in interpreting the significance of 3' mismatches (ala the TaqMAMA method of Li et al 2004).
	- Returned to the simple Makefile format
	- Changed the alignment output from "5' " -> "5'-" and " 3'" -> "-3'"
	- Added support for database sequences with degenerate nucleotides. While non-degenerate bases are still needed 
	  to initiate an alignment, we now align query assay oligos with degnerate base target sequences by making
	  the optimistic assumption of complementary base-pairing (when possible).
- Version 2.11 (April 15, 2020)
	- When seeding alignments between assay oligos and target sequences, only include the "real" bases (A, T, G, and C).
	- Don't accept alignments between assay oligos and target sequences that include more than 50% degenerate bases.
