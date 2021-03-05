- Version 2.18 (March 1, 2021)
	- Fixed multiplex assay enumeration code to correctly handle groups of assays that share oligos (i.e., the same
	forward primer may appear in multiple assays, or a rever primer in one assay may be the forward primer in a
	different assay).
	- Removed redundant calls to set the salt concentration in the DNA melting engine. These calls triggered an expensive
	recalculation of dynamic programming parameters.
	- Made the NucCruc::strand() function call consistent between tntblast_local and tntblast_worker.
	- Added code to cache melting temperature calculation results for each target sequence. This significantly accelerates searching
	assays in multiplex mode (i.e. `--plex T`) by avoiding redundant searching of the same oligo (when the oligo appears in multiple assays).
	against the same template sequence more than once.
- Version 2.17 (October 26, 2020)
	- Simplified Makefile to remove the need to modify *both* the `USE_BLAST_DB` and `BLAST_DIR`
	variables. Now, the `-DUSE_BLAST_DB` option is automatically added to the `FLAG` Makefile variable
	when the user has set the `BLAST_DIR` Makfile variable.
- Version 2.16 (September 28, 2020)
	- Fixed bug in seq.h that broke the relationship between "Seq_strand_plus", "Seq_strand_minus" and "Seq_strand_both". This bug prevented the search of probe and padlock queries.
- Version 2.15 (July 17, 2020)
	- Fixed sorting crtieria in bind_oligo.cpp -> sort_by_bound_match() to ensure that we preferentially display the more informative assay-template alignments.
- Version 2.14 (June 24, 2020)
	- Added the accession to results when search BLAST-formatted databases.
	- Fixed the "5'-" -> "5' " and "-3'" -> " 3'" that was supposed to have taken place in version 2.1
- Version 2.13 (May 27, 2020)
	- Cleaned up code to remove some historical complexity.
	- Removed dependence on the old NCBI C toolkit for reading BLAST databases.
		- Removed support for reading ASN.1 annotation files.
	- Added support for reading the new (version 5) BLAST database produced by the NCBI BLAST+ application.
		- The current version of [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) requires a modern C++ compiler to build.
		- The current Clang compiler on OS X works fine. Older (< 5) versions of GCC don't work.
		- For older Centos Linux systems with GCC version < 5, the [following steps](https://linuxize.com/post/how-to-install-gcc-compiler-on-centos-7/) will install a newer GCC compiler:
			- `sudo yum install centos-release-scl`
			- `sudo yum install devtoolset-7`
			- `scl enable devtoolset-7 bash`
		- For those who need to work with older x86 hardware, you may need to compile the NCBI BLAST+ program using `--without-sse42` to disable the use of potentially unsupported SSE instructions.
- Version 2.12 (May 3, 2020)
	- Modified the alignment style (again!) to indicate complementary, but unaligned, base pairs with a ':' symbol (thanks to Paul Li for suggesting this). The complementary (but unaligned) base pairs are usually found at the ends of oligos and occur
	when the thermodynamic penatly of incorporating a mismatch is larger than the benefit of incorporating terminal matching bases.
- Version 2.11 (April 15, 2020)
	- When seeding alignments between assay oligos and target sequences, only include the "real" bases (A, T, G, and C).
	- Don't accept alignments between assay oligos and target sequences that include more than 50% degenerate bases.
- Version 2.1 (March 28, 2020)
	- Modified alignment style in output to show unaligned bases in the assay oligo and target sequences. This is
	 intended to assist in interpreting the significance of 3' mismatches (ala the TaqMAMA method of Li et al 2004).
	- Returned to the simple Makefile format
	- Changed the alignment output from "5' " -> "5'-" and " 3'" -> "-3'"
- Version 2.03
	- Fixed bug introduced when databases were no longer searched by length.
	- Updated annotation_gbk.cpp to enable parsing of GBK files that contain the "CONTIG" field.
- Version 2.02
	- Stop sorting the database sequences by length. This stratgey was initially used to provide better load balancing
	  (by searching longer sequences first), but now only slows things down! It is much faster to read sequences in the
	  order in which they are stored in the database, rather than jumping around the database to read the larger sequences
	  first.
	- When tntblast is run without *any* command line arguments, it now prints a list of command line arguments and exists
	  (rather than complaining about a missing input file).
- Version 2.01
	- Fixed compiler errors by adding "#include <string.h>" and "#include <limits.h>" to a number of source files
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
- Version 1.06
	- Added GFF3 sequence format
	- Fixed bug in fasta parser
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
- Version 1.04
	- First public release
	- Added the NCBI toolkit to the Visual C++ project and exe (Windows)
