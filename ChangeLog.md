# Version 2.61 (September 1, 2023)
- Fixed a multi-threading race condition when reading BLAST-formatted databases (thanks to Steven Higgins for pointing this out!).

# Version 2.6 (August 4, 2023)
- Updated command line options output to indicate optional arguments
- Fixed algorithm that was intended to remove spurious matches to target sequences with a large fraction of degerate bases.
  There was a bug that had the side effect of removing matches to primers with more than 50% unalignable bases.
  This test has been fixed and the `num_read_bases()` member function in `nuc_cruc.h` has been replaced by the `fraction_aligned_real_base_pairs()` member function. The test for fraction of aligned real bases in applied in the
  `bind_oligo_to_*()` functions contained in `bind_oligo.cpp`.
- Fixed bug in the display of unaligned bases in primer alignments (in `nuc_cruc_output.cpp`). 
  This bug only appears when the primer sequence dangled beyond the ends of the target sequence.
- Fixed a bug in masking primer binding sites (see `mask_primer_5()` and `mask_primer_3()`) that generated incorrect masking
  when primer sequences dangled off the end of target sequences.
- Fixed a bug in the display of amplicon sequences for primers that overhang the target sequence (see `amplicon_search.cpp`).
  Only effected amplicons that needed to be displayed as the reverse complement. 
- Fixed buffer overflow error in `file_type()` (see `annotation.cpp`).
- Fixed bug when searching fasta and fastq formatted databases using the MPI version of the code. This bug was caused by
  not broadcasting the total size of the fasta or fastq sequence file to the workers. 

# Version 2.5 (Dec 11, 2021)
- Fixed run-time error when only a single MPI rank is specified (tntblast.cpp).
- Fixed the temperature units in the program usage output for `--temperature`. Should read 'K' for Kelvin (options.cpp).

# Version 2.4 (October 25, 2021)
- Fixed parsing error for annotation files (GBK and EMBL).

# Version 2.3 (October 22, 2021)
- Cleaned up `seq_hash.h` to address compiler warnings about clearing/setting an object of non-trivial type.
- Removed WIN32 support for reading fasta sequences.
- Removed support for memory mapping fasta and fastq sequences (to make the code easier to maintain -- if you need speed, please convert the sequences to a BLAST+ database).
- Removed unused variable from `annotation_util.cpp`.
- Removed support for sequence fragmentation.
- Added support for reading gzip-compressed fasta, fastq, GBK and EMBL files. Please note that this is convience feature only! Reading these gzip-compressed files is slower than reading uncompressed sequence files.
- Improved the reading of sequence annotation files by reading the sequence data in a single pass (as opposed to the previous two-pass scheme, which did not work well for compressed files).
- Unified the previously chaotic encoding of ascii bases to binary bases. 

# Version 2.2 (May 27, 2021)
- Enabled accession and NCBI TaxId limits on BLAST database searching. As with the current version of NCBI BLAST+, TaxId values must be at the species-level (or below). Specifying a higher-level TaxId will generate an error.
- Changed the default TaqMan probe concentration (for melting temperature calculations) to be 2.5e-7 M (used to be the same as the primer concentration).
- Updated the user parameter implementation.
- Updated MPI transport layer.
- Fixed compilation errors when *not* using MPI (i.e., when `USE_MPI` is not defined in the Makefile).

# Version 2.192 (May 15, 2021)
- Updated the Makefile to specify '-std=c++14', which is now required when using the ncbi-blast-2.11 code for reading BLAST+
formatted databases.
- Fixed a crash caused by an unhandled exception in sequence_data::read_bio_seq_ncbi(). There appears to be at
least one sequence in the current (May 10, 2021) BLAST 'nt' database (OID #3877023) that causes an exception in
ncbi::GetBioseqNoData(). While the cause of this exception is currently unknown, the tntblast code has been modified to catch
ncbi exceptions and treat the corresponding sequence as having zero bytes.
	- This particular OID also causes an issue with the `blastdbcmd` utility that is provided by NCBI. In particular, `./blastdbcmd -db nt -entry all -outfmt "%o %a %t"| awk '{if( ($1 > 3877020) && ($1 < 3877025) ){print $0}}'` generates an error when `blastdbcmd` encounters OID 3877023:
```
3877021 NM_202183.2 Arabidopsis thaliana AP2 domain-containing transcription factor family protein (CRF12), mRNA
3877022 NM_202187.2 Arabidopsis thaliana Mediator complex, subunit Med10 (AT1G26665), mRNA
Error: [blastdbcmd] error while reading seqid
```
- Removed old NCBI C-toolkit related code.

# Version 2.191 (April 2, 2021)
- Fixed the bug in tntblast.cpp that caused the number of OpenMP threads to always be reported as 1 (thanks to Dylan Skola for pointing this out).
- Updated the Makefile to simplify building tntblast with OpenMP support on on OS X using the Clang c++
compiler.

# Version 2.19 (March 15, 2021)
- Free multiplex-associated cache memory in `plus_strand_melt_cache` and `minus_strand_melt_cache` by swapping with an empty `unordered_map`, rather than relying on `unordered_map::clear()` to free memory.

# Version 2.18 (March 1, 2021)
- Fixed multiplex assay enumeration code to correctly handle groups of assays that share oligos (i.e. the same
forward primer may appear in multiple assays, or a rever primer in one assay may be the forward primer in a
different assay).
- Removed redundant calls to set the salt concentration in the DNA melting engine. These calls triggered an expensive
recalculation of dynamic programming parameters.
- Made the NucCruc::strand() function call consistent between tntblast_local and tntblast_worker.
- Cache melting temperature calculation results for each target sequence. This significantly accelerates searching
multiplex assays (or singleplex assays searched in multiplex mode, i.e. `--plex T`) by avoiding redundant searching of the same oligo against the same template sequence more than once.

# Version 2.17 (October 26, 2020)
- Simplified Makefile to remove the need to modify *both* the `USE_BLAST_DB` and `BLAST_DIR`
variables. Now, the `-DUSE_BLAST_DB` option is automatically added to the `FLAG` Makefile variable
when the user has set the `BLAST_DIR` Makfile variable.

# Version 2.16 (September 28, 2020)
- Fixed bug in seq.h that broke the relationship between "Seq_strand_plus", "Seq_strand_minus" and "Seq_strand_both". This bug prevented the search of probe and padlock queries.

# Version 2.15 (July 17, 2020)
- Fixed sorting crtieria in bind_oligo.cpp -> sort_by_bound_match() to ensure that we preferentially display the more informative primer template alignments.

# Version 2.14 (June 24, 2020)
- Added the accession to results when search BLAST-formatted databases.
- Fixed the "5'-" -> "5' " and "-3'" -> " 3'" that was supposed to have taken place in version 2.1

# Version 2.13 (May 27, 2020)
- Cleaned up code to remove some historical complexity.
- Removed dependence on the old NCBI C toolkit for reading BLAST databases.
	- Removed support for reading ASN.1 annotation files.
- Added support for reading the new (version 5) BLAST database produced by the NCBI C++ BLAST application.
	- The current version of the NCBI C++ toolkit requires a modern C++ compiler to build.
	- The current Clang compiler on OS X works fine. Older (< 5) versions of GCC don't work.
	- For older Centos Linux systems with GCC version < 5, the following steps will install a newer GCC compiler (from https://linuxize.com/post/how-to-install-gcc-compiler-on-centos-7/):
		- sudo yum install centos-release-scl
		- sudo yum install devtoolset-7
		- scl enable devtoolset-7 bash
	- For those who need to work with older x86 hardware, you may need to compile the NCBI C++ toolkit using "--without-sse42" to disable the use of potentially unsupported SSE instructions.

# Version 2.12 (May 3, 2020)
- Modified the alignment style (again!) to indicate complementary, but unaligned, base pairs with a ':' symbol (thanks to Paul Li for suggesting this). The complementary (but unaligned) base pairs are usually found at the ends of oligos and occur
when the thermodynamic penatly of incorporating a mismatch is larger than the benefit of incorporating terminal matching bases.

# Version 2.11 (April 15, 2020)
- When seeding alignments between assay oligos and target sequences, only include the "real" bases (A, T, G, and C).
- Don't accept alignments between assay oligos and target sequences that include more than 50% degenerate bases.

# Version 2.1 (March 28, 2020)
- Modified alignment style in output to show unaligned bases in the assay oligo and target sequences. This is
	intended to assist in interpreting the significance of 3' mismatches (ala the TaqMAMA method of Li et al 2004).
- Returned to the simple Makefile format
- Changed the alignment output from "5' " -> "5'-" and " 3'" -> "-3'"

# Version 2.03
- Fixed bug introduced when databases were no longer searched by length.
- Updated annotation_gbk.cpp to enable parsing of GBK files that contain the "CONTIG" field.

# Version 2.02
- Stop sorting the database sequences by length. This stratgey was initially used to provide better load balancing
	(by searching longer sequences first), but now only slows things down! It is much faster to read sequences in the
	order in which they are stored in the database, rather than jumping around the database to read the larger sequences
	first.
- When tntblast is run without *any* command line arguments, it now prints a list of command line arguments and exists
	(rather than complaining about a missing input file).

# Version 2.01
- Fixed compiler errors by adding "#include <string.h>" and "#include <limits.h>" to a number of source files

# Version 2.0
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

# Version 1.06
- Added GFF3 sequence format
- Fixed bug in fasta parser

# Version 1.05
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

# Version 1.04
- First public release
- Added the NCBI toolkit to the Visual C++ project and exe (Windows)
