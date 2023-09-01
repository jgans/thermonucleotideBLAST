![tntblast logo](https://github.com/jgans/thermonucleotideBLAST/blob/master/logo.png)

# Current version is 2.61 (September 1, 2023)
# Overview
ThermonucleotideBLAST is a software program for searching a target database of nucleic acid sequences using an assay-specific query. ThermonucleotideBLAST queries are based on biochemical assays (i.e. a pair of oligonucleotide sequences representing PCR primers or Padlock probes, a triplet of oligos representing PCR primers and a TaqMan probe or a single oligo representing a hybridization probe). Unlike existing programs (e.g. BLAST) which use heuristic measures of sequence similarity for identifying matches between a query and target sequence, ThermonucleotideBLAST uses physically relevant measures of sequence similarity -- free energy and melting temperature. For example, given a pair of PCR primers, a database of DNA targets and an annealing temperature, ThermonucleotideBLAST will return a list of predicted amplicons that will (ideally) match experimental PCR results. To enable searching of very large sequence databases (i.e. all of Genbank), ThermonucleotideBLAST can use run-time database and query segmentation to distribute the computational load across multiple CPUs.

The reference for ThermonucleotideBLAST is: "[Improved assay-dependent searching of nucleic acid sequence databases.](https://www.ncbi.nlm.nih.gov/pubmed/18515842?ordinalpos=1&itool=EntrezSystem2.PEntrez.Pubmed.Pubmed_ResultsPanel.Pubmed_DefaultReportPanel.Pubmed_RVDocSum)" by J. D. Gans and M. Wolinsky, Nucleic Acids Res. 2008 Jul;36(12):e74. doi: 10.1093/nar/gkn301.

Please send bug-reports/suggestions/comments/rants to jgans@lanl.gov

## Features
- Search by
  - PCR primer pair
  - PCR primer pair + hybridization probe (i.e. TaqMan PCR)
  - Padlock primer pair
  - Hybridization probe

- Match criteria specified by
  - Melting temperature (Tm), computed by the nearest neighbor method
  - Free energy change upon binding, computed by the nearest neighbor method
  - Match criteria can be supplemented with experiment specific heuristics (e.g. maximum amplicon length,  3' exact matches , etc.) 
- Supported file formats
  - FASTA and gzip-compressed FASTA
  - GBK (GenBank files containing sequence and annotation information) and gzip-compressed GBK
  - EMBL (EMBL files containing sequence and annotation information) and gzip-compressed EMBL
     - Please note that reading compressed FASTA/GBK/EMBL files is quite slow (due to the IO pattern in tntblast, not due to the cost of decompressing). If you need high performance searching, please search sequences stored in a BLAST-formatted database.
  - BLAST-formatted databases (version 5 -- requires the BLAST+ source code)
- Adaptive query and database segmentation allows efficient parallel execution with low memory overhead (ideal for diskless clusters)
- Annotation of genomic regions matches by query assays (requires a target database that contains annotation information)
- Different PCR primer strand concentrations can be specified for the forward and reverse primer to model asymmetric PCR conditions
- Non-ATGC bases are allowed
  - Inosine is a valid base
  - Assay oligos containing degenerate bases are handled by enumerating all possible nucleotide combinations. The melting temperature calculation accounts for the reduced concentration of each oligo.
  - Matches to database sequences with degenerate nucleotides are now supported. However, matches will only be found if there
  is at least one, short perfect match "seed" of real (i.e. A, T, G, C, I) target bases to initiate a pairwise alignment.

- Supported platforms are Linux and OS X
  - Single CPU computers
  - Multi-core and multi-CPU shared memory machines (using OpenMP)
  - Clusters (using MPI)
  
# Getting started

ThermonucleotideBLAST is written in C++ and has been tested on both shared memory and cluster based computers and compiles under Linux and OS X. 

ThermonucleotideBLAST reads target databases in a number of common formats, including FASTA, BLAST-databases and Genbank flat files (GBK). To handle the task of distributing database searches across multiple CPUs, ThermonucleotideBLAST can use either MPI or OpenMP.

Assay query files are tab-delimited text files with two or more columns:
- PCR assays:
  - Column 1: Assay name (no spaces allowed)
  - Column 2: forward primer sequence
  - Column 3: reverse primer sequence
- TaqManPCR assays:
  - Column 1: Assay name (no spaces allowed)
  - Column 2: forward primer sequence
  - Column 3: reverse primer sequence
  - Column 4: probe sequence
- Simple hybridization probe:
  - Column 1: Assay name (no spaces allowed)
  - Column 2: probe sequence
- Padlock probe
  - Column 1: Assay name (no spaces allowed)
  - Column 2: 5'-probe sequence
  - Column 3: 3'-probe sequence
  
## Optional prerequisites are:

- MPI (for parallel execution on a cluster computer)
  - Development and testing uses [OpenMPI](https://www.open-mpi.org/)
- A compiler that supports OpenMP (for multithreaded execution on a single, multi-core computer)
  - gcc/g++ versions 4.2 and higher
  - clang (however, OS X users will not have OpenMP support "out-of-the-box". I recommend following [these](https://iscinumpy.gitlab.io/post/omp-on-high-sierra/) instructions to get OpenMP working on OS X).
- The compiled NCBI BLAST+ source code (for searching BLAST-formatted databases)
  - Download, build and install the [NCBI BLAST+ program](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
  - The current version of the NCBI BLAST+ program requires a modern C++ compiler to build
    - The current Clang compiler on OS X works fine. Older (< 5) versions of GCC don't work.
    - For older Centos Linux systems with GCC version < 5, the [following steps](https://linuxize.com/post/how-to-install-gcc-compiler-on-centos-7/) will install a newer GCC compiler:
      - `sudo yum install centos-release-scl`
      - `sudo yum install devtoolset-7`
      - `scl enable devtoolset-7 bash`
  - For those who need to work with older x86 hardware, you may need to compile the NCBI BLAST+ code using `--without-sse42` to disable the use of potentially unsupported SSE instructions.
  - After running `make` and `make install` to build and install the NCBI BLAST+ program, note the BLAST+ installation location. This location (which includes the `include/` and `lib/` sub-directories created by running the BLAST+ `make install`) must be used to set the `BLAST_DIR` variable that appears in the thermonucleotideBLAST `Makefile`.
  
ThermonucleotideBLAST uses the DNA hybridization parameters published by the SantaLucia lab to compute duplex stability. However, the parameters for computing the delta G contribution of terminal mismatches (for internal loops) have not been published. Since we have not obtained permission to include these parameters in ThermonucleotideBLAST, they have not been included. The parameters are available as part of the excellent UNAFold suite of programs. To use these terminal mismatch hybridization parameters in ThermonucleotideBLAST, first download the UNAFold package. Then, run the parse_tstacki.pl script with the path to UNAFold as the only argument and pipe the script output to the file "nuc_cruc_santa_lucia_tstacki.cpp".  For example: "parse_tstacki.pl /path/to/UNAFold > nuc_cruc_santa_lucia_tstacki.cpp". After this file has been created, ThermonucleotideBLAST can be compiled as described below. Note that including these parameters is optional: you do not need to include these terminal mismatch parameters to compile and run ThermonucleotideBLAST (although calculation accuracy will be improved by including them).

Building ThermonucleotideBLAST without any of these prerequisites will generate a single CPU version that will not read BLAST-formated databases. 

# How to build and install ThermonucleotideBLAST on Unix-like systems
For now, you will need to directly edit the very simple Makefile.
- The `USE_MPI` macro (defined by adding `-DUSE_MPI` to the `FLAGS` variable) enables MPI-based, multi-computer parallel searching.
- The `USE_BLAST_DB` macro (defined by adding `-DUSE_BLAST_DB` to the `FLAGS` variable) enables the searchig of sequences stored in BLAST-formatted databases.
  - This option also requires a compiled and installed NCBI C++ toolkit. 
  - Indicate the path of the local NCBI BLAST+ installation with the `BLAST_DIR` variable.
- When compiling *without* the NCBI BLAST+ package, removed the `-DUSE_BLAST_DB` flag and comment out the `BLAST_DIR` variable (i.e. `#BLAST_DIR`).
- Adding `-fopenmp` to the `OPENMP` variable enables OpenMP-based, multiprocessor parallel searching.

Currently, thermonucleotideBLAST will only use OpenMP (multiple cores on a single computer) or MPI (single-threaded task on multiple computers), but not both.

# How to run ThermonucleotideBLAST

On a single CPU machine:

`./tntblast <args>`

On multiple CPUs using OpenMP:

`./tntblast <args>`

Note that the number of parallel workers used by OpenMP is determined by the environment variable `OMP_NUM_THREADS`. Try setting this variable to the number of CPU's times the number of cores per CPU on your machine.

On multiple CPU's using MPI:

The command line required to invoke ThermonucleotideBLAST varies depending on the version and configuration of MPI. For most systems, the command line will be similar to:

`mpirun -np N [--hostfile my_hosts] ./tntblast <args>`

where N is the number of CPUs to run on, my_hosts is an optional list of machine names and <args> are a list of command arguments that are passed to ThermonucleotideBLAST. ThermonucleotideBLAST employs a master/worker architecture that requires N >= 2.

## BLAST database searching

If compiled with the NCBI BLAST+ libraries, ThermonucleotideBLAST can read BLAST-formatted nucleotide sequence databases. By default, the user-specified queries are searched against *all* sequences in the provided BLAST database. However, thanks to the NCBI BLAST+ implementation, it is now possible to [limit database searches by sequence accession or species-level NCBI taxonomic id](https://www.ncbi.nlm.nih.gov/books/NBK569846/). 
- The `--blast-include` command-line flag accepts a single accession or NCBI TaxID value. All corresponding sequence entires in the database will be searched. This command-line argument can be reapeated multiple times.
- The `--blast-exclude` command-line flag accepts a single accession or NCBI TaxID value. All corresponding sequence entires in the database will *not* be searched. This command-line argument can be reapeated multiple times.

By default, the set of BLAST database sequences to select for searching is entire database.
If one or more sequences are specified with the `--blast-include` command, then the search set is limited to these sequences. Then, if one more sequences are specified with the `--blast-exclude` command, then these excluded sequences are removed from the search set.

Specifing an accession or taxonomic id that does not exist in the database will generate an error. Please note that taxonomic id values are currently required to be at, or below, species-level. Providing a higher level taxonomic id (e.g., genus-level) will generate an error.

## Example 1: Testing viral hemorrhagic fever virus assays from the literature

In this example, we will computationally test the specificity of PCR and TaqMan PCR based hemorrhagic fever virus (HFV) detection assays obtained from the literature.

First, we will need a database of sequences to screen. Downloading all [filoviruses](https://www.ncbi.nlm.nih.gov/nuccore/?term=txid11266[Organism:exp]) from NCBI will allow us to test the specificity of assays that target a number of HFVs (including Ebola, Marburg, Junin and Lassa viruses) as well as their near neighbors. When downloading these sequences, choose the "FASTA" file format (under "Send to") and "File" (under "Destination"). When the data has finished downloading, rename the resulting file "filo.fna".

Second, using a text editor, copy the following queries into a text file called "HFV_query.txt": 
```
gibb-marburg    TTCCCCTTTGGAGGCATC    GGAGGATCCAACAGCAAGG    CGATGGGCTTTCAGGACAGGTGT
towner-ebola    GAAAGAGCGGCTGGCCAAA    AACGATCTCCAACCTTGATCTTT    TGACCGAAGCCATCACGACTGCAT
zhai-filo    TATTCTCYCTACAAAAGCATTGGG    GCTTCTGCGAGTGTTTGGACATT   
```

Assay "gibb-marburg" is from Gibb et. al., Molecular and Cellular Probes. 15, 259-266, 2001

Assay "towner-ebola" is from Towner et. al., J. Virology 78, 4330-4341, 2004

Assay "zhai-filo" is from Zhai et. al., Clin. Microbiol. 45, 224-226, 2007

Please note that:
- PCR and TaqMan PCR assays can be mixed in the same query file
- Query oligos can have degenerate bases (as in zhai-filo).

Third, in the same directory that contains the query and target sequence files, run ThermonucleotideBLAST using the following command line:

`./tntblast -i HFV_query.txt -d filo.fna -e 40 -E 45 -o HFV.out`

These command arguments instruct ThermonucleotideBLAST to search the target sequences (contained in `filo.fna`) using the assays (contained in query file HFV_query.txt) and return all matches where the primer oligos bind at 40 °C (or higher; -e 40) and the probe sequences bind at 45 °C (or higher; -E 45). Since we haven't specified as assay format, the queries are assumed to contain PCR primers with an optional TaqMan probe. Unless specified with the -l flag, the maximum amplicon length that will be tested is 2000 bases.

Fourth, after running the command, the output will be in the file HFV.out. A truncated version of the output (showing only a single match) is displayed below:
```
#####################################################################################
name = gibb-marburg
forward primer = 5'-TTCCCCTTTGGAGGCATC-3'
reverse primer = 5'-GGAGGATCCAACAGCAAGG-3'
forward primer tm = 54.9758
reverse primer tm = 56.0109
forward primer hairpin tm = 21.2137
reverse primer hairpin tm = 45.2609
forward primer homodimer tm = 0
reverse primer homodimer tm = 0
heterodimer tm = 0
forward primer dG[-16.8574] = dH[-135.5] - T*dS[-0.382533]
reverse primer dG[-17.8955] = dH[-146.5] - T*dS[-0.414653]
min 3' clamp = 18
forward primer %GC = 55.5556
reverse primer %GC = 57.8947
forward primer heuristics = MULTI_5_GC, NO_POLY_RUNS
reverse primer heuristics = MULTI_5_GC, NO_POLY_RUNS
amplicon range = 6121 .. 6267
amplicon length = 147
probe = 5'-CGATGGGCTTTCAGGACAGGTGT-3'
probe tm = 62.2047
probe hairpin tm = 26.5141
probe homodimer tm = 0
probe dG[-22.9778] = dH[-180.2] - T*dS[-0.506923]
probe %GC = 56.5217
probe range = 6143 .. 6165
probe contained in forward strand (+)
forward align 5' TTCCCCTTTGGAGGCATC 3'
forward align    ||||||||||||||||||
forward align 3' AAGGGGAAACCTCCGTAG 5'
forward align dimer alignment size = 18
reverse align 5' GGAGGATCCAACAGCAAGG 3'
reverse align    |||||||||||||||||||
reverse align 3' CCTCCTAGGTTGTCGTTCC 5'
reverse align dimer alignment size = 19
probe align 5' CGATGGGCTTTCAGGACAGGTGT 3'
probe align    |||||||||||||||||||||||
probe align 3' GCTACCCGAAAGTCCTGTCCACA 5'
probe align dimer alignment size = 23
>gi|13489275|ref|NC_001608.2| Lake Victoria marburgvirus, complete genome
TTCCCCTTTGGAGGCATCCAAGCGATGGGCTTTCAGGACAGGTGTACCTCCCAAGAATGTTGAGTATACAGAAGGGGAGGAAGCCAAAACATGCTACAATATAAGTGTAACGGATCCCTCTGGAAAATCCTTGCTGTTGGATCCTCC
```

Output description:

- All temperatures are in °C and energy is in Kcal/Mol.
- Heterodimer tm refers to interaction between the forward and reverse primer.
- Min 3' clamp is the number of 3' oligo bases that exactly complement the target sequence.
- %GC is the percentage of bases in the oligo that are G or C (as opposed to A or T).
- The forward and reverse primer heuristics show how each primer satisfies a list of empirical primer design rules;
  - POLY_3_GC: Avoid runs of 3 or more G's or C's at the 3' prime end.
  - MULTI_5_GC: The 5 bases at the 5' end should contain no more than 3 G's or C's if no two pyrimidines (T, C) are adjacent, 2 G's or C's otherwise.
  - NO_POLY_RUNS: Polyprimidine (T, C) and polypurine (A, G) runs should be avoided.
  - NO_3_T: Avoid 3' terminal T's.
- The alignment range is the zero based location of amplicon in the target sequence.
- When a '$' character appears in the sequence alignments, it is a "virtual" base that indicates a dangling base contribution to the free energy calculation (enabled by the `--dangle5` or `--dangle3` command line arguments that were not used in this example).
- The sequence shown at the bottom of the output is the predicted full length amplicon (including primer binding sites).

## Example 2: Extracting the annotations of targeted genes

In this example we will extract the annotations of the genes targeted by PCR assays. Using the hemorrhagic fever virus assays from example 1 (shown above), we will use ThermonucleotideBLAST to extract all of the matching annotations from the genome of the Lake Victoria marburgvirus.

First, since FASTA files only contain very limited annotation information, we first need to download annotations for the [Lake Victoria marburgvirus from NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NC_001608). Save the file to disk in "GenBank (full)" display format and call it "NC_001608.gbk". This file contains both annotation and the complete genome sequence for the Lake Victoria marburgvirus.

Second, in the same directory that contains the query (the same file from example 1; HFV_query.txt) and target sequence (NC_001608.gbk) files, run ThermonucleotideBLAST using the following command line:

`./tntblast -i HFV_query.txt -d NC_001608.gbk -e 40 -E 45 -o HFV.annot.out`

Third, after running the command, the output will be in the file HFV.annot.out. A truncated version of the output (showing only a single match) is displayed below:

```
name = gibb-marburg
forward primer = 5'-TTCCCCTTTGGAGGCATC-3'
reverse primer = 5'-GGAGGATCCAACAGCAAGG-3'
forward primer tm = 54.9758
reverse primer tm = 56.0109
forward primer hairpin tm = 21.2137
reverse primer hairpin tm = 45.2609
forward primer homodimer tm = 0
reverse primer homodimer tm = 0
heterodimer tm = 0
forward primer dG[-16.8574] = dH[-135.5] - T*dS[-0.382533]
reverse primer dG[-17.8955] = dH[-146.5] - T*dS[-0.414653]
min 3' clamp = 18
forward primer %GC = 55.5556
reverse primer %GC = 57.8947
forward primer heuristics = MULTI_5_GC, NO_POLY_RUNS
reverse primer heuristics = MULTI_5_GC, NO_POLY_RUNS
amplicon range = 6121 .. 6267
amplicon length = 147
probe = 5'-CGATGGGCTTTCAGGACAGGTGT-3'
probe tm = 62.2047
probe hairpin tm = 26.5141
probe homodimer tm = 0
probe dG[-22.9778] = dH[-180.2] - T*dS[-0.506923]
probe %GC = 56.5217
probe range = 6143 .. 6165
probe contained in forward strand (+)
forward align 5' TTCCCCTTTGGAGGCATC 3'
forward align    ||||||||||||||||||
forward align 3' AAGGGGAAACCTCCGTAG 5'
forward align dimer alignment size = 18
reverse align 5' GGAGGATCCAACAGCAAGG 3'
reverse align    |||||||||||||||||||
reverse align 3' CCTCCTAGGTTGTCGTTCC 5'
reverse align dimer alignment size = 19
probe align 5' CGATGGGCTTTCAGGACAGGTGT 3'
probe align    |||||||||||||||||||||||
probe align 3' GCTACCCGAAAGTCCTGTCCACA 5'
probe align dimer alignment size = 23
annotation(gene)[5939..7984]+: gi|13489279|gb|NP_042028.1| GP GP protein
>gi|13489275|ref|NC_001608.2| Lake Victoria marburgvirus, complete genome
TTCCCCTTTGGAGGCATCCAAGCGATGGGCTTTCAGGACAGGTGTACCTCCCAAGAATGTTGAGTATACAGAAGGGGAGGAAGCCAAAACATGCTACAATATAAGTGTAACGGATCCCTCTGGAAAATCCTTGCTGTTGGATCCTCC
```

Note that this output is almost identical to the output shown above for example 1. The only change is the addition of annotation information for all of the genome features that overlap the amplicon produced by the primer pair gibb-marburg. In this case, the amplicon is contained within the GP protein. The location of the GP protein in the genome is given by the zero based gene coordinates ([5939..7984]), the gene is encoded on the positive strand (as indicated by the "+", genes on the negative strand will be indicated by a "-"). The gene identifier (gi|13489279) and accession (gb|NP_042028.1) are also shown.

