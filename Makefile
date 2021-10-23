OBJS = tntblast.o degenerate_na.o primer.o input.o nuc_cruc_anchor.o \
	nuc_cruc.o nuc_cruc_output.o nuc_cruc_santa_lucia.o bind_oligo.o \
	amplicon_search.o probe_search.o padlock_search.o tntblast_master.o \
	tntblast_worker.o tntblast_util.o options.o tntblast_local.o bitmask.o \
	sequence_data.o sequence_data_annot.o annotation.o annotation_embl.o \
	annotation_gbk.o annotation_util.o \
	sequence_data_fastx.o tntblast_timer.o mpi_util.o

CC = mpic++

PROFILE = #-g -pg

# If you're using gcc and would like to enable single computer multi-threading,
# uncomment this OPENMP variable
OPENMP = -fopenmp

# If you're using the clang C++ compiler on OS X and would like to enable single 
# computer multi-threading, please uncomment this CLANG_OPENMP variable.
#	- Before attemping to compile tntblast, please read the very helpfull
# 		blog post at https://iscinumpy.gitlab.io/post/omp-on-high-sierra/
# 	- Follow the process documented in the above blog post to install the required
#		OpenMP libraries that are needed by the Clang compiler.
#	- When running tntblast on OS X, please note that DYLD_LIBRARY_PATH must be set 
#		to the directory that contains the OpenMP libomp.dylib file, i.e., :
#	
#		export DYLD_LIBRARY_PATH=/$HOME/llvm-project/build-openmp/runtime/src
# With the introduction of ncbi-blast-2.11, we must specify -std=c++14 (as opposed to the
# -std=c++11 that was used for ncbi-blast-2.10

#CLANG_OPENMP = -Xpreprocessor -fopenmp

# Define USE_MPI to enable MPI
FLAGS = $(PROFILE) -O3 -Wall $(OPENMP) -std=c++14 -DUSE_MPI

INC = -I. 
LIBS = -lm -lz

ifdef CLANG_OPENMP
	OPENMP = $(CLANG_OPENMP)
	INC += -I$(HOME)/llvm-project/build-openmp/runtime/src
	LIBS += -L$(HOME)/llvm-project/build-openmp/runtime/src -lomp
endif

# The BLAST_DIR variable should only be defined if you wish to be able to
# read NCBI BLAST-formatted database files. This functionality is optional, and
# can be disabled by commenting out the BLAST_DIR variable by adding the '#' symbol
# to the start of the line below (i.e. "#BLAST_DIR").
BLAST_DIR = $(HOME)/ncbi-blast-2.12.0+-src

ifdef BLAST_DIR

	FLAGS += -DUSE_BLAST_DB

	# Compile with the NCBI C++ toolkit (tntblast will be able to read BLAST-formatted databases)
	INC += -I$(BLAST_DIR)/include/ncbi-tools++
	LIBS += -L $(BLAST_DIR)/lib \
		-lseqdb -lxobjmgr -lblastdb -lgeneral -lgenome_collection -llmdb \
		-lseq -lpub -lmedline -lseqcode -lseqset -lsequtil -lxser -lxutil -lxncbi -lsubmit -lbiblio -ldl

endif

.SUFFIXES : .o .cpp .c
.cpp.o:	
	$(CC) $(FLAGS) $(INC) -c $<

all: tntblast

tntblast : $(OBJS)
	$(CC) $(PROFILE) -o tntblast $(OBJS) $(LIBS) $(OPENMP)

clean:
	-rm -f *.o

