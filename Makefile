OBJS = tntblast.o degenerate_na.o primer.o input.o nuc_cruc_anchor.o \
	nuc_cruc.o nuc_cruc_output.o nuc_cruc_santa_lucia.o bind_oligo.o \
	amplicon_search.o probe_search.o padlock_search.o tntblast_master.o \
	tntblast_worker.o tntblast_util.o options.o tntblast_local.o bitmask.o \
	sequence_data.o sequence_data_annot.o annotation.o annotation_embl.o \
	annotation_gbk.o annotation_ptt.o annotation_util.o annotation_gff3.o gff3.o \
	sequence_data_fastx.o tntblast_timer.o

CC = mpic++

PROFILE = #-pg
OPENMP = #-fopenmp

# Define USE_MPI to enable MPI
FLAGS = $(PROFILE) -O3 -Wall $(OPENMP) -std=c++11 -DUSE_MPI

# The BLAST_DIR variable should only be defined if you wish to be able to
# read NCBI BLAST-formatted database files. This functionality is optional, and
# can be disabled by commenting out the BLAST_DIR variable by adding the '#' symbol
# to the start of the line below (i.e. "#BLAST_DIR"). To *enable* the functionality,
# the set the BLAST_DIR variable to the install path of the NCBI BLAST+ program.
BLAST_DIR = ncbi-blast+-dir

ifdef BLAST_DIR

FLAGS += -DUSE_BLAST_DB

# Compile with the NCBI C++ toolkit (tntblast will be able to read BLAST-formatted databases)
BLAST_INC = $(BLAST_DIR)/include/ncbi-tools++
BLAST_LIBS = -L $(BLAST_DIR)/lib \
	-ldl -lseqdb -lxobjmgr -lblastdb -lgeneral -lgenome_collection -llmdb \
	-lseq -lpub -lmedline -lseqcode -lseqset -lsequtil -lxser -lxutil -lxncbi -lsubmit -lbiblio

INC = -I. -I$(BLAST_INC)
LIBS = -lm $(BLAST_LIBS)
else

# Compile without the NCBI C++ toolkit (tntblast will not read BLAST-formatted databases)
BLAST_INC =
BLAST_LIBS =
INC = -I.
LIBS = -lm
endif

.SUFFIXES : .o .cpp .c
.cpp.o:
	$(CC) $(FLAGS) $(INC) -c $<

all: tntblast

tntblast : $(OBJS)
	$(CC) $(PROFILE) -o tntblast $(OBJS) $(LIBS) $(OPENMP)

clean:
	-rm -f *.o

