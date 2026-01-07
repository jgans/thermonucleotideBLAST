OBJS = tntblast.o degenerate_na.o primer.o input.o nuc_cruc_anchor.o \
	nuc_cruc.o nuc_cruc_output.o nuc_cruc_santa_lucia.o bind_oligo.o \
	amplicon_search.o probe_search.o padlock_search.o tntblast_master.o \
	tntblast_worker.o tntblast_util.o options.o tntblast_local.o bitmask.o \
	sequence_data.o sequence_data_annot.o annotation.o annotation_embl.o \
	annotation_gbk.o annotation_util.o \
	sequence_data_fastx.o tntblast_timer.o mpi_util.o

# -----------------------
# Toolchain defaults (override from CLI / CI)
# -----------------------
CC ?= g++
BIN ?= tntblast

PROFILE ?= # -g -pg
OPENMP ?=
FLAGS_EXTRA ?=

# Set to 1 to compile with MPI support (will switch compiler to mpic++ unless CC is explicitly set)
USE_MPI ?= 0

# Optional: set BLAST_DIR to enable BLAST DB support; default is empty (disabled)
BLAST_DIR ?=

# -----------------------
# Flags
# -----------------------
BASE_FLAGS = $(PROFILE) -O3 -Wall -std=c++14
FLAGS = $(BASE_FLAGS) $(OPENMP) $(FLAGS_EXTRA)

INC = -I.
LIBS = -lm -lz

# If MPI requested and user didn't override CC, use mpic++
ifeq ($(USE_MPI),1)
  ifeq ($(origin CC), default)
    CC = mpic++
  endif
  FLAGS += -DUSE_MPI
endif

# -----------------------
# Clang OpenMP (optional)
# -----------------------
ifdef CLANG_OPENMP
	OPENMP = $(CLANG_OPENMP)
	INC += -I$(HOME)/llvm-project/build-openmp/runtime/src
	LIBS += -L$(HOME)/llvm-project/build-openmp/runtime/src -lomp
endif

# -----------------------
# Optional BLAST DB support
# -----------------------
ifneq ($(strip $(BLAST_DIR)),)
	FLAGS += -DUSE_BLAST_DB
	INC += -I$(BLAST_DIR)/include/ncbi-tools++
	LIBS += -L $(BLAST_DIR)/lib \
		-lseqdb -lxobjmgr -lblastdb -lgeneral -lgenome_collection -llmdb \
		-lseq -lpub -lmedline -lseqcode -lseqset -lsequtil -lxser -lxutil -lxncbi -lsubmit -lbiblio -ldl
endif

.SUFFIXES : .o .cpp .c
.cpp.o:
	$(CC) $(FLAGS) $(INC) -c $<

all: $(BIN)

$(BIN): $(OBJS)
	$(CC) $(PROFILE) -o $(BIN) $(OBJS) $(LIBS) $(OPENMP)

clean:
	-rm -f *.o $(BIN)