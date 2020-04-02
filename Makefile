OBJS = tntblast.o degenerate_na.o primer.o input.o nuc_cruc_anchor.o \
	nuc_cruc.o nuc_cruc_output.o nuc_cruc_santa_lucia.o bind_oligo.o \
	amplicon_search.o probe_search.o padlock_search.o tntblast_master.o \
	tntblast_worker.o tntblast_util.o options.o tntblast_local.o bitmask.o \
	sequence_data.o sequence_data_annot.o annotation.o annotation_embl.o \
	annotation_gbk.o annotation_ptt.o annotation_util.o annotation_gff3.o gff3.o \
	sequence_data_fastx.o ncbi_util.o tntblast_timer.o

CC = mpic++

PROFILE = #-pg
OPENMP = #-Xpreprocessor -fopenmp

# Define USE_MPI to enable MPI
FLAGS = $(PROFILE) -O3 -Wall $(OPENMP) -std=c++11 -DUSE_MPI

INC = -I. 

LIBS = -lm 

.SUFFIXES : .o .cpp .c
.cpp.o:
	$(CC) $(FLAGS) $(INC) -c $<

all: tntblast

tntblast : $(OBJS)
	$(CC) $(PROFILE) -o tntblast $(OBJS) $(LIBS) $(OPENMP)

clean:
	-rm -f *.o


