CC=gcc
CFLAGS=-O3 -fopenmp

# choose one from below
#CFLAGS += -D__3_PNT_APPROX__
CFLAGS += -D__5_PNT_APPROX__
#CFLAGS += -D__7_PNT_APPROX__

#CFLAGS += -D__SELF_GRAV__
#CFLAGS += -D__HARMONIC__
CFLAGS += -D__GRAV__

SCH_1D_OBJ = sch_1d.o init_run.o setup_IC.o io.o evolve.o moment.o df.o potential.o
SCH_1D_DEP := $(SCH_1D_OBJ)

POTENTIAL_TEST_OBJ = potential_test.o init_run.o setup_IC.o moment.o potential.o
POTENTIAL_TEST_DEP := $(POTENTIAL_TEST_OBJ)

all: sch_1d

sch_1d : $(SCH_1D_DEP)
	$(CC) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lm -lfftw3

setup_IC_multi: multi_stream_IC.o init_run.o setup_IC.o io.o df.o
	$(CC) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lm

potential_test: $(POTENTIAL_TEST_DEP)
	$(CC) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lm -lfftw3

clean:
	-rm -rf *.o

distclean: clean
	-rm -rf *~
	-rm -rf sch_1d
