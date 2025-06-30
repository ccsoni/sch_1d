CC=gcc
CFLAGS=-O3 -fopenmp

SCH_1D_OBJ = sch_1d.o init_run.o setup_IC.o io.o evolve.o moment.o
SCH_1D_DEP := $(SCH_1D_OBJ)

all: sch_1d

sch_1d : $(SCH_1D_DEP)
	$(CC) $(CFLAGS) -o $@ $^ -lgsl -lgslcblas -lm

clean:
	-rm -rf *.o

distclean: clean
	-rm -rf *~
	-rm -rf sch_1d
