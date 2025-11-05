#CXX=icpx
CXX=g++
CXXFLAGS=-O3 -ffast-math

#CXXFLAGS += -D__PERIODIC__
CXXFLAGS += -D__GRAV__

OPENMP_FLAG = -fopenmp

SCH_1D_OBJ = main.o setup_IC.o moment.o evolve.o io.o df.o potential.o
SCH_1D_DEP := $(SCH_1D_OBJ)

POTENTIAL_TEST_OBJ = potential_test.o setup_IC.o df.o potential.o moment.o
POTENTIAL_TEST_DEP := $(POTENTIAL_TEST_OBJ)

sch_1d: $(SCH_1D_DEP)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm $(OPENMP_FLAG)

potential_test: $(POTENTIAL_TEST_DEP)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lfftw3 -lm $(OPENMP_FLAG)

sch_1d.o: sch_1d.h

clean:
	-rm -rf *.o

distclean: clean
	-rm -rf *~
	-rm -rf sch_1d
