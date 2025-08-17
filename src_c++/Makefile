#CXX=icpx
CXX=g++
CXXFLAGS=-O3 -ffast-math

CXXFLAGS += -D__PERIODIC__

SCH_1D_OBJ = main.o setup_IC.o moment.o evolve.o io.o df.o
SCH_1D_DEP = $(SCH_1D_OBJ)

sch_1d: $(SCH_1D_DEP)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm

sch_1d.o: sch_1d.h

clean:
	-rm -rf *.o

distclean: clean
	-rm -rf *~
	-rm -rf sch_1d
