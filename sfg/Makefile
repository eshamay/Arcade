include ../Makefile

TEST	= $(MDSYSTEM) moritah2o.o test.o
SFG = $(ANALYZER) moritah2o.o cp2k-morita2002.o

test : $(TEST)
	$(CXX) $(TEST) -o test

sfg : $(SFG)
	$(CXX) $(SFG) $(LIBS) -lmkl_lapack -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lgsl -lgslcblas -o ../bin/structure-analyzer

cleansfg :
	rm -f *.o test
