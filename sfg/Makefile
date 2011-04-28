include ../Makefile

TEST	= $(MDSYSTEM) moritah2o.o test.o 
SFG = $(ANALYZER) moritah2o.o cp2k-morita2002.o

test : $(TEST)
	$(CXX) $(TEST) -o test

sfg : $(SFG)
	$(CXX) $(SFG) $(LIBS) -lgsl -lgslcblas -o ../bin/morita-sfg

cleansfg :
	rm -f *.o test
