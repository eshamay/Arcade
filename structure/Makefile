include ../Makefile

STRUCTURE				= $(ANALYZER)
SYSTEMANALYSES	= $(STRUCTURE) neighbor-analysis.o structure-analyzer.o

structure-analyzer : $(SYSTEMANALYSES)
	$(CXX) $(SYSTEMANALYSES) $(LIBS) -lboost_thread -o ../bin/structure-analyzer

cleanstructure :
	rm -f *.o 
