include ../Makefile

STRUCTURE				= $(ANALYZER)
SYSTEMANALYSES	= $(STRUCTURE) manipulators.o structure-analyzer.o

structure-analyzer : $(SYSTEMANALYSES)
	$(CXX) $(SYSTEMANALYSES) $(LIBS) -lpthread -lgsl -lgslcblas -o ../bin/structure-analyzer

cleanstructure :
	rm -f *.o 
