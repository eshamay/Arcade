include ../Makefile

ANALYSES = h2o-analysis.o so2-system-analysis.o

STRUCTURE				= $(ANALYZER)
SYSTEMANALYSES	= $(STRUCTURE) manipulators.o $(ANALYSES) structure-analyzer.o

structure-analyzer : $(SYSTEMANALYSES)
	$(CXX) $(SYSTEMANALYSES) $(LIBS) -lpthread -lgsl -lgslcblas -o ../bin/structure-analyzer

cleanstructure :
	rm -f *.o 
