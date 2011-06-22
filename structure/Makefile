include ../Makefile

ANALYSES = histogram-analysis.o h2o-analysis.o so2-system-analysis.o rdf-analysis.o bond-analysis.o neighbor-analysis.o angle-analysis.o 

STRUCTURE				= $(ANALYZER)
SYSTEMANALYSES	= $(STRUCTURE) manipulators.o $(ANALYSES) structure-analyzer.o

structure-analyzer : $(SYSTEMANALYSES)
	$(CXX) $(SYSTEMANALYSES) $(LIBS) -lpthread -lgsl -lgslcblas -o ../bin/structure-analyzer

cleanstructure :
	rm -f *.o 
