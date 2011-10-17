include ../Makefile

ANALYSES = histogram-analysis.o so2-system-analysis.o rdf-analysis.o bond-analysis.o cycle-analysis.o neighbor-analysis.o angle-analysis.o angle-bond-analysis.o malonic-analysis.o succinic-analysis.o h2o-analysis.o so2-analysis.o so2-angle-analysis.o

STRUCTURE				= $(ANALYZER)
SYSTEMANALYSES	= $(STRUCTURE) manipulators.o $(ANALYSES) structure-analyzer.o 

structure-analyzer : $(SYSTEMANALYSES)
	$(CXX) $(SYSTEMANALYSES) $(LIBS) -lpthread -lgsl -lgslcblas -o ../bin/structure-analyzer

cleanstructure :
	rm -f *.o 
