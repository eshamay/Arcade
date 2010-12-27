include ../Makefile

STRUCTURE				= $(ANALYZER)
SYSTEMANALYSES	= $(STRUCTURE) structure-analyzer.o

structure-analyzer : $(SYSTEMANALYSES)
	$(CXX) $(SYSTEMANALYSES) $(LIBS) -lpthread -o ../bin/structure-analyzer

cleanstructure :
	rm -f *.o 
