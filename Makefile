CXX			= icpc  -I$(BOOST) -wd981,444,383,177,1418,1782,869,1572,47
OPTIMIZE 	= -finline-functions -finline -funroll-all-loops -O3 -DNDEBUG -m64 -fast -restrict
#DEBUG		= -O0 -g3 -ggdb -D_GLIBCXX_DEBUG -Wno-deprecated -debug #-wd981,1599,1572,383
DEBUG		= -O0 -g3 -ggdb -D_GLIBCXX_DEBUG -Wno-deprecated -DNDEBUG #-wd981,1599,1572,383

CPPFLAGS    = -Wall -ftemplate-depth-100 $(OPTIMIZE)
#CPPFLAGS    = -Wall -ftemplate-depth-100 $(DEBUG)

LIBS		= -lconfig++

MOLECULES = $(MDSRC)/h2o.o $(MDSRC)/oh.o $(MDSRC)/h.o $(MDSRC)/h3o.o $(MDSRC)/hno3.o $(MDSRC)/so2.o $(MDSRC)/ctc.o $(MDSRC)/alkane.o
MDSYSTEM = $(MDSRC)/utility.o $(MDSRC)/atom.o $(MDSRC)/molecule.o $(MOLECULES) $(MDSRC)/moleculefactory.o $(MDSRC)/mdsystem.o $(MDSRC)/bondgraph.o
AMBERSYSTEM = $(MDSRC)/crdfile.o $(MDSRC)/topfile.o $(MDSRC)/ambersystem.o
XYZSYSTEM = $(MDSRC)/xyzfile.o $(MDSRC)/wannier.o $(MDSRC)/xyzsystem.o $(MDSRC)/molgraph.o $(MDSRC)/molgraphfactory.o
ANALYZER = $(MDSYSTEM) $(AMBERSYSTEM) $(XYZSYSTEM) $(MDSRC)/dataoutput.o

%.o: %.cpp %.h
	$(CXX) $(CPPFLAGS) -c -o $@ $<
