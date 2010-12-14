CXX			= icpc  -I$(EIGEN) -I$(BOOST) -wd981,444,383,177,1418,1782,869,1572
OPTIMIZE 	= -finline-functions -finline -funroll-all-loops -O3 -DNDEBUG -m64 -fast -restrict
DEBUG		= -O0 -g3 -ggdb -D_GLIBCXX_DEBUG -Wno-deprecated -DNDEBUG -debug #-wd981,1599,1572,383

CPPFLAGS    = -Wall -ftemplate-depth-100 $(OPTIMIZE)
#CPPFLAGS    = -Wall -ftemplate-depth-100 $(DEBUG)

LIBS		= -L$(HOME)/share/lib -lconfig++ #-L$(MPI_HOME)/lib -L$(ATLAS)/lib 

MOLECULES = h2o.o oh.o h.o h3o.o hno3.o so2.o ctc.o alkane.o decane.o
MDSYSTEM = utility.o atom.o molecule.o $(MOLECULES) moleculefactory.o mdsystem.o
AMBERSYSTEM = $(MDSYSTEM) crdfile.o topfile.o ambersystem.o
TEST = $(AMBERSYSTEM) test.o

test: $(TEST)
	$(CXX) $(TEST) $(LIBS)

%.o: %.cpp %.h
	$(CXX) $(CPPFLAGS) -c -o $@ $<
