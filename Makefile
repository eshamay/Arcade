MOLECULES = h2o.o oh.o h.o h3o.o hno3.o so2.o ctc.o alkane.o decane.o
MDSYSTEM = utility.o atom.o molecule.o moleculefactory.o mdsystem.o wannier.o
AMBERSYSTEM = $(MDSYSTEM) $(MOLECULES) crdfile.o topfile.o ambersystem.o
TEST = $(AMBERSYSTEM)

test: $(TEST)
	g++ -g $(TEST) test.cpp

%.o: %.cpp %.h
	$(CXX) -g $(CPPFLAGS) -c -o $@ $<
