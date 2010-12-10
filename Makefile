TEST = crdfile.o atom.o

test: $(TEST)
	g++ $(TEST) test.cpp

%.o: %.cpp %.h
	$(CXX) $(CPPFLAGS) -c -o $@ $<
