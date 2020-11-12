include ../Makefile.rules
HEADERS = $(shell ls ../*.h)
LIB_NAMES = ../libfluid.a ../libmonte_carlo.a
LINK_ARGS = -L.. -lfluid -lmonte_carlo
main: main.cpp $(LIB_NAMES) $(HEADERS) Makefile ../Makefile.rules
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $< $(LINK_ARGS)
clean:
	rm -f main
	rm -f *.a *.o
