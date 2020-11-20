include Makefile.rules
HEADERS = $(shell ls *.h)

EXAMPLE_SUBDIRS = $(shell ls -d example_* | ls -d solution_* | ls -d heatrod | ls -d convectiondiffusion | sed "s/example_template//")


.PHONY: example_subdirs $(EXAMPLE_SUBDIRS)
example_subdirs: $(EXAMPLE_SUBDIRS)


$(EXAMPLE_SUBDIRS): $(FLUID_LIB) $(MC_LIB)
	$(MAKE) -C $@

all: $(FLUID_LIB) $(MC_LIB) example_subdirs


clean: 
	for d in $(EXAMPLE_SUBDIRS); do make -C $$d clean ; done
	rm -f *.a *.o 

MC_OBJS = monte_carlo.o
FLUID_OBJS = 3diagsys.o bndcond.o disc.o grid.o lutable.o phivar.o poisson.o

%.o: %.cpp $(HEADERS) Makefile Makefile.rules
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

$(FLUID_LIB): $(FLUID_OBJS) $(HEADERS) Makefile Makefile.rules
	ar r $(FLUID_LIB) $(FLUID_OBJS)

$(MC_LIB): $(MC_OBJS) $(HEADERS) Makefile Makefile.rules
	ar r $(MC_LIB) $(MC_OBJS)

