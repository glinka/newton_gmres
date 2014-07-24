ES_SRCS=esmain.cc eigen_solvers.cc
ES_OBJECTS=$(ES_SRCS:.cc=.o)
SRC=main.cc gmres.cc newton.cc test_fns.cc
OBJECTS=$(SRCS:.cc=.o)

CXX = g++

CXXFLAGS = -g -Wall -std=c++0x #-O3

all: newtongmres eigsolver

%.o: %.c
	$(CXX) $(CXXFLAGS) -c $<

newtongmres: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

eigsolver: $(ES_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^ > ./.depend

clean:
	$(RM) *.o 

include .depend
