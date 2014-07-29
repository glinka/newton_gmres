ES_SRCS=test_es.cc eigen_solvers.cc
ES_OBJECTS=$(ES_SRCS:.cc=.o)
SRCS=main.cc gmres.cc newton.cc test_fns.cc
OBJECTS=$(SRCS:.cc=.o)

CXX = g++

CXXFLAGS = -g -Wall -std=c++0x -Wno-sign-compare -O3 -fPIC

all: newtongmres test_es

%.o: %.c
	$(CXX) $(CXXFLAGS) -c $<

newtongmres: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

test_es: $(ES_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^ > ./.depend

clean:
	$(RM) *.o 

include .depend
