CXX := g++
CXXFLAGS := -g -std=c++1z -O3 -Wall
vpath %.cc src
vpath %.hh src

all: digraph

digraph: main.o graph.o util.o
	$(CXX) $^ -o $@

clean:
	rm digraph *.o

include Makefile.dep
