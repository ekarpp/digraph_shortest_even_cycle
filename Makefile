CXX := g++
CXXFLAGS := -g -std=c++1z -O3 -Wall -Wextra -march=native
vpath %.cc src
vpath %.hh src

all: digraph

digraph: main.o graph.o util.o gf.o
	$(CXX) $^ -o $@

clean:
	rm digraph *.o

include Makefile.dep
