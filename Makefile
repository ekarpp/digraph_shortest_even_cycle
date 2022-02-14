CXX := g++
CXXFLAGS := -g -std=c++1z -O3 -Wall -Wextra -march=native
VPATH = src:tests
BIN := digraph digraph-tests

all: $(BIN)

digraph: main.o graph.o util.o gf.o extension.o fmatrix.o
	$(CXX) $^ -o $@

digraph-tests: tests.o util.o gf.o gf_test.o extension.o extension_test.o fmatrix.o matrix_test.o
	$(CXX) $^ -o $@

clean:
	rm $(BIN) *.o

include Makefile.dep
