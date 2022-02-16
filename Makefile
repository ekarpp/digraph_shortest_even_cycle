CXX := g++
CXXFLAGS := -g -std=c++1z -O3 -Wall -Wextra -march=native
VPATH = src:tests
BIN := digraph digraph-tests

all: $(BIN)

digraph: main.o graph.o util.o gf.o extension.o fmatrix.o ematrix.o polynomial.o
	$(CXX) $^ -o $@

digraph-tests: tests.o util.o gf.o gf_test.o extension.o extension_test.o fmatrix.o matrix_test.o ematrix.o fmatrix_test.o polynomial.o
	$(CXX) $^ -o $@

clean:
	rm -f $(BIN) *.o

test: digraph-tests
	./digraph-tests -egfm -d20 -n15 -t1000

include Makefile.dep
