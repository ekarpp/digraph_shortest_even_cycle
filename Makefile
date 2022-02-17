CXX := g++
CXXFLAGS := -g -std=c++1z -O3 -Wall -Wextra -march=native
VPATH = src:tests/unit
BIN := digraph digraph-tests
OBJ := graph.o util.o gf.o extension.o fmatrix.o ematrix.o polynomial.o solver.o

all: $(BIN)

digraph: main.o $(OBJ)
	$(CXX) $^ -o $@

digraph-tests: tests.o $(OBJ) gf_test.o extension_test.o matrix_test.o fmatrix_test.o util_test.o
	$(CXX) $^ -o $@

clean:
	rm -f $(BIN) *.o

test: digraph-tests
	./digraph-tests -egfmu -d20 -n15 -t1000

include Makefile.dep
