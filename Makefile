CXX := g++
CXXFLAGS := -g -std=c++1z -O3 -Wall -Wextra -march=native
VPATH = src:tests/unit:tests/perf
BIN := digraph digraph-tests digraph-scale
OBJ := graph.o util.o gf.o extension.o fmatrix.o ematrix.o polynomial.o solver.o

all: $(BIN)

digraph: main.o $(OBJ)
	$(CXX) $^ -o $@

digraph-tests: tests.o $(OBJ) gf_test.o extension_test.o matrix_test.o fmatrix_test.o util_test.o solver_test.o ematrix_test.o
	$(CXX) $^ -o $@

digraph-scale: scale.o $(OBJ)
	$(CXX) $^ -o $@

clean:
	rm -f $(BIN) *.o *.s *.asm1 *.asm2

test: digraph-tests
	./digraph-tests -egfmux -d20 -n15 -t1000

test-solver: digraph-tests
	./digraph-tests -s -n10 -t100

%.s: %.cc
	$(CXX) -S $(CXXFLAGS) -fverbose-asm $^

%.asm1: %.s
	c++filt < $^ > $@

%.asm2: %.o
	objdump -d -S $^ > $@

include Makefile.dep
