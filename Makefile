CXX := g++
CXXFLAGS := -g -std=c++1z -O3 -Wall -Wextra -march=native
VPATH = src:tests/unit:tests/perf
OBJ := graph.o util.o gf.o extension.o fmatrix.o ematrix.o polynomial.o solver.o
BIN := digraph digraph-tests digraph-scale gf-perf extension-perf

ifdef 32bits
	CXXFLAGS += -D GF2_bits=32
	BITS := 32
else
	CXXFLAGS += -D GF2_bits=16
	BITS := 16
endif

all: $(BIN)

digraph: main.o $(OBJ)
	$(CXX) $^ -o $@
	mv digraph digraph$(BITS)

digraph-tests: tests.o $(OBJ) gf_test.o extension_test.o fmatrix_test.o util_test.o solver_test.o ematrix_test.o geng_test.o
	$(CXX) $^ -o $@
	mv digraph-tests digraph-tests$(BITS)

digraph-scale: scale.o $(OBJ)
	$(CXX) $^ -o $@
	mv digraph-scale digraph-scale$(BITS)

gf-perf: gf_perf.o $(OBJ)
	$(CXX) $^ -o $@
	mv gf-perf gf-perf$(BITS)

extension-perf: extension_perf.o $(OBJ)
	$(CXX) $^ -o $@
	mv extension-perf extension-perf$(BITS)

clean:
	rm -f $(addsuffix 16, $(BIN)) $(addsuffix 32, $(BIN)) *.o *.s *.asm1 *.asm2

nauty:
	cd nauty && ./configure && make geng && make listg && make directg && cd ..

test: digraph-tests
	./digraph-tests$(BITS) -egfux -d20 -t1000
	./digraph-tests$(BITS) -s -d10 -t100

geng-test: digraph-tests
	nauty/geng -q $(vert) | nauty/directg -q | nauty/listg -aq | ./digraph-tests$(BITS) -c

%.s: %.cc
	$(CXX) -S $(CXXFLAGS) -fverbose-asm $^

%.asm1: %.s
	c++filt < $^ > $@

%.asm2: %.o
	objdump -d -S $^ > $@

include Makefile.dep
