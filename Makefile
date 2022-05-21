# CLEAN THIS UP! compilation takes a looooooong time

CXX := g++
CXXFLAGS := -g -std=c++1z -O3 -Wall -Wextra -march=native
VPATH = src:tests/unit:tests/perf
BIN := digraph digraph-tests extension-perf gf-perf
OBJ := graph.o util.o gf.o extension.o fmatrix.o ematrix.o polynomial.o solver.o
TEST_OBJ := gf_test.o extension_test.o fmatrix_test.o util_test.o solver_test.o ematrix_test.o geng_test.o

ifeq ($(bits), 32)
	CXXFLAGS += -D GF2_bits=32
else ifeq ($(bits), 16)
	CXXFLAGS += -D GF2_bits=16
else
	CXXFLAGS += -D GF2_bits=0
endif

all: $(BIN) nauty/geng nauty/directg nauty/listg

###################
# SOLVER BINARIES #
###################

digraph32:
	$(MAKE) digraphX bits=32
	rm -f *.o

digraph16:
	$(MAKE) digraphX bits=16
	rm -f *.o

digraph0:
	$(MAKE) digraphX bits=0
	rm -f *.o

digraphX: main.o $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@
	mv $@ digraph$(bits)

digraph: digraph32 digraph16 digraph0

#################
# TEST BINARIES #
#################

digraph-tests32:
	$(MAKE) digraph-testsX bits=32
	rm -f *.o

digraph-tests16:
	$(MAKE) digraph-testsX bits=16
	rm -f *.o

digraph-tests0:
	$(MAKE) digraph-testsX bits=0
	rm -f *.o

digraph-testsX: tests.o $(OBJ) $(TEST_OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@
	mv $@ digraph-tests$(bits)

digraph-tests: digraph-tests32 digraph-tests16 digraph-tests0

####################
# GF PERF BINARIES #
####################

gf-perf32:
	$(MAKE) gf-perfX bits=32
	rm -f *.o

gf-perf16:
	$(MAKE) gf-perfX bits=16
	rm -f *.o

gf-perf0:
	$(MAKE) gf-perfX bits=0
	rm -f *.o

gf-perfX: gf_perf.o $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@
	mv $@ gf-perf$(bits)

gf-perf: gf-perf32 gf-perf16 gf-perf0

#####################
# EXT PERF BINARIES #
#####################

extension-perf32:
	$(MAKE) extension-perfX bits=32
	rm -f *.o

extension-perf16:
	$(MAKE) extension-perfX bits=16
	rm -f *.o

extension-perf0:
	$(MAKE) extension-perfX bits=0
	rm -f *.o

extension-perfX: extension_perf.o $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@
	mv $@ extension-perf$(bits)

extension-perf: extension-perf32 extension-perf16 extension-perf0


.PHONY: clean
clean:
	rm -f $(addsuffix 0, $(BIN)) $(addsuffix 16, $(BIN)) $(addsuffix 32, $(BIN)) *.o *.s *.asm1 *.asm2 && cd nauty && git clean -xf && git checkout .

#########
# NAUTY #
#########

nauty/geng:
	cd nauty && [ -f config.log ] || ./configure && make geng

nauty/listg:
	cd nauty && [ -f config.log ] || ./configure && make listg

nauty/directg:
	cd nauty && [ -f config.log ] || ./configure && make directg

#################
# TEST COMMANDS #
#################

test32: digraph-tests32
	./digraph-tests32 -efux -d20 -t1000
	./digraph-tests32 -s -d10 -t100

test16: digraph-tests16
	./digraph-tests16 -egfux -d20 -t1000
	./digraph-tests16 -s -d10 -t100

test0: digraph-tests0
	./digraph-tests0 -egfux -d20 -n20 -t1000
	./digraph-tests0 -s -d10 -n20 -t100

test: test32 test16 test0

#geng-test: digraph-tests nauty/geng nauty/directg nauty/listg
#	mkdir -p geng-fail/$(vert)
#	nauty/geng -q $(vert) | nauty/directg -q | nauty/listg -aq | ./digraph-tests$(BITS) -c

#############
# ASM STUFF #
#############

%.s: %.cc
	$(CXX) -S $(CXXFLAGS) -fverbose-asm $^

%.asm1: %.s
	c++filt < $^ > $@

%.asm2: %.o
	objdump -d -S $^ > $@

include Makefile.dep
