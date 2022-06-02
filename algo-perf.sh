#!/bin/bash

uname -a
cat /etc/*release
g++ --version

make clean
make digraph-tests16
make digraph16

make test16

REPEATS=50
DEG=6

for v in 8 16 24 32 40 48 56 64
do
	cd graphs
	./k_gen.py $v
	./c_gen.py $v
	./config_model.py $v $DEG $REPEATS
	./erdos_renyi.py $v "edges" $REPEATS
	cd ..

	echo "COMPLETE + CYCLE"
	for (( i=1; i<=$REPEATS; i++ ))
	do
		./digraph16 -qtf graphs/complete/k$v -s $RANDOM >> "results/complete_$v.out"
		./digraph16 -qtf graphs/cycle/c$v -s $RANDOM >> "results/cycle_$v.out"
	done

	echo "CONFIG"
	for file in graphs/config/*
	do
		./digraph16 -qtf $file -s $RANDOM >> "results/config_$v.out"
	done
	rm graphs/config/cm*

	echo "ERDOS"
	for file in graphs/erdos_renyi/*
	do
		./digraph16 -qtf $file -s $RANDOM >> "results/erdos_$v.out"
	done
	rm graphs/erdos_renyi/er*
done
