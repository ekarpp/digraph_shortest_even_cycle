#!/usr/bin/bash

make digraph

for f in geng-fail/$1/*
do
	echo "algo"
	./digraph16 -qf $f
	echo "brute"
	./digraph16 -qbf $f
	echo
done
