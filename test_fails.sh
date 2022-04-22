#!/usr/bin/bash

make digraph

total=0
fail=0

for f in geng-fail/$1/*
do
	algo=$(./digraph16 -qf $f | grep -v seed)
	brute=$(./digraph16 -qbf $f | grep -v seed)

	if [ "$algo" != "$brute" ]; then
		fail=$((fail+1))
	fi

	total=$((total+1))
done

echo "$fail/$total failed"
