#!/usr/bin/bash

make digraph

total=0
fail=0

echo

for d in geng-fail/*
do
	echo $d
	for f in $d/*
	do
		algo=$(./digraph16 -qf $f | grep -v seed)
		brute=$(./digraph16 -qbf $f | grep -v seed)

		if [ "$algo" != "$brute" ]; then
			fail=$((fail+1))
		fi

		total=$((total+1))
	done
	echo "$fail/$total failed"
	echo
done

