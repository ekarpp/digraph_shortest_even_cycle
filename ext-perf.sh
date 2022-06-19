#!/bin/bash

hostname
uname -a
cat /etc/*release
g++ --version

make clean
make objects -j 20
make digraph-tests
make extension-perf

make test

ARGS="-s123 -t 30"

for b in 0 16 32
do
    echo "$b bits"
    ./extension-perf${b} $ARGS
    ./extension-perf${b}-PAR $ARGS
done
