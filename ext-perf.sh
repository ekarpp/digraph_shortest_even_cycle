#!/bin/bash

uname -a
cat /etc/*release
g++ --version

make clean
make digraph-tests
make extension-perf

make test

ARGS="-s123 -t10000000"

echo "32 bit"
./extension-perf32 $ARGS
echo "16 bit"
./extension-perf16 $ARGS
echo "0 bit"
./extension-perf0 $ARGS
