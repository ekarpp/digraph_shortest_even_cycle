#!/bin/bash

#SBATCH --time=00:30:00
#SBATCH --mem=35G
#SBATCH --partition=batch-hsw
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=ext-perf-seq.out

hostname
uname -a
cat /etc/*release
g++ --version

make clean
make objects
make digraph-tests
make extension-perf

make test

ARGS="-s123 -t 30"

for b in 0 16 32
do
    echo "$b bits"
    ./extension-perf${b} $ARGS
done
