#!/bin/bash

#SBATCH --time=00:35:00
#SBATCH --mem=52G
#SBATCH --partition=batch-hsw
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=ext-perf-seq.out
module load gcc/11.2.0

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
