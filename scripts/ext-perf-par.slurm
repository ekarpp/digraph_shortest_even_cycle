#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --mem=52G
#SBATCH --partition=batch-skl
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1
#SBATCH --output=ext-perf-par.out
module load gcc/11.2.0

hostname
uname -a
cat /etc/*release
g++ --version

make clean
make objects -j 24
make digraph-tests
make extension-perf

make test

ARGS="-s123 -t 30"

for b in 0 16 32
do
    echo "$b bits"
    ./extension-perf${b}-PAR $ARGS
done
