#!/bin/bash

#SBATCH --time=00:02:30
#SBATCH --mem=10G
#SBATCH --partition=batch-skl
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1
#SBATCH --output=mem-bench.out
module load gcc/11.2.0

hostname
uname -a
cat /etc/*release
g++ --version

make clean
make mem-bench

./mem-bench
