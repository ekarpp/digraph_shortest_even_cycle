#!/bin/bash

#SBATCH --time=00:01:00
#SBATCH --nodelist=pe[1-48,65-81]
#SBATCH --output=extension-perf.log

uname -a
cat /etc/*release
g++ --version

make clean
make extension-perf

ARGS="-s123 -t10000000"

echo "32 bit"
./extension-perf32 $ARGS
echo "16 bit"
./extension-perf16 $ARGS
echo "0 bit"
./extension-perf0 $ARGS
