#!/bin/bash

# what do here?
#SBATCH --time=04:00:00
#SBATCH --mem=2G
#SBATCH --output=ScriptOutput.log

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
