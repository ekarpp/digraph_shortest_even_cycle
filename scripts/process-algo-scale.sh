#!/bin/bash

REPEATS=6
graphs=("CONFIG" "COMPLETE" "CYCLE" "ERDOS")

for g in ${graphs[@]}; do
    grep -A $(($REPEATS*3)) $g $1 | grep computed > $g_$1
done
