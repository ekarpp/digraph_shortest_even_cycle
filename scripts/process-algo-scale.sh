#!/bin/bash

REPEATS=6
graphs=("CONFIG" "COMPLETE" "CYCLE" "ERDOS")

for g in ${graphs[@]}; do
    grep -A $(($REPEATS*3)) $g $1 | grep --line-buffered computed \
        | awk 'NR%6 !=1 { printf ("%s %s\n", $4, $7) }' > ${g}_${1}.plot
done
