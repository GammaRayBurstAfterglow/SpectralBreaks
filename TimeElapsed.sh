#!/bin/bash

# Sorts through all directories recursively, which searching all of the output files for the elapsed time, then outputting into elapsed.txt
find ~/GRBA -mindepth 2 -maxdepth 2 -type d | while read -r dir; do
        find "$dir" -type f -name "*_out_*.dat" -exec grep "Emission calculated!" {} + | sed 's/.*time = //; s/ sec//'
done > elapsed.txt

# Sorts so head and tail show any outliars in the data
sort elapsed.txt -o elapsed.txt
