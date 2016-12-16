#!/bin/bash -l

# Small script to conctatenate all of the output files from pabeam.py


BASENAME=$1 # base name of the output files (everything before the .{rank}.dat)
FNAME=$2 # name of output file
NPROCESSES=$3 # total number of processes used = total number of data files

# cat the data from each file into
for i in {0..${NPROCESSES}};do cat ${BASENAME}.${i}.dat >> $FNAME;done


