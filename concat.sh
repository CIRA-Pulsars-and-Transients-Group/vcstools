#!/bin/bash -l

# Small script to conctatenate all of the output files from pabeam.py


BASENAME=$1 # base name of the output files (everything before the .{rank}.dat)
FNAME=$2 # name of output file
NPROCESSES=$3 # total number of processes used = total number of data files

# cat the data from each file into
n=`bc -l <<< "${NPROCESSES}-1"`
if [ $n -eq 0 ]; then
	echo "only 1 file created - just renaming..."
	echo "mv ${BASENAME}.0.dat ${FNAME}"
	mv ${BASENAME}.0.dat ${FNAME}
else
	echo "concatenating ${n} files together..."
	for i in $(eval echo "{0..$n}");do echo "concatenating file $i"; cat ${BASENAME}.${i}.dat >> $FNAME;done
fi

