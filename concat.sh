#!/bin/bash -l

# Small script to conctatenate all of the output files from pabeam.py


BASENAME=$1 # base name of the output files (everything before the .dat)
FNAME=$2 # name of output file

# write the header part of the file
{
echo "##File Type: Far field"
echo "##File Format: 3"
echo "##Source: mwa_tiedarray"
echo "##Date: 2016-11-14 15:14:00"
echo "** File exported by FEKO kernel version 7.0.1-482"
echo
echo "#Request Name: FarField"
echo "#Frequency:   184960000.0"
echo "#Coordinate System: Spherical"
echo "#No. of Theta Samples: 9000"
echo "#No. of Phi Samples: 36000"
echo "#Result Type: Gain"
echo "#No. of Header Lines: 1"
echo '#       "Theta"             "Phi"           "Re(Etheta)"       "Im(Etheta)"        "Re(Ephi)"         "Im(Ephi)"       "Gain(Theta)"       "Gain(Phi)"       "Gain(Total)"'
} >> $FNAME

# cat the data from each file into
for i in {${START}..${END}};do cat ${BASENAME}.${i}.dat >> $FNAME;done


