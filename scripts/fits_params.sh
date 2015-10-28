#!/bin/bash

len=`ls /scratch/mwaops/vcs/$1/combined/*ics.dat | wc -l`
echo "Length of Obs: $len"

#OLDIFS=$IFS
#IFS='.'
first_file=(`ls /scratch/mwaops/vcs/$1/combined/*ics.dat | head -n 1`)
first_time=${first_file:(-18):10}
#IFS=$' \t\n'

time_out=(`timeconvert.py --gps=$first_time 2>/dev/null`)

start=${time_out[2]}T${time_out[3]}
echo "Obs. Date: $start"
#IFS=$OLDIFS


fits_db.py $1
get_meta.py $1


#IFS=':'
lst=${time_out[12]}
lst_s=`echo "(${lst:0:2}*3600) + (${lst:3:2}*60) + ${lst:6:2}" | bc -l`
echo "LST Seconds: $lst_s"

utc=${time_out[3]}
utc_s=`echo "(${utc:0:2}*3600) + (${utc:3:2}*60) + ${utc:6:2}" | bc -l`
echo "UTC Seconds: $utc_s"

#IFS='.'
mjd=(${time_out[9]})
mjd_whole=${mjd:0:5}
echo "Truncated MJD: $mjd_whole"



#IFS=$OLDIFS
