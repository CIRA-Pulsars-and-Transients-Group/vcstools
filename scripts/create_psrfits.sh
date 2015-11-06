#!/bin/bash

if [ $# -lt 1 ] ; then
    echo "Usage: $0 obsID"  1>&2
    echo "    script will produce psrfits files in folder fits under obsID"  1>&2
    exit 1
fi

if [ -f $FILE /scratch/mwaops/vcs/$1/fits/mk_psrfits_in ];
then
    rm /scratch/mwaops/vcs/$1/fits/mk_psrfits_in
fi

export PATH=${PATH}:/home/fkirsten/software/galaxy-scripts/scripts/
export PYTHONPATH=${PYTHONPATH}:/home/fkirsten/software/galaxy-scripts/scripts/
export PYTHONPATH=${PYTHONPATH}:/group/mwaops/stremblay/MandC_Core/
source ~/.modules-gnu
module load python/2.6.8 psycopg2 pyephem matplotlib

length=`ls /scratch/mwaops/vcs/$1/combined/*ics.dat | wc -l`
#echo "Length of Obs: $len"

#OLDIFS=$IFS
#IFS='.'
first_file=(`ls /scratch/mwaops/vcs/$1/combined/*ics.dat | head -n 1`)
first_time=${first_file:(-18):10}
#IFS=$' \t\n'

time_out=(`timeconvert.py --gps=$first_time 2>/dev/null`)

start=${time_out[2]}T${time_out[3]}
#echo "Obs. Date: $start"

#IFS=$OLDIFS

# get the centre frequency
freq_out=( `fits_db.py $1` )
freq=${freq_out[2]}

# get RA, Dec, Azimuth, Zenith
coords=( `get_meta.py $1` )
RA=${coords[1]}
Dec=${coords[3]}
Azimuth=${coords[5]}
Zenith=${coords[7]}
#IFS=':'
lst=${time_out[12]}
lst_s=`echo "(${lst:0:2}*3600) + (${lst:3:2}*60) + ${lst:6:2}" | bc -l`
#echo "LST Seconds: $lst_s"

utc=${time_out[3]}
utc_s=`echo "(${utc:0:2}*3600) + (${utc:3:2}*60) + ${utc:6:2}" | bc -l`
#echo "UTC Seconds: $utc_s"

#IFS='.'
mjd=(${time_out[9]})
mjd_whole=${mjd:0:5}
#echo "Truncated MJD: $mjd_whole"

mkdir /scratch/mwaops/vcs/$1/fits
cd /scratch/mwaops/vcs/$1/fits
echo -e "/scratch/mwaops/vcs/$1/fits/mk_psrfits_in\n\n$1\n\n${length}\n${USER}\n\n$1\n\n\nMWA-G0024\n${start}\n\n${freq}\n30.72\n${RA}\n${Dec}\n${Azimuth}\n${Zenith}\n\n${lst_s}\n${utc_s}\n${mjd_whole}\n\n\n\n\n\n\n\n\n" | ~/software/galaxy-scripts/bin/make_psrfits &
cat /scratch/mwaops/vcs/$1/combined/*ics.dat > /scratch/mwaops/vcs/$1/fits/mk_psrfits_in &
wait

#IFS=$OLDIFS
