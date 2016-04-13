#!/bin/bash

if [ $# -lt 1 ] ; then
    echo "Usage: $0 obsID [base_directory]"  1>&2
    echo "    script will produce psrfits files in folder fits under obsID"  1>&2
    echo "    base_directory is optional, default is /scratch2/mwaops/vcs/obsID"  1>&2
    echo "    if base_directory is provided script expects the combined folder to exist, will create fits folder."  1>&2
    exit 1
fi

obsID=$1
basedir=${2:-/scratch2/mwaops/vcs/${obsID}}

# check for old pipe in target directory, delete if exists
if [ -p $FILE ${basedir}/fits/mk_psrfits_in ];then
    rm -rf ${basedir}/fits/mk_psrfits_in
    if [[ ! $? -eq 0 ]];then
	echo "Cannot remove old data pipe. Aborting..."
	exit $?
    fi
fi
# check for fits files in target directory that have the same obsID
# abort if found, will not delete anything.
if [ -f $FILE ${basedir}/fits/${obsID}*.fits ];then
    echo "There already exist fits files with obsID ${obsID} in ${basedir}/fits."
    echo "Please delete those before running this script. Aborting here..."
    exit 1
fi
length=`ls ${basedir}/combined/*ics.dat | wc -l`
#echo "Length of Obs: $length"

#OLDIFS=$IFS
#IFS='.'
first_file=(`ls ${basedir}/combined/*ics.dat | head -n 1`)
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

mkdir ${basedir}/fits
cd ${basedir}/fits
echo -e "${basedir}/fits/mk_psrfits_in\n\n$1\n\n${length}\n${USER}\n\n$1\n\n\nMWA-G0024\n${start}\n\n${freq}\n30.72\n${RA}\n${Dec}\n${Azimuth}\n${Zenith}\n\n${lst_s}\n${utc_s}\n${mjd_whole}\n\n\n\n\n\n\n\n\n" | make_psrfits & sleep 1.0
cat ${basedir}/combined/*ics.dat > ${basedir}/fits/mk_psrfits_in
wait
rm -rf ${basedir}/fits/mk_psrfits_in
