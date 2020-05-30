#!/bin/bash

usage()
{
cat <<EOF
Usage: $0 obsID [base_directory] [lowest channel] [# of channels]
    script will produce psrfits files in folder fits under obsID
    base_directory is optional, default is /group/mwavcs/vcs/obsID
    if base_directory is provided script expects the combined folder
    to exist, will create fits folder.
    For picket fence observations you have to supply base_directory,
    the lowest coarse channel of the desired subband, and the number
    of contiguous channels in the subband. Each subband needs to be
    dealt with individually
EOF
}
if [ $# -lt 1 ] ; then
    usage
    exit 1
fi
if [ $1 == '-h' ];then
    usage
    exit 0
fi

obsID=$1
basedir=${2:-/group/mwavcs/vcs/${obsID}}
lowchan=${3:-0}
n_chans=${4:-24}


# check for old pipe in target directory, delete if exists
if [ -p $FILE ${basedir}/ics/mk_psrfits_in ];then
    rm -rf ${basedir}/ics/mk_psrfits_in
    if [[ ! $? -eq 0 ]];then
	echo "Cannot remove old data pipe. Aborting..."
	exit $?
    fi
fi
# check for fits files in target directory that have the same obsID
# abort if found, will not delete anything.
if [ -f $FILE ${basedir}/ics/${obsID}*.fits ];then
    echo "There already exist fits files with obsID ${obsID} in ${basedir}/ics."
    echo "Please delete those before running this script. Aborting here..."
    exit 1
fi
length=`ls ${basedir}/combined/*ics.dat | wc -l`
#echo "Length of Obs: $length"

#OLDIFS=$IFS
#IFS='.'
first_file=(`ls ${basedir}/combined/*ics.dat | head -n 1`)
first_time=${first_file:-18:10}
#IFS=$' \t\n'

time_out=(`timeconvert.py --gps=$first_time 2>/dev/null`)

start=${time_out[2]}T${time_out[3]}
#echo "Obs. Date: $start"

#IFS=$OLDIFS

# determine the bandwidth
bandwidth=`bc -l <<< "1.28 * ${n_chans}"`

# get the centre frequency
if [ ${lowchan} -eq 0 ];then
    freq=`obs_query.py $1 -cf`
    #freq=${freq_out}
else
    freq=`bc -l <<< "(1.28 * ${lowchan} - 0.64) + ${bandwidth} / 2."`
fi

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

mkdir ${basedir}/ics
cd ${basedir}/ics
echo -e "${basedir}/ics/mk_psrfits_in\n\n$1\n\n${length}\n${USER}\n\n$1\n\n\nMWA-G0024\n${start}\n\n${freq}\n${bandwidth}\n${RA}\n${Dec}\n${Azimuth}\n${Zenith}\n\n${lst_s}\n${utc_s}\n${mjd_whole}\n\n\n\n\n\n\n\n\n" | make_psrfits & sleep 1.0
cat ${basedir}/combined/*ics.dat > ${basedir}/ics/mk_psrfits_in
wait
rm -rf ${basedir}/ics/mk_psrfits_in
