#!/usr/bin/env bash

pwait() {
    while [ $(jobs -p | wc -l) -ge $1 ]; do
        sleep 0.33
    done
}

message="
\nYou have to at least supply the obsID via -o \n
By default 2 untar jobs run in parralel \n
Option are:\n
\t -o \t obsID \n
\t -j \t number of parallel jobs (default 2)\n
\t\t if run outside mpirun/aprun should at least be 2 to make sure \n
\t\t script works (script itself counts as job outside mpirun/aprun)\n
\t -w \t working directory that contains the tarballs \n
\t\t (default: /astro/mwavcs/vcs/\${obsID}/combined)\n
\t -b \t gpstime of first tarball\n
\t -e \t gpstime of last tarball (by default all tarballs in ${workdir} will be untarred.)\n
\t -k \t if supplied: tarballs will be kept, else deleted after unpacking\n
"
while getopts ":j:w:b:e:o:k" flag;do
    case $flag in
	o)
	    if [[ $OPTARG = -* ]];then
		echo -e $message
		exit 1
	    fi
	    obsid=$OPTARG
	    ;;
	j)
	    maxjobs=$OPTARG
	    ;;
	w)
	    workdir=$OPTARG
	    ;;
	b)
	    start=$OPTARG
	    ;;
	e)
	    stop=$OPTARG
	    ;;
	k)
	    keep=1
	    ;;
    esac
done

if [[ -z $obsid ]];then
    echo -e $message
    exit 1
fi
if [[ -z $maxjobs ]];then
    maxjobs=1
fi
if [[ -z $workdir ]];then
    workdir=/astro/mwavcs/vcs/${obsid}/combined/
fi
if [ ! -d "${workdir}" ];then
    echo "working directory ${workdir} does not exist."
    exit 1
fi
if [[ -z $start ]];then
    if [[ -z $stop ]];then
	# as no times were supplied we'll untar all tarballs
	# to that end we extract the gpstimes of the first
	# and last files
	cd ${workdir}
	first_file=`ls ${obsid}_*_combined.tar | head -1`
	last_file=`ls ${obsid}_*_combined.tar | tail -1`
	start=${first_file:11:10}
	stop=${last_file:11:10}
    else
	echo -e "\nIf you supply a stop time you should also supply a start time\n"
	echo -e $message
	exit 1
    fi
else
    if [[ -z $stop ]];then
	echo -e "\nIf you supply a start time you should also supply a stop time\n"
	echo -e $message
	exit 1
    fi
fi
if [[ $stop -lt $start ]];then
    echo "starttime $start is later than stoptime $stop. Aborting"
    exit 1
fi
if [[ -z $keep ]];then
    keep=0
fi


# okay, we should have all the info we need from here
# and can just loop through all times
cd $workdir
for gpstime in `seq $start $stop`;do
    file=${obsid}_${gpstime}_combined.tar
    if [ ! -f ${file} ];then
	echo "${file} does not exist."
	continue
    fi
    if [[ $keep -eq 0 ]];then
        (tar xvf ${file} && rm -rfv ${file}) &
        pwait $maxjobs
    else
        tar xvf $file &
        pwait $maxjobs
    fi
done
wait
echo "Finished untarring this set of tar balls."

# if -k is NOT supplied the tar balls will be deleted:
#if [[ $keep -eq 0 ]];then
#    for gpstime in `seq $start $stop`;do
#	rm -rf ${obsid}_${gpstime}_combined.tar
#    done
#fi

