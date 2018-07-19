#!/bin/bash -l
#SBATCH --time=00:20:00
#SBATCH --nodes=25
 
module load cudatoolkit
module load cfitsio
export RTSBIN=/group/mwaops/PULSAR/bin/
 
if [ $# -lt 1 ] ; then
    echo "Usage: $0 working_dir rts_infile"  1>&2
    echo "    script will CD to working dir then execute RTS with rts_infile"  1>&2
    exit 1
fi
 
if [ "$RTSBIN" == "" ] ; then
    echo "RTSBIN is not set. Exiting."
    exit 1
fi
 
echo "changing to $1"
cd $1
res=$?
if [ $res -ne 0 ] ; then
    echo "Failed to cd to working dir. exiting..."
    exit 1
fi
 
echo "Executing RTS from $RTSBIN on $HOSTNAME"
 
aprun -n 25 -N 1  $RTSBIN/rts_gpu $2
exit $?
#####
# end of script
#####
