#!/usr/bin/env python

import numpy as np
import argparse
import sys


def write_batch_files(tmin,tmax,step,thetares,phires,nnodes,ra,dec,obsid,freq,eff,flaggedtiles,write,odir):

    times = np.arange(tmin,tmax,step=step)
    times = np.append(times,tmax)

    nprocesses = 20*nnodes
    flags = " ".join(flaggedtiles)


    for i in range(len(times)-1):
        fname = "make_pabeam_{0}_{1:.2f}MHz.batch".format(times[i],freq/1e6)
        onamebase = "{0}_{1}_{2:.2f}MHz_{3}_{4}".format(obsid,float(times[i]),freq/1e6,ra,dec)
        with open(fname,'w') as f:
            f.write("#!/bin/bash -l\n\n")
            f.write("#SBATCH --account=mwaops\n#SBATCH --partition=workq\n#SBATCH --nodes={0}\n".format(nnodes))
            f.write("#SBATCH --time=12:00:00\n#SBATCH --output={0}\n\n".format(fname.replace(".batch",".out")))
            f.write("nprocesses={0}\nobsid={1}\nfreq={2}\neff={3}\n".format(nprocesses,obsid,freq,eff))
            f.write('ra=\'"{0}"\'\ndec=\'"{1}"\'\nflags="{2}"\ntres={3}\npres={4}\n'.format(ra,dec,flags,thetares,phires))
            f.write('obstime={0}\nodir="{1}"\n\n'.format(times[i],odir))
            # submit the next script with a dependency on this one
            f.write("sbatch --depend=afterany:{0} {1}\n".format("${SLURM_JOB_ID}","make_pabeam_{0}_{1}MHz.batch".format(times[i+1],freq/1e6)))  
            if write:
                f.write('echo "aprun -n ${nprocesses} pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir} --write"\n')
                f.write("aprun -n ${nprocesses} pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir} --write\n\n")
            else:
                f.write('echo "aprun -n ${nprocesses} pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir}"\n')
                f.write("aprun -n ${nprocesses} pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir}\n\n")

            # write the concatenation bash commands
            concatstr= """
BASENAME={0} # base name of the output files (everything before the .[rank].dat)
FNAME={1} # name of output file

if [ -f $FNAME ]; then
    echo "found existing file ${{FNAME}}. Removing..."
    rm -f $FNAME
fi

# cat the data from each file into 1 beam patten file
n=$((nprocesses-1))
if [ $n -eq 0 ]; then
    echo "only 1 file created - just renaming..."
    echo "mv ${{BASENAME}}.0.dat ${{FNAME}}"
    mv ${{BASENAME}}.0.dat ${{FNAME}}
else
    echo "concatenating $nprocesses files together..."
    for i in $(eval echo "{{0..$n}}");do 
        echo "concatenating file $i" 
        cat ${{BASENAME}}.${{i}}.dat >> $FNAME
    done
fi\n""".format(onamebase, onamebase+".dat")

            #f.write("bash /group/mwaops/bmeyers/code/pabeam/concat.sh {0} {1} {2}\n".format(onamebase,onamebase+".dat",nprocesses))
            f.write(concatstr)
            
            # remove the partial beam pattern files written by processes
            f.write("rm {0}\n".format(onamebase+".*.dat"))


    #write the last file without submitting the next one
    fname = "make_pabeam_{0}_{1:.2f}MHz.batch".format(times[-1],freq/1e6)
    onamebase = "{0}_{1}_{2:.2f}MHz_{3}_{4}".format(obsid,float(times[-1]),freq/1e6,ra,dec)
    with open(fname,'w') as f:
        f.write("#!/bin/bash -l\n\n")
        f.write("#SBATCH --account=mwaops\n#SBATCH --partition=workq\n#SBATCH --nodes={0}\n".format(nnodes))
        f.write("#SBATCH --time=12:00:00\n#SBATCH --output={0}\n\n".format(fname.replace(".batch",".out")))
        f.write("nprocesses={0}\nobsid={1}\nfreq={2}\neff={3}\n".format(nprocesses,obsid,freq,eff))
        f.write('ra=\'"{0}"\'\ndec=\'"{1}"\'\nflags="{2}"\ntres={3}\npres={4}\n'.format(ra,dec,flags,thetares,phires))
        f.write('obstime={0}\nodir="{1}"\n\n'.format(times[-1],odir))
        if write:
            f.write('echo "aprun -n ${nprocesses} pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir} --write"\n')
            f.write("aprun -n ${nprocesses} pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir} --write\n\n")
        else:
            f.write('echo "aprun -n ${nprocesses} pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir}"\n')
            f.write("aprun -n ${nprocesses} pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir}\n\n")
        
        # write the concatenation bash commands
        concatstr= """
BASENAME={0} # base name of the output files (everything before the .[rank].dat)
FNAME={1} # name of output file

if [ -f $FNAME ]; then
    echo "found existing file ${{FNAME}}. Removing..."
    rm -f $FNAME
fi

# cat the data from each file into 1 beam patten file
n=$((nprocesses-1))
if [ $n -eq 0 ]; then
    echo "only 1 file created - just renaming..."
    echo "mv ${{BASENAME}}.0.dat ${{FNAME}}"
    mv ${{BASENAME}}.0.dat ${{FNAME}}
else
    echo "concatenating $nprocesses files together..."
    for i in $(eval echo "{{0..$n}}");do 
        echo "concatenating file $i" 
        cat ${{BASENAME}}.${{i}}.dat >> $FNAME
    done
fi\n""".format(onamebase, onamebase+".dat")
        f.write(concatstr)
        
        # remove the partial beam pattern files written by processes
        f.write("rm {0}\n".format(onamebase+".*.dat"))





parser = argparse.ArgumentParser(description="Simple script to help write the batch scripts required for running the tied-array beam simulations over multiple epochs")

parser.add_argument("--tmin",type=int,help="GPS time for first evaluation")
parser.add_argument("--tmax",type=int,help="GPS time for last evaluation")
parser.add_argument("--step",type=int,help="Time step between each evaluation (in seconds) [default: 500]",default=500)
parser.add_argument("--thetares",type=float,help="Resolution of theta (zenith angle) grid, in degrees [default: 0.01]",default=0.01)
parser.add_argument("--phires",type=float,help="Resolution of phi (azimuth) grid, in degrees [default: 0.01]",default=0.01)
parser.add_argument("--nodes",type=int,help="Number of nodes to use per run (more than 10 is overkill and will actually be slower...) [default: 1]",default=1)
parser.add_argument("--obsid",type=int,help="Observation ID")
parser.add_argument("--freq",type=float,help="Observing frequency in Hz")
parser.add_argument("--eff",type=float,help="Radiation efficiency [default: 1.0]",default=1.0)
parser.add_argument("--flagged",nargs='+',help="Flagged tiles (as in RTS flagged_tiles.txt)")

parser.add_argument("--ra",type=str,help="RAJ2000 of target")
parser.add_argument("--dec",type=str,help="DECJ2000 of target (use = to assign option)") # only because of the old version of argparse Galaxy has can't handle negative arguments...

parser.add_argument("--write",type=bool,help="Write output (True) or just report statistics (False)  [default: True]",default=True)
parser.add_argument("--odir",type=str,help="Output directory")

args = parser.parse_args()


ra = args.ra
dec = args.dec
write_batch_files(args.tmin,args.tmax,args.step,args.thetares,args.phires,args.nodes,ra,dec,args.obsid,args.freq,args.eff,args.flagged,args.write,args.odir)
