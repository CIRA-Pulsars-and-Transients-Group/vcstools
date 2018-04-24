#!/usr/bin/env python

import numpy as np
from astropy.time import Time
import argparse
import sys
import os

def write_batch_files(tmin,tmax,step,thetares,phires,nnodes,ra,dec,obsid,freq,eff,flaggedtiles,maploc,write,odir):

    times = np.arange(tmin,tmax,step=step)
    times = np.append(times,tmax)

    nprocesses = 20*nnodes
    flags = " ".join(flaggedtiles)


    for i in range(len(times)-1):
        fname = "make_pabeam_{0}_{1:.2f}MHz.batch".format(times[i],freq/1e6)
        onamebase = "{0}_{1}_{2:.2f}MHz_{3}_{4}".format(obsid,float(times[i]),freq/1e6,ra,dec)
        with open(fname,'w') as f:
            f.write("#!/bin/bash -l\n\n")
            f.write("#SBATCH --account=mwaops\n#SBATCH --gid=mwaops\n#SBATCH --cluster=galaxy\n#SBATCH --partition=workq\n")
            f.write("#SBATCH --nodes={0}\n".format(nnodes))
            f.write("#SBATCH --time=12:00:00\n#SBATCH --output={0}\n\n".format(fname.replace(".batch",".out")))
            f.write("nprocesses={0}\nobsid={1}\nfreq={2}\neff={3}\n".format(nprocesses,obsid,freq,eff))
            f.write('ra=\'"{0}"\'\ndec=\'"{1}"\'\nflags="{2}"\ntres={3}\npres={4}\n'.format(ra,dec,flags,thetares,phires))
            f.write('obstime={0}\nodir="{1}"\n\n'.format(times[i],odir))
            # submit the next script with a dependency on this one
            f.write("sbatch --depend=afterany:{0} {1}\n".format("${SLURM_JOB_ID}","make_pabeam_{0}_{1}MHz.batch".format(times[i+1],freq/1e6)))  
            if write:
                f.write('echo "srun -u -n ${nprocesses} pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir} --write"\n')
                f.write("srun -u -n ${nprocesses} pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir} --write\n\n")
            else:
                f.write('echo "srun -u -n ${nprocesses} pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir}"\n')
                f.write("srun -u -n ${nprocesses} pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir}\n\n")

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

            
            # Now write the showspec batch for this time
            write_showspec_batch(times[i],obsid,ra,dec,freq,(90/thetares)+1,360/phires,onamebase+".dat")

    #write the last file without submitting the next one
    fname = "make_pabeam_{0}_{1:.2f}MHz.batch".format(times[-1],freq/1e6)
    onamebase = "{0}_{1}_{2:.2f}MHz_{3}_{4}".format(obsid,float(times[-1]),freq/1e6,ra,dec)
    with open(fname,'w') as f:
        f.write("#!/bin/bash -l\n\n")
        f.write("#SBATCH --account=mwaops\n#SBATCH --gid=mwaops\n#SBATCH --cluster=galaxy\n#SBATCH --partition=workq\n")
        f.write("#SBATCH --nodes={0}\n".format(nnodes))
        f.write("#SBATCH --time=12:00:00\n#SBATCH --output={0}\n\n".format(fname.replace(".batch",".out")))
        f.write("nprocesses={0}\nobsid={1}\nfreq={2}\neff={3}\n".format(nprocesses,obsid,freq,eff))
        f.write('ra=\'"{0}"\'\ndec=\'"{1}"\'\nflags="{2}"\ntres={3}\npres={4}\n'.format(ra,dec,flags,thetares,phires))
        f.write('obstime={0}\nodir="{1}"\n\n'.format(times[-1],odir))
        if write:
            f.write('echo "srun -u -n ${nprocesses} pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir} --write"\n')
            f.write("srun -u -n ${nprocesses} pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir} --write\n\n")
        else:
            f.write('echo "srun -u -n ${nprocesses} pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir}"\n')
            f.write("srun -u -n ${nprocesses} pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir}\n\n")
        
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

        # Now write the showspec batch for this time
        write_showspec_batch(times[-1],obsid,ra,dec,freq,(90/thetares)+1,360/phires,onamebase+".dat",maploc)


def write_showspec_batch(time,obsid,ra,dec,freq,ntheta,nphi,out,maploc):
    fname = "showspec_{0}_{1:.2f}MHz.batch".format(time,freq/1e6)
    oname = fname.replace(".batch",".out")
    
    unix = int(Time(int(time), format='gps').unix)
    fstart = int(np.floor(freq/1e6))
    fend = int(np.ceil(freq/1e6))

    azcol = 1
    zacol = 0
    gaincol = 8

    with open(fname,'w') as f:
        f.write("#!/bin/bash -l\n\n#SBATCH --account=mwaops\n#SBATCH --cluster=galaxy\n#SBATCH --partition=workq\n#SBATCH --nodes=1\n#SBATCH --time=3:00:00\n")
        f.write("#SBATCH --output={0}\n\n".format(oname))
        f.write("maploc={0}\nshowspec=`which showspec`\ncreatemap=`which create_skymap.sh`\n\n".format(maploc))
        
        infostr = """
obsid={0}    # observation ID
ux={1}       # unix time for actual beam evalutation 
gps={2}      # GPS time for the actual beam evaluation
ra="{3}"   # RA of source
dec="{4}"  # DEC of source

freq={5:.2f}        # actual frequency (in MHz)
freq_start={6}      # floor of $freq
freq_end={7}        # ceiling of $freq
file={8}            # output beam pattern file name\n\n""".format(obsid, unix, time, ra, dec, freq/1e6, fstart, fend, out)
        
        f.write(infostr)

        checkstr = """
# if the skymap doesn't already exist for the required frequency, make it...
if [ ! -f ${maploc}/skymaps/gsm_${freq}MHz_ni${freq}.out ]; then
    $createmap $freq $maploc
fi\n\n"""

        f.write(checkstr)

        runstr = """
# now run showspec
echo "$showspec ni_list -s ${{ux}} -i 0 -f ${{gps}}.spec -c 0 -b $file -q sun=0 -q save_map_fits=0 -q map_file_base=${{ux}}_ -q save_map_at_freq=${{freq}} -q freq_start=${{freq_start}} -q freq_end=${{freq_end}} -q cache_on=0 -q ant_cache_on=1 -q site=mwa -q ant_eff=0 -q binary_input=1 -q ant_rotation_deg=0.00 -q save_pattern_map_db=0 -p $maploc -q max_feko_theta=90 -q bradley=1 -q sources=1 -q MAX_THETA_COUNT={0} -q MAX_PHI_COUNT={1} -q phi_column_index={2} -q theta_column_index={3} -q gain_column_index={4}"

$showspec ni_list -s ${{ux}} -i 0 -f ${{gps}}.spec -c 0 -b $file -q sun=0 -q save_map_fits=0 -q map_file_base=${{ux}}_ -q save_map_at_freq=${{freq}} -q freq_start=${{freq_start}} -q freq_end=${{freq_end}} -q cache_on=0 -q ant_cache_on=1 -q site=mwa -q ant_eff=0 -q binary_input=1 -q ant_rotation_deg=0.00 -q save_pattern_map_db=0 -p $maploc -q max_feko_theta=90 -q bradley=1 -q sources=1 -q MAX_THETA_COUNT={0} -q MAX_PHI_COUNT={1} -q phi_column_index={2} -q theta_column_index={3} -q gain_column_index={4}""".format(int(ntheta), int(nphi), azcol, zacol, gaincol)

        f.write(runstr)





# Argument parsing

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

parser.add_argument("--maploc",type=str,help="Path to directory where skymaps/ exists, within which there are GSM temp. maps.", default="$PWD")

parser.add_argument("--write",action='store_true',help="Write output files to disk")
parser.add_argument("--odir",type=str,help="Output directory")

args = parser.parse_args()


ra = args.ra
dec = args.dec
write_batch_files(args.tmin,args.tmax,args.step,args.thetares,args.phires,args.nodes,ra,dec,args.obsid,args.freq,args.eff,args.flagged,args.maploc,args.write,args.odir)
