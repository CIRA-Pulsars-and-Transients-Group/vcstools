#!/usr/bin/env python3

import numpy as np
from astropy.time import Time
import argparse
import sys
import os
import logging

logger = logging.getLogger(__name__)

pabeam_sbatch_header = """#!/bin/bash -l
#SBATCH --export=NONE
#SBATCH --account=mwavcs
#SBATCH --cluster=galaxy
#SBATCH --partition=gpuq
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={nprocesses}
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=make_pabeam_{obsid}_{time}_{freq:.2f}MHz_%j.out

module load python/3.6.3
module load argparse
module load numpy/1.13.3
module load astropy
module load mpi4py
module use /group/mwa/software/modulefiles
module load vcstools/master
"""

pabeam_params = """nprocesses={nprocesses}
obsid={obsid}
freq={freq}
eff={eff}
ra='"{ra}"'
dec='"{dec}"'
flags="{flags}"
tres={tres}
pres={pres}
obstime={time}
odir="{odir}"
pabeam=/group/mwa/software/vcstools/vcstools/beam_sim/python/pabeam.py
"""

pabeam_base_cmd = """srun --export=all -u -n ${nprocesses} python3 ${pabeam} -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir}"""

pabeam_concat_cmd= """BASENAME={0} # base name of the output files (everything before the .[rank].dat)
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
fi
"""


showspec_sbatch_header = """#!/bin/bash -l
#SBATCH --export=NONE
#SBATCH --account=mwavcs
#SBATCH --cluster=galaxy
#SBATCH --partition=gpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --output={outfile}

module use /group/mwa/software/modulefiles
module load gsm
module load showspec
"""

showspec_params = """maploc={maploc}
showspec=`which showspec`
createmap=`which create_skymap.sh`

obsid={obsid}  # observation ID
ux={unixtime}  # unix time for actual beam evalutation 
gps={gpstime}  # GPS time for the actual beam evaluation
ra="{ra}"  # RA of source
dec="{dec}"  # DEC of source

freq={freq:.2f}  # actual frequency (in MHz)
freq_start={freq_start}  # floor of $freq
freq_end={freq_end}  # ceiling of $freq
file={infile}  # input beam pattern file name

# if the skymap doesn't already exist for the required frequency, make it...
if [ ! -f ${{maploc}}/skymaps/gsm_${{freq}}MHz_ni${{freq}}.out ]; then
    $createmap $freq $maploc
fi
"""

showspec_base_cmd = """
echo "$showspec ni_list -s ${{ux}} -i 0 -f ${{gps}}.spec -c 0 -b $file -q sun=0 -q save_map_fits=0 -q map_file_base=${{ux}}_ -q save_map_at_freq=${{freq}} -q freq_start=${{freq_start}} -q freq_end=${{freq_end}} -q cache_on=0 -q ant_cache_on=1 -q site=mwa -q ant_eff=0 -q binary_input=1 -q ant_rotation_deg=0.00 -q save_pattern_map_db=0 -p $maploc/skymaps -q max_feko_theta=90 -q bradley=1 -q sources=1 -q MAX_THETA_COUNT={ntheta} -q MAX_PHI_COUNT={nphi} -q phi_column_index={azcol} -q theta_column_index={zacol} -q gain_column_index={gaincol}"

srun --export=all -u $showspec ni_list -s ${{ux}} -i 0 -f ${{gps}}.spec -c 0 -b $file -q sun=0 -q save_map_fits=0 -q map_file_base=${{ux}}_ -q save_map_at_freq=${{freq}} -q freq_start=${{freq_start}} -q freq_end=${{freq_end}} -q cache_on=0 -q ant_cache_on=1 -q site=mwa -q ant_eff=0 -q binary_input=1 -q ant_rotation_deg=0.00 -q save_pattern_map_db=0 -p $maploc/skymaps -q max_feko_theta=90 -q bradley=1 -q sources=1 -q MAX_THETA_COUNT={ntheta} -q MAX_PHI_COUNT={nphi} -q phi_column_index={azcol} -q theta_column_index={zacol} -q gain_column_index={gaincol}
"""




def write_batch_files(tmin, tmax, step, thetares, phires, nnodes, ra, dec, obsid, freq, eff, flaggedtiles, maploc, write, odir):

    times = np.arange(tmin, tmax, step=step)
    times = np.append(times, tmax)

    nprocesses = 8 * nnodes
    flags = " ".join(flaggedtiles)

    if write:
        pabeam_run_cmd = pabeam_base_cmd + " --write"
    else:
        pabeam_run_cmd = pabeam_base_cmd


    # Loop over all times
    for i in range(len(times)):
        fname = "make_pabeam_{0}_{1:.2f}MHz.batch".format(times[i], freq/1e6)
        onamebase = "{0}_{1}_{2:.2f}MHz_{3}_{4}".format(obsid, float(times[i]), freq/1e6, ra, dec)

        header_str = pabeam_sbatch_header.format(nodes=nnodes, nprocesses=nprocesses, obsid=obsid, time=float(times[i]), freq=freq/1e6)
        params_str = pabeam_params.format(nprocesses=nprocesses, obsid=obsid, freq=freq, eff=eff, ra=ra, dec=dec, flags=flags, tres=thetares, pres=phires, time=times[i], odir=odir)

        with open(fname,'w') as f:
            f.write(header_str)
            f.write("\n\n")
            f.write(params_str)
            f.write("\n")
            f.write("echo {0}\n".format(pabeam_run_cmd))
            f.write("{0}\n\n".format(pabeam_run_cmd))

            # Write the concatenation bash commands
            f.write(pabeam_concat_cmd.format(onamebase, onamebase + ".dat"))
            
            # Remove the partial beam pattern files written by processes
            f.write("rm {0}\n".format(onamebase + ".*.dat"))

        # Now write the showspec batch for this time
        write_showspec_batch(times[i], obsid, ra, dec, freq, (90/thetares)+1, 360/phires, onamebase+".dat", maploc)



def write_showspec_batch(time, obsid, ra, dec, freq, ntheta, nphi, infile, maploc):
    fname = "showspec_{0}_{1:.2f}MHz.batch".format(time, freq/1e6)
    oname = fname.replace(".batch", ".out")
    
    unix = int(Time(int(time), format='gps').unix)
    fstart = int(np.floor(freq/1e6))
    fend = int(np.ceil(freq/1e6))

    # by default from the Python simulation code, the azimtuh, zenith angle and beam powers 
    # are in the following columns
    azcol = 1
    zacol = 0
    gaincol = 8

    with open(fname,'w') as f:

        header_str = showspec_sbatch_header.format(outfile=oname)
        params_str = showspec_params.format(maploc=maploc, obsid=obsid, unixtime=unix, gpstime=time,
                                             ra=ra, dec=dec, freq=freq/1e6, freq_start=fstart, freq_end=fend,
                                             infile=infile)
        run_str = showspec_base_cmd.format(ntheta=int(ntheta), nphi=int(nphi), azcol=azcol, zacol=zacol, gaincol=gaincol)
            
        f.write(header_str)
        f.write("\n\n")
        f.write(params_str)
        f.write("\n")
        f.write(run_str)



if __name__ == "__main__":

    loglevels = dict(DEBUG=logging.DEBUG,
                    INFO=logging.INFO,
                    WARNING=logging.WARNING,
                    ERROR = logging.ERROR)


    # Argument parsing
    parser = argparse.ArgumentParser(description="Simple script to help write the batch scripts required for running the tied-array beam simulations over multiple epochs")

    parser.add_argument("--tmin", type=int, help="GPS time for first evaluation")
    parser.add_argument("--tmax", type=int, help="GPS time for last evaluation")
    parser.add_argument("--step", type=int, help="Time step between each evaluation (in seconds) [default: 500]", default=500)
    parser.add_argument("--thetares", type=float, help="Resolution of theta (zenith angle) grid, in degrees [default: 0.05]", default=0.05)
    parser.add_argument("--phires", type=float, help="Resolution of phi (azimuth) grid, in degrees [default: 0.05]", default=0.05)
    parser.add_argument("--nodes", type=int, help="Number of nodes to use per run (more than 10 is overkill and will actually be slower...) [default: 1]", default=1)
    parser.add_argument("--obsid", type=int, help="Observation ID")
    parser.add_argument("--freq", type=float, help="Observing frequency in Hz")
    parser.add_argument("--eff", type=float, help="Radiation efficiency [default: 1.0]", default=1.0)
    parser.add_argument("--flagged", nargs='+', help="Flagged tiles (as in RTS flagged_tiles.txt)")

    parser.add_argument("--ra", type=str, help="RAJ2000 of target in hh:mm:ss.ss (use = to assign option)")
    parser.add_argument("--dec", type=str, help="DECJ2000 of target in dd:mm:ss.ss (use = to assign option)") # only because argparse can't handle negative arguments...

    parser.add_argument("--maploc", type=str, help="Path to directory where skymaps/ exists (or is to be created). The GSM temperature maps exist there or will be created.", default="$PWD")

    parser.add_argument("--write", action='store_true', help="Write output files to disk")
    parser.add_argument("--odir", type=str, help="Output directory")
    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO", choices=loglevels.keys(), default="INFO")



    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    logger.warn("FULL_EE BEAM MODEL NOT AVAILABLE IN PYTHON 3 YET. ANALYTIC BEAM MODEL WILL BE USED. FOR GREATER ACCURACY PLEASE USE PYTHON 2 VRESION.")

    write_batch_files(args.tmin, args.tmax, args.step, args.thetares, args.phires, 
                      args.nodes, args.ra, args.dec, args.obsid, args.freq, args.eff, args.flagged,
                      args.maploc, args.write, args.odir)
