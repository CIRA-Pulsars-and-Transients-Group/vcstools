#!/usr/bin/env python3

import numpy as np
from astropy.time import Time
import argparse
import sys
import os
import logging

from vcstools.metadb_utils import obs_max_min, get_common_obs_metadata, ensure_metafits
from vcstools.job_submit import submit_slurm
from vcstools.general_utils import mdir
from vcstools.config import load_config_file

logger = logging.getLogger(__name__)

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



def write_batch_files(obsid, begin, end,
                      ra, dec, freq, flaggedtiles,
                      step=500, thetares=0.05, phires=0.05,
                      nnodes=1, eff=1, beam_model='hyperbeam',
                      maploc="$PWD", odir=None, delays=[0] * 16,
                      write=True, write_showspec=False,
                      vcstools_version='master', metafits_loc=None):

    comp_config = load_config_file()

    times = np.arange(begin, end, step=step)
    times = np.append(times, end)

    nprocesses = 32 * nnodes
    flags = " ".join(flaggedtiles)

    if odir is None:
        # Make default directories
        product_dir = os.path.join(comp_config['base_data_dir'], obsid, 'pabeam', '{}_{}'.format(ra, dec))
        batch_dir   = os.path.join(comp_config['base_data_dir'], obsid, 'batch')
    else:
        product_dir = batch_dir = odir
    mdir(product_dir, 'Product Dir', gid=comp_config['gid'])
    mdir(batch_dir,   'Batch Dir',   gid=comp_config['gid'])

    # Loop over all times
    for i in range(len(times)):
        fname = "make_pabeam_{0}_{1}_{2}_{3:.2f}MHz".format(ra, dec, times[i], freq/1e6)
        onamebase = "{0}_{1}_{2:.2f}MHz_tres{3}_pres{4}_{5}_{6}".format(obsid, float(times[i]), freq/1e6, thetares, phires, ra, dec)

        commands = []
        # Write out params
        commands.append("nprocesses={}".format(nprocesses))
        commands.append("obsid={}".format(obsid))
        commands.append("""ra='"{}"'""".format(ra))
        commands.append("""dec='"{}"'""".format(dec))
        commands.append("freq={}".format(freq))
        commands.append("eff={}".format(eff))
        commands.append('flags="{}"'.format(flags))
        commands.append('delays="{}"'.format(delays))
        commands.append("tres={}".format(thetares))
        commands.append("pres={}".format(phires))
        commands.append("obstime={}".format(times[i]))
        commands.append("odir={}".format(product_dir))
        commands.append("metafits_loc={}".format(metafits_loc))
        commands.append('beam="{}"'.format(beam_model))
        # TODO remove this once hyperbeam is installed with python
        commands.append("export PYTHONPATH=$PYTHONPATH:/pawsey/mwa/software/python3/hyperbeam/v0.3.0/lib/python3.8/site-packages")

        # Main command
        pabeam_command = "srun --export=all -u -n ${nprocesses} pabeam.py " +\
                         "-o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --metafits ${metafits_loc} " +\
                         "--flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir} --beam_model ${beam} --delays ${delays}"
        if write:
            pabeam_command = pabeam_command + " --write"

        commands.append('cd {}'.format(product_dir))
        commands.append('echo "{}"'.format(pabeam_command))
        commands.append(pabeam_command)

        # Combine the output files into one
        commands.append(pabeam_concat_cmd.format(onamebase, onamebase + ".dat"))

        # Remove the partial beam pattern files written by processes
        commands.append("rm {0}\n".format(onamebase + ".*.dat"))
        module_list = ['mpi4py', 'hyperbeam/v0.3.0']
        submit_slurm(fname, commands, batch_dir=batch_dir,
                     module_list=module_list,
                     slurm_kwargs={"time": "12:00:00",
                                   "nodes": nnodes,
                                   "ntasks-per-node" : nprocesses},
                     vcstools_version=vcstools_version,
                     queue='cpuq')

        if write_showspec:
            # Now write the showspec batch for this time
            write_showspec_batch(times[i], obsid, ra, dec, freq, (90/thetares)+1, 360/phires, onamebase+".dat", maploc)



def write_showspec_batch(time, obsid, ra, dec, freq, ntheta, nphi, infile, maploc):
    # Nick: I'm not sure what this is or what it's for so I'm not going to touch it
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

    beam_models = ['analytic', 'advanced', 'full_EE', 'hyperbeam']

    # Argument parsing
    parser = argparse.ArgumentParser(description="Simple script to help write the batch scripts required for running the tied-array beam simulations over multiple epochs")

    required_options = parser.add_argument_group('Required Options')
    required_options.add_argument("-o", "--obsid", type=int, help="Observation ID")
    required_options.add_argument("-b", "--begin", type=int, help="GPS time for first evaluation")
    required_options.add_argument("-e", "--end", type=int, help="GPS time for last evaluation")
    required_options.add_argument("-a", "--all", action="store_true", default=False, help="Perform on entire observation span. Use instead of -b & -e. ")
    required_options.add_argument("--ra", type=str, help="RAJ2000 of target in hh:mm:ss.ss (use = to assign option)")
    required_options.add_argument("--dec", type=str, help="DECJ2000 of target in dd:mm:ss.ss (use = to assign option)") # only because argparse can't handle negative arguments...

    optional_options = parser.add_argument_group('Optional Options')
    optional_options.add_argument('-O', '--cal_obs', type=int, default=None,
                                  help="Observation ID of calibrator you want to process. "
                                       "Only required in to work out the default location of flagged_tiles. [no default]")
    optional_options.add_argument("--freq", type=float, help="Observing frequency in Hz [default: The observations centre frequency]", default=None)
    optional_options.add_argument("--flagged_tiles", type=str, default=None,
                                  help="Path (including file name) to file containing the flagged tiles as used in the RTS, will be used by get_delays. ")
    optional_options.add_argument("--beam_model", type=str, default='hyperbeam',
                                  help='Decides the beam approximation that will be used. Options: "analytic" the analytic beam model (2012 model, fast and reasonably accurate), "advanced" the advanced beam model (2014 model, fast and slighty more accurate) or "full_EE" the full EE model (2016 model, slow but accurate). " Default: "analytic"')
    optional_options.add_argument("--step", type=int, help="Time step between each evaluation (in seconds) [default: 500]", default=500)
    optional_options.add_argument("--thetares", type=float, help="Resolution of theta (zenith angle) grid, in degrees [default: 0.05]", default=0.05)
    optional_options.add_argument("--phires", type=float, help="Resolution of phi (azimuth) grid, in degrees [default: 0.05]", default=0.05)
    optional_options.add_argument("--nodes", type=int, help="Number of nodes to use per run (more than 10 is overkill and will actually be slower...) [default: 1]", default=1)
    optional_options.add_argument("--eff", type=float, help="Radiation efficiency [default: 1.0]", default=1.0)
    optional_options.add_argument("--maploc", type=str, default="$PWD",
                                  help="Path to directory where skymaps/ exists (or is to be created). The GSM temperature maps exist there or will be created.")
    optional_options.add_argument("--dont_write", action='store_true', help="Don't write output files to disk")
    optional_options.add_argument("--write_showspec", action='store_true', help="Write out the showspec batch script")
    optional_options.add_argument("--odir", type=str, default=None, help="Output directory")
    optional_options.add_argument("--vcstools_version", type=str, default="master", help="VCSTools version to load in jobs (i.e. on the queues) ")
    optional_options.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO", choices=loglevels.keys(), default="INFO")

    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # get metadata
    common_metadata = get_common_obs_metadata(args.obsid)

    # Option parsing
    if not args.obsid:
        logger.error("Observation ID required, please put in with -o or --obsid")
        sys.exit(0)
    if args.all and (args.begin or args.end):
        logger.error("Please specify EITHER (-b,-e) OR -a")
        sys.exit(0)
    elif args.all:
        args.begin, args.end = obs_max_min(args.obsid)
    if args.begin and args.end:
        if args.begin > args.end:
            logger.error("Starting time is after end time")
            sys.exit(0)
    if args.beam_model not in beam_models:
        logger.error("Unknown beam model. Please use one of {0}. Exiting.".format(beam_models))
        quit()

    # Default parsing
    if args.freq is None:
        args.freq = common_metadata[5] * 1e6 #Hz
        logger.info("Using the observations centre frequency: {} MHz".format(args.freq / 1e6))

    # Flagged tile parsing
    comp_config = load_config_file()
    if args.flagged_tiles:
        flagged_tiles_file = os.path.abspath(args.flagged_tiles)
        if not os.path.isfile(args.flagged_tiles):
            logger.error("Your are not pointing at a file with your input "
                            "to --flagged_tiles. Aborting here as the "
                            "beamformer will not run...")
            sys.exit(0)
    elif args.cal_obs:
        DI_dir = os.path.join(comp_config['base_data_dir'], str(args.obsid), 'cal', str(args.cal_obs), 'rts')
        if os.path.isfile("{0}/flagged_tiles.txt".format(DI_dir)):
            flagged_tiles_file = "{0}/flagged_tiles.txt".format(DI_dir)
            logger.info("Found tiles flags in {0}/flagged_tiles.txt. "
                        "Using it by default".format(DI_dir))
        else:
            logger.error("No file in {0}/flagged_tiles.txt. Use --flagged_tiles")
            sys.exit(0)
    else:
        logger.error("Please either use --cal_id to work out the default location or "
                     "--flagged_tiles for non-standard locations")
        sys.exit(0)
    # Converting the file into a list
    with open(flagged_tiles_file) as f:
        flagged_tiles= f.readlines()
    flagged_tiles = [x.strip() for x in flagged_tiles]

    # Ensure metafits file
    metafits_file = "{0}_metafits_ppds.fits".format(args.obsid)
    data_dir = os.path.join(comp_config['base_data_dir'], str(args.obsid))
    metafits_file_loc = os.path.join(data_dir, metafits_file)
    ensure_metafits(data_dir, args.obsid, metafits_file)

    write_batch_files(str(args.obsid), args.begin, args.end,
                      args.ra, args.dec, args.freq, flagged_tiles,
                      step=args.step, thetares=args.thetares, phires=args.phires,
                      nnodes=args.nodes, eff=args.eff, beam_model=args.beam_model,
                      maploc=args.maploc, odir=args.odir, delays=common_metadata[4][0],
                      write=(not args.dont_write), write_showspec=args.write_showspec,
                      vcstools_version=args.vcstools_version, metafits_loc=metafits_file_loc)
