#!/usr/bin/env python3

import subprocess
import os
import sys
import glob
from mpi4py import MPI
import logging
import argparse

from vcstools.general_utils import setup_logger

logger = logging.getLogger(__name__)

if __name__ == '__main__':
    loglevels = dict(DEBUG=logging.DEBUG,
                    INFO=logging.INFO,
                    WARNING=logging.WARNING,
                    ERROR = logging.ERROR)

    #Arguments
    parser = argparse.ArgumentParser(description="""Combines raw VCS data""")
    parser.add_argument("-o", "--obsid", type=int, help="Input Observation ID.", required=True)
    parser.add_argument("-s", "--start", type=int, help="First GPS second of data to recombine.", required=True)
    parser.add_argument("-d", "--duration", type=int, help="Seconds of data to recomine", required=True)
    parser.add_argument("-w", "--data_dir", type=str, help="Directory containing the raw data")
    parser.add_argument("-p", "--output_base_dir", type=str, help="The base directory to put the recombined data. Eg /astro/mwavcs/vcs/<obsid>.")
    parser.add_argument("-e", "--recombine_command", type=str, help="the filename of the recombine function. Default: recombine", default="recombine")
    parser.add_argument("-L", "--loglvl", type=str, choices=loglevels.keys(), default="INFO", help="Logger verbosity level. Default: INFO")
    parser.add_argument("-V", "--version", action="store_true", help="Print version and quit")
    args = parser.parse_args()

    # set up the logger for stand-alone execution
    logger = setup_logger(logger, log_level=loglevels[args.loglvl])

    if args.version:
        try:
            import version
            logger.info(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            logger.error("Couldn't import version.py - have you installed vcstools?")
            logger.error("ImportError: {0}".format(ie))
            sys.exit(0)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if size != args.duration:
        logger.error("Wrong number of CPU cores (cores != --duration). Exiting.")
        sys.exit(1)

    # Core dependent time to recombine
    time_to_combine = args.start + rank

    # Check if files are already succesfully combined
    combined_files_glob = f"{args.data_dir}/combined/{args.obsid}_{time_to_combine}_ch*.dat"
    # Check if
    left_to_combined = 24
    for to_check in sorted(glob.glob(combined_files_glob)):
        file_statinfo = os.stat(to_check)

        if (file_statinfo.st_size == 327680000):
            # File with correct size
            left_to_combined -= 1
        else:
            # File with incorrect size so delete it
            os.remove(to_check)

    logger.debug(f"Thread {rank} :: Final left to combined check: {left_to_combined}")

    if (left_to_combined > 0):
        logger.info(f"Thread {rank} :: Combining raw files")

        # Create file names based on vcs box names
        files_to_combine = []
        for vcs in range(1,17):
            for stream in [0,1]:
                files_to_combine.append(f"{args.obsid}_{time_to_combine}_vcs{vcs:0=2d}_{stream}.dat ")

        # Create the recombine commad
        recombine_line = f"{args.recombine_command} -o {args.obsid} -t {time_to_combine} -m {args.output_base_dir}/{args.obsid}_metafits_ppds.fits -i {args.output_base_dir}/combined -f "
        for raw_file in files_to_combine:
            # Add each file
            recombine_line += f"{args.data_dir}/{raw_file}"

        subprocess.call(recombine_line, shell=True)
    else:
        logger.warn("Combined files already present. Exiting.")
    comm.Barrier()





