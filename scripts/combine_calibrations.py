#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import filecmp
from shutil import copyfile
import numpy as np

from vcstools.config import load_config_file
from vcstools.general_utils import mdir, setup_logger
from vcstools.metadb_utils import get_common_obs_metadata

logger = logging.getLogger(__name__)

def order_chans(channels):
    # reogranises the channels to take into account the flip above channel id 128
    channels = np.array(channels, dtype=np.int)
    hichans = [c for c in channels if c>128]
    lochans = [c for c in channels if c<=128]
    lochans.extend(list(reversed(hichans)))
    ordered_channels = lochans
    return ordered_channels

if __name__ == '__main__':
    # Dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING)

    try:
        comp_config = load_config_file()
    except Exception:
        # No computer found so making a default for the argparse help to work
        comp_config = {'base_data_dir' : f"/astro/mwavcs/{os.environ['USER']}/"}
    parser = argparse.ArgumentParser(description="If there is no calibration observation with the correct channels, " +\
                                     "use this script to combine the result of two calibration solutions to make a new " +\
                                     "calibration solution with all the required frequency channels.")
    parser.add_argument("-o", "--obsid", type=int,
                        help="Observation ID you want to process")
    parser.add_argument("-c1", "--cal1", type=int,
                        help="The first calibration observation ID you want to combine")
    parser.add_argument("-c2", "--cal2", type=int,
                        help="The second calibration observation ID you want to combine")
    parser.add_argument('--out_dir', type=str,\
                        help="The output directory of the combined solution. " +\
                             "Default is {0}[obsid]/cal/[cal1]_[cal2]/rts".format(comp_config['base_data_dir']))
    parser.add_argument("-V", "--version", action="store_true", help="Print version and quit")
    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO",
                                    choices=loglevels.keys(), default="INFO")
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

    if args.out_dir:
        out_dir = args.out_dir
    else:
        out_dir = os.path.join(comp_config['base_data_dir'], str(args.obsid), "cal", f"{args.cal1}_{args.cal2}", "rts")
    cal_base_dir = os.path.join(comp_config['base_data_dir'], str(args.obsid), "cal")
    logger.debug(cal_base_dir)

    obs_chans = get_common_obs_metadata(args.obsid)[6]
    obs_chans_reordered = order_chans(obs_chans)
    cal1_chans = get_common_obs_metadata(args.cal1)[6]
    cal1_chans_reordered = order_chans(cal1_chans)
    cal2_chans = get_common_obs_metadata(args.cal2)[6]
    cal2_chans_reordered = order_chans(cal2_chans)

    logger.debug(obs_chans)
    logger.debug(obs_chans_reordered)
    logger.debug(cal1_chans)
    logger.debug(cal1_chans_reordered)
    logger.debug(cal2_chans)
    logger.debug(cal2_chans_reordered)


    # Check input calbration obs have all the frequency channels we need
    for fc in obs_chans:
        if not (fc in cal1_chans or fc in cal2_chans):
            logger.error(f"Frequency channel {fc} is not found in either of the calibration obs. Exiting")
            sys.exit(1)

    # Make output dir
    mdir(out_dir, 'Combined Calibration', gid=comp_config['gid'])

    # Check flagged*.txt files are the same and move them to new directory
    if not filecmp.cmp(f"{cal_base_dir}/{args.cal1}/rts/flagged_tiles.txt",
                       f"{cal_base_dir}/{args.cal1}/rts/flagged_tiles.txt"):
        logger.warn(f"{cal_base_dir}/{args.cal1}/rts/flagged_tiles.txt and {cal_base_dir}/{args.cal1}"+\
                    "/rts/flagged_tiles.txt are different. It is safer to recalibrate with identical flags")
    copyfile(f"{cal_base_dir}/{args.cal1}/rts/flagged_tiles.txt", f"{out_dir}/flagged_tiles.txt")
        
    if not filecmp.cmp(f"{cal_base_dir}/{args.cal1}/rts/flagged_channels.txt",
                       f"{cal_base_dir}/{args.cal1}/rts/flagged_channels.txt"):
        logger.warn(f"{cal_base_dir}/{args.cal1}/rts/flagged_channels.txt and {cal_base_dir}/{args.cal1}"+\
                    "/rts/flagged_channels.txt are different. It is safer to recalibrate with identical flags")
    copyfile(f"{cal_base_dir}/{args.cal1}/rts/flagged_channels.txt", f"{out_dir}/flagged_channels.txt")


    # Loop over observation channels and work out the matching calibration solution IDs so we can rename the files correctly
    # Check first calibrator first
    logger.info(f"Using the following files from {args.cal1}")
    for ci, cal_chan in enumerate(cal1_chans_reordered, 1):
        for oi, obs_chan in enumerate(obs_chans_reordered, 1):
            if obs_chan == cal_chan:
                logger.info(f"node{ci:03d} is becoming node{oi:03d}")
                copyfile(f"{cal_base_dir}/{args.cal1}/rts/DI_JonesMatrices_node{ci:03d}.dat",
                         f"{out_dir}/DI_JonesMatrices_node{oi:03d}.dat")
                copyfile(f"{cal_base_dir}/{args.cal1}/rts/BandpassCalibration_node{ci:03d}.dat",
                         f"{out_dir}/BandpassCalibration_node{oi:03d}.dat")
    # Then second calibrator
    logger.info(f"Using the following files from {args.cal2}")
    for ci, cal_chan in enumerate(cal2_chans_reordered, 1):
        for oi, obs_chan in enumerate(obs_chans_reordered, 1):
            if obs_chan == cal_chan:
                logger.info(f"node{ci:03d} is becoming node{oi:03d}")
                copyfile(f"{cal_base_dir}/{args.cal2}/rts/DI_JonesMatrices_node{ci:03d}.dat",
                         f"{out_dir}/DI_JonesMatrices_node{oi:03d}.dat")
                copyfile(f"{cal_base_dir}/{args.cal2}/rts/BandpassCalibration_node{ci:03d}.dat",
                         f"{out_dir}/BandpassCalibration_node{oi:03d}.dat")

    

    

