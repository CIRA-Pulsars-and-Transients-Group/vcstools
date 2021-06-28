#!/usr/bin/env python3

import subprocess
import sys
import argparse
import glob
import re
import logging
logger = logging.getLogger(__name__)

from vcstools.config import load_config_file
from vcstools.general_utils import setup_logger

def opt_parser():
    # Dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING)
    beam_models = ['analytic', 'advanced', 'full_EE']
    parser = argparse.ArgumentParser(description="Tool for cleaning up after processing MWA Pulsar Data on Galaxy")
    parser.add_argument("-o", "--obs", type=int, default=None, help="Observation ID you want to process [no default]")
    parser.add_argument("--raw", action="store_true", default=False,
                        help="Add this option if you want to remove the raw and/or recombined files for an observation")
    parser.add_argument("--beamformed", action="store_true", default=False,
                        help="Add this option if you want to remove the recombined files for an observation")
    parser.add_argument("-V", "--version", action="store_true", help="Print version and quit")
    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO", 
                                    choices=loglevels.keys(), default="INFO")
    return parser.parse_args(), loglevels

def munlink_files(folder, file_type):
    while True:
        question = 'Preparing to delete {0} files from {1}. Proceed [y/n]? '.format(file_type, folder)
        authority = raw_input(question)
        if (authority == "Y") or (authority =="y"):
            for file in glob.iglob(folder+"/*"):
                if re.search(file_type,file):
                    command = "munlink {0}".format(file)
                    subprocess.call(command, shell=True)
            break
        elif (authority == "N") or (authority ==  "n"):
            break
        else:
            logger.error("Unrecognized input option, please try again.")


def remove_raw(obs):
    comp_config = load_config_file()
    raw_folder = os.path.join(comp_config['base_data_dir'], str(obs), "raw")
    combined_folder = os.path.join(comp_config['base_data_dir'], str(obs), "combined")
    
    raw_files = False
    tar_files = False
    combined_files = False
    ics_files = False
    for file in glob.iglob("{0}/*".format(raw_folder)):
        if re.search('dat', file):
            raw_files = True
    for file in glob.iglob("{0}/*".format(combined_folder)):
        if re.search('tar', file):
            tar_files = True
        if re.search('ch\d{3}', file):
            combined_files = True
        if re.search('ics', file):
            ics_files = True
    
    if raw_files:
        munlink_files(raw_folder, "dat")
    if tar_files:
        munlink_files(combined_folder, "tar")
    if combined_files:
        munlink_files(combined_folder, "ch\d{3}")
    if ics_files:
        munlink_files(combined_folder, "ics")
        
def remove_beamformed(obs, pointing=None):
    comp_config = load_config_file()
    pointing_folder = os.path.join(comp_config['base_data_dir'], str(obs), "pointings")
    if not pointing:
        authority = ('No pointing specified, would you like to remove all pointings for this observation?')
        if (authority == "Y") or (authority == "y"):
            pointings = glob.glob("{0}/*:*:*:*:*")
            if not pointings:
                logger.error("No valid pointings in {0}. Exiting...")
                sys.exit(0)
            else:
                for pointing in pointings:
                    logger.info("Checking if pointing {0} has been uploaded to MWA Pulsar Database...".format(pointing))
                    
            # Upload to MWA Pulsar Database if not there already
            # Remove each pointing
    return




if __name__ == '__main__':
    args, loglevels = opt_parser()
    if args.version:
        try:
            import version
            logger.info(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            logger.error("Couldn't import version.py - have you installed vcstools?")
            logger.error("ImportError: {0}".format(ie))
            sys.exit(0)

    # set up the logger for stand-alone execution
    logger = setup_logger(logger, log_level=loglevels[args.loglvl])

    if args.raw:
        remove_raw(args.obs)
    if args.beamformed:
        logger.error("Sorry, this option is not implemented yet :(")




