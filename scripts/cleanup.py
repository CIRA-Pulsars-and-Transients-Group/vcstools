#!/usr/bin/env python

import subprocess
import sys
import argparse
import glob
import re


def opt_parser():
    parser = argparse.ArgumentParser(description="Tool for cleaning up after processing MWA Pulsar Data on Galaxy")
    parser.add_argument("-o", "--obs", type=int, default=None, help="Observation ID you want to process [no default]")
    parser.add_argument("--raw", action="store_true", default=False,
                        help="Add this option if you want to remove the raw and/or recombined files for an observation")
    parser.add_argument("--beamformed", action="store_true", default=False,
                        help="Add this option if you want to remove the recombined files for an observation")
    parser.add_argument("-V", "--version", action="store_true", help="Print version and quit")
    
    return parser.parse_args()

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
            print "Unrecognized input option, please try again."


def remove_raw(obs):
    raw_folder = "/astro/mwaops/vcs/{0}/raw".format(obs) #This is a terrible thing to do and needs to be replaced with a config file
    combined_folder = "/astro/mwaops/vcs/{0}/combined".format(obs) #This is a terrible thing to do and needs to be replaced with a config file
    
    raw_files = False
    tar_files = False
    combined_files = False
    ics_files = False
    for file in glob.iglob("{0}/*".format(raw_folder)):
        if re.search('dat', file):
            raw_files = True;
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
        
def remove_beamformed(obs,pointing=None):
    pointing_folder = "/group/mwaops/vcs/{0}/pointings".format(obs) # TODO: Replace this with proper config file
    if not pointing:
        authority = ('No pointing specified, would you like to remove all pointings for this observation?')
        if (authority == "Y") or (authority == "y"):
            pointings = glob.glob("{0}/*:*:*:*:*")
            if not pointings:
                print "No valid pointings in {0}. Exiting..."
                sys.exit(0)
            else:
                for pointing in pointings:
                    print "Checking if pointing {0} has been uploaded to MWA Pulsar Database...".format(pointing)
                    
            # Upload to MWA Pulsar Database if not there already
            # Remove each pointing
    return




if __name__ == '__main__':
    args = opt_parser()
    if args.version:
        try:
            import version
            print(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            print("Couldn't import version.py - have you installed vcstools?")
            print("ImportError: {0}".format(ie))
            sys.exit(0)
    
    if args.raw:
        remove_raw(args.obs)
    if args.beamformed:
        print "Sorry, this option is not implemented yet :("




