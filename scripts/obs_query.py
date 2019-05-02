#!/usr/bin/env python

import logging
import argparse

logger = logging.getLogger(__name__)
from mwa_metadb_utils import getmeta 

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def print_info(obs_id):
    """
    Snippet function to show how to get the first and last files for a given observation
    """

    obsinfo = getmeta(service='obs', params={'obs_id':str(obs_id)})

    logger.info("Obs ID: {0}".format(obs_id))
    logger.info("Observation Name: {0}".format(obsinfo['obsname']))
    logger.info("Channels: {0}".format(obsinfo['rfstreams']['0']['frequencies']))
    logger.info("Duration: {0}".format(obsinfo['stoptime'] - obsinfo['starttime'], "seconds"))
   
if __name__ == '__main__':
    from sys import argv
    
    #dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                    INFO=logging.INFO,
                    WARNING=logging.WARNING)
    
    #Arguments
    parser = argparse.ArgumentParser(description="""Returns information on a given OBS ID""")
    parser.add_argument('-o','--obsid',type=str, help='Input OBS IDs in the format " -o 1099414416"')
    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO",
                                    choices=loglevels.keys(), default="INFO")
    parser.add_argument("-V", "--version", action="store_true", help="Print version and quit")
    args = parser.parse_args()
    

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)


    if args.version:
        import sys
        try:
            import version
            print(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            print("Couldn't import version.py - have you. installed vcstools?")
            print("ImportError: {0}".format(ie))
            sys.exit(0)


    print_info(args.obsid)

