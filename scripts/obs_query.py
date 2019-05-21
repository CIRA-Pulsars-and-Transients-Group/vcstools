#!/usr/bin/env python3

import logging
import argparse

logger = logging.getLogger(__name__)
from mwa_metadb_utils import getmeta, is_number 


def print_info(obs_id):
    """
    Snippet function to show how to get the first and last files for a given observation
    """

    obsinfo = getmeta(service='obs', params={'obs_id':str(obs_id)})
    
    logger.info("Obs ID: {0}".format(obs_id))
    #See if OBS ID is valid. If not, exit
    try:
        logger.info("Observation Name: {0}".format(obsinfo['obsname']))
    except:
        import sys
        logger.error("Could not get obs name - is the obs id correct?: {0}".format(obs_id))
        sys.exit(0)   
    logger.info("Channels: {0}".format(obsinfo['rfstreams']['0']['frequencies']))
    logger.info("Start Time: {0}".format(obsinfo['starttime']))
    logger.info("Stop Time: {0}".format(obsinfo['stoptime']))
    logger.info("Duration: {0}".format(obsinfo['stoptime'] - obsinfo['starttime'], "seconds"))
    
   
if __name__ == '__main__':
    from sys import argv
    
    #dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                    INFO=logging.INFO,
                    WARNING=logging.WARNING,
                    ERROR = logging.ERROR)
    
    #Arguments
    parser = argparse.ArgumentParser(description="""Returns information on a given OBS ID""")
    parser.add_argument("obsid", type=int, help="Input Observation ID")
    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO",
                                    choices=loglevels.keys(), default="INFO")
    parser.add_argument("-V", "--version", action="store_true", help="Print version and quit. Currently this requires an obsid to work. Any number will suffice.")
    #TODO: Fix the above line so that obsid is not required for the user to see the version  

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
            logger.info("Version: {0}".format(version.__version__))
            sys.exit(0)
        except ImportError as ie:
            logger.error("Couldn't import version.py - have you. installed vcstools?")
            logger.error(("ImportError: {0}".format(ie)))

            sys.exit(0)


    print_info(args.obsid)

