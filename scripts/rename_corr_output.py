#!/usr/bin/env python3
from astropy.time import Time
from datetime import datetime
import argparse
import os

from vcstools.general_utils import setup_logger

import logging
logger = logging.getLogger(__name__)

def splitTimestring(timestring):
    yr=int(timestring[0:4])
    month=int(timestring[4:6])
    day=int(timestring[6:8])
    hr=int(timestring[8:10])
    minute=int(timestring[10:12])
    sec=int(timestring[12:14])
    return yr,month,day,hr,minute,sec

def computeTimes(startsec, n_secs=100):
    yr,month,day,hr,minute,sec = splitTimestring(startsec)
    startgps=int(Time(datetime(yr,month,day,hr,minute,sec),scale='utc').gps)
    times = Time(range(startgps + 1, startgps + n_secs), scale='utc',format='gps').datetime
    return " ".join([time.strftime(format='%Y%m%d%H%M%S') for time in times])

def get_nsecs(startsec, endsec):
    yr,month,day,hr,minute,sec = splitTimestring(startsec)
    startgps=int(Time(datetime(yr,month,day,hr,minute,sec),scale='utc').gps)
    yr,month,day,hr,minute,sec = splitTimestring(endsec)
    endgps=int(Time(datetime(yr,month,day,hr,minute,sec),scale='utc').gps)
    return endgps - startgps + 1

def rename(obsID, times, workdir=None):
    if not workdir:
        workdir='/astro/mwavcs/vcs/{0}/vis'.format(obsID)
    n_secs = len(times.split())
    cmd="cd {0};k=1;times=({1});for time in \"${{times[@]}}\";do for box in {2}_${{time}}_gpubox*_00.fits; do mv ${{box}} ${{box%00.fits}}$(printf %02d $k).fits;done;let k=$k+1;if [[ $k -gt {3} ]];then break;fi;done".format(workdir, times, obsID, n_secs)
    return os.system(cmd)

def opt_parser():
    parser=argparse.ArgumentParser(description="Script to rename the files that are produced by the offline correlator to the convention assumed by Andr Offringa's tools.")
    parser.add_argument("-o", "--obs", metavar="OBS ID", type=int, dest='obsID',\
                            help="Observation ID you want to process [no default]",\
                            required=True)
    parser.add_argument("-b", "--begin", metavar="start", type=str, dest='begin',\
                            help="Time string of first file that needs renaming. Format=YYYYMMDDHHMMSS. [default=%(default)s]",\
                            required=True, default=None)
    parser.add_argument("-e", "--end", metavar="stop", type=str, dest='end',\
                            help="Time string of last file that needs renaming. Format=YYYYMMDDHHMMSS. [default=%(default)s]",\
                            required=False, default=None)
    parser.add_argument("-n", "--nsec", metavar="n_secs", type=int, \
                            dest='n_secs', default=100,\
                            help="Number of seconds that need renaming. Maximal 100." +\
                            "[default=%(default)s]",\
                            required=False)
    parser.add_argument('-w', '--workdir', type=str, dest='workdir',\
                            help="Directory " + \
                            "that contains the output files of the offline Correlator. Default is /astro/mwavcs/vcs/<obsID/vis.")
    parser.add_argument("-V", "--version", action="store_true", help="Print version and quit")
    return parser.parse_args()

if __name__ == '__main__':
    args = opt_parser()

    # set up the logger for stand-alone execution
    logger = setup_logger(logger, log_level=loglevels[args.loglvl])

    if args.version:
        import sys
        try:
            import version
            logger.info(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            logger.error("Couldn't import version.py - have you installed vcstools?")
            logger.error("ImportError: {0}".format(ie))
            sys.exit(0)

    if args.end:
        n_secs = get_nsecs(args.begin, args.end)
    else:
        n_secs = args.n_secs
    if n_secs > 100:
        logger.warning("n_secs is too large ({0}). Cannot have more than 100 time steps -- reducing n_secs to 100.\n".format(n_secs))
        #print "n_secs is too large ({0}). Cannot have more than 100 time steps -- reducing n_secs to 100.".format(n_secs)
        n_secs = 100
    quit(rename(args.obsID, computeTimes(args.begin, n_secs), workdir=args.workdir))
