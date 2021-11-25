#!/usr/bin/env python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import os
import sys
import argparse
import logging

from vcstools.config import load_config_file
from vcstools.check_files import check_download, check_recombine
from vcstools.general_utils import setup_logger

logger = logging.getLogger(__name__)


def opt_parser(loglevels):
    try:
        comp_config = load_config_file()
    except Exception:
        # No computer found so making a default for the argparse help to work
        comp_config = {'base_data_dir' : "/astro/mwavcs/vcs/"}
    parser = argparse.ArgumentParser(description="scripts to check sanity of downloads and recombine.")
    parser.add_argument("-m", "--mode", type=str, choices=['download','recombine'],\
                          help="Mode you want to run: download, recombine", dest='mode', default=None)
    parser.add_argument("-d", "--data_type", type=str, choices=['11','15','16', 'raw','ics','tar_ics'],\
                          help="Only necessary when checking downloads. Types refer to those as definded " + \
                            "in voltdownload.py: 11 = Raw, 15 = ICS only, 16 = ICS and tarballs of recombined data.", \
                            dest='data_type', default=None)
    parser.add_argument("-o", "--obs", metavar="OBS ID", type=int, dest='obsID',\
                            help="Observation ID you want to process [no default]", default=None)
    parser.add_argument("-b", "--begin", metavar="start", type=int, dest='begin',\
                            help="gps time of first file to check on [default=%(default)s]",\
                            default=None)
    parser.add_argument("-e", "--end", metavar="stop", type=int, dest='end',\
                            help="gps time of last file to check on [default=%(default)s]",\
                            default=None)
    parser.add_argument("-a", "--all", action="store_true", default=False, help="Perform on entire observation span. Use instead of -b & -e. [default=%(default)s]")
    parser.add_argument("-i", "--increment", metavar="time increment", type=int, \
                            dest='increment',\
                            help="Effectively the number of seconds to ckeck for " +\
                            "starting at start time [default=%(default)s]",\
                            default=None)
    parser.add_argument("-s", "--size", type=int, dest='size',\
                          help="The files size in bytes that you expect all files" +\
                          " to have. Per default will figure this out from files on he archive" +\
                            "We expect 253440000 (download raw), 327680000" +\
                          " (recombined, not ics), 7865368576 (tarballs)", default=None)
    parser.add_argument("-S", "--size_ics", type=int, help='Size in bytes that' +\
                            "you expect the ics files to have. Default = %(default)s", \
                            dest='size_ics', default=30720000)
    parser.add_argument('-w', '--work_dir', type=str, dest='work_dir',\
                            help="Directory to check the files in. " +\
                                 "Default is {0}[obsID]/[raw,combined]".format(comp_config['base_data_dir']))
    parser.add_argument("-V", "--version", action="store_true", help="Print version and quit")
    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO",
                                    choices=loglevels.keys(), default="INFO")
    return parser.parse_args()

if __name__ == '__main__':
    # Dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING)

    args = opt_parser(loglevels)
    comp_config = load_config_file()
    work_dir_base = os.path.join(comp_config['base_data_dir'], str(args.obsID))

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

    if (args.mode is None) or (args.obsID is None):
        logger.error("You must specify BOTH a mode and observation ID")
        sys.exit(1)

    if args.all:
        from vcstools.metadb_utils import obs_max_min
        args.begin, args.end = obs_max_min(args.obsID)
    if args.end:
        if not args.begin:
            logger.error("If you supply and end time you also *have* to supply a begin time.")
            sys.exit(1)
        args.increment = args.end - args.begin + 1
        logger.info("Checking {0} seconds.".format(args.increment))
    if not args.all and not args.begin:
        logger.error("You have to either set the -a flag to process the whole obs or povide a start and stop time with -b and -e")
        sys.exit(1)
    if args.begin:
        if not args.end and not args.increment:
            logger.error("If you specify a begin time you also have to provide either an end time (-e end) or the number of seconds to check via the increment flag (-i increment)")
            sys.exit(1)
    if args.mode == 'download':
        if not args.data_type:
            logger.error("In download mode you need to specify the data type you downloaded.")
            sys.exit(1)
        if args.data_type == '11':
            data_type = 'raw'
        elif args.data_type == '15':
            data_type = 'ics'
        elif args.data_type == '16':
            data_type = 'tar_ics'
        else:
            data_type = args.data_type
        if data_type == 'raw':
            work_dir = work_dir_base + '/raw/'
        else:
            work_dir = work_dir_base + '/combined/'
        if args.work_dir:
            work_dir = args.work_dir
        if data_type == 'raw' or data_type == 'tar_ics' or data_type == 'ics':
            sys.exit(check_download(args.obsID, directory=work_dir,
                                    startsec=args.begin, n_secs=args.increment, data_type=data_type))
    elif args.mode == 'recombine':
        required_size = 327680000
        if args.size:
            required_size = args.size
        required_size_ics = args.size_ics
        work_dir = work_dir_base + '/combined/'
        if args.work_dir:
            work_dir = args.work_dir
        sys.exit(check_recombine(args.obsID, directory=work_dir, required_size=required_size, \
                            required_size_ics=required_size_ics, startsec=args.begin, n_secs=args.increment))
    else:
        logger.error("No idea what you want to do. This mode is not supported. Ckeck the help.")
        sys.exit(1)
