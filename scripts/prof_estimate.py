#!/usr/bin/env python3

import sys
import os
import logging
import argparse

from vcstools.prof_utils import get_from_bestprof, get_from_ascii, subprocess_pdv
from vcstools.gfit import gfit
from vcstools.general_utils import setup_logger

logger = logging.getLogger(__name__)

if __name__ == '__main__':

    loglevels = dict(DEBUG=logging.DEBUG,\
                    INFO=logging.INFO,\
                    WARNING=logging.WARNING,\
                    ERROR=logging.ERROR)

    parser = argparse.ArgumentParser(description="""A utility file for calculating a variety of pulse profile properties""",\
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    file_parse = parser.add_argument_group("File Input", description="Choose only one option")
    file_parse.add_argument("--bestprof", type=str, help="The pathname of the file containing the pulse profile. Use if in .pfd.bestprof format")
    file_parse.add_argument("--ascii", type=str, help="The pathname of the file containing the pulse profile. Use if in ascii text format")
    file_parse.add_argument("--archive", type=str, help="The pathname of the file containing the pulse profile. Use if in archive format")

    g_inputs = parser.add_argument_group("Gaussian Inputs")
    g_inputs.add_argument("--max_N", type=int, default=6, help="The maximum number of gaussian components to attempt to fit")
    g_inputs.add_argument("--min_comp_len", type=int, default=None,\
                        help="Minimum length of a component to be considered real. Measured in bins. If none, will use 1 percent of total profile length")
    g_inputs.add_argument("--clip_type", type=str, default="regular", choices=["regular", "noisy", "verbose"], help="The verbosity of clipping used\
                        by the Gaussian fitter.")

    other_inputs = parser.add_argument_group("Other Inputs")
    other_inputs.add_argument("--plot_name", type=str, help="The name of the output plot file. If none, will not plot anything")
    other_inputs.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbostity level")
    args = parser.parse_args()

    # set up the logger for stand-alone execution
    logger = setup_logger(logger, log_level=loglevels[args.loglvl])

    if args.bestprof:
        _, _, _, _, _, _, _, _, profile, _, _  = get_from_bestprof(args.bestprof)
    elif args.ascii:
        profile = get_from_ascii(args.ascii)[0]
    elif args.archive:
        subprocess_pdv(args.archive, outfile="archive.txt")
        profile = get_from_ascii("archive.txt")[0]
        os.remove("archive.txt")
    else:
        logger.error("Please supply either an ascii or bestprof profile")
        sys.exit(1)

    g_fitter = gfit(profile, max_N=args.max_N, min_comp_len=args.min_comp_len,
                    plot_name=args.plot_name, clip_type=args.clip_type)
    g_fitter.auto_gfit()
    if args.plot_name:
        g_fitter.plot_fit()