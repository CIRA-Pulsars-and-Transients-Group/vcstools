#!/usr/bin/env python3

import sys
import os
import logging
import argparse

from vcstools.prof_utils import get_from_bestprof, get_from_ascii, auto_gfit, prof_eval_gfit, subprocess_pdv

logger = logging.getLogger(__name__)

if __name__ == '__main__':

    loglevels = dict(DEBUG=logging.DEBUG,\
                    INFO=logging.INFO,\
                    WARNING=logging.WARNING,\
                    ERROR=logging.ERROR)

    parser = argparse.ArgumentParser(description="""A utility file for calculating a variety of pulse profile properties""",\
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    inputs = parser.add_argument_group("Inputs")
    inputs.add_argument("--bestprof", type=str, help="The pathname of the file containing the pulse profile. Use if in .pfd.bestprof format")
    inputs.add_argument("--ascii", type=str, help="The pathname of the file containing the pulse profile. Use if in ascii text format")
    inputs.add_argument("--archive", type=str, help="The pathname of the file containing the pulse profile. Use if in archive format")
    inputs.add_argument("--ignore_threshold", type=float, default=0.02, help="Maxima with values below this fraction of the profile maximum will be ignored.")
    inputs.add_argument("--min_comp_len", type=int, default=None,\
                        help="Minimum length of a component to be considered real. Measured in bins. If none, will use 1 percent of total profile length")
    inputs.add_argument("--alpha", type=float, default=2, help="Used by the clipping function to determine the noise level. A lower value indicates\
                        a higher verbosity level in the noise clipping function.")
    inputs.add_argument("--auto", action="store_true", help="Used to automatically find the best alpha value to clip this profile")
    inputs.add_argument("--period", type=float, help="The period of the puslar in ms. Found automatically if bestprof supplied.\
                        Used in S/N calculation. Not required")

    g_inputs = parser.add_argument_group("Gaussian Inputs")
    g_inputs.add_argument("--max_N", type=int, default=6, help="The maximum number of gaussian components to attempt to fit")

    other_inputs = parser.add_argument_group("Other Inputs")
    other_inputs.add_argument("--plot_name", type=str, help="The name of the output plot file. If none, will not plot anything")
    other_inputs.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbostity level")
    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if args.bestprof:
        _, _, _, period, _, _, _, profile, _  = get_from_bestprof(args.bestprof)
    elif args.ascii:
        profile = get_from_ascii(args.ascii)[0]
        period = args.period
    elif args.archive:
        subprocess_pdv(args.archive, outfile="archive.txt")
        profile = get_from_ascii("archive.txt")[0]
        os.remove("archive.txt")
        period = args.period
    else:
        logger.error("Please supply either an ascii or bestprof profile")
        sys.exit(1)

    if args.auto:
        auto_gfit(profile, max_N=args.max_N, ignore_threshold=args.ignore_threshold,\
                        plot_name=args.plot_name, min_comp_len=args.min_comp_len, period=period)
    else :
        prof_eval_gfit(profile, max_N=args.max_N, ignore_threshold=args.ignore_threshold,\
                        plot_name=args.plot_name, min_comp_len=args.min_comp_len, alpha=args.alpha, period=period)