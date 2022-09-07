#!/usr/bin/env python

import argparse
import numpy as np

from vcstools.client import upload_beam


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Upload rows to the Beams table in the SMART webapp""")
    parser.add_argument('file', type=str, help="Pointing grid file")
    parser.add_argument('-o', '--obsid', type=int,
            help="The observation ID of the MWA observation")
    parser.add_argument('-b', '--begin', type=int,
            help="Begin time gps")
    parser.add_argument('-e', '--end', type=int,
            help="End time gps")
    parser.add_argument('-c', '--command', type=str, default="",
            help="mwa_search_pipeline.nf command")
    parser.add_argument('-s', '--super_computer_id', type=int, default=1,
            help="Super computer ID. 1=Ozstar, 2=Garrawarla, 3=SHAO, 4=Galaxy, 5=Magnus (default=1)")
    parser.add_argument('-v', '--mwa_search_version', type=str, default="v3.0",
            help="The version of mwa_search used to create these beams")
    args = parser.parse_args()

    print("Uploading beams")
    pointing_list = np.loadtxt(args.file, dtype=str)
    if pointing_list.ndim == 0:
        pointing_list.shape = (1,)

    upload_beam(
            pointing_list,
            args.obsid,
            args.begin,
            args.end,
            mwa_search_command=args.command,
            mwa_search_version=args.mwa_search_version,
            supercomputer_id=args.super_computer_id
            )

