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
    parser.add_argument('-c', '--command', type=str,
            help="mwa_search_pipeline.nf command")
    parser.add_argument('-s', '--super_computer_id', type=int, default=1,
            help="Super computer ID. Ozstar: 1. Garrawarla: 2. SHAO: 3. Galaxy: 4. Magnus: 5. Default 1.")
    args = parser.parse_args()

    print("Uploading beams")
    pointing_list = np.loadtxt(args.file, dtype=str)
    if pointing_list.ndim == 0:
        pointing_list.shape = (1,)
    upload_beam(pointing_list, args.obsid, args.begin, args.end,
                mwa_search_command=args.command, mwa_search_version='v3.0', supercomputer_id=args.super_computer_id)

