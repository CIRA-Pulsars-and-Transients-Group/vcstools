#!/usr/bin/env python

import argparse
import glob
import concurrent.futures
import sys
import numpy as np
import requests

from vcstools.client import upload_beam, upload_obsid, upload_cand
from vcstools.progress_bar import progress_bar
from vcstools.pointing_utils import sex2deg


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Upload the candidates on prometheus""")
    parser.add_argument('-o', '--obsid', type=int,
            help='The observation ID of the MWA observation')
    parser.add_argument('-b', '--begin', type=int,
            help="Begin time gps")
    parser.add_argument('-e', '--end', type=int,
            help="End time gps")
    parser.add_argument('-f', '--file', type=str,
            help="Pointing grid file")
    parser.add_argument('-c', '--command', type=str,
            help="mwa_search_pipeline.nf command")
    parser.add_argument('-m', '--max_threads', type=int, default=8,
            help="Maximum number of threads to use for concurrent submission.")
    parser.add_argument('-s', '--super_computer_id', type=int, default=1,
            help="Super computer ID. Ozstar: 1. Garrawarla: 2. SHAO: 3. Galaxy: 4. Magnus: 5. Default 1.")
    args = parser.parse_args()

    print("Uploading obsid")
    try:
        upload_obsid(args.obsid)
    except requests.exceptions.HTTPError:
        print("Obsid already uploaded")

    if args.file:
        print("Uploading beams")
        pointing_list = np.loadtxt(args.file, dtype=str)
        upload_beam(pointing_list, args.obsid, args.begin, args.end,
                    mwa_search_command=args.command, mwa_search_version='v3.0', supercomputer_id=args.super_computer_id)
    else:
        print("WARNING: Not uploading beams")

    print("Uploading cands")
    pfd_files = glob.glob('B*pfd')
    upload_cand(pfd_files, obsid=args.obsid,
                search_type='SMART_10min', search_params_id=1)

