#!/usr/bin/env python

import argparse
import glob
from vcstools.client import upload_cand


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Upload candidates to the SMART webapp""")
    parser.add_argument('-o', '--obsid', type=int,
            help="The observation ID of the MWA observation")
    parser.add_argument('-a', '--all', action='store_true',
            help="Load all PFD files in the current directory (overrides -p)")
    parser.add_argument('-p', '--pfd_file', type=str,
            help="The PFD file to be uploaded. This can be overridden by -a")
    args = parser.parse_args()

    if args.pfd_file is None and args.all == False:
        raise RuntimeError("Either -p or -a must be supplied")

    print("Uploading cands")
    if args.all:
        pfd_files = glob.glob('*pfd')
    else:
        pfd_files = [args.pfd_file]

    upload_cand(pfd_files, obsid=args.obsid,
                search_type='SMART_10min', search_params_id=1)

