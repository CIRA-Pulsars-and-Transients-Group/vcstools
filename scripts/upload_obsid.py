#!/usr/bin/env python

import argparse
import glob
import concurrent.futures
import sys
import numpy as np
import requests

from vcstools.client import upload_obsid
from vcstools.progress_bar import progress_bar
from vcstools.pointing_utils import sex2deg


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Upload the candidates on prometheus""")
    parser.add_argument('-o', '--obsid', type=int,
            help='The observation ID of the MWA observation')
    args = parser.parse_args()

    print("Uploading obsid")
    try:
        upload_obsid(args.obsid)
    except requests.exceptions.HTTPError:
        print("Obsid already uploaded")

