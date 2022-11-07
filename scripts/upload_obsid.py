#!/usr/bin/env python

import argparse
import requests
from vcstools.client import upload_obsid

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

