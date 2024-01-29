#!/usr/bin/python3
"""
Contains a single function that reads in a list of integers (channel numbers),
and returns a list of the same integers whose order depends on whether or not
any of the channels are above 128.
"""
from astropy.io import fits as pyfits
import argparse
import os
from vcstools.general_utils import sfreq
from vcstools.config import load_config_file

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Print the coarse channel sequence as in gpubox/subband order. Make sure you have *_metafits_ppds.fits")
    parser.add_argument("-o", "--obsid", type=int, help="Observation ID of target", default=None, required=True)
    parser.add_argument("-m", "--mfits", type=str, help="Please provide the path of the *_metafits_ppds.fits if it is not /scratch/mwavcs/$USER/[OBSID]/*_metafits_ppds.fits")
    parser.add_argument("-V", "--version", action="store_true", help="Print version and quit")
    args = parser.parse_args()

    if args.version:
        from vcstools.general_utils import print_version
        print_version()

    if args.mfits:
        metafits = args.mfits
    else:
        config = load_config_file()
        metafits = config["base_data_dir"] + '/' + str(args.obsid) + '/' + str(args.obsid) + '_metafits_ppds.fits'

    hdulist    = pyfits.open(metafits)
    freq_str   = hdulist[0].header['CHANNELS']
    freq_array = [int(f) for f in freq_str.split(',')]

    print(sfreq(freq_array))
