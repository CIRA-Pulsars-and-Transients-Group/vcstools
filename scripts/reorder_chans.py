#!/usr/bin/python3
"""
Contains a single function that reads in a list of integers (channel numbers),
and returns a list of the same integers whose order depends on whether or not
any of the channels are above 128.
"""
from astropy.io import fits as pyfits
import argparse
from vcstools.general_utils import sfreq

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Print the coarse channel sequence as in gpubox/subband order. Make sure you have *_metafits_ppds.fits")
    parser.add_argument("-o", "--obsid", type=int, help="Observation ID of target", default=None, required=True)
    parser.add_argument("-m", "--mfits", type=str, help="Please provide the path of the *_metafits_ppds.fits if it is not /astro/mwavcs/vcs/[OBID]/*_metafits_ppds.fits")
    parser.add_argument("-V", "--version", action="store_true", help="Print version and quit")
    args = parser.parse_args()

    if args.version:
        import sys
        try:
            import version
            print(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            print("Couldn't import version.py - have you installed vcstools?")
            print("ImportError: {0}".format(ie))
            sys.exit(0)

    if args.mfits:
        metafits = args.mfits
    else:
        metafits = '/astro/mwavcs/vcs/' + str(args.obsid) + '/' + str(args.obsid) + '_metafits_ppds.fits'

    hdulist    = pyfits.open(metafits)
    freq_str   = hdulist[0].header['CHANNELS']
    freq_array = [int(f) for f in freq_str.split(',')]

    print(sfreq(freq_array))
