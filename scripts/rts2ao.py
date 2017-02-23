#!/usr/bin/env python

import numpy as np
from aocal import AOCal
from astropy.io import fits
import glob
from reorder_chans import *
from optparse import OptionParser

def real2cmplx2x2mat(riririri):
    """
    Converts a list of 8 floats (real, imag, real, imag...)
    to a 2x2 complex matrix:

      (ra, ia, rb, ib, rc, ic, rd, id) -->

          [ ra+ia*I  rb+ib*I ]
          [ rc+ic*I  rd+id*I ]

    """
    result = np.empty((2, 2,),dtype=np.complex128)
    for i in range(4):
        result[i/2][i%2] = riririri[2*i] + riririri[2*i+1]*1j
    return result

def rtsfile(metafits, rts_filename_pattern="DI_JonesMatrices_node[0-9]*.dat"):
    """
    Read DI Jones matrices from RTS output files and convert to "aocal" format.
    Assumes RTS solutions are one per coarse channel.
    """
    # (Relative) comparison of RTS and OFFRINGA polarisation ordering:
    # OFFRINGA:  XX-R  XX-I  XY-R  XY-I  YX-R  YX-I  YY-R  YY-I
    # RTS:       YY-R  YY-I  YX-R  YX-I  XY-R  XY-I  XX-R  XX-R
    pol_map = [3, 2, 1, 0]

    # Get antenna reording from metafits file:
    f = fits.open(metafits)
    d = f[1].data
    nants = f[1].header[4] / 2 # I THINK this is the number of antennas
    ant_map = list(d.field('Antenna'))
    if "TileName" in f[1].header.values():
        tilenames = list(d.field('TileName'))
        ant_names = tilenames[0:-1:2]
    elif "Tile" in f[1].header.values():
        tilenames = list(d.field('Tile'))
        ant_names = ["Tile{0:03d}".format(tilenames[i*2]) for i in range(len(tilenames)/2)]
    else:
        print("Error: Antenna Name field not found")
        exit()
    ao_order  = [0 for i in range(len(ant_map)/2)]
    for i in range(len(ant_map)/2):
        ao_order[ant_map[i*2]] = i
    chans = [int(f) for f in f[0].header['CHANNELS'].split(',')]
    ch_order = np.argsort(sfreq(chans))

    # Assumptions:
    nintervals = 1
    npols = 4

    # Get file names
    rts_filenames = sorted(glob.glob(rts_filename_pattern)) # <-- Assumes file names are appropriately ordered by channel
    #rts_filenames.reverse()
    nchannels = len(rts_filenames)

    firsttime = True
    for chan in ch_order:
        rts_filename = rts_filenames[chan]
        with open(rts_filename, "r") as rts_file:
            # Common factor of all gains is a single number in the first line of the file
            amp = float(rts_file.readline())

            # The second line contains the model primary beam Jones matrix (in the direction of the calibrator)
            Jrefline = rts_file.readline()
            Jref = real2cmplx2x2mat([float(i) for i in Jrefline.split(",")])
            invJref = np.linalg.inv(Jref)

            # Read in the remainder of the Jones matrices (one per antenna)
            lines = rts_file.readlines()

            # If first time through, get number of antennas and set up data array for solution
            if firsttime:
                nantennas = len(lines)
                # Create numpy array structure
                data = np.empty((nintervals, nantennas, nchannels, npols,),dtype=np.complex128)
                data[:] = np.nan
                firsttime = False
            else:
                assert len(lines) == nantennas, "Files contain different numbers of antennas"

            # Parse each line
            for ant in range(len(lines)):
                line = lines[ant]
                jonesline = line.split(",")
                assert len(jonesline) == 2*npols, "Incorrect number of elements in Jones matrix"

                # Convert into a complex matrix
                M = real2cmplx2x2mat([float(i) for i in jonesline])
                G = np.matmul(M, invJref)

                ant_idx = ant_map[ant*2]
                for pol in range(len(pol_map)):
                    px = pol_map[pol]/2
                    py = pol_map[pol]%2
                    # Put in the Jones matrix TIMES the inverse of the reference Jones matrix
                    data[0,ant_idx,chan,pol] = G[py][px]

    new_aocal = AOCal(data, 0, 0)
    return new_aocal

if __name__ == "__main__":

    # Parse command line
    parser = OptionParser(description="rts2ao.py is a tool for converting an RTS solution to an Offringa-style solution.")

    # Add valid command line options
    parser.add_option("-d", "--dir", metavar="DIR", default="./", help="Directory containing RTS solution files (DI_Jones...) [default=%default]")
    parser.add_option("-m", "--metafits", type="string", default=None, help="Metafits file from which to determine antenna and channel order")
    parser.add_option("-r", "--rts_glob", type="string", default="DI_JonesMatrices_node[0-9]*.dat", help="Filename globbing pattern for the RTS direction independent solution files [default=%default]")
    parser.add_option("-o", "--outfile", type="string", default="calibration_solution.bin", help="Name of the file to write [default=%default]")

    # Parse!
    (opts, args) = parser.parse_args()

    if not opts.metafits:
        print "Metafits file required. [-m]"
        quit()

    ao = rtsfile(opts.metafits, opts.dir + '/' + opts.rts_glob)
    ao.tofile(opts.outfile)
