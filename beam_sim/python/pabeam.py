#!/usr/bin/env python


"""
Script to calculate the array factor and determine the tied-array beam of the MWA.

Author: Bradley Meyers
Date: 2016-12-8
"""



# numerical and maths modules
import numpy as np
from astropy.time import Time
from astropy.io import fits
import astropy.units as u
from astropy.constants import c, k_B
# use updated astropy ephemerides

import gc
gc.enable()

#utility and processing modules
import os
from mpi4py import MPI
import argparse
import logging

#from mwapy import ephem_utils,metadata
from vcstools.metadb_utils import getmeta
from vcstools.pointing_utils import getTargetAZZA
from vcstools.progress_bar import progress_bar
from mwa_pb import primary_beam as pb
from mwa_pb import config
import mwa_hyperbeam
beam = mwa_hyperbeam.FEEBeam(config.h5file)

logger = logging.getLogger(__name__)

def getTileLocations(obsid, flags=[], metafits=None):
    """
    Function grab the MWA tile locations for any given observation ID. Downloads the relevant metafits file from the database, saves it as <obdID>_metafits_ppds.fits.

    Input:
      obsid - the GPS observation ID
      flags - RTS tile flags (i.e. the first entry in the metafits correspond to "tile 0", irrespective of what the antenna name is)

    Return:
      a list of lists containing the following:
        list[0] = a list of tile positions East of the array centre
        list[1] = a list of tile positions North of the array centre
        list[2] = a list of tile heights about sea-level
    """

    f = fits.open(metafits)

    east = f[1].data['East'][::2]
    north = f[1].data['North'][::2]
    height = f[1].data['Height'][::2] # height above sea-level

    # MWA array centre height above sea-level
    mwacentre_h = 377.827
    height = height - mwacentre_h

    # flag the tiles from the x,y,z positions
    east = np.delete(east,flags)
    north = np.delete(north,flags)
    height = np.delete(height,flags)

    return east, north, height


def get_obstime_duration(obsid, metafits):
    """
    Funciton to grab the recorded start-time and duration of the observation

    Input:
      obsid - the GPS observation ID

    Return:
      a list containing the folloing two items (in order):
        list[0] = observation starting time in UTC
        list[1] = observation duration in seconds
    """

    # metafits file will already have been downloaded
    f = fits.open(metafits)

    return [f[0].header['DATE-OBS'], f[0].header['EXPOSURE']]



def calcWaveNumbers(wl, p, t):
    """
    Function to calculate the 3D wavenumbers for a given wavelength and az/za grid.

    Input:
      wl - central wavelength for the observation
      p - azimuth/phi (either a scalar or an array)
      t - zenith angle/theta (either a scalar or an array)

    Return:
      [kx, ky, kz] - the 3D wavenumbers
    """
    # the standard equations are:
    #   a = 2 * pi / lambda
    #   kx = a * sin(theta) * cos(phi)
    #   ky = a * sin(theta) * sin(phi)
    #   kz = a * cos(theta)
    # this is assuming that theta,phi are in the convention from Sutinjo et al. 2015
    #   i.e. phi = pi/2 - az
    kx = (2*np.pi/wl) * np.sin(t) * np.cos(p)
    ky = (2*np.pi/wl) * np.sin(t) * np.sin(p)
    kz = (2*np.pi/wl) * np.cos(t)

    return [kx,ky,kz]


# Make generator functions for Azimuth and Zenith so we don't have to create the whole
# meshgrid at once and then store it.
# This is far more memory efficient that storing the entire AZ and ZA planes in memory,
# especially with the required resolution.
def genAZZA(start, stop, step, end=False):
    """
    Generator function to use for iterating over angles (both ZA and Azimuth).

    Input:
      start - angle to start iteration from
      stop - angle to finish iteration before
      step - step size between iterations
      end - return the "stop" parameter to make range [start,stop] rather than [start,stop)

    Return:
      None - is a generator
    """
    i = 0
    num = int((stop - start) / step)
    while i < num:
        yield start + i * step
        i += 1
    if end:
        yield stop
    return


def createArrayFactor(data, rank):
    """
    Primary function to calculate the array factor with the given information.

    Input:
      delays - the beamformer delays required for point tile beam (should be a set of 16 numbers)
      time -  the time at which to evaluate the target position (as we need Azimuth and Zenith angle, which are time dependent)
      obsfreq - the centre observing frequency for the observation
      eff - the array efficiency (frequency and pointing dependent, currently require engineers to calculate for us...)
      flagged_tiles - the flagged tiles from the calibration solution (in the RTS format)
      xpos - x position of the tiles
      ypos - y position of the tiles
      zpos - z position of the tiles
      theta_res - the zenith angle resolution in degrees
      phi_res - the azimuth resolution in degrees
      za_chunk - list of ZA to compute array factor for
      write -  whether to actually write a file to disk
      rank - mpi rank (thread index)
    Return:
      results -  a list of lists cotaining [ZA, Az, beam power], for all ZA and Az in the given band
    """
    time = data['time']
    delays = data['delays']
    obsfreq = data['freqs']
    obswl = data['obswl']
    xpos = data['x']
    ypos = data['y']
    zpos = data['z']
    target_kx = data['tkx']
    target_ky = data['tky']
    target_kz = data['tkz']
    theta_res = data['tres']
    phi_res = data['pres']
    eff = data['eff']
    beam_model = data['beam']
    za = data['za_arr'][rank]
    az = data['az_arr'][rank]
    oname = data['oname']
    write = data['write']

    # calculate the relevent wavenumber for (theta,phi)
    kx, ky, kz = calcWaveNumbers(obswl, (np.pi / 2) - az, za)
    logger.debug("kx {} ky {} kz {}".format(kx[0], ky[0], kz[0]))
    logger.debug("rank {:2d} Wavenumbers shapes: {} {} {}".format(rank, kx.shape, ky.shape, kz.shape))
    #array_factor = 0 + 0.j
    array_factor = np.zeros(za.shape, dtype=np.complex_)

    # determine the interference pattern seen for each tile
    #for x, y, z in zip(xpos, ypos, zpos):
    for xyz in progress_bar(zip(xpos, ypos, zpos), "Performing tile calculations: "):
        x, y, z = xyz
        #ph = kx * x + ky * y + kz * z
        #logger.debug("kz {} z {} kz * z {}".format(kz[0], z, np.multiply(kz, z)[0]))
        ph = np.add(np.multiply(kx, x), np.multiply(ky, y))
        ph = np.add(ph, np.multiply(kz, z))
        ph_target = target_kx * x + target_ky * y + target_kz * z
        #array_factor += np.cos(ph - ph_target) + 1.j * np.sin(ph - ph_target)
        #logger.debug("ph {} ph_target {}".format(ph[0], ph_target))
        array_factor = np.add(array_factor, np.cos(ph - ph_target) + 1.j * np.sin(ph - ph_target))
        # With garbage collection (gc) enabeled this should clear the memory
        ph = ph_target = None
    kx = ky = kz = None
    logger.debug("array_factor {}".format(array_factor[0]))
    logger.debug("rank {:2d} array_factor shapes: {}".format(rank, array_factor.shape))

    # normalise to unity at pointing position
    array_factor /= len(xpos)
    array_factor_power = np.abs(array_factor)**2
    logger.debug("array_factor_power {}".format(array_factor_power[0]))

    array_factor_max = np.amax(array_factor_power)

    logger.info("rank {:2d} array factor maximum = {}".format(rank, array_factor_max))

    # calculate the tile beam at the given Az,ZA pixel
    if beam_model == 'hyperbeam':
        logger.debug("rank {:2d} begin hyperbeam".format(rank))
        jones = beam.calc_jones_array(az, za, obsfreq, delays, [1.0] * 16, True)
        logger.debug("rank {:2d} end hyperbeam".format(rank))
        jones = jones.reshape(za.shape[0], 1, 2, 2)
        vis = pb.mwa_tile.makeUnpolInstrumentalResponse(jones, jones)
        tile_xpol, tile_ypol = (vis[:, :, 0, 0].real, vis[:, :, 1, 1].real)
    elif beam_model == 'analytic':
        tile_xpol, tile_ypol = pb.MWA_Tile_analytic(za, az,
                                                freq=obsfreq, delays=delays,
                                                zenithnorm=True,
                                                power=True)
    elif beam_model == 'advanced':
        tile_xpol, tile_ypol = pb.MWA_Tile_advanced(za, az,
                                                freq=obsfreq, delays=delays,
                                                zenithnorm=True,
                                                power=True)
    elif beam_model == 'full_EE':
        tile_xpol, tile_ypol = pb.MWA_Tile_full_EE(za, az,
                                                freq=obsfreq, delays=delays,
                                                zenithnorm=True,
                                                power=True)
    logger.debug("rank {:2d} primary beam done".format(rank))
    tile_pattern = np.divide(np.add(tile_xpol, tile_ypol), 2.0)
    tile_pattern = tile_pattern.flatten()
    logger.debug("tile_pattern {}".format(tile_pattern[0]))
    logger.debug("rank {:2d} tile_pattern done".format(rank))

    # calculate the phased array power pattern
    logger.debug("rank {:2d} tile_pattern.shape {} array_factor_power.shape {}".format(rank, tile_pattern.shape, array_factor_power.shape))
    logger.debug("rank {:2d} tile_pattern {}".format(rank, tile_pattern[0]))
    logger.debug("rank {:2d} array_factor_power {}".format(rank, array_factor_power[0]))
    phased_array_pattern = np.multiply(tile_pattern, array_factor_power)
    #phased_array_pattern = tile_pattern[0][0] * np.abs(array_factor)**2 # indexing due to tile_pattern now being a 2-D array
    logger.debug("rank {:2d} phased_array_pattern done".format(rank))
    # add this contribution to the beam solid angle
    #omega_A_array = np.sin(za) * array_factor_power * np.radians(theta_res) * np.radians(phi_res)
    omega_A_array = np.sin(za)

    for A in [array_factor_power, np.radians(theta_res), np.radians(phi_res)]:
        omega_A_array = np.multiply(omega_A_array, A)
    logger.debug("rank {:2d} omega_A_array shapes: {}".format(rank, omega_A_array.shape))
    omega_A = np.sum(omega_A_array)
    logger.debug("omega_A 1 vals: za {}, array_factor_power {}, theta_res {}, phi_res {}".format(za[1], array_factor_power[1], theta_res, phi_res))
    logger.debug("omega_A 1 {}".format(omega_A_array[1]))

    # write a file based on rank of the process being used
    if write:
        logger.info("writing file from worker {0}".format(rank))
        with open(oname.replace(".dat",".{0}.dat".format(rank)), 'w') as f:
            if rank == 0:
                # if master process, write the header information first and then the data
                f.write("##File Type: Far field\n##File Format: 3\n##Source: mwa_tiedarray\n##Date: {0}\n".format(time.iso))
                f.write("** File exported by FEKO kernel version 7.0.1-482\n\n")
                f.write("#Request Name: FarField\n#Frequency: {0}\n".format(obsfreq))
                f.write("#Coordinate System: Spherical\n#No. of Theta Samples: {0}\n#No. of Phi Samples: {1}\n".format(ntheta, nphi))
                f.write("#Result Type: Gain\n#No. of Header Lines: 1\n")
                #f.write('#\t"Theta"\t"Phi"\t"Re(Etheta)"\t"Im(Etheta)"\t"Re(Ephi)"\t"Im(Ephi)"\t"Gain(Theta)"\t"Gain(Phi)"\t"Gain(Total)"\n')
                f.write('#"Theta"\t"Phi"\t"Gain(Total)"\n')

            #for res in results:
            for zad, azd, pap in zip(np.degrees(za), np.degrees(az), phased_array_pattern):
                # write each line of the data
                # we actually need to rotate out phi values by: phi = pi/2 - az because that's what FEKO expects.
                # values are calculated using that convetion, so we need to represent that here
                #f.write("{0:.5f}\t{1:.5f}\t0\t0\t0\t0\t0\t0\t{2}\n".format(res[0], res[1], res[2]))
                f.write("{0:.5f}\t{1:.5f}\t{2}\n".format(zad, azd, pap))

    else:
        logger.warn("worker {0} not writing".format(rank))

    return omega_A

def parse_options():#comm):
    #####################
    ##  OPTION PARSING ##
    #####################
    beam_models = ['analytic', 'advanced', 'full_EE', 'hyperbeam']
    parser = argparse.ArgumentParser(description="""Script to calculate the array factor required to model the tied-array beam for the MWA.
                            This is an MPI-based simulation code and will use all available processes when run
                            (i.e. there is no user choice in how many to use)""",\
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-o", "--obsid", type=str, action='store', metavar="obsID",
                        help="""Observation ID (used to figure out spatial configuration of array).
                             Also used to retrieve observation start-time and duration.""",default=None)

    parser.add_argument("--flagged_tiles", type=str, nargs='+', action='store', metavar="tile",
                        help="The tiles flagged as in when running the RTS. Must be a list of space separated tile numbers, e.g. 0 1 2 5",default=None)

    parser.add_argument("--delays", type=int, nargs='+', action='store', metavar="tile",
                        help="The tile delays", default=[0] * 16)

    parser.add_argument("-p","--target",nargs=2,metavar=("ra","dec"),
                help="The RA and DEC of the target pointing (i.e the desired phase-centre). Should be formtted like: hh:mm:ss.sss dd\"mm\'ss.ssss",
                default=("00:00:00.0000","00:00:00.0000"))

    parser.add_argument("-t","--time",type=float,action='store',metavar="time",\
                help="""The GPS time desired for beam evaluation. This will override whatever is read from the metafits.""",default=None)

    parser.add_argument("-f","--freq",type=float,action='store',metavar="frequency",help="The centre observing frequency for the observation (in Hz!)",default=184.96e6)

    parser.add_argument("-e","--efficiency",type=float,action='store',metavar="eta",help="Frequency and pointing dependent array efficiency",default=1)

    parser.add_argument("--grid_res",type=float,action='store',nargs=2,metavar=("theta_res","phi_res"),
                help="""Resolution of the Azimuth (Az) and Zenith Angle (ZA) grid to be created in degrees.
                    Be warned: setting these too small will result in a MemoryError and the job will die.""",default=(0.1,0.1))

    parser.add_argument("--coplanar",action='store_true',help="Assume the array is co-planar (i.e. height above array centre is 0 for all tiles)")

    parser.add_argument("--zenith",action='store_true',help="Assume zenith pointing (i.e  delays are 0), ZA = 0 and AZ = 0")

    parser.add_argument("--out_dir",type=str,action='store',help="Location (full path) to write the output data files",default=".")

    parser.add_argument("--write",action='store_true',
                help="""Write the beam pattern to disk when done calculating.
                    If this option is not passed, you will just get a '.stats' files containing basic information about simulation parameters and the calculated gain.""")

    parser.add_argument("--metafits", type=str, help="Metafits file location")

    parser.add_argument("--beam_model", type=str, default='hyperbeam', help='Decides the beam approximation that will be used. Options: "analytic" the analytic beam model (2012 model, fast and reasonably accurate), "advanced" the advanced beam model (2014 model, fast and slighty more accurate) or "full_EE" the full EE model (2016 model, slow but accurate). " Default: "analytic"')

    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbositylevel. Default: INFO", choices=loglevels.keys(), default="INFO")

    """
    args = None
    try:
        if comm.Get_rank() == 0:
            # parse the arguments
            args = parser.parse_args()
    finally:
        args = comm.bcast(args, root=0)

    if args is None:
        comm.Abort()
    """
    args = parser.parse_args()

    return args

###############
## Setup MPI ##
###############
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if __name__ == "__main__":
    #dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                    INFO=logging.INFO,
                    WARNING=logging.WARNING,
                    ERROR = logging.ERROR)
    args = parse_options()#comm)

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.propagate = False
    if rank == 0:
        logger.info("will use {0} processes".format(size))

        # small calculations and data gathering from arguments is fine and won't run into trouble by having multiple processes do it simultaneously
        ra = str(args.target[0]).replace('"', '').replace("'", "")
        dec = args.target[1].replace('"', '').replace("'", "")
        tres, pres = args.grid_res
        ntheta, nphi = 90 // tres, 360 // pres

        if args.flagged_tiles is None:
            flags = []
        else:
            flags = args.flagged_tiles

        if args.zenith:
            delays = [0] * 16
        else:
            delays = args.delays

        # same for obs time, master reads and then distributes
        logger.info("getting times")
        if args.time is None:
            time = Time(get_obstime_duration(args.obsid, args.metafits)[0], format='gps', fdir=args.out_dir)
        else:
            time = Time(args.time, format='gps')

        # get the target azimuth and zenith angle in radians and degrees
        # these are defined in the normal sense: za = 90 - elevation, az = angle east of North (i.e. E=90)
        logger.info("getting source position")
        if args.zenith:
            srcAz, srcZA, srcAz_deg, srcZA_deg = 0, 0, 0, 0
        else:
            srcAz, srcZA, srcAz_deg, srcZA_deg = getTargetAZZA(ra, dec, time)

        # get the tile locations from the metafits
        logger.info("getting tile locations from metafits file")
        xpos, ypos, zpos = getTileLocations(args.obsid, flags, metafits=args.metafits)
        if args.coplanar:
            zpos = np.zeros_like(xpos)

        # convert frequency to wavelength
        obswl =  args.freq / c.value

        # calculate the target (kx,ky,kz)
        target_kx, target_ky, target_kz = calcWaveNumbers(obswl, (np.pi / 2) - srcAz, srcZA)

        # set the base output file name (will be editted based on the worker rank)
        oname = "{0}/{1}_{2}_{3:.2f}MHz_tres{4}_pres{5}_{6}_{7}.dat".format(args.out_dir, args.obsid, time.gps, args.freq / 1e6, args.grid_res[0], args.grid_res[1], ra, dec)


        # figure out how many chunks to split up ZA into
        totalZAevals = ((np.pi / 2) / np.radians(tres)) + 1 # total number of ZA evaluations required
        assert totalZAevals >= size, "Total number of ZA evalutions must be >= the number of processors available"

        za = np.radians(np.linspace(0, 90, int(ntheta) + 1))
        az = np.radians(np.linspace(0, 360, int(nphi) + 1))
        logger.debug("za.shape {} az.shape {}".format(za.shape, az.shape))
        #azv, zav = np.meshgrid(az, za)
        zav, azv = np.meshgrid(za, az)
        logger.debug("zav.shape {} azv.shape {}".format(zav.shape, azv.shape))
        azv = azv.flatten()
        zav = zav.flatten()
        logger.debug("flattened zav.shape {} azv.shape {}".format(zav.shape, azv.shape))
        logger.debug("zav[0] {} zav[-1] {} azv[0] {} azv[-1] {}".format(zav[0], zav[-1], azv[0], azv[-1]))
        # split the sky into ZA chunks based on the number of available processes
        za_array = np.array_split(zav, size)
        az_array = np.array_split(azv, size)



        # create the object that contains all the data from master
        data = {'obsid' : args.obsid,
                'time':time,
                'delays':delays,
                'flags':flags,
                'freqs': args.freq,
                'obswl': obswl,
                'x':xpos,
                'y':ypos,
                'z':zpos,
                'tkx':target_kx,
                'tky':target_ky,
                'tkz':target_kz,
                'tres':tres,
                'pres':pres,
                'eff':args.efficiency,
                'beam':args.beam_model,
                'za_arr':za_array,
                'az_arr':az_array,
                'oname':oname,
                'write':args.write,
                }
    else:
        # create the data object, but leave it empty (will be filled by bcast call)
        data = None

    # now broadcast the data from the master to the slaves
    data = comm.bcast(data, root=0)
    if data:
        logger.info("broadcast received by worker {0} successfully".format(rank))
    else:
        logger.error("broadcast failed to worker {0}".format(rank))
        logger.error("!!! ABORTING !!!")
        comm.Abort(errorcode=1)

    # wait for all processes to have recieved the data
    comm.barrier()

    # create array factor for given ZA band and write to file
    beam_area = createArrayFactor(data, rank)
    logger.debug("rank {:2d} beam_area: {}".format(rank, beam_area))

    # collect results for the beam area calculation
    if rank != 0:
        comm.send(beam_area, dest=0)
    elif rank == 0:
        result = beam_area
        for i in range(1, size):
            result += comm.recv(source=i)

    # wait for everything to be collected (not really sure if this is necessary...)
    comm.barrier()

    # calculate the gain for that pointing and frequency and write a "stats" file
    if rank == 0:
        eff_area = args.efficiency * (c.value / args.freq)**2 * (4 * np.pi / result)
        gain = (1e-26) * eff_area / (2 * k_B.value)

        logger.info("==== Summary ====")
        logger.info("** Pointing **")
        logger.info("(RA, Dec)      : {0} {1}".format(ra, dec))
        logger.info("(Az, ZA) [rad] : {0:.3f} {1:.3f}".format(srcAz, srcZA))
        logger.info("(Az, ZA) [deg] : {0:.3f} {1:.3f}".format(srcAz_deg, srcZA_deg))
        logger.info("------------------------------------------")
        logger.info("** Time **")
        logger.info("GPS : {0}".format(time.gps))
        logger.info("UTC : {0}".format(time.iso))
        logger.info("MJD : {0}".format(time.mjd))
        logger.info("------------------------------------------")
        logger.info("** Telescope parameters **")
        logger.info("Observation ID       : {0}".format(args.obsid))
        logger.info("Frequency            : {0} MHz".format(args.freq/1e6))
        logger.info("Radiation efficiency : {0}".format(args.efficiency))
        logger.info("Flagged tiles        : {0}".format(" ".join(flags)))
        logger.info("------------------------------------------")
        logger.info("** Simulation resolution **")
        logger.info("Theta (ZA) resolution    [deg] : {0}".format(tres))
        logger.info("                      [arcmin] : {0}".format(tres * 60))
        logger.info("Phi   (Az) resolution    [deg] : {0}".format(pres))
        logger.info("                      [arcmin] : {0}".format(pres * 60))
        logger.info("------------------------------------------")
        logger.info("** Calculated quantities **")
        logger.info("Beam solid angle [steradians] : {0}".format(result))
        logger.info("Beam solid angle    [sq. deg] : {0}".format(np.degrees(np.degrees(result))))
        logger.info("Effective area          [m^2] : {0:.3f}".format(eff_area))
        logger.info("Effective gain         [K/Jy] : {0:.6f}".format(gain))

        # write the stats file
        with open(oname.replace(".dat", ".stats"), "w") as f:
            f.write("==== Summary ====\n")
            f.write("** Pointing **\n")
            f.write("(RA, Dec)      : {0} {1}\n".format(ra, dec))
            f.write("(Az, ZA) [rad] : {0:.3f} {1:.3f}\n".format(srcAz, srcZA))
            f.write("(Az, ZA) [deg] : {0:.3f} {1:.3f}\n".format(srcAz_deg, srcZA_deg))
            f.write("\n")
            f.write("** Time **\n")
            f.write("GPS : {0}\n".format(time.gps))
            f.write("UTC : {0}\n".format(time.iso))
            f.write("MJD : {0}\n".format(time.mjd))
            f.write("\n")
            f.write("** Telescope parameters **\n")
            f.write("Observation ID       : {0}\n".format(args.obsid))
            f.write("Frequency            : {0} MHz\n".format(args.freq/1e6))
            f.write("Radiation efficiency : {0}\n".format(args.efficiency))
            f.write("Flagged tiles        : {0}\n".format(" ".join(flags)))
            f.write("\n")
            f.write("** Simulation resolution **\n")
            f.write("Theta (ZA) resolution    [deg] : {0}\n".format(tres))
            f.write("                      [arcmin] : {0}\n".format(tres * 60))
            f.write("Phi   (Az) resolution    [deg] : {0}\n".format(pres))
            f.write("                      [arcmin] : {0}\n".format(pres * 60))
            f.write("\n")
            f.write("** Calculated quantities **\n")
            f.write("Beam solid angle [steradians] : {0}\n".format(result))
            f.write("Beam solid angle    [sq. deg] : {0}\n".format(np.degrees(np.degrees(result))))
            f.write("Effective area          [m^2] : {0:.3f}\n".format(eff_area))
            f.write("Effective gain         [K/Jy] : {0:.6f}".format(gain))


