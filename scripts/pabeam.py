#!/usr/bin/env python


"""
Script to calculate the array factor and determine the tied-array beam of the MWA.

Author: Bradley Meyers
Date: 2016-12-8
"""

# numerical and maths modules
import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.constants import c, k_B
# use updated astropy ephemerides

#utility and processing modules
import os
import sys
from mpi4py import MPI
import argparse
import logging
import time as timing
import requests
import random
import gc
gc.enable()

#from mwapy import ephem_utils,metadata
from vcstools.metadb_utils import getmeta, get_common_obs_metadata
from vcstools.pointing_utils import getTargetAZZA
from vcstools.general_utils import setup_logger
from mwa_pb import primary_beam as pb
from mwa_pb import config
import mwa_hyperbeam
beam = mwa_hyperbeam.FEEBeam(config.h5file)
from vcstools.beam_sim import getTileLocations, get_obstime_duration, partial_convolve_sky_map,\
                              calcWaveNumbers, calcSkyPhase, calcArrayFactor

logger = logging.getLogger(__name__)


def createArrayFactor(za, az, data):
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
    obsid = data['obsid']
    ra = data['ra']
    dec = data['dec']
    times = data['time']
    delays = data['delays']
    xpos = data['x']
    ypos = data['y']
    zpos = data['z']
    theta_res = data['tres']
    phi_res = data['pres']
    eff = data['eff']
    beam_model = data['beam']
    #za = data['za_arr']
    #az = data['az_arr']
    out_dir = data['out_dir']
    write = data['write']
    coplanar = data['coplanar']
    plot_jobid = data['plot_jobid']

    # work out core depent part of data to process
    nfreq = len(data['freqs'])
    ifreq = rank % nfreq
    ichunk = rank // nfreq
    obsfreq = data['freqs'][ifreq]
    #za = zav_chunks[ichunk]
    #az = azv_chunks[ichunk]


    # convert frequency to wavelength
    obswl =  obsfreq / c.value
    # calculate the relevent wavenumber for (theta,phi)
    logger.info( "rank {:3d} Calculating wavenumbers".format(rank))
    kx, ky, kz = calcWaveNumbers(obswl, (np.pi / 2) - az, za)
    #logger.debug("kx[0] {} ky[0] {} kz[0] {}".format(kx[0], ky[0], kz[0]))
    logger.debug("rank {:3d} Wavenumbers shapes: {} {} {}".format(rank, kx.shape, ky.shape, kz.shape))

    # phase of each sky position for each tile
    logger.info( "rank {:3d} Calculating phase of each sky position for each tile".format(rank))
    ph_tile = calcSkyPhase(xpos, ypos, zpos, kx, ky, kz, coplanar=coplanar)

    requests.get('https://ws.mwatelescope.org/progress/update',
                 params={'jobid':plot_jobid,
                         'workid':rank,
                         'current':2,
                         'total':5,
                         'desc':'WaveN_phase'})

    # calculate the tile beam at the given Az,ZA pixel
    logger.info( "rank {:3d} Calculating tile beam".format(rank))
    if beam_model == 'hyperbeam':
        # This method is no longer needed as mwa_pb uses hyperbeam
        jones = beam.calc_jones_array(az, za, obsfreq, delays, [1.0] * 16, True)
        jones = jones.reshape(za.shape[0], 1, 2, 2)
        vis = pb.mwa_tile.makeUnpolInstrumentalResponse(jones, jones)
        tile_xpol, tile_ypol = (vis[:, :, 0, 0].real, vis[:, :, 1, 1].real)
    elif beam_model == 'analytic':
        tile_xpol, tile_ypol = pb.MWA_Tile_analytic(za, az,
                                                freq=obsfreq, delays=[delays, delays],
                                                zenithnorm=True,
                                                power=True)
    elif beam_model == 'advanced':
        tile_xpol, tile_ypol = pb.MWA_Tile_advanced(za, az,
                                                freq=obsfreq, delays=[delays, delays],
                                                zenithnorm=True,
                                                power=True)
    elif beam_model == 'full_EE':
        tile_xpol, tile_ypol = pb.MWA_Tile_full_EE(za, az,
                                                freq=obsfreq, delays=np.array([delays, delays]),
                                                zenithnorm=True,
                                                power=True,
                                                interp=False)
    #logger.info("rank {:3d} Combining tile pattern".format(rank))
    tile_pattern = np.divide(np.add(tile_xpol, tile_ypol), 2.0)
    tile_pattern = tile_pattern.flatten()
    logger.debug("max(tile_pattern) {}".format(max(tile_pattern)))

    requests.get('https://ws.mwatelescope.org/progress/update',
                 params={'jobid':plot_jobid,
                         'workid':rank,
                         'current':3,
                         'total':5,
                         'desc':'tilebeam'})

    omega_A_times = []
    sum_B_T_times = []
    sum_B_times   = []
    for time in times:
        # set the base output file name (will be editted based on the worker rank)
        oname = "{0}/{1}_{2}_{3:.2f}MHz_tres{4}_pres{5}_{6}_{7}.dat".format(out_dir, obsid, time, obsfreq / 1e6, theta_res, phi_res, ra, dec)

        # get the target azimuth and zenith angle in radians and degrees
        # these are defined in the normal sense: za = 90 - elevation, az = angle east of North (i.e. E=90)
        #logger.info("getting source position")
        srcAz, srcZA, srcAz_deg, srcZA_deg = getTargetAZZA(ra, dec, time)# calculate the target (kx,ky,kz)
        target_kx, target_ky, target_kz = calcWaveNumbers(obswl, (np.pi / 2) - srcAz, srcZA)

        # Get phase of target for each tile
        # phase of each sky position for each tile
        logger.info( "rank {:3d} Calculating phase of each sky position for the target at time {}".format(rank, time))
        ph_target = calcSkyPhase(xpos, ypos, zpos, target_kx, target_ky, target_kz)

        # determine the interference pattern seen for each tile
        logger.info( "rank {:3d} Calculating array_factor".format(rank))
        array_factor, array_factor_power = calcArrayFactor(ph_tile, ph_target)

        #logger.debug("array_factor_power[0] {}".format(array_factor_power[0]))
        logger.debug("rank {:3d} array_factor[0] {}".format(rank, array_factor[0]))
        logger.debug("rank {:3d} array_factor shapes: {}".format(rank, array_factor.shape))
        logger.debug("rank {:3d} array factor maximum = {}".format(rank, np.amax(array_factor_power)))

        requests.get('https://ws.mwatelescope.org/progress/update',
                    params={'jobid':plot_jobid,
                            'workid':rank,
                            'current':4,
                            'total':5,
                            'desc':'array_factor'})

        # calculate the phased array power pattern
        logger.info( "rank {:3d} Calculating phased array pattern".format(rank))
        logger.debug("rank {:3d} tile_pattern.shape {} array_factor_power.shape {}".format(rank, tile_pattern.shape, array_factor_power.shape))
        phased_array_pattern = np.multiply(tile_pattern, array_factor_power)
        #phased_array_pattern = tile_pattern[0][0] * np.abs(array_factor)**2 # indexing due to tile_pattern now being a 2-D array

        # add this contribution to the beam solid angle
        logger.info( "rank {:3d} Calculating omega_A".format(rank))
        #omega_A_array = np.sin(za) * array_factor_power * np.radians(theta_res) * np.radians(phi_res)
        omega_A_array = np.multiply(np.multiply(np.multiply(np.sin(za), array_factor_power), np.radians(theta_res)), np.radians(phi_res))
        logger.debug("rank {:3d} omega_A_array shapes: {}".format(rank, omega_A_array.shape))
        omega_A = np.sum(omega_A_array)
        logger.debug( "rank {:3d} freq {:.2f}MHz beam_area: {}".format(rank, obsfreq/1e6, omega_A))
        omega_A_times.append(omega_A)
        #logger.debug("omega_A[1] vals: za[1] {}, array_factor_power[1] {}, theta_res {}, phi_res {}".format(za[1], array_factor_power[1], theta_res, phi_res))
        #logger.debug("omega_A[1] {}".format(omega_A_array[1]))

        # Do a partial sum over the sky and frequency to be finished outside of the function
        sum_B_T, sum_B = partial_convolve_sky_map(az, za, phased_array_pattern, obsfreq, time)
        sum_B_T_times.append(sum_B_T)
        sum_B_times.append(sum_B)


        # write a file based on rank of the process being used
        if write:
            logger.info("writing file from worker {0}".format(rank))
            with open(oname, 'w') as f:
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
            logger.info("rank {:3d} done and not writing".format(rank))

    # Average over time
    omega_A = np.average(omega_A_times)
    sum_B_T = np.average(sum_B_T_times)
    sum_B   = np.average(sum_B_times)

    # return lists of results for each time step
    return [omega_A, sum_B_T, sum_B]

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
    #####################
    ##  OPTION PARSING ##
    #####################
    beam_models = ['analytic', 'advanced', 'full_EE', 'hyperbeam']
    parser = argparse.ArgumentParser(description="""Script to calculate the array factor required to model the tied-array beam for the MWA.
                                     This is an MPI-based simulation code and will use all available processes when run
                                     (i.e. there is no user choice in how many to use)""",\
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-o", "--obsid", type=str,
                        help="""Observation ID (used to figure out spatial configuration of array).
                             Also used to retrieve observation start-time and duration.""")
    parser.add_argument("--flagged_tiles", type=str, nargs='+',
                        help="The tiles flagged as in when running the RTS. Must be a list of space separated tile numbers, e.g. 0 1 2 5")
    parser.add_argument("--delays", type=int, nargs='+', default=[0] * 16,
                        help="The tile delays")
    parser.add_argument("-f", "--freq", type=float, nargs='+', default=None,
                        help="The centre observing frequency for the observation (in Hz!)")

    parser.add_argument("-p", "--target", type=str,
                        help="The RA and DEC of the target pointing (i.e the desired phase-centre). Should be formtted like: hh:mm:ss.ss_dd:mm:ss.ss")
    parser.add_argument("-b", "--begin", type=int,
                        help="""The GPS time to begin the simulation from. This will override whatever is read from the metafits.""")
    parser.add_argument("-d", "--duration", type=int,
                        help="""The duration of the simulation in seconds from the begin time. This will override whatever is read from the metafits.""")
    parser.add_argument("-s", "--step", type=int, default=1800,
                        help="""The step between simulations in seconds. Default: 1800 (30 mins).""")

    parser.add_argument("-e", "--efficiency", type=float, default=1,
                        help="Frequency and pointing dependent array efficiency")
    parser.add_argument("--grid_res", type=float, action='store', nargs=2, metavar=("theta_res","phi_res"), default=(0.1,0.1),
                        help="""Resolution of the Azimuth (Az) and Zenith Angle (ZA) grid to be created in degrees.
                             Be warned: setting these too small will result in a MemoryError and the job will die.""")
    parser.add_argument("--coplanar", action='store_true',
                        help="Assume the array is co-planar (i.e. height above array centre is 0 for all tiles)")
    parser.add_argument("--out_dir", type=str, action='store', default=".",
                        help="Location (full path) to write the output data files")
    parser.add_argument("--write", action='store_true',
                        help="""Write the beam pattern to disk when done calculating.
                            If this option is not passed, you will just get a '.stats' files containing basic information about simulation parameters and the calculated gain.""")
    parser.add_argument("--metafits", type=str,
                        help="Metafits file location")
    parser.add_argument("--beam_model", type=str, default='full_EE',
                        help='Decides the beam approximation that will be used. Options: "analytic" the analytic beam model (2012 model, fast and reasonably accurate), "advanced" the advanced beam model (2014 model, fast and slighty more accurate) or "full_EE" the full EE model (2016 model, slow but accurate). " Default: "analytic"')
    parser.add_argument("-L", "--loglvl", type=str, choices=loglevels.keys(), default="INFO",
                        help="Logger verbositylevel. Default: INFO")
    args = parser.parse_args()

    # set up the logger for stand-alone execution
    logger = setup_logger(logger, log_level=loglevels[args.loglvl])

    if rank == 0:
        logger.info("will use {0} processes".format(size))

        # small calculations and data gathering from arguments is fine and won't run into trouble by having multiple processes do it simultaneously
        ra, dec = args.target.split("_")
        tres, pres = args.grid_res
        ntheta, nphi = 90 // tres, 360 // pres

        if args.flagged_tiles is None:
            flags = []
        else:
            flags = args.flagged_tiles
        delays = args.delays

        # same for obs time, master reads and then distributes
        logger.info("getting times")
        obs_begin, obs_duration = get_obstime_duration(args.obsid, args.metafits)
        if not args.begin:
            # Use observation begin time from metafits
            args.begin = Time(obs_begin, format='isot', scale='utc').gps
        if not args.duration:
            # Use observation duration time from metafits
            args.duration = int(obs_duration)
        sim_end = args.begin + args.duration -1
        # Make a range of times to test
        times_range = np.arange(args.begin, sim_end, args.step)
        # Centre the range
        times_range = times_range + (sim_end - times_range[-1]) // 2
        logger.info("running for {} time steps: {}".format(len(times_range), times_range))
        # Convert to astropy format
        time = Time(times_range, format='gps')

        if args.freq is None:
            # By default use the all coarse channel centre frequencies
            args.freq = np.array(get_common_obs_metadata(1275751816)[-1])*1.28*1e6

        # get the tile locations from the metafits
        logger.info("getting tile locations from metafits file")
        xpos, ypos, zpos = getTileLocations(args.obsid, flags, metafits=args.metafits)
        logger.debug("xpos: {}".format(xpos))
        logger.debug("ypos: {}".format(ypos))
        logger.debug("zpos: {}".format(zpos))
        #if args.coplanar:
        #    zpos = np.zeros_like(xpos)

        assert size % len(args.freq) == 0, "Frequencies must be a factor of the number of processors available"
        # number of chunks to split the positions grid into
        nchunks = size // len(args.freq)

        za = np.radians(np.linspace(0, 90, int(ntheta) + 1))
        #za_s = np.array_split(za, nchunks)
        az = np.radians(np.linspace(0, 360, int(nphi) + 1))
        logger.debug("za.shape {} az.shape {}".format(za.shape, az.shape))
        zav, azv = np.meshgrid(za, az)
        logger.debug("zav.shape {} azv.shape {}".format(zav.shape, azv.shape))
        azv = azv.flatten()
        zav = zav.flatten()
        logger.debug("flattened zav.shape {} azv.shape {}".format(zav.shape, azv.shape))
        logger.debug("zav.dtype {} azv.dtype {}".format(zav.dtype, azv.dtype))
        azv_chunks = np.array_split(azv, nchunks)
        zav_chunks = np.array_split(zav, nchunks)

        npositions = [len(chunk) for chunk in azv_chunks]

        # Create a progress plot
        plot_jobid = random.randint(0, 1000) # random job id
        requests.get('https://ws.mwatelescope.org/progress/create',
                     params={'jobid':plot_jobid,
                             'desc':'pabeam simultation',
                             'maxage':43200}) # 12 hours
        logger.info('plot_jobid: {}'.format(plot_jobid))
        logger.info('View the plot at: https://ws.mwatelescope.org/progress/show?jobid={}'.format(plot_jobid))


        # create the object that contains all the data from master
        data = {'obsid' : args.obsid,
                'ra':ra,
                'dec':dec,
                'time':time,
                'delays':delays,
                'flags':flags,
                'freqs': args.freq,
                'x':xpos,
                'y':ypos,
                'z':zpos,
                'tres':tres,
                'pres':pres,
                'eff':args.efficiency,
                'beam':args.beam_model,
                #'za_arr':za_s,
                #'az_arr':az,
                'out_dir':args.out_dir,
                'write':args.write,
                'nchunks':nchunks,
                'coplanar':args.coplanar,
                'npositions':npositions,
                'plot_jobid':plot_jobid
                }
        logger.debug("data: {}".format(data))

        # Start timing
        start_time = timing.perf_counter()
    else:
        # create the data object, but leave it empty (will be filled by bcast call)
        data = None

    # now broadcast the data from the master to the slaves
    data = comm.bcast(data, root=0)

    # broadcasting az and za arrays seperately as it can cause overflow errors
    nfreq = len(data['freqs'])
    ichunk = rank // nfreq
    zav_chunk = np.empty(data['npositions'][ichunk], dtype=np.float64)
    azv_chunk = np.empty(data['npositions'][ichunk], dtype=np.float64)

    for i in range(1, size):
        # send data to each rank
        if rank == 0:
            logger.debug("send rank {}".format(i))
            ichunk = i // nfreq
            zav_chunk = comm.Send([zav_chunks[ichunk], MPI.DOUBLE], dest=i)
        # receive the data for each rank
        if rank == i:
            ichunk = i // nfreq
            comm.Recv([zav_chunk, MPI.DOUBLE], source=0)
        comm.barrier()
    if rank == 0:
        zav_chunk = zav_chunks[ichunk]
        azv_chunk = azv_chunks[ichunk]
        logger.debug("rank {:3d} zav_chunk {}".format(rank, zav_chunk))
        logger.debug("rank {:3d} azv_chunk {}".format(rank, azv_chunk))
    #if rank != 0:
    #    zav_chunks = [None] * data['nchunks']
    #    azv_chunks = [None] * data['nchunks']
    #zav = comm.bcast([zav, MPI.DOUBLE], root=0)
    #for i in range(data['nchunks']):
    #    zav_chunks[i] = comm.bcast(zav_chunks[i], root=0)
    #    azv_chunks[i] = comm.bcast(azv_chunks[i], root=0)
    #any(x is None for x in list_1)
    #if data:
    logger.debug("rank {:3d} zav_chunk {}".format(rank, zav_chunk))
    logger.debug("rank {:3d} azv_chunk {}".format(rank, azv_chunk))
    if data and (zav_chunk is not None) and (azv_chunk is not None):
        logger.info("broadcast received by worker {0} successfully".format(rank))
    else:
        logger.error("broadcast failed to worker {0}".format(rank))
        logger.error("!!! ABORTING !!!")
        comm.Abort(errorcode=1)


    # wait for all processes to have recieved the data
    comm.barrier()

    # Show progress began
    requests.get('https://ws.mwatelescope.org/progress/update',
                 params={'jobid':data['plot_jobid'],
                         'workid':rank,
                         'current':1,
                         'total':5,
                         'desc':'start'})

    # Perform cacluation
    results_array = createArrayFactor(zav_chunk, azv_chunk, data)

    requests.get('https://ws.mwatelescope.org/progress/update',
                 params={'jobid':data['plot_jobid'],
                         'workid':rank,
                         'current':5,
                         'total':5,
                         'desc':'Done'})

    # collect results for the beam area calculation
    if rank != 0:
        comm.send(results_array, dest=0)
    elif rank == 0:
        omega_A_sky_freq = np.empty(size)
        sum_B_T_sky_freq = np.empty(size)
        sum_B_sky_freq   = np.empty(size)
        omega_A, sum_B_T, sum_B = results_array
        omega_A_sky_freq[0] = omega_A
        sum_B_T_sky_freq[0] = sum_B_T
        sum_B_sky_freq[0]   = sum_B
        for i in range(1, size):
            results_array = comm.recv(source=i)
            omega_A, sum_B_T, sum_B = results_array
            omega_A_sky_freq[i] = omega_A
            sum_B_T_sky_freq[i] = sum_B_T
            sum_B_sky_freq[i]   = sum_B

    # wait for everything to be collected (not really sure if this is necessary...)
    comm.barrier()

    # calculate the gain for that pointing and frequency and write a "stats" file
    if rank == 0:
        logger.info("pabeam benchmark: {} s".format(timing.perf_counter()-start_time))
        # split the omega_A values into the sky chuncks
        omega_A_sky_freq = np.array_split(omega_A_sky_freq, nchunks)
        # Sum over the sky chuncks
        omega_A_freq = np.sum(omega_A_sky_freq, axis=0)
        # Average over frequency
        omega_A = np.average(omega_A_freq)
        # Opposite method
        #omega_A_sky = np.average(omega_A_sky_freq, axis=1)
        #omega_A = np.sum(omega_A_sky)

        # Convert to eff area and gain (work out for each frequency then average)
        eff_area = args.efficiency * 4 * np.pi * np.average( (c.value / np.array(args.freq))**2 / omega_A_freq)
        gain = (1e-26) * eff_area / (2 * k_B.value)

        # Sum all the partial sums then do the final division
        t_ant = np.sum(sum_B_T_sky_freq) / np.sum(sum_B_sky_freq)

        logger.info("==== Summary ====")
        logger.info("** Pointing **")
        logger.info("(RA, Dec)      : {0} {1}".format(ra, dec))
        #logger.info("(Az, ZA) [rad] : {0:.3f} {1:.3f}".format(srcAz, srcZA))
        #logger.info("(Az, ZA) [deg] : {0:.3f} {1:.3f}".format(srcAz_deg, srcZA_deg))
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
        logger.info("Delays               : {0}".format(' '.join(np.array(delays, dtype=str))))
        logger.info("------------------------------------------")
        logger.info("** Simulation resolution **")
        logger.info("Theta (ZA) resolution    [deg] : {0}".format(tres))
        logger.info("                      [arcmin] : {0}".format(tres * 60))
        logger.info("Phi   (Az) resolution    [deg] : {0}".format(pres))
        logger.info("                      [arcmin] : {0}".format(pres * 60))
        logger.info("------------------------------------------")
        logger.info("** Calculated quantities **")
        logger.info("Beam solid angle [steradians] : {0}".format(omega_A))
        logger.info("Beam solid angle    [sq. deg] : {0}".format(np.degrees(np.degrees(omega_A))))
        logger.info("Effective area          [m^2] : {0:.3f}".format(eff_area))
        logger.info("Effective gain         [K/Jy] : {0:.6f}".format(gain))
        logger.info("Antena Temperature        [K] : {0:6.1f}".format(t_ant))

        # write the stats file
        oname = "{0}/{1}_{2:.2f}-{3:.2f}MHz_tres{4}_pres{5}_{6}_{7}.stats".format(args.out_dir, args.obsid, args.freq[0]/1e6, args.freq[-1]/1e6, tres, pres, ra, dec)
        with open(oname, "w") as f:
            f.write("==== Summary ====\n")
            f.write("** Pointing **\n")
            f.write("(RA, Dec)      : {0} {1}\n".format(ra, dec))
            #f.write("(Az, ZA) [rad] : {0:.3f} {1:.3f}\n".format(srcAz, srcZA))
            #f.write("(Az, ZA) [deg] : {0:.3f} {1:.3f}\n".format(srcAz_deg, srcZA_deg))
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
            f.write("Delays               : {0}".format(' '.join(np.array(delays, dtype=str))))
            f.write("\n")
            f.write("** Simulation resolution **\n")
            f.write("Theta (ZA) resolution    [deg] : {0}\n".format(tres))
            f.write("                      [arcmin] : {0}\n".format(tres * 60))
            f.write("Phi   (Az) resolution    [deg] : {0}\n".format(pres))
            f.write("                      [arcmin] : {0}\n".format(pres * 60))
            f.write("\n")
            f.write("** Calculated quantities **\n")
            f.write("Beam solid angle [steradians] : {0}\n".format(omega_A))
            f.write("Beam solid angle    [sq. deg] : {0}\n".format(np.degrees(np.degrees(omega_A))))
            f.write("Effective area          [m^2] : {0:.3f}\n".format(eff_area))
            f.write("Effective gain         [K/Jy] : {0:.6f}\n".format(gain))
            f.write("Antena Temperature        [K] : {0:6.1f}\n".format(t_ant))