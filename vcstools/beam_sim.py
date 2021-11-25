"""
The functions required to simulate the tied-array beam response of the MWA. 
All equations can be found in https://ui.adsabs.harvard.edu/abs/2018IAUS..337..378M/abstract
"""
import numpy as np
import sys
import os
from itertools import chain

from astropy.constants import c
from astropy.io import fits

#from mwapy import ephem_utils,metadata
from mwa_pb.primarybeammap_tant import get_Haslam, map_sky

import logging
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
        list[3] = a list of tile cable delays
    """

    f = fits.open(metafits)

    east = f[1].data['East'][::2]
    north = f[1].data['North'][::2]
    height = f[1].data['Height'][::2] # height above sea-level
    cable_delays_raw = f[1].data['Length']

    # Format Cable delays to floats
    cable_delays = []
    for cab in cable_delays_raw:
        cable_delays.append(float(cab[3:]))
    cable_delays = np.array(cable_delays)

    # MWA array centre height above sea-level
    mwacentre_h = 377.827
    height = height - mwacentre_h

    # flag the tiles from the x,y,z positions
    east = np.delete(east,flags)
    north = np.delete(north,flags)
    height = np.delete(height,flags)
    cable_delays = np.delete(cable_delays,flags)

    return east, north, height, cable_delays


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


def calc_pixel_area(za, az_res, za_res):
    """
    Calculate the area of a pixel on the sky from their height, width and zenith angle

    Parameters:
    -----------
    za: float
        The zenith angle of the pixel in radians
    az_res: float
        The azuimuth resolution of the pixel in degrees
    za_res: float
        The zenith   resolution of the pixel in degrees

    Returns:
    --------
    area: float
        Area of the pixel in square radians
    """
    return np.radians(az_res) * np.radians(za_res) * np.sin(za)


def calcWaveNumbers(freq, p, t):
    """
    Function to calculate the 3D wavenumbers for a given wavelength and az/za grid. 
    This is the part of equation 7 within the square brackets not includeing x_n, y_n and z_n

    Parameters:
    -----------
    freq: float
        Central frequency for the observation in Hz
    p: float
        azimuth/phi (either a scalar or an array)
        this is assuming that theta,phi are in the convention from Sutinjo et al. 2015
    t: float
        zenith angle/theta (either a scalar or an array)

    Returns:
    --------
      [kx, ky, kz] - the 3D wavenumbers
    """
    C =  2 * np.pi * freq / c.value
    kx = C * np.multiply(np.sin(t), np.cos(p))
    ky = C * np.multiply(np.sin(t), np.sin(p))
    kz = C * np.cos(t)

    return [kx,ky,kz]


def calcSkyPhase(xpos, ypos, zpos, kx, ky, kz, coplanar=False):
    """
    Completes the calculation of equation 7 to get the phase of the tiles for each position on the sky

    Parameters:
    -----------
    xpos[N]:
        A list of tile positions East of the array centre
    ypos[N]:
        A list of tile positions North of the array centre
    zpos[N]:
        A list of tile heights about sea-level
    kx[az/za]:
        The x 3D wavenumbers for a given wavelength and az/za grid
    ky[az/za]:
        The y 3D wavenumbers for a given wavelength and az/za grid
    kz[az/za]:
        The z 3D wavenumbers for a given wavelength and az/za grid

    Returns:
    --------
    ph_tile[N][az/za]:
        A list of the phases for each tile
    """
    #ph_tile = []
    #for x, y, z in zip(xpos, ypos, zpos):
        #ph = kx * x + ky * y + kz * z
    logger.debug("np.multiply(np.tile(kx, (len(xpos), 1))).shape : {}".format(np.tile(kx, (len(xpos), 1)).shape))
    if coplanar:
        ph_tile = list(chain(np.add(np.multiply(kx, x), np.multiply(ky, y)) for x, y in zip(xpos, ypos)))
    else:
        ph_tile = list(chain(np.add(np.add(np.multiply(kx, x), np.multiply(ky, y)), np.multiply(kz, z)) for x, y, z in zip(xpos, ypos, zpos)))
    logger.debug("ph_tile.shape[0] : {}".format(ph_tile[0].shape))
    return ph_tile


def calc_geometric_delay_distance(p, t):
    """
    Equation 1 of Ord 2019. Changed cos(el) to sin(za)
    """
    gx = np.multiply(np.sin(t), np.sin(p))
    gy = np.multiply(np.sin(t), np.cos(p))
    #gx = np.multiply(np.sin(t), np.cos(p))
    #gy = np.multiply(np.sin(t), np.sin(p))
    gz = np.cos(t)
    return gx, gy, gz


def cal_phase_ord(xpos, ypos, zpos, delays, gx, gy, gz, freq, coplanar=False, no_delays=False):
    """
    Equation 2 and 3 of Ord 2019.

    Parameters:
    -----------
    xpos[N]:
        A list of tile positions East of the array centre
    ypos[N]:
        A list of tile positions North of the array centre
    zpos[N]:
        A list of tile heights about sea-level
    gx[az/za]:
        The x 3D wavenumbers for a given wavelength and az/za grid
    gy[az/za]:
        The y 3D wavenumbers for a given wavelength and az/za grid
    gz[az/za]:
        The z 3D wavenumbers for a given wavelength and az/za grid

    Returns:
    --------
    ph_tile[N][az/za]:
        A list of the phases for each tile
    """
    if no_delays:
        # Set cable delays to zeros so that they're not included
        delays = np.zeros_like(delays)
    #ph_tile = []
    #for x, y, z in zip(xpos, ypos, zpos):
        #ph = kx * x + ky * y + kz * z
    # Equation 1 and 2 (TODO there is a minus in the ord paper so check if this is required)
    # /delta t * c = gx * x + gy * y + gz * z + L
    if coplanar:
        wl = list(chain(np.subtract(np.add(np.multiply(gx, x), np.multiply(gy, y)), l) for x, y, l in zip(xpos, ypos, delays)))
    else:
        wl = list(chain(np.subtract(np.add(np.add(np.multiply(gx, x), np.multiply(gy, y)), np.multiply(gz, z)), l) for x, y, z, l in zip(xpos, ypos, zpos, delays)))
    # Equation 3
    ph_tile = 2 * np.pi * np.array(wl) * freq / c.value
    logger.debug("ph_tile.shape[0] : {}".format(ph_tile[0].shape))
    return ph_tile


def calcArrayFactor(ph_tiles, ph_targets):
    """
    Calculates array factor pointed at some target zenith angle (za) and azimuth (az) (equation 11)

    Parameters:
    -----------
    ph_tiles[N][az/za]:
        List of the sky phases for the tiles
    ph_targets[N]:
        List of the sky phases for the target

    Returns:
    --------
    array_factor[az/za]:
        The array factor for each za and az
    array_factor_power[az/za]:
        The array factor power for each za and az
    """
    #array_factor = np.zeros(za.shape, dtype=np.complex_)
    #for i, _ in enumerate(xpos):
        #array_factor += np.cos(ph - ph_target) + 1.j * np.sin(ph - ph_target)
    #array_factor_tiles = list(chain(np.cos(ph_tile - ph_target) + 1.j * np.sin(ph_tile - ph_target) for ph_tile, ph_target in zip(ph_tiles, ph_targets)))
    array_factor_tiles = list(chain( np.multiply( np.cos(ph_tile) + 1.j * np.sin(ph_tile),  np.cos(ph_target) - 1.j * np.sin(ph_target) ) for ph_tile, ph_target in zip(ph_tiles, ph_targets)))
    array_factor = np.sum(array_factor_tiles, axis=0)

    # normalise to unity at pointing position
    array_factor = np.divide(array_factor, len(ph_tiles))
    array_factor_power = np.abs(array_factor)**2

    return array_factor, array_factor_power


def partial_convolve_sky_map(az_grid, za_grid, pixel_area, phased_array_pattern, freq, time):
    # Get Haslam and interpolate onto grid
    # taken from https://github.com/MWATelescope/mwa_pb/blob/master/mwa_pb/primarybeammap_tant.py
    # then edited to include pixel size
    logger.debug("freq: {}".format(freq))
    my_map = get_Haslam(freq)
    #mask = numpy.isnan(za_grid)
    #za_grid[numpy.isnan(za_grid)] = 90.0  # Replace nans as they break the interpolation

    #Supress print statements of the primary beam model functions
    sys.stdout = open(os.devnull, 'w')
    sky_grid = map_sky(my_map['skymap'], my_map['RA'], my_map['dec'], time, az_grid, za_grid)
    logger.debug("np.max(sky_grid): {}".format(np.max(sky_grid.flatten())))
    logger.debug("sky_grid.shape: {}".format(sky_grid.flatten().shape))
    logger.debug("phased_array_pattern.shape: {}".format(phased_array_pattern.shape))
    sys.stdout = sys.__stdout__

    logger.debug("np.average(sky_grid.flatten()): {}".format(np.average(sky_grid.flatten())))
    sum_B_T = np.sum(phased_array_pattern * sky_grid.flatten() * pixel_area)
    sum_B   = np.sum(phased_array_pattern * pixel_area)
    return sum_B_T, sum_B


def read_sefd_file(sefd_file, all_data=False):
    """
    Read in the output sefd file from a pabeam.py simulation.

    Parameters:
    -----------
    sefd_file: str
        The location of the sefd file to be read in.
    all_data: boolean
        If True will return the freq, sefd, t_sys, t_ant, gain, effective_area.
        If False will only return the sefd. Default False
    """
    input_array = np.loadtxt(sefd_file, dtype=float)
    if all_data:
        freq  = input_array[:,0]
        time  = input_array[:,1]
        sefd  = input_array[:,2]
        t_sys = input_array[:,3]
        t_ant = input_array[:,4]
        gain  = input_array[:,5]
        effective_area = input_array[:,6]
        beam_solid_angle = input_array[:,7]
        return freq, time, sefd, t_sys, t_ant, gain, effective_area, beam_solid_angle
    else:
        time  = input_array[:,1]
        sefd_freq_time  = np.array(input_array[:,2])
        # Work out the number of time steps
        ntime = list(time).count(list(time)[0])
        # Reshape into [freq][time] array
        sefd_freq_time = sefd_freq_time.reshape((sefd_freq_time.shape[0]//ntime, ntime))
        sefd_time = np.average(sefd_freq_time, axis=0)
        # Calc mean and std
        sefd_mean = np.mean(sefd_time)
        sefd_std  = np.std(sefd_time)
        return sefd_freq_time, sefd_mean, sefd_std