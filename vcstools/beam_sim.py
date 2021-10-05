"""
The functions required to simulate the tied-array beam response of the MWA. 
All equations can be found in https://ui.adsabs.harvard.edu/abs/2018IAUS..337..378M/abstract
"""
# numerical and maths modules
import numpy as np
from itertools import chain

#utility and processing modules
import sys
import os
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
    This is the part of equation 7 within the square brackets not includeing x_n, y_n and z_n

    Parameters:
    -----------
      wl - central wavelength for the observation
      p - azimuth/phi (either a scalar or an array)
      t - zenith angle/theta (either a scalar or an array)

    Returns:
    --------
      [kx, ky, kz] - the 3D wavenumbers
    """
    # the standard equations are:
    #   a = 2 * pi / lambda
    #   kx = a * sin(theta) * cos(phi)
    #   ky = a * sin(theta) * sin(phi)
    #   kz = a * cos(theta)
    # this is assuming that theta,phi are in the convention from Sutinjo et al. 2015
    #   i.e. phi = pi/2 - az
    kx = (2*np.pi/wl) * np.multiply(np.sin(t), np.cos(p))
    ky = (2*np.pi/wl) * np.multiply(np.sin(t), np.sin(p))
    kz = (2*np.pi/wl) * np.cos(t)

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
    array_factor_tiles = list(chain(np.cos(ph_tile - ph_target) + 1.j * np.sin(ph_tile - ph_target) for ph_tile, ph_target in zip(ph_tiles, ph_targets)))
    array_factor = np.sum(array_factor_tiles, axis=0)

    # normalise to unity at pointing position
    array_factor /= len(ph_tiles)
    array_factor_power = np.abs(array_factor)**2

    return array_factor, array_factor_power


def convolve_sky_map(az_grid, za_grid, phased_array_pattern, freq, time):
    # Get Haslam and interpolate onto grid
    # taken from https://github.com/MWATelescope/mwa_pb/blob/master/mwa_pb/primarybeammap_tant.py
    # then edited to include pixel size
    my_map = get_Haslam(freq)
    #mask = numpy.isnan(za_grid)
    #za_grid[numpy.isnan(za_grid)] = 90.0  # Replace nans as they break the interpolation

    #Supress print statements of the primary beam model functions
    sys.stdout = open(os.devnull, 'w')
    sky_grid = map_sky(my_map['skymap'], my_map['RA'], my_map['dec'], time, az_grid, za_grid)
    sys.stdout = sys.__stdout__

    t_ant = np.sum(phased_array_pattern * sky_grid * np.sin(za_grid)) / \
            np.sum(phased_array_pattern * np.sin(za_grid))
    return t_ant


def partial_convolve_sky_map(az_grid, za_grid, phased_array_pattern, freq, time):
    # Get Haslam and interpolate onto grid
    # taken from https://github.com/MWATelescope/mwa_pb/blob/master/mwa_pb/primarybeammap_tant.py
    # then edited to include pixel size
    my_map = get_Haslam(freq)
    #mask = numpy.isnan(za_grid)
    #za_grid[numpy.isnan(za_grid)] = 90.0  # Replace nans as they break the interpolation

    #Supress print statements of the primary beam model functions
    sys.stdout = open(os.devnull, 'w')
    sky_grid = map_sky(my_map['skymap'], my_map['RA'], my_map['dec'], time, az_grid, za_grid)
    sys.stdout = sys.__stdout__

    sum_B_T = np.sum(phased_array_pattern * sky_grid * np.sin(za_grid))
    sum_B   = np.sum(phased_array_pattern * np.sin(za_grid))
    return sum_B_T, sum_B