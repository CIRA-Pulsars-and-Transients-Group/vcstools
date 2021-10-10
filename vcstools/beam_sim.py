"""
The functions required to simulate the tied-array beam response of the MWA. 
All equations can be found in https://ui.adsabs.harvard.edu/abs/2018IAUS..337..378M/abstract
"""
# numerical and maths modules
import numpy as np
from astropy.constants import c

#utility and processing modules
import sys
import os
from itertools import chain
from astropy.io import fits
import matplotlib.pyplot as plt

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
    return np.radians(az_res * za_res) * np.sin(za)


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
    wl =  freq / c.value
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

    sum_B_T = np.sum(phased_array_pattern * sky_grid.flatten() * pixel_area)
    sum_B   = np.sum(phased_array_pattern * pixel_area)
    return sum_B_T, sum_B


def plot_vcsbeam_psf(psf_file, output_name="vcsbeam_psf.png"):
    input_array = np.loadtxt(psf_file, dtype=float)
    print(input_array.shape[0])
    ra    = input_array[:,0]
    dec   = input_array[:,1]
    power = input_array[:,2]
    ra.shape = dec.shape = power.shape = (int(np.sqrt(input_array.shape[0])),
                                          int(np.sqrt(input_array.shape[0])))
    fig, ax = plt.subplots()
    im = ax.pcolormesh(ra, dec, power)

    plt.xlabel(r"Right Ascension ($^{\circ}$)")
    plt.ylabel(r"Declination ($^{\circ}$)")
    plt.colorbar(im,#spacing='uniform', shrink = 0.65, #ticks=[2., 10., 20., 30., 40., 50.],
                 label="Normalised array factor power")
    plt.savefig(output_name)


def plot_track_beam_response(response_file, output_name="vcsbeam_response.png", time_max=-1):
    input_array = np.loadtxt(response_file, dtype=float)
    sec          = input_array[:,0]
    freq         = input_array[:,1]
    stokes_I     = input_array[:,5]
    array_factor = input_array[:,9]
    response = stokes_I * array_factor**2

    # Split into freq chunks
    nsec = max(sec)
    nfreq = input_array.shape[0] // nsec
    sec_map   = np.array(np.array_split(sec,      nfreq))[:,:time_max]
    freq_map  = np.array(np.array_split(freq,     nfreq))[:,:time_max]
    power_map = np.array(np.array_split(response, nfreq))[:,:time_max]

    # Sum over frequency
    power_time = np.average(power_map, axis=0)
    sec_time = sec_map[0]

    # Sum over time
    power_freq = np.average(power_map, axis=1)
    freq_freq = freq_map[:,0]

    size = 3
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(2*size,size))

    # Plot rack over time
    ax1.plot(sec_time, power_time)
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Array Factor")

    # Plot rack over time
    ax2.plot(freq_freq, power_freq)
    ax2.set_xlabel("Frequency (MHz)")
    ax2.set_ylabel("Array Factor")
    plt.tight_layout()
    plt.savefig(output_name)


def plot_pabeam(dat_file, output_name="pabeam_psf.png"):
    with open(dat_file) as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]
        nza = int(lines[7].split(" ")[-1])
        naz = int(lines[8].split(" ")[-1])
    input_array = np.loadtxt(dat_file, dtype=float)
    print(input_array.shape)
    za    = input_array[:,0]
    az    = input_array[:,1]
    power = input_array[:,2]
    print(za[0:10])
    print(az[0:10])
    #za.shape = az.shape = power.shape = (nza, naz)
    za.shape = az.shape = power.shape = (naz, nza)
    ax = plt.subplot(1, 1, 1, projection='polar')
    #ax = plt.subplot(1, 1, 1)
    #ax.set_rlim(1, 100)
    #ax.set_rscale('log')
    print(za.shape)
    print(za[0])
    #print(za[nza-1])
    print(az.shape)
    print(az[0])
    #print(az[nza-1])

    im = plt.pcolormesh(np.radians(az), za, power)
    plt.xlabel(r"Azimuth ($^{\circ}$)")
    plt.ylabel(r"Zenith ($^{\circ}$)")
    plt.colorbar(im, label="Normalised array factor power")
    plt.savefig(output_name, dpi=500)