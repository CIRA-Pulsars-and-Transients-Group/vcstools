#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess

import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.table import Table
from astropy.time import Time

import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches

from mwa_pb import primary_beam
from vcstools.metadb_utils import getmeta, get_common_obs_metadata
from find_pulsar_in_obs import get_beam_power_over_time

import logging
logger = logging.getLogger(__name__)

def get_beam_power(obsid_data, sources, beam_model="analytic", centeronly=True):

    obsid, _, _, duration, delays, centrefreq, channels = obsid_data

    # Figure out GPS times to evaluate the beam in order to map the drift scan
    midtimes = np.array([float(obsid) + 0.5 * duration]) # although only includes one element, could in future be extended
    Ntimes = len(midtimes)


    # Make an EarthLocation for the MWA
    MWA_LAT = -26.7033   # degrees
    MWA_LON = 116.671    # degrees
    MWA_HEIGHT = 377.827 # metres
    mwa_location = EarthLocation(lat=MWA_LAT*u.deg, lon=MWA_LON*u.deg, height=MWA_HEIGHT*u.m)

    # Pre-allocate arrays for the beam patterns
    if not centeronly:
        PowersX = np.zeros((len(sources), Ntimes, len(channels)))
        PowersY = np.zeros((len(sources), Ntimes, len(channels)))
        frequencies = np.array(channels) * 1.28e6
    else:
        PowersX = np.zeros((len(sources), Ntimes, 1))
        PowersY = np.zeros((len(sources), Ntimes, 1))
        frequencies = [centrefreq]

    # Create arrays of the RAs and DECs to be sampled (already in degrees)
    #print "Constructing RA and DEC arrays for beam model..."
    logger.info("Constructing RA and DEC arrays for beam model...")
    RAs = np.array([x[0] for x in sources])
    Decs = np.array([x[1] for x in sources])
    if len(RAs) == 0:
        sys.stderr.write('Must supply >=1 source positions\n')
        return None
    if not len(RAs) == len(Decs):
        sys.stderr.write('Must supply equal numbers of RAs and Decs\n')
        return None

    # Turn RAs and DECs into SkyCoord objects that are to be fed to Alt/Az conversion
    coords = SkyCoord(RAs, Decs, unit=(u.deg, u.deg))

    # Make times to feed to Alt/Az conversion
    obstimes = Time(midtimes, format='gps', scale='utc')

    #print "Converting RA and DEC to Alt/Az and computing beam pattern..."
    logger.info("Converting RA and DEC to Alt/Az and computing beam pattern...")
    for t, time in enumerate(obstimes):
        # Convert to Alt/Az given MWA position and observing time
        altaz = coords.transform_to(AltAz(obstime=time, location=mwa_location))

        # Change Altitude to Zenith Angle
        theta = np.pi / 2 - altaz.alt.rad
        phi = altaz.az.rad

        # Calculate beam pattern for each frequency, and store the results in a ndarray
        for f, freq in enumerate(frequencies):
            if beam_model == "analytic":
                rX, rY = primary_beam.MWA_Tile_analytic(theta, phi, freq=freq, delays=delays, zenithnorm=True, power=True)
            elif beam_model == "FEE":
                rX, rY = primary_beam.MWA_Tile_full_EE(theta, phi, freq=freq, delays=delays, zenithnorm=True, power=True)
            else:
                #print "Unrecognised beam model '{0}'. Defaulting to 'analytic'."
                logger.warning("Unrecognised beam model '{0}'. Defaulting to 'analytic'.")
                rX, rY = primary_beam.MWA_Tile_analytic(theta, phi, freq=freq, delays=delays, zenithnorm=True, power=True)

            PowersX[:, t, f] = rX
            PowersY[:, t, f] = rY

    # Convert X and Y powers into total intensity
    Powers = 0.5 * (PowersX + PowersY)

    return Powers


def read_data(fname, delim=",", read_format="csv", coords=True):
    #print "reading from", fname
    logger.info("reading from {0}".format(fname))
    # Read the data file as an Astropy Table
    tab = Table.read(fname, delimiter=delim, format=read_format)

    if coords:
        # If we just want the coordinates, create a list of SkyCoords and return it
        tab_coords = SkyCoord(tab["RAJ"], tab["DECJ"], unit=(u.hourangle, u.deg))
        return tab_coords
    else:
        # Otherwise, return the whole table
        return tab



def plotSkyMap(obsfile, targetfile, oname, show_psrcat=False, show_mwa_sky=False, show_mwa_unique=False):


    fig = plt.figure()
    ax = fig.add_subplot(111, projection="mollweide")

    mwa_dec_lim = 30
    mwa_only_dec_lim = -50

    if show_mwa_sky:
        # Make a patch that is transformable through axis projection
        # and highlight the visible sky for the MWA
        Path = mpath.Path
        ra_start = -np.pi
        ra_end = np.pi
        dec_start = -np.pi / 2
        dec_end = np.radians(mwa_dec_lim)
        path_data = [
                    (Path.MOVETO, (ra_start, dec_start)),
                    (Path.LINETO, (ra_start, dec_end)),
                    (Path.LINETO, (ra_end, dec_end)),
                    (Path.LINETO, (ra_end, dec_start)),
                    (Path.CLOSEPOLY, (ra_end, dec_start)),
                    ]
        codes, verts = zip(*path_data)
        path = mpath.Path(verts, codes)
        patch = mpatches.PathPatch(path, lw=0, ec="lightskyblue", fc="lightskyblue", label="All MWA sky")
        ax.add_patch(patch)

    if show_mwa_unique:
        # Make a patch that is transformable through axis projection
        # and highlight the part of the sky ONLY visible to the MWA
        ra_start = -np.pi
        ra_end = np.pi
        dec_start = -np.pi / 2
        dec_end = np.radians(mwa_only_dec_lim)
        path_data = [
                    (Path.MOVETO, (ra_start, dec_start)),
                    (Path.LINETO, (ra_start, dec_end)),
                    (Path.LINETO, (ra_end, dec_end)),
                    (Path.LINETO, (ra_end, dec_start)),
                    (Path.CLOSEPOLY, (ra_end, dec_start)),
                    ]
        codes, verts = zip(*path_data)
        path = mpath.Path(verts, codes)
        patch = mpatches.PathPatch(path, lw=0, ec="lightgreen", fc="lightgreen", label="Unique MWA sky")
        ax.add_patch(patch)

    if show_psrcat:
        # Retrieve the local installed PSRCAT catalogue and plot those pulsar positions
        cmd = 'psrcat -o short_csv -nocand -nonumber -c "Jname RAJ DECJ" | sed "2d" > psrcat.csv'
        subprocess.call(cmd, shell=True)
        psrcat_coords = read_data("psrcat.csv", delim=";")
        os.remove("psrcat.csv")

        # Create a mask for pulsar in and out of the declination limit of the MWA
        maskGood = psrcat_coords.dec.wrap_at(180*u.deg).deg < mwa_dec_lim
        maskBad = psrcat_coords.dec.wrap_at(180*u.deg).deg >= mwa_dec_lim

        psrcat_ra_good = -psrcat_coords.ra.wrap_at(180*u.deg).rad[maskGood] # negative because RA increases to the West
        psrcat_dec_good = psrcat_coords.dec.wrap_at(180*u.deg).rad[maskGood]

        psrcat_ra_bad = -psrcat_coords.ra.wrap_at(180*u.deg).rad[maskBad] # negative because RA increases to the West
        psrcat_dec_bad = psrcat_coords.dec.wrap_at(180*u.deg).rad[maskBad]

        # Now plot the pulsar locations
        ax.scatter(psrcat_ra_good, psrcat_dec_good, 0.01, marker="x", color="0.4", zorder=1.4)
        ax.scatter(psrcat_ra_bad, psrcat_dec_bad, 0.01, marker="x", color="0.8", zorder=1.4)


    # Calculate beam patterns and plot contours
    levels = np.arange(0.25, 1., 0.05)
    cmap = plt.get_cmap("cubehelix_r")

    obsids = read_data(obsfile, coords=False)["OBSID"]
    for obsid in obsids:
        #print "Accessing database for observation: {0}".format(obsid)
        logger.info("Accessing database for observation: {0}".format(obsid))

        # TODO: I think this process is now implemented in a function in mwa_metadb_utils, need to double check
        beam_meta_data = getmeta(service='obs', params={'obs_id': obsid})

        ra = beam_meta_data[u'metadata'][u'ra_pointing']
        dec = beam_meta_data[u'metadata'][u'dec_pointing']
        duration = beam_meta_data[u'stoptime'] - beam_meta_data[u'starttime'] #gps time
        delays = beam_meta_data[u'rfstreams'][u'0'][u'xdelays']

        minfreq = float(min(beam_meta_data[u'rfstreams'][u"0"][u'frequencies']))
        maxfreq = float(max(beam_meta_data[u'rfstreams'][u"0"][u'frequencies']))
        centrefreq = 1.28e6 * (minfreq + (maxfreq-minfreq)/2) #in MHz
        channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']

        beam_meta_data = get_common_obs_metadata(obsid)

        # Create a meshgrid over which to iterate
        Dec = [] ; RA = []
        map_dec_range = np.arange(-89, 90, 3)
        map_ra_range = np.arange(0, 360, 3)
        for i in map_dec_range:
            for j in map_ra_range:
                Dec.append(i)
                RA.append(j)
        RADecIndex = np.arange(len(RA))
        names_ra_dec = np.column_stack((RADecIndex, RA, Dec))

        #print "Creating beam patterns..."
        time_intervals = 600 # seconds

        powout = get_beam_power_over_time(names_ra_dec, common_metadata=beam_meta_data, dt=time_intervals, degrees=True)
        z=[] ; x=[] ; y=[]
        for c in range(len(RA)):
            temppower = 0.
            for t in range(powout.shape[1]):
                power_ra = powout[c,t,0]
                if power_ra > temppower:
                    temppower = power_ra
            z.append(temppower)
            if RA[c] > 180:
                x.append(-RA[c]/180.*np.pi+2*np.pi)
            else:
                x.append(-RA[c]/180.*np.pi)

            y.append(Dec[c]/180.*np.pi)

        nx=np.array(x) ; ny=np.array(y); nz=np.array(z)

        #print "Plotting beam pattern contours..."
        logger.info("Plotting beam pattern contours...")
        # Set vmin and vmax to ensure the beam patterns are on the same color scale
        # and plot the beam pattern contours on the map

        #c = ax.tricontour(wrapped_ra, wrapped_dec, final_pattern_power, levels=levels, cmap=cmap, vmin=levels.min(), vmax=levels.max(), zorder=1.3)
        c = ax.tricontour(nx, ny, nz, levels=levels, cmap=cmap, vmin=levels.min(), vmax=levels.max(), zorder=1.3)

    # Make a figure color scale based on the contour sets (which all have the same max/min values
    fig.colorbar(c, fraction=0.02, pad=0.03, label="Zenith normalised power")


    # Plot the target positions, using Astropy to wrap at correct RA
    target_coords = read_data(targetfile)
    wrapped_target_ra = -target_coords.ra.wrap_at(180*u.deg).rad # negative because RA increases to the West
    wrapped_target_dec = target_coords.dec.wrap_at(180*u.deg).rad
    #print "Plotting target source positions..."
    logger.info("Plotting target source positions...")
    ax.scatter(wrapped_target_ra, wrapped_target_dec, 10, marker="x", color="red", zorder=1.6, label="Target sources")

    xtick_labels = ["10h", "8h", "6h", "4h", "2h", "0h", "22h", "20h", "18h", "16h", "14h"]
    ax.set_xticklabels(xtick_labels, color="0.2")
    ax.set_xlabel("Right Ascension")
    ax.set_ylabel("Declination")
    ax.grid(True, color="k", lw=0.5, ls=":")

    # Place upper-right corner of legend at specified Axis coordinates
    ax.legend(loc="upper right", bbox_to_anchor=(1.02, 0.08), numpoints=1, borderaxespad=0., fontsize=6)


    plt.savefig(oname, format="eps", bbox_inches="tight")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--obsfile", type=str, help="File containing a list of observation IDs, one ID per line, column labelled 'OBSID'")

    parser.add_argument("-t", "--targetfile", type=str, help="File containing a list of target sources, with positions in columns labelled 'RAJ' and 'DECJ'")

    parser.add_argument("--oname", type=str, help="Output EPS file name [default: skymap.eps]", default="skymap.eps")

    parser.add_argument("--show_psrcat", action="store_true", help="Whether to show background pulsars on map (assumes psrcat is on PATH)")
    parser.add_argument("--show_mwa_sky", action="store_true", help="Whether to split the sky based on MWA visibility")
    parser.add_argument("--show_mwa_unique", action="store_true", help="Show the portion of sky ONLY accessible by the MWA")

    args = parser.parse_args()

    plotSkyMap(args.obsfile, args.targetfile, args.oname, args.show_psrcat, args.show_mwa_sky, args.show_mwa_unique)
