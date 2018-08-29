#!/usr/bin/env python 

import os
import sys
import argparse
import subprocess
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.table import Table
from astropy.time import Time

from mwapy.pb import primary_beam
#import ephem
#from mwapy import ephem_utils,metadata

import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
#import matplotlib.tri as tri
#import matplotlib.cm as cm

from mwa_metadb_utils import getmeta


def get_beam_power(obsid_data, sources, dt=296, beam_model="analytic", centeronly=True):
    
    obsid, ra, dec, duration, delays, centrefreq, channels = obsid_data
    
    # Figure out GPS times to evaluate the beam in order to map the drift scan
    starttimes = np.arange(0, duration, duration)
    stoptimes = starttimes + dt
    stoptimes[stoptimes > duration] = duration
    Ntimes = len(starttimes)
    midtimes = float(obsid) + 0.5 * (starttimes + stoptimes)


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
    print "Constructing RA and DEC arrays for beam model..."
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

    print "Converting RA and DEC to Alt/Az and computing beam pattern..."
    for t, time in enumerate(obstimes):
        # Convert to Alt/Az given MWA position and observing time
        altaz = coords.transform_to(AltAz(obstime=time, location=mwa_location))

        # Change Altitude to Zenith Angle
        theta = np.pi / 2 - altaz.alt.rad
        phi = altaz.az.rad
        
        # Calcualte beam pattern for each frequency, and store the results in a ndarray
        for f, freq in enumerate(frequencies):
            if beam_model == "analytic":
                rX, rY = primary_beam.MWA_Tile_analytic(theta, phi, freq=freq, delays=delays, zenithnorm=True, power=True)
            elif beam_model == "FEE":
                rX, rY = primary_beam.MWA_Tile_full_EE(theta, phi, freq=freq, delays=delays, zenithnorm=True, power=True)
            else:
                print "Unrecognised beam model '{0}'. Defaulting to 'analytic'."
                rX, rY = primary_beam.MWA_Tile_analytic(theta, phi, freq=freq, delays=delays, zenithnorm=True, power=True)
                
            PowersX[:, t, f] = rX
            PowersY[:, t, f] = rY
    
    # Convert X and Y powers into total intensity
    Powers = 0.5 * (PowersX + PowersY)
                 
    return Powers


def read_data(fname, delim=",", format="csv", coords=True):
    print "reading from", fname
    # Read the data file as an Astropy Table
    tab = Table.read(fname, delimiter=delim, format=format)

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
        
        psrcat_ra_good = psrcat_coords.ra.wrap_at(180*u.deg).rad[maskGood]
        psrcat_dec_good = psrcat_coords.dec.wrap_at(180*u.deg).rad[maskGood]
        
        psrcat_ra_bad = psrcat_coords.ra.wrap_at(180*u.deg).rad[maskBad]
        psrcat_dec_bad = psrcat_coords.dec.wrap_at(180*u.deg).rad[maskBad]

        # Now plot the pulsar locations
        ax.scatter(-psrcat_ra_good, psrcat_dec_good, 0.01, marker="x", color="0.4", zorder=1.4)
        ax.scatter(-psrcat_ra_bad, psrcat_dec_bad, 0.01, marker="x", color="0.8", zorder=1.4)




    
    # Calculate beam patterns and plot contours
    levels = np.arange(0.25, 1., 0.05)
    cmap = plt.get_cmap("cubehelix_r")

    obsids = read_data(obsfile, coords=False)["OBSID"]
    for obsid in obsids:
        print "Accessing database for observation: {0}".format(obsid)

        # TODO: I think this process is now implemented in a function in mwa_metadb_utils, need to double check
        beam_meta_data = getmeta(service='obs', params={'obs_id': obsid})
            
        ra = beam_meta_data[u'metadata'][u'ra_pointing']
        dec = beam_meta_data[u'metadata'][u'dec_pointing']
        time = beam_meta_data[u'stoptime'] - beam_meta_data[u'starttime'] #gps time
        delays = beam_meta_data[u'rfstreams'][u'0'][u'xdelays']

        minfreq = float(min(beam_meta_data[u'rfstreams'][u"0"][u'frequencies']))
        maxfreq = float(max(beam_meta_data[u'rfstreams'][u"0"][u'frequencies']))
        centrefreq = 1.28e6 * (minfreq + (maxfreq-minfreq)/2) #in MHz
        channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']

        cord = [obsid, ra, dec, time, delays, centrefreq, channels]


        RA = []
        Dec = []
        z = []
        x = []
        y = []
        
        # TODO: This next section, until the tricontour command, needs documentation...
        for i in range(-87, 88, 3):
            for j in range(0, 361, 3):
                Dec.append(i)
                RA.append(j)

        print "Creating beam patterns..."
        powout = get_beam_power(cord, zip(RA,Dec), dt=600)

        # TODO: this currently just stretches the beam pattern as calculated in the middle of the observation from the start to finish times
        # We need to double check that this is doing the right thing...
        for i in range(len(RA)):
            temppower = powout[i, 0, 0]
            for t in range(0, (time + 361)/720):
                if i % 121 >= t:
                    power_ra = powout[i-t, 0, 0] #ra+t*15./3600 3deg
                else:
                    power_ra = powout[i+121-t, 0, 0]
                if power_ra > temppower:
                    temppower = power_ra
            z.append(temppower)
            if RA[i] > 180:
                x.append(np.radians(-RA[i]) + 2 * np.pi)
            else:
                x.append(np.radians(-RA[i]))

            y.append(np.radians(Dec[i]))


        print "Plotting beam pattern contours..."
        # Set vmin and vmax to ensure the beam patterns are on the same color scale
        # and plot the beam pattern contours on the map
        c = ax.tricontour(x, y, z, levels=levels, cmap=cmap, vmin=levels.min(), vmax=levels.max(), zorder=1.3)

    # Make a figure color scale based on the contour sets (which all have the same max/min values
    fig.colorbar(c, fraction=0.02, pad=0.03, label="Zenith normalised power")


    # Plot the target positions
    target_coords = read_data(targetfile)
    print "Plotting target source positions..."
    ax.scatter(-target_coords.ra.wrap_at(180*u.deg).rad, target_coords.dec.wrap_at(180*u.deg).rad, 10, marker="x", color="red", zorder=1.6, label="Target sources")
    
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
