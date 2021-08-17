#!/usr/bin/env python3
"""
This script will accept a source position and a UTC day in which to calculate the ephemerides for 24 hours.
This can be done for multiple observatories simultaneously.

Date: 2018/02/13
Author: Bradley Meyers
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, DateFormatter
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, get_body
from astropy import units as u
from astropy.time import Time, TimezoneInfo
from datetime import datetime
import logging

from vcstools.general_utils import setup_logger

logger = logging.getLogger(__name__)

# Radio telescope dictionary, formatted as:
#       Name:(latitude, longitude, elevation)
site_dict = dict(MWA=(-26.7033, 116.671, 377.827),
                 PKS=(-32.999944, 148.262306, 372),
                 ATCA=(-30.312778, 149.550278, 209),
                 ASKAP=(-26.696, 116.637, 372),
                 MOST=(-35.370707, 149.424658, 737),
                 GMRT=(19.096517, 74.049742, 656),
                 LWA=(34.07, -107.63, 2124),
                 GBT=(38.4331, -79.8397, 807),
                 VLA=(34.078749, -107.618283, 2124),
                 Arecibo=(18.34417, -66.75278, 323),
                 LOFAR=(52.90889, 6.86889, 6),
                 JodrellBank=(53.23625, -2.307139, 77),
                 Effelsberg=(50.5247, 6.8828, 346),
                 MeerKAT=(-30.721, 21.411, 1300))


# Nominal telescope elevation limits (in degrees)
elev_limits_dict = dict(MWA=10, PKS=30, GBT=5,
                        ATCA=12, ASKAP=15, GMRT=17,
                        VLA=20, Arecibo=70, Effelsberg=8,
                        MeerKAT=15)


class Observatory(object):

    """Class to hold information regarding the observatory and the ephemeris it observes for the given target."""

    def __init__(self, name, latitude=None, longitude=None, elevation=None):
        # check if the requested telescope is in the known site dictionary
        if name in site_dict.keys():
            logger.info("Found observatory coordinates for '{0}': {1}".format(name, site_dict[name]))
            self.name = name
            self.latitude = site_dict[name][0]
            self.longitude = site_dict[name][1]
            self.elevation = site_dict[name][2]
        else:
            # if not, check to make sure lat, long and elev were provided, otherwise abort
            logger.warning("Observatory '{0}' not found in internal list".format(name))
            if (latitude is None) or (longitude is None) or (elevation is None):
                logger.critical("Aborting here, there was no observatory found and one or more coordinates were "
                                "not provided.")
                sys.exit(1)
            else:
                logger.info("Using provided name and location")
                logger.debug("    telescope = {0}".format(name))
                logger.debug("    lon, lat, elev = ({0}, {1}, {2})".format(longitude, latitude, elevation))
                self.name = name
                self.latitude = latitude
                self.longitude = longitude
                self.elevation = elevation

        # determine the location of the site for later use
        self.location = EarthLocation(lat=self.latitude * u.deg,
                                      lon=self.longitude * u.deg,
                                      height=self.elevation * u.m)

        # set default values
        self.target = None
        self.altaz = None
        self.alt = None
        self.altmaxidx = None
        self.altmax = None
        self.tz = None
        self.maxtimeUTC = None
        self.maxtimeLST = None
        self.maxtimeLocalStr = None
        self.utcoffsetStr = None
        self.sun = None

        # check to see if there is a known elevation limit for the given site name
        if name in elev_limits_dict.keys():
            self.elev_limit = elev_limits_dict[name]
        else:
            self.elev_limit = None

    def compute_target_position(self, coords, times, tz):
        # compute the ephemeris of a given set of coordinates and set relevant attributes
        if coords and times:
            self.target = coords
            self.altaz = self.target.transform_to(AltAz(obstime=times, location=self.location))
            self.alt = self.altaz.alt.deg

            self.altmaxidx = np.argmax(self.alt)
            self.altmax = self.alt[self.altmaxidx]
            self.maxtimeUTC = times[self.altmaxidx]
            self.maxtimeLST = self.maxtimeUTC.sidereal_time('apparent', "{0}d".format(self.longitude))
            self.maxtimeLocalStr = str(self.maxtimeUTC.to_datetime(timezone=tz))[:-6]
            self.utcoffsetStr = str(self.maxtimeUTC.to_datetime(timezone=tz))[-6:]
        else:
            logger.critical("The RA, DEC or time was not provided and so we cannot calculate the altaz of the target. "
                            "Aborting.")
            sys.exit(1)

    def compute_sun_position(self, times):
        # compute the Sun's sky position over time
        sun_position = get_body('sun', times, self.location)
        self.sun = sun_position.transform_to(AltAz(obstime=times, location=self.location))


def plot_ephem(ax, times, obs, plot_sun=False, draw_peaks=False):
    """Given an axis object, the times (x) and observed ephemeris (y), plot the source track."""
    # plot the ephemeris over 24 hours
    eph = ax.plot_date(times, obs.alt, xdate=True, ls="-", lw=2, marker='', alpha=0.5, label="{0}".format(obs.name))

    # plot the horizon line
    ax.axhline(0, ls="-", lw=3, color='k')

    # if the telescope elevation limit is known, plot it too
    if obs.elev_limit is not None:
        ax.axhline(obs.elev_limit, lw=3, color=eph[0].get_color(), label="{0} elev. limit".format(obs.name))

    # plot a line marking the time of maximum elevation
    if draw_peaks:
        ax.axvline(date2num(obs.maxtimeUTC.datetime), ls="--", lw=2, color=eph[0].get_color())

    # if request, plot the Sun's ephemeris for the site
    if plot_sun:
        ax.plot_date(times[::30], obs.sun.alt, xdate=True, ls=":", lw=2, marker='', alpha=0.4, color=eph[0].get_color())

    return "{0}: {1} UTC{2}".format(obs.name, obs.maxtimeLocalStr, obs.utcoffsetStr)
    

def calculate_ephem(ra, dec, date, tzoffset=0, site=["MWA"], show_sun=False, draw_peaks=False):
    """Compute the ephemeris for a target on a given day at a target position."""
    # first, let's set up the times we want to evaluate the source position for
    year, month, day = date.split("/")
    hours = np.arange(0, 24, 0.01)
    hr = [int(i) for i in hours]
    mi = [int((hours[i]-hr[i]) * 60) for i in range(len(hours))]
    se = [0]*len(hours)
    t = [datetime(int(year), int(month), int(day), h, m, s) for h, m, s in zip(hr, mi, se)]

    times = Time(t, scale='utc', format='datetime')  # list of Time objects at which the ephemeris is computed
    tz = TimezoneInfo(utc_offset=tzoffset * u.hour)  # timezone object to convert times
    plttimes = date2num([t for t in times.datetime])  # times in plottable format

    # set up figure for plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)

    # site can be a list, so we need to create an Observatory for each
    # and calculate the ephemeris for each location
    site_max = []
    coords = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

    logger.info("Target: {0}".format(coords.to_string('hmsdms')))
    for s in site:
        o = Observatory(s)
        o.compute_target_position(coords, times, tz)
        logger.info("    max. target elevation  = {0}".format(o.altmax))
        logger.info("    time of max. elevation = {0} UTC{1}".format(o.maxtimeLocalStr, o.utcoffsetStr))

        if show_sun:
            o.compute_sun_position(times[::30])
        site_max.append(plot_ephem(ax, plttimes, o, plot_sun=show_sun, draw_peaks=draw_peaks))

    # plotting aesthetics
    ax.set_xlim(min(plttimes), max(plttimes))
    if ax.get_ylim()[1] > 90:
        ax.set_ylim(-10, 91)
    else:
        ax.set_ylim(-10, None)

    ys = [ax.get_ylim()[0]] * len(plttimes)
    ax.fill_between(plttimes, ys, interpolate=True, color='gray')  # fill from axis lower limit to 0

    title_str = "Target :: {0}\n{1}".format(coords.to_string('hmsdms'), "\n".join(site_max))
    ax.set_title(title_str)
    ax.set_xlabel("Time (UTC)", fontsize=14)
    ax.set_ylabel("Elevation  [deg]", fontsize=14)

    ax.xaxis.set_major_formatter(DateFormatter("%m/%d %H:%M"))
    plt.xticks(rotation=30, ha="right")
    ax.legend()
    ax.grid()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Determine ephemeris for source object over 24hrs.")
    parser.add_argument('-r', '--ra', type=str, metavar="RAJ2000",
                        help="RAJ2000 coordinate of target source (hh:mm:ss.ss)", default=None)
    parser.add_argument('-d', '--dec', type=str, metavar="DECJ2000",
                        help="DECJ2000 coordinate of target source (dd:mm:ss.ss)", default=None)
    parser.add_argument('-u', '--utcdate', type=str,   metavar="date", help="Desired ephemeris UTC date (YYYY/MM/DD)",
                        default=None)
    parser.add_argument('-o', '--utcoff',  type=float, metavar="offset", help="Hour offset from UTC [default = 0]",
                        default=0)
    parser.add_argument('-s', '--site',    type=str,   metavar="name", nargs='+', choices=site_dict.keys(),
                        help="Common radio telescope sites to use as observer position. "
                             "Choose from: {0} [default = MWA]".format(tuple(site_dict.keys())), default=["MWA"])
    parser.add_argument('--sun', action='store_true', help="Plot the Sun's ephemeris for each location")
    parser.add_argument('--draw_peak', action='store_true',
                        help="Draw a vertical line representing the peak elevation time")
    # TODO: need to allow any number of manually defined observing positions, but for now just use those in the
    #  site_dict
    #parser.add_argument('--observer', type=float, nargs=3, metavar=("lat", "lon", "elev"),
    # help="Latitude (deg), longitude (deg) and elevation (m) of observer. No default.",default=(None,None,None))
    # TODO: is centering really necessary?
    parser.add_argument('-V', '--version', action='store_true', help="Print version and quit")
    args = parser.parse_args()

    if args.version:
        try:
            import version
            logger.debug(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            logger.error("Couldn't import version.py - have you installed vcstools?")
            logger.error("ImportError: {0}".format(ie))
            sys.exit(0)

    # set up the logger for stand-alone execution
    logger = setup_logger(logger, log_level=loglevels[args.loglvl])

    if (args.ra is None) or (args.dec is None) or (args.utcdate is None):
        logger.error("You must specify a RA, Dec and UTC date to plot the ephemeris")

    if args.site is None:
        logger.warning("You didn't provide a site name, assuming MWA...")
        args.site = ["MWA"]

    calculate_ephem(args.ra, args.dec, args.utcdate,
                    tzoffset=args.utcoff, site=args.site, show_sun=args.sun, draw_peaks=args.draw_peak)
