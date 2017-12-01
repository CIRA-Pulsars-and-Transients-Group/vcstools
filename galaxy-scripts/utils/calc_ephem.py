#!/usr/bin/env python

"""
Updated to be a little more general and not require LSL libraries...
I suppose the next step would be to move over to using astropy given that pyephem is depricated now...

Date: 2016/06/13
Author: Bradley Meyers
"""

#import ephem
#from lsl.common.stations import lwa1
#from psr_constants import *
#from math import pi
import matplotlib.pyplot as plt
from matplotlib import dates
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord, AltAz
from astropy import units as u
from astropy.time import Time, TimezoneInfo
from astropy.utils import iers
from datetime import datetime, timedelta
import sys
import argparse

# get updated IERS tables
iers.IERS.iers_table = iers.IERS_A.open(iers.IERS_A_URL)


#radio telescope sites: Name, latitude, longitude, elevation
site_dict = {'MWA':(-26.7033,116.671,377.827),\
		'PKS':(-32.999944,148.262306,372),\
		'ATCA':(-30.312778,149.550278,209),\
		'ASKAP':(-26.696,116.637,372),\
		'MOST':(-35.370707,149.424658,737),\
		'GMRT':(19.096517,74.049742,656),\
		'LWA':(34.07,-107.63,2124),\
		'GBT':(38.4331,-79.8397,807),\
		'VLA':(34.078749,-107.618283,2124),\
		'Arecibo':(18.34417,-66.75278,323),\
		'LOFAR':(52.90889,6.86889,6),\
		'JodrellBank':(53.23625,-2.307139,77),\
		'Effelsberg':(50.5247,6.8828,346)}


class Observatory(object):

    def __init__(self, name, latitude=None, longitude=None, elevation=None):
        if name in site_dict.keys():
            print "Found observatory coordinate"
            self.name = name
            self.latitude = site_dict[name][0]
            self.longitude = site_dict[name][1]
            self.elevation = site_dict[name][2]
        else:
            print "Observatory not found in internal list"
            if (latitude is None) or (longitude is None) or (elevation is None):
                print "Aborting here, there was no observatory found and one or more coordinates were not provided."
                sys.exit(1)
            else:
                self.name = name
                self.latitude = latitude
                self.longitude = longitude
                self.elevation = elevation    
            
        self.location = EarthLocation(lat=self.latitude * u.deg, \
                                      lon=self.longitude * u.deg, \
                                      height=self.elevation * u.m)
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


    def compute_target_position(self, ra, dec, times, tzoffset):
        if ra and dec and times:
            self.target = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            self.altaz = self.target.transform_to(AltAz(obstime=times, location=self.location))
            self.alt = self.altaz.alt.deg

            self.altmaxidx = np.argmax(self.alt)
            self.altmax = self.alt[self.altmaxidx]
            self.maxtimeUTC = times[self.altmaxidx]
            self.maxtimeLST = self.maxtimeUTC.sidereal_time(('apparent'), "{0}d".format(self.longitude))
            tmp = str(self.maxtimeUTC.to_datetime(timezone=tzoffset))
            self.maxtimeLocalStr = tmp[:-6] 
            self.utcoffsetStr = tmp[-6:]
        else:
            print "The RA, DEC or time was not provided and so we cannot calculate the altaz of the target. Aborting."
            sys.exit(1)
        

def plot_ephem(ax, times, obs):
    p = ax.plot_date(times, obs.alt, marker='', ls="-", lw=2, alpha=0.6, label=obs.name)
    ax.axhline(0, ls="--", color='k')
    ax.axvline(dates.date2num(obs.maxtimeUTC.datetime), ls="--", lw=2, color=p[0].get_color())
    if ax.get_ylim()[1] > 90:
	    ax.set_ylim(None, 90)
	
    ylims = ax.get_ylim()
    
    #zero = np.zeros(len(obs.alt))
    ax.fill_between(times, [ax.get_ylim()[0]]*len(times), interpolate=True, color='gray')

    coords = obs.target.to_string('hmsdms')
    ax.set_title("Source: {0}".format(coords))
    ax.set_xlabel("Time (UTC)")
    ax.set_ylabel("Elevation  [deg]")


def calculate_ephem(ra, dec, date, tzoffset, site, center):
    
    # first, let's set up the times we want to evaluate the source position for
    year, month, day = date.split("/")
    hours = np.arange(0,24,0.01)
    hr = [int(i) for i in hours]
    mi = [int((hours[i]-hr[i])*60) for i in range(len(hours))]
    se = [0]*len(hours)
    t = [datetime(int(year), int(month), int(day), h, m, s) for h,m,s in zip(hr,mi,se)]

    times = Time(t, scale='utc', format='datetime') # list of Time objects at which the ephemeris is computed
    tz = TimezoneInfo(utc_offset=tzoffset * u.hour) # timezone object to convert times
    plttimes = dates.date2num([t for t in times.datetime]) # times in plottable format

#    # TODO: need to re-implement this in a consistent way...   
#    if center:
#        # shifting things to the centre of the plot
#        dt = (times[-1]-times[0]) / len(times)
#        times = (maxtime - 12*u.hour) + dt * np.arange(len(times))
#        altaz = target.transform_to(AltAz(obstime=times,location=location))
#        alt = altaz.alt.deg
#        #times += tzoffset*u.hour
#        # for things to be unaltered downstream
#        maxidx = 1200


    # set up figure for plot
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111)

    # site can be a list, so we need to create an Observatory for each
    for s in site:
        o = Observatory(s)
        o.compute_target_position(ra, dec, times, tz)
        print "For {0}:".format(o.name)
        print "Maximum elevation = {0:.2} deg".format(o.altmax)
        print "Time of max. elevation:"
        print "  UTC = {0}".format(o.maxtimeUTC)
        print "local = {0}".format(o.maxtimeLocalStr)
        print "  LST = {0}".format(o.maxtimeLST)
        print "-------------------------------"
        plot_ephem(ax, plttimes, o)
    
    plt.legend(loc=2)
    plt.tight_layout()
    plt.show()
    
       








parser = argparse.ArgumentParser(description="Determine ephemeris for source object over 24hrs.")
parser.add_argument('--ra', type=str, metavar="RAJ2000", help="RAJ2000 coordinate of target source (hh:mm:ss.ss)",required=True)
parser.add_argument('--dec', type=str, metavar="DECJ2000", help="DECJ2000 coordinate of target source (dd:mm:ss.ss)",required=True)
parser.add_argument('--utcdate', type=str, metavar="date", help="Desired ephemeris UTC date (YYYY/MM/DD)",required=True)
parser.add_argument('--utcoff', type=float, metavar="offset", help="Hour offset from UTC [default = 0]",default=0)
parser.add_argument('--site', type=str, metavar="name", nargs='+', choices=site_dict.keys(), help="Common radio telescope sites to use as observer position. Choose from: {0}. Or use --observer option. No default.".format(site_dict.keys()),default=None)
# TODO: need to allow any number of manually defined observing positions, but for now just use those in the site_dict
#parser.add_argument('--observer', type=float, nargs=3, metavar=("lat", "lon", "elev"), help="Latitude (deg), longitude (deg) and elevation (m) of observer. No default.",default=(None,None,None))
parser.add_argument('-c', '--center', action='store_true', help='Center the time of maximum elevation on the plot')
args = parser.parse_args()

calculate_ephem(args.ra, args.dec, args.utcdate, args.utcoff, args.site, args.center)
