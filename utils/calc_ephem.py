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
from astropy.coordinates import EarthLocation,SkyCoord,AltAz
from astropy import units as u
from astropy.time import Time,TimezoneInfo
from datetime import datetime, timedelta
import sys
import argparse

def calculate_ephem(ra,dec,date,tzoffset,center, lat,lon,elev):
    
    location = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=elev*u.m)

    target = SkyCoord(ra,dec,unit=(u.hourangle,u.deg))
    
    year,month,day = date.split("/")
    hours = np.arange(0,24,0.01)
    hr = [int(i) for i in hours]
    mi = [int((hours[i]-hr[i])*60) for i in range(len(hours))]
    se = [0]*len(hours)
    t = [datetime(int(year),int(month),int(day),h,m,s) for h,m,s in zip(hr,mi,se)]
 
    times = Time(t,scale='utc',format='datetime')
    altaz = target.transform_to(AltAz(obstime=times,location=location))
    alt = altaz.alt.deg

    maxidx = alt.argmax() #np.where(alt==alt.max())[0][0]
    maxtime = times[maxidx]
    if center:
        # shifting things to the centre of the plot
        dt = (times[-1]-times[0]) / len(times)
        times = (maxtime - 12*u.hour) + dt * np.arange(len(times))
        altaz = target.transform_to(AltAz(obstime=times,location=location))
        alt = altaz.alt.deg
        times += tzoffset*u.hour
        # for things to be unaltered downstream
        maxidx = 1200

    # converting the times to something plottable
    times = dates.date2num([t for t in times.datetime])
    tz = TimezoneInfo(utc_offset=tzoffset*u.hour)
    lst = str(maxtime.to_datetime(timezone=tz))
    utcoff = lst[-6:]
    lst = lst[:-6]
    print "Local time (UTC{0}) of maximum elevation: {1}".format(utcoff,lst)

    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111)
    ax.plot_date(times,alt,color='r',lw=2,alpha=0.6)
    ax.axhline(0,ls="--",color='k')
    ax.axvline(times[maxidx],ls="--",color='r',lw=2)
    #ax.set_xlim(times[0], times[-1])
    ylims = ax.get_ylim()
    if ylims[1] > 90:
	ax.set_ylim(ylims[0],90)
	
    zero = np.zeros(len(alt))
    ax.fill_between(times,[ax.get_ylim()[0]]*len(hours),interpolate=True,color='gray')
    ax.set_title("Source: {0} {1}\n site coords: lon={2:.3f}d lat={3:.3f}d elev.={4:.2f}m\n max. elev: {5:.2f}d @  {6} UTC{7}".format(ra,dec,lat,lon,elev,alt.max(),lst,utcoff))
    ax.set_xlabel("Time (UTC{0})  ".format(utcoff))
    ax.set_ylabel("Elevation  [deg]")
    plt.show()




#radio telescope sites
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



parser = argparse.ArgumentParser(description="Determine ephemeris for source object over 24hrs.")
parser.add_argument('--ra',type=str,metavar="RAJ2000",help="RAJ2000 coordinate of target source (hh:mm:ss.ss)",required=True)
parser.add_argument('--dec',type=str,metavar="DECJ2000",help="DECJ2000 coordinate of target source (dd:mm:ss.ss)",required=True)
parser.add_argument('--utcdate',type=str,metavar="date",help="Desired ephemeris UTC date (YYYY/MM/DD)",required=True)
parser.add_argument('--utcoff',type=float,metavar="offset",help="Hour offset from UTC [default = 0]",default=0)
parser.add_argument('--site',type=str,metavar="name",nargs=1,choices=site_dict.keys(),help="Common radio telescope sites to use as observer position. Choose from: {0}. No default.".format(site_dict.keys()),default=None)
parser.add_argument('--observer',type=float,nargs=3,metavar=("lat", "lon", "elev"),help="Latitude (deg), longitude (deg) and elevation (m) of observer. No default.",default=(None,None,None))
parser.add_argument('-c', '--center', action='store_true', help='if set will center the plot the time of maximum elevation')
args = parser.parse_args()

# chekc to see if multiple sites were given
if args.site:
	if args.observer:
		print "--site takes priority over --observer argument"
	observatory = site_dict[args.site[0]]
elif args.observer:
	observatory = args.observer
else:
	print "Somehow managed to get by with providing an observer location?"

calculate_ephem(args.ra,args.dec,args.utcdate,args.utcoff, args.center, *observatory)
