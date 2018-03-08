#!/usr/bin/env python

import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord, AltAz
from astropy import units as u
from astropy.time import Time

import matplotlib as mpl
#mpl.rc("text", usetex=True)
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from mwapy.pb import primary_beam as pb

import urllib
import urllib2
import json
import sys
import argparse


def getmeta(service='obs', params=None):
    """
    getmeta(service='obs', params=None)
    
    Given a JSON web service ('obs', find, or 'con') and a set of parameters as
    a Python dictionary, return the RA and Dec in degrees from the Python dictionary.
    """
    BASEURL = 'http://mwa-metadata01.pawsey.org.au/metadata/'
    if params:
        data = urllib.urlencode(params)  # Turn the dictionary into a string with encoded 'name=value' pairs
    else:
        data = ''
    #Validate the service name
    if service.strip().lower() in ['obs', 'find', 'con']:
        service = service.strip().lower()
    else:
        print "invalid service name: %s" % service
        return
    #Get the data
    try:
        result = json.load(urllib2.urlopen(BASEURL + service + '?' + data))
    except urllib2.HTTPError as error:
        print "HTTP error from server: code=%d, response:\n %s" % (error.code, error.read())
        return
    except urllib2.URLError as error:
        print "URL or network error: %s" % error.reason
        return
    #Return the result dictionary
    return result



def get_obs_metadata(obs):
    from mwapy.pb import mwa_db_query as mwa_dbQ

    beam_meta_data = getmeta(service='obs', params={'obs_id':obs})
    channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
    freqs = [float(c)*1.28 for c in channels]
    xdelays = beam_meta_data[u'rfstreams'][u"0"][u'xdelays']
    pointing_AZ, pointing_EL, pointing_ZA = mwa_dbQ.get_beam_pointing(obs)
    
    return {"channels":channels,
            "frequencies":freqs,
            "delays":xdelays,
            "az":pointing_AZ,
            "za":pointing_ZA
            }


def compute_target_position(ra, dec, time):
    MWA_LAT  = -26.7033
    MWA_LON  = 116.671
    MWA_ELEV = 377.827
    MWA_LOCATION = EarthLocation(lat=MWA_LAT * u.deg,
                                 lon=MWA_LON * u.deg,
                                 height=MWA_ELEV * u.m)
    
    coords = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    
    return coords.transform_to(AltAz(obstime=time, location=MWA_LOCATION))


def log_normalise(data, vmin, vmax):
    """ This is almost directly copied from matplotlib's color.py """
    result = np.ma.masked_less_equal(data, 0, copy=False)
    if vmin > vmax:
        raise ValueError("minvalue must be less than or equal to maxvalue")
    elif vmin <= 0:
        raise ValueError("values must all be positive")
    elif vmin == vmax:
        result.fill(0)
    else:
        mask = np.ma.getmask(result)
        result = np.ma.array(np.clip(result.filled(vmax), vmin, vmax), mask=mask)

        resdat = result.data
        mask = result.mask
        if mask is np.ma.nomask:
            mask = (resdat <= 0)
        else:
            mask |= (resdat <= 0)
        
        np.log(resdat, resdat)
        resdat -= np.log(vmin)
        resdat /= (np.log(vmax) - np.log(vmin))
        result = np.ma.array(resdat, mask=mask, copy=False)

    return result


def plot_beam(obs, target, cal, freq, time):

    metadata = get_obs_metadata(obs)

    phi = np.linspace(0,360,3600)
    theta = np.linspace(0,90,900)

    # make coordinate grid
    az, za = np.meshgrid(np.radians(phi), np.radians(theta))
        
    # compute beam and plot
    gx, gy = pb.MWA_Tile_full_EE(za, az, freq=freq*1e6, delays=metadata["delays"], power=True, zenithnorm=True)
    beam = (gx + gy) / 2.0

    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111, polar=True, aspect='auto')
    
    # filled contour setup
    lower_contour = 7e-3
    upper_contour = beam.max()
    fill_min = 1e-2 # 1% of zenith power
    fill_max = 0.95 * beam.max() # 95% max beam power ( != zenith power)

    Z = np.copy(beam)
    Z[Z <= fill_min] = 0
    Z[Z >= fill_max] = fill_max

    cc_levels = np.logspace(np.log10(lower_contour), np.log10(upper_contour), num=20)
    cc_cmap = plt.get_cmap('gray_r')
    
   
    cc_norm = LogNorm(vmin=cc_levels.min(), vmax=beam.max())
    cc = ax.contourf(az, za, Z, levels=cc_levels, cmap=cc_cmap, norm=cc_norm, zorder=1000)
    cc.cmap.set_over('black')  # anything over max is black
    cc.cmap.set_under('white') # anything under min is white

    # contour lines steup
    cs_levels = np.logspace(np.log10(fill_min), np.log10(fill_max), num=5)
    cs_cmap = plt.get_cmap('gist_heat')
    cs_norm = LogNorm(vmin=cs_levels.min(), vmax=cs_levels.max())
    cs = ax.contour(az, za, beam, levels=cs_levels, cmap=cs_cmap, norm=cs_norm, zorder=1001)

    # color bar setup
    cbarCC = plt.colorbar(cc, pad=0.08, shrink=0.9)
    cbarCC.set_label(label="zenith normalised power", size=20)
    cbarCC.set_ticks(cs_levels)
    cbarCC.set_ticklabels([r"${0:.2f}$".format(x) for x in cs_levels])
    cbarCC.ax.tick_params(labelsize=18)
    #add lines from the contours to the filled color map
    cbarCC.add_lines(cs)

    # plot the pointing centre of the tile beam
    ax.plot(np.radians(metadata["az"]), np.radians(metadata["za"]), ls="", marker="+", ms=8, color='cyan', zorder=1002, label="pointing centre")

    if target is not None:
        target_az = target.altaz.az.rad
        target_za = np.pi/2 - target.altaz.alt.rad

        # get beam power for target
        bpt_x, bpt_y = pb.MWA_Tile_full_EE([[target_za]], [[target_az]], freq=freq*1e6, delays=metadata["delays"], power=True, zenithnorm=True)
        bpt = (bpt_x + bpt_y) / 2.0
        lognormbpt = log_normalise(bpt, cc_levels.min(), beam.max())
        print "Beam power @ source:",bpt[0][0]
        print "log-normalised:",lognormbpt[0][0]
        
        # plot the target position on sky
        ax.plot(target_az, target_za, ls="", marker="o", color='C3', zorder=1002, label="target ({0:.2f})".format(bpt[0][0]))
        
        # plot the target on the color bar
        cbarCC.ax.plot(0.5, lognormbpt, color='C3', marker="o")

    if cal is not None:
        cal_az = cal.altaz.az.rad
        cal_za = np.pi/2 - cal.altaz.alt.rad
        
        # get the beam power for calibrator
        bpc_x, bpc_y = pb.MWA_Tile_full_EE([[cal_za]], [[cal_az]], freq=freq*1e6, delays=metadata["delays"], power=True, zenithnorm=True)
        bpc = (bpc_x + bpc_y) / 2.0
        lognormbpc = log_normalise(bpc, cc_levels.min(), beam.max())
        print "Beam power @ calibrtor:",bpc[0][0]
        print "log-normalised:",lognormbpc[0][0]
        
        # plot the calibrator position on sky
        ax.plot(cal_az, cal_za, ls="", marker="o", color='C2', zorder=1002, label="calibrator ({0:.2f})".format(bpc[0][0]))

        # plot the calibrator on the color bar
        cbarCC.ax.plot(0.5, lognormbpc, color='C2', marker="o")

    # draw grid
    ax.grid(color='k', ls="--", lw=0.5)

    # azimuth labels
    ax.set_xticks(np.radians([0, 45, 90, 135, 180, 225, 270, 315]))
    ax.set_xticklabels([r"${0:d}^\circ$".format(int(np.ceil(x))) for x in np.degrees(ax.get_xticks())], color='k')

    #Zenith angle labels
    ax.set_rlabel_position(250)
    ax.set_yticks(np.radians([20, 40, 60, 80]))
    ax.set_yticklabels([r"${0:d}^\circ$".format(int(np.ceil(x))) for x in np.degrees(ax.get_yticks())], color='k')

    # Title
    ax.set_title("MWA Tile beam (FEE)\naz = {0:.2f}, za = {1:.2f}, freq = {2:.2f}MHz\n{3}".format(metadata["az"], metadata["za"], freq, time.iso))

    plt.legend(bbox_to_anchor=(1.05,1.05))
    #plt.savefig("{0}_{1:.2f}MHz_tile.eps".format(obs, freq), bbox_inches="tight", format="eps")
    plt.savefig("{0}_{1:.2f}MHz_tile.png".format(obs, freq), bbox_inches="tight")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plotting tool for tile beam pattern")

    parser.add_argument("-o", "--obsid", type=int, help="Observation ID ", required=True)
    parser.add_argument("-f", "--freq", type=float, help="Frequency to evaluate beam at (MHz)", required=True)
    parser.add_argument("--ra", type=str, help="Target RA (J2000)")
    parser.add_argument("--dec", type=str, help="Target DEC (J2000)")
    parser.add_argument("--cal_ra", type=str, help="Calibrator RA (J2000)")
    parser.add_argument("--cal_dec", type=str, help="Calibrator DEC (J2000)")
    parser.add_argument("--gps", type=int, help="Time at which to evaluate target/calibrator position (GPS seconds)")
    parser.add_argument("--utc", type=str, help="Time at which to evaluate target/calibrator position (YYYY-MM-DDThh:mm:ss.ss UTC)")
    args = parser.parse_args()
    
    if args.gps and not args.utc:
        time = Time(args.gps, format="gps")
    elif args.utc and not args.gps:
        time = Time(args.utc, scale="utc", format="isot")
    elif args.gps and arg.utc:
        print "You supplied GPS and UTC times: only using GPS time"
        time = Time(args.gps, format="gps")
    else:
        time = None

    if args.ra and args.dec and time:
        target = compute_target_position(args.ra, args.dec, time) 
    else:
        print "No target RA, Dec or Time given. Not plotting target position."
        target = None
       
    if args.cal_ra and args.cal_dec and time:
        calibrator = compute_target_position(args.cal_ra, args.cal_dec, time)
    else:
        print "No calibrator RA, Dec or Time given. Not plotting calibrator position."
        calibrator = None
     

        
    plot_beam(args.obsid, target, calibrator, args.freq, time)
