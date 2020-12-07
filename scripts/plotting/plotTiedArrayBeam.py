#!/usr/bin/env python3

import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord, AltAz
from astropy import units as u
from astropy.time import Time

import matplotlib.pyplot as plt

from mwa_pb import primary_beam as pb
from vcstools.metadb_utils import getmeta, mwa_alt_az_za

import argparse

import logging
logger = logging.getLogger(__name__)

def get_obs_metadata(obs):
    beam_meta_data = getmeta(service='obs', params={'obs_id':obs})
    channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
    freqs = [float(c)*1.28 for c in channels]
    xdelays = beam_meta_data[u'rfstreams'][u"0"][u'xdelays']
    #pythodelays = beam_meta_data[u'rfstreams'][u"0"][u'xdelays']
    ydelays = beam_meta_data[u'rfstreams'][u"0"][u'ydelays']
    _, pointing_AZ, pointing_ZA = mwa_alt_az_za(obs)

    return {"channels":channels,
            "frequencies":freqs,
            "xdelays":xdelays,
            "ydelays":ydelays,
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
    """This is almost directly copied from matplotlib's color.py """
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


def load_data(fname):
    # TODO: the column numbers have changed, so will need to make sure simulation code + this is consistent
    theta, phi, beam = np.genfromtxt(fname, usecols=(0,1,8), skip_header=14, unpack=True)
    return theta, phi, beam


def plot_beam(obs, fname, target, freq, time):

    logger.info("Getting observation metadata")
    metadata = get_obs_metadata(obs)

    logger.info("Loading data from file (this may take a while...)")
    theta, phi, beam = load_data(fname)

    logger.info("Re-normalising beam-pattern")
    az, za = np.meshgrid(np.radians(sorted(set(phi))), np.radians(sorted(set(theta))))
    delays = [metadata['xdelays'], metadata['ydelays']]
    gx, gy = pb.MWA_Tile_analytic(za, az, freq=freq*1e6, delays=delays, power=True, zenithnorm=True)
    tile_beam = (gx + gy) / 2.0

    # re-normalise tied-array beam to 1 at maximum and then apply tile beam
    # this acts to reduce the effects of simulation resolution (for plotting purposes)
    beam /= beam.max()
    beam *= np.ravel(tile_beam)



    # start setup for plotting
    lower_contour = 1e-2
    upper_contour = beam.max()

    fill_min = 7e-3
    fill_max = 0.95 * beam.max()
    beam[beam <= fill_min] = 0
    beam[beam >= fill_max] = fill_max


    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111, polar=True, aspect='equal')

    # plot the beam pattern
    logger.info("Plotting beam pattern (this may take a while, too...)")
    cf_levels = np.linspace(lower_contour, upper_contour, num=20)
    #cf_levels_log = np.logspace(np.log10(lower_contour), np.log10(upper_contour), num=20)
    logger.info("     contour levels: max,min = {0}, {1}".format(cf_levels.max(), cf_levels.min()))
    cf_cmap = plt.get_cmap('gray_r')
    logger.info("     beam levels: max,min = {0},{1}".format( beam.max(), beam.min()))

    #sys.exit()
    logger.info("     plotting tied-array beam pattern contours")
    #cf_norm = LogNorm(vmin=cf_levels.min(), vmax=beam.max())
    #cf = ax.tricontourf(np.radians(phi), np.radians(theta), beam, cmap=cf_cmap, norm=cf_norm, levels=cf_levels)
    cf = ax.tricontourf(np.radians(phi), np.radians(theta), beam, cmap=cf_cmap, levels=cf_levels)
    cf.cmap.set_over('k')
    cf.cmap.set_under('white')

    # color bar setup
    logger.info("     assigning colorbar")
    cbar_levels = np.linspace(fill_min, fill_max, num=6)
    #cbar_levels_log = np.logspace(np.log10(fill_min), np.log10(fill_max), num=6)
    cbar = plt.colorbar(cf, shrink=0.9, pad=0.08)
    cbar.set_label(label="zenith normalised power", size=20, labelpad=20)
    cbar.set_ticks(cbar_levels)
    cbar.set_ticklabels([r"${0:.3f}$".format(x) for x in cbar_levels])
    cbar.ax.tick_params(labelsize=18)

    # plot the pointing centre
    logger.info("     plotting observation pointing centre")
    ax.plot(np.radians(metadata["az"]), np.radians(metadata["za"]), ls="", marker="+", ms=8, color='cyan', zorder=1002, label="pointing centre")

    # plot the target position on sky
    if target is not None:
        logger.info("     plotting target position")
        target_az = target.altaz.az.rad
        target_za = np.pi/2 - target.altaz.alt.rad

        # get beam power for target
        bpt_x, bpt_y = pb.MWA_Tile_analytic(target_za, target_az, freq=freq*1e6, delays=delays, power=True, zenithnorm=True)
        bpt = (bpt_x + bpt_y) / 2.0
        lognormbpt = log_normalise(bpt, cf_levels.min(), beam.max())
        logger.info("Beam power @ source: {0}".format(bpt))
        logger.info("   log-normalised: {0}".format(lognormbpt))

        # plot the target position on sky
        ax.plot(target_az, target_za, ls="", marker="o", color='C3', zorder=1002, label="   target:{} ".format(bpt))

        # plot the target on the color bar
        cbar.ax.plot(0.5, lognormbpt, color='C3', marker="o")

    # draw grid
    ax.grid(color='k', ls="--", lw=0.5)

    # azimuth labels
    ax.set_xticks(np.radians([0, 45, 90, 135, 180, 225, 270, 315]))
    ax.set_xticklabels([r"${0:d}^\circ$".format(int(np.ceil(x))) for x in np.degrees(ax.get_xticks())], color='k')

    # zenith angle labels
    ax.set_rlabel_position(250)
    ax.set_ylim(0, np.pi/2)
    ax.set_yticks(np.radians([20, 40, 60, 80]))
    ax.set_yticklabels([r"${0:d}^\circ$".format(int(np.ceil(x))) for x in np.degrees(ax.get_yticks())], color='k')

    # title
    ax.set_title("MWA tied-array beam\naz = {0:.2f}, za = {1:.2f}, freq = {2:.2f}MHz\n{3}".format(metadata["az"], metadata["za"], freq, time.iso))

    logger.info("Saving figure")
    #plt.savefig("{0}_{1:.2f}MHz_tabeam.eps".format(obs, freq), bbox_inches="tight", format="eps")
    plt.savefig("{0}_{1:.2f}MHz_tabeam.png".format(obs, freq), bbox_inches="tight")




if __name__ == "__main__":

    #dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                    INFO=logging.INFO,
                    WARNING=logging.WARNING,
                    ERROR = logging.ERROR)

    parser = argparse.ArgumentParser(description="Plotting tool for the tied-array beam simulations")

    parser.add_argument("-o", "--obsid", type=int, help="Observation ID ", required=True)
    parser.add_argument("-f", "--freq", type=float, help="Frequency to evaluate beam at (MHz)", required=True)
    parser.add_argument("--file", type=str, help="File name storing the tied-array beam data", required=True)
    parser.add_argument("--ra", type=str, help="Target RA (J2000)")
    parser.add_argument("--dec", type=str, help="Target DEC (J2000)")
    parser.add_argument("--gps", type=int, help="Time at which to evaluate target position (GPS seconds)")
    parser.add_argument("--utc", type=str, help="Time at which to evaluate target position (YYYY-MM-DDThh:mm:ss.ss UTC)")
    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO", choices=loglevels.keys(), default="INFO")
    args = parser.parse_args()

    #set log levels
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    logger.propagate = False


    if args.gps and not args.utc:
        time = Time(args.gps, format="gps")
    elif args.utc and not args.gps:
        time = Time(args.utc, scale="utc", format="isot")
    elif args.gps and args.utc:
        logger.warn("You supplied GPS and UTC times: only using GPS time")
        time = Time(args.gps, format="gps")
    else:
        time = None

    if args.ra and args.dec and time:
        target = compute_target_position(args.ra, args.dec, time)
    else:
        logger.warn("No target RA, Dec or Time given. Not plotting target position.")
        target = None


    plot_beam(args.obsid, args.file, target, args.freq, time)
