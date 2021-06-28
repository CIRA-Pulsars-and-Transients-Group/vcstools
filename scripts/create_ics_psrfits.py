#!/usr/bin/env python3

import logging
import argparse
import sys
import os
import subprocess
import glob
import getpass

from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy import units as u

from vcstools.metadb_utils import get_common_obs_metadata, mwa_alt_az_za
from vcstools.config import load_config_file
from vcstools.general_utils import setup_logger

logger = logging.getLogger(__name__)


if __name__ == "__main__":

    # Dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING)

    parser = argparse.ArgumentParser()

    parser.add_argument("obsID", type=int, help="Observation ID for which to create the incoherent sum PSRFITS")
    parser.add_argument("-b", "--base_dir", type=str, help="Base directory to look and then create the PSRFITS. "
                                                           "In the directory provided, there MUST exist a 'combined' "
                                                           "subfolder. An additional subfolder named 'ics' will be "
                                                           "created as a result of the incoherent sum creation. "
                                                           "Default is /group/mwavcs/vcs/[obsID]", default=None)
    parser.add_argument("-l", "--loglvl", type=str, choices=loglevels.keys(), help="Desired logging verbosity level",
                        default="INFO")
    parser.add_argument("-V", "--version", action="store_true", help="Print version and quit")
    args = parser.parse_args()

    # set up the logger for stand-alone execution
    logger = setup_logger(logger, log_level=loglevels[args.loglvl])

    if args.version:
        try:
            import version
            logger.info(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            logger.error("Couldn't import version.py - have you installed vcstools?")
            logger.error("ImportError: {0}".format(ie))
            sys.exit(0)

    comp_config = load_config_file()
    if args.base_dir is None:
        data_dir = os.path.join(comp_config['base_data_dir'], str(args.obsID), "combined")
        ics_dir = os.path.join(comp_config['base_data_dir'], str(args.obsID), "ics")
    else:
        data_dir = "{base_dir}/{obsid}/combined".format(base_dir=args.base_dir, obsid=args.obsID)
        ics_dir = "{base_dir}/{obsid}/ics".format(base_dir=args.base_dir, obsid=args.obsID)

    logger.info("Data directory: {data_dir}".format(data_dir=data_dir))
    logger.info("ICS output directory: {ics_dir}".format(ics_dir=ics_dir))

    data_files = sorted(glob.glob("{data_dir}/*ics*.dat".format(data_dir=data_dir)))
    nfiles = len(data_files)

    ics_files = glob.glob("{ics_dir}/*.fits".format(ics_dir=ics_dir))
    if len(ics_files) > 0:
        logger.error("There are already PSRFITS files in the ics directory. Please delete or move them elsewhere")
        sys.exit(1)

    if os.path.isfile("{ics_dir}/mk_psrfits_in".format(ics_dir=ics_dir)):
        logger.error("Pipe file already exists in ICS directory. Please delete it")
        sys.exit(1)

    first_gps_second = int(data_files[0].split("/")[-1].split("_")[1])
    logger.info("First GPS second present: {0}".format(first_gps_second))
    mwa_loc = EarthLocation.from_geodetic(lon="116:40:14.93", lat="-26:42:11.95", height=377.8)
    first_gps_time = Time(first_gps_second, format="gps", scale='utc', location=mwa_loc)

    utc_isot = first_gps_time.isot
    utc_time = utc_isot.split("T")[-1]
    utc_dt = first_gps_time.to_datetime()
    utc_sec = float(utc_dt.hour) * 3600 + float(utc_dt.minute) * 60 + float(utc_dt.second)
    logger.debug("UTC = {0}".format(utc_isot))
    logger.debug("UTC seconds = {0}".format(utc_sec))

    lst = str(first_gps_time.sidereal_time('apparent').to_string(unit=u.hour, sep=":"))
    lst_sec = float(first_gps_time.sidereal_time('apparent').to_string(unit=u.hour, decimal=True)) * 3600.0
    logger.debug("LST = {0}".format(lst))
    logger.debug("LST seconds = {0}".format(lst_sec))

    mjd = first_gps_time.mjd
    mjd_trunc = int(mjd)
    logger.debug("MJD = {0}".format(mjd))

    logger.info("Getting observation metadata")
    metadata = get_common_obs_metadata(args.obsID)

    logger.info("Organising channels")
    channels = metadata[-1]
    nchans = len(channels)

    ch_offset = channels[-1] - channels[0]
    if ch_offset != nchans - 1:
        logger.error("Picket fence observation - cannot combine picket fence incoherent sum data (with this script)")
        sys.exit(1)

    bw = nchans * 1.28
    cfreq = (1.28 * min(channels) - 0.64) + bw / 2.0
    logger.info("Centre frequency: {0} MHz".format(cfreq))
    logger.info("Bandwidth: {0} MHz".format(bw))

    logger.info("Acquiring pointing position information")
    ra, dec = metadata[1:3]
    alt, az, za = mwa_alt_az_za(args.obsID, ra=ra, dec=dec)

    user = getpass.getuser()
    os.system("mkdir -p {ics_dir}".format(ics_dir=ics_dir))

    # TODO: this is bloody awful, surely there's a better way to do this?
    make_command = 'cd {ics_dir} && echo -e "{ics_dir}/mk_psrfits_in\n' \
                   '\n' \
                   '{obsid}\n' \
                   '\n' \
                   '{nfiles}\n' \
                   '{user}\n' \
                   '\n' \
                   '{obsid}\n' \
                   '\n' \
                   '\n' \
                   'MWA-G0024\n' \
                   '{utc_isot}\n' \
                   '\n' \
                   '{cfreq}\n' \
                   '{bw}\n' \
                   '{ra}\n' \
                   '{dec}\n' \
                   '{az}\n' \
                   '{za}\n' \
                   '\n' \
                   '{lst_sec}\n' \
                   '{utc_sec}\n' \
                   '{mjd_trunc}\n' \
                   '\n' \
                   '\n' \
                   '\n' \
                   '\n' \
                   '\n' \
                   '\n' \
                   '\n' \
                   '\n" | make_psrfits & sleep 1.0'.format(ics_dir=ics_dir, obsid=args.obsID, nfiles=nfiles,
                                                           user=user, utc_isot=utc_isot, cfreq=cfreq, bw=bw,
                                                           ra=ra, dec=dec, az=az, za=za, lst_sec=lst_sec,
                                                           utc_sec=utc_sec, mjd_trunc=mjd_trunc)

    cat_command = 'cat {data_dir}/*ics.dat > {ics_dir}/mk_psrfits_in'.format(data_dir=data_dir, ics_dir=ics_dir)

    logfile = "{ics_dir}/make_psrfits.log".format(ics_dir=ics_dir)
    logger.info("Now piping data to make_psrfits, writing stdout/stderr to:")
    logger.info("   {0}".format(logfile))
    with open(logfile, "w") as f:
        p1 = subprocess.Popen(make_command, shell=True, stdout=f, stderr=f)
        p2 = subprocess.Popen(cat_command, shell=True, stdout=f, stderr=f)
        p2.wait()
    os.system("rm -rf {ics_dir}/mk_psrfits_in".format(ics_dir=ics_dir))


