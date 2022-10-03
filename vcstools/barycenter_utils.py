#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import astropy.units as u
import numpy as np
from astropy.constants import c
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.coordinates import solar_system_ephemeris
from astropy.time import Time

# define MWA's Earth location
TEL_LAT = -26.703319
TEL_LON = 116.67081
TEL_ELEV = 377.827
TEL_LOCATION = EarthLocation(
    lat=TEL_LAT * u.deg, lon=TEL_LON * u.deg, height=TEL_ELEV * u.m
)
SEC_PER_DAY = 86400.0

# use the JPL planetary ephemeris
solar_system_ephemeris.set("jpl")

LOGGING_CONFIG = {}
logging_format = "[%(asctime)s] %(levelname)s::%(module)s: %(message)s"
logging.basicConfig(format=logging_format, level=logging.INFO)
logger = logging.getLogger()


def get_barycentric_correction(
    ra: str,
    dec: str,
    mjd: float,
    convention: str = "radio",
    for_presto: bool = False,
):
    """
    Compute the instantaneous barycentric velocity correction towards a given
    direction at a given time.

    Parameters
    ----------
    ra: str
        Right ascension formatted as hh:mm:ss.ss
    dec: str
        Declination formatted as dd:mm:ss.ss
    mjd: float
        Time in MJD
    convention: str, optional
        One of "optical", "radio", or "relativistic" which specifies the
        velocity convention used. Default is "radio"
    for_presto: bool, optional
        If the desired factor is for use with PRESTO tools, we have to
        negate the value.

    Returns
    -------
    float
        Barycentric velocity correction as a fraction of the speed of light.

    NOTES: This code produces values in the opposite sense to PRESTO
           (e.g., this barycentric value is the negative of what should
           be provided to PRESTO's utilities).
    """
    coord = SkyCoord(ra, dec, frame="icrs", unit=(u.hourangle, u.deg))
    t = Time(mjd, format="mjd", scale="utc")

    logger.info(f"RA (J2000) : {coord.ra.to_string(u.hour)}")
    logger.info(f"Dec (J2000) : {coord.dec.to_string(u.degree)}")
    logger.info(f"Topocentric MJD : {t.mjd}")

    vb_optical = coord.radial_velocity_correction(
        kind="barycentric", obstime=t, location=TEL_LOCATION
    )

    if convention == "optical":
        vb = vb_optical
    elif convention == "radio":
        vb = vb_optical.to(u.Hz, u.doppler_optical(1 * u.Hz)).to(
            vb_optical.unit, u.doppler_radio(1 * u.Hz)
        )
    elif convention == "relativistic":
        vb = vb_optical.to(u.Hz, u.doppler_optical(1 * u.Hz)).to(
            vb_optical.unit, u.doppler_relativistic(1 * u.Hz)
        )
    else:
        logger.warning(
            "convention not recognised, must be one of: optical, radio, relativistic"
        )
        logger.warning("Defaulting to RADIO convention")
        vb = vb_optical.to(u.Hz, u.doppler_optical(1 * u.Hz)).to(
            vb_optical.unit, u.doppler_radio(1 * u.Hz)
        )

    logger.info(
        "barycentric velocity (fraction of speed of light) = "
        "{0}".format(vb.value / c.value)
    )

    if for_presto:
        logger.debug(
            f"Negating computed value for PRESTO use (sign-convention differences)"
        )
        return -vb.value / c.value
    else:
        return vb.value / c.value


def get_mean_barycentric_correction(
    ra: str,
    dec: str,
    start_mjd: float,
    duration: float,
    nsteps: int = 100,
    convention: str = "radio",
    for_presto: bool = False,
):
    """
    Compute the mean barycentric velocity correction towards a given direction
    over a certain time span.

    Parameters
    ----------
    ra: str
        Right ascension formatted as hh:mm:ss.ss
    dec: str
        Declination formatted as dd:mm:ss.ss
    start_mjd: float
        Start time in MJD
    duration: float
        Duration over which to calculate the mean correction, in seconds
    nsteps: int, optional
        Number of linearly-spaced time steps to evaluate the correction (default: 100)
    convention: str, optional
        One of "optical", "radio", or "relativistic" which specifies the
        velocity convention used. Default is "radio"
    for_presto: bool, optional
        If the desired factor is for use with PRESTO tools, we have to
        negate the value.

    Returns
    -------
    float
        The mean Barycentric velocity correction as a fraction of the speed of
        light

    NOTES: This code produces values in the opposite sense to PRESTO
           (e.g., this barycentric value is the negative of what should
           be provided to PRESTO's utilities).
    """
    mjds = np.linspace(start_mjd, start_mjd + duration / SEC_PER_DAY, nsteps)
    coord = SkyCoord(ra, dec, frame="icrs", unit=(u.hourangle, u.deg))
    times = Time(mjds, format="mjd", scale="utc")

    logger.info(f"RA (J2000) : {coord.ra.to_string(u.hour)}")
    logger.info(f"Dec (J2000) : {coord.dec.to_string(u.degree)}")
    logger.info(f"Topocentric MJD range : {np.min(mjds)} - {np.max(mjds)}")

    vb_corr = []
    for t in times:
        vb_optical = coord.radial_velocity_correction(
            kind="barycentric", obstime=t, location=TEL_LOCATION
        )

        if convention == "optical":
            vb = vb_optical
        elif convention == "radio":
            vb = vb_optical.to(u.Hz, u.doppler_optical(1 * u.Hz)).to(
                vb_optical.unit, u.doppler_radio(1 * u.Hz)
            )
        elif convention == "relativistic":
            vb = vb_optical.to(u.Hz, u.doppler_optical(1 * u.Hz)).to(
                vb_optical.unit, u.doppler_relativistic(1 * u.Hz)
            )
        else:
            logger.warning(
                "convention not recognised, must be one of: optical, radio, "
                "relativistic"
            )
            logger.warning("Defaulting to RADIO convention")
            vb = vb_optical.to(u.Hz, u.doppler_optical(1 * u.Hz)).to(
                vb_optical.unit, u.doppler_radio(1 * u.Hz)
            )

        vb_corr.append(vb.value)

    mean_vb = np.mean(vb_corr)

    logger.info(
        "mean barycentric velocity (fraction of speed of light) = "
        "{0}".format(mean_vb / c.value)
    )

    if for_presto:
        logger.debug(
            f"Negating computed value for PRESTO use (sign-convention differences)"
        )
        return -mean_vb / c.value
    else:
        return mean_vb / c.value
