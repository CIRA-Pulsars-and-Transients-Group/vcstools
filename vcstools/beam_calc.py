"""
Functions used to calculate or simulate the MWA tile beam.
"""

import os
import sys
import numpy as np
from astropy.table import Table
import psrqpy

#MWA scripts
from mwa_pb import primary_beam
from mwa_pb import config

from vcstools import data_load
from vcstools.pointing_utils import sex2deg
from vcstools.metadb_utils import mwa_alt_az_za, get_common_obs_metadata,\
                                  getmeta, obs_max_min

try:
    import mwa_hyperbeam
except ImportError:
    logger.warning('Could not import mwa_hyperbeam; using pure Python implementation')


import logging
logger = logging.getLogger(__name__)

def pixel_area(ra_min, ra_max, dec_min, dec_max):
    """Calculate the area of a pixel on the sky from the pixel borders

    Parameters
    ----------
    ra_min : `float`
        The Right Acension minimum in degrees.
    ra_max : `float`
        The Right Acension maximum in degrees.
    dec_min : `float`
        The Declination minimum in degrees.
    dec_max : `float`
        The Declination maximum in degrees.

    Returns
    -------
    area : `float`
        Area of the pixel in square degrees.
    """
    return (ra_max - ra_min) * np.degrees(np.sin(np.radians(dec_max)) - np.sin(np.radians(dec_min)))


def field_of_view(obsid,
                  common_metadata=None, dur=None):
    """Will find the field-of-view of the observation (including the drift) in square degrees.

    Parameters
    ----------
    obsid : `int`
        The MWA observation ID.
    common_metadata : `list`, optional
        The list of common metadata generated from :py:meth:`vcstools.metadb_utils.get_common_obs_metadata`
    dur : `int`, optional
        Duration of observation to calculate for in seconds.
        By default will use the entire observation duration.

    Returns
    -------
    area : `float`
        The field-of-view of the observation in square degrees.
    """
    if common_metadata is None:
        common_metadata = get_common_obs_metadata(obsid)

    if dur is None:
        dt = 296
    else:
        dt = 100
        # Change the dur to the inpur dur
        obsid, ra, dec, _, delays, centrefreq, channels = common_metadata
        common_metadata = [obsid, ra, dec, dur, delays, centrefreq, channels]

    # Make a pixel for each degree on the sky
    names_ra_dec = []
    for ra in range(0,360):
        for dec in range(-90,90):
            names_ra_dec.append(["sky_pos", ra+0.5, dec+0.5])

    # Get tile beam power for all pixels
    sky_powers = get_beam_power_over_time(names_ra_dec,
                                          common_metadata=common_metadata,
                                          degrees=True, dt=dt)

    # Find the maximum power over all time
    max_sky_powers = []
    for pixel_power in sky_powers:
        temp_power = 0.
        for time_power in pixel_power:
            if time_power[0] > temp_power:
                temp_power = time_power[0]
        max_sky_powers.append(temp_power)

    # Find all pixels greater than the half power point and sum their area
    half_power_point = max(max_sky_powers) / 2
    i = 0
    area_sum = 0
    for ra in range(0,360):
        for dec in range(-90,90):
            if max_sky_powers[i] > half_power_point:
                area_sum = area_sum + pixel_area(ra, ra+1, dec, dec+1)
            i = i + 1
    return area_sum


def from_power_to_gain(power, cfreq, n, coh=True):
    """Estimate the gain from the tile beam power.

    Parameters
    ----------
    power : `float`
        The tile beam power at the source position.
    cfreq : `float`
        The centre frequency of the observation in Hz
    n : `int`
        The number of non-flagged MWA tiles.
    coh : `boolean`, optional
        `True` if the observation is coherent (tied-array beam) or `False` if it's incoherent. |br| Default: `True`.

    Returns
    -------
    gain : `float`
        Gain in K/Jy.
    """

    from astropy.constants import c,k_B
    from numpy import sqrt

    obswl = c.value/cfreq
    #for coherent
    if coh:
        coeff = obswl**2*16*n/(4*np.pi*k_B.value)
    else:
        coeff = obswl**2*16*sqrt(n)/(4*np.pi*k_B.value)
    logger.debug("Wavelength {} m".format(obswl))
    logger.debug("Gain coefficient: {}".format(coeff))
    SI_to_Jy = 1e-26
    return (power*coeff)*SI_to_Jy


def get_Trec(obsfreq, trcvr_file=None):
    """Get receiver temperature from the temperature receiver file.

    Parameters
    ----------
    obsfreq : `float`
        The observing frequency in MHz.
    trcvr_file : `str`, optional
        The Trec file location to read in. If none is supplied, will use the Trec file in the data directory.

    Returns
    -------
    Trec : `float`
        The receiver temperature in K.
    """
    if trcvr_file is None:
        # Use trcvr file from vcstools data location
        trcvr_file = data_load.TRCVR_FILE
    tab = Table.read(data_load.TRCVR_FILE, format="csv")
    Trec = 0.0
    for r in range(len(tab)-1):
        if tab[r][0]==obsfreq:
            Trec = tab[r][1]
        elif tab[r][0] < obsfreq < tab[r+1][0]:
            Trec = ((tab[r][1] + tab[r+1][1])/2)
    if Trec == 0.0:
        logger.debug("ERROR getting Trec")
    return Trec


def beam_enter_exit(powers, duration, dt=296, min_z_power=0.3):
    """Calculates when the source enters and exits the beam

    Parameters
    ----------
    powers : `list`, (ntimes, nfreqs)
        Powers for the duration every dt and freq.
    duration : `int`
        Duration of the observation according to the metadata in seconds.
    dt : `int`, optional
        The time interval of how often powers are calculated. Default: 296.
    min_z_power : `float`, optional
        Zenith normalised power cut off. |br| Default: 0.3.

    Returns
    -------
    dect_beg_norm, dect_end_norm : `float`
        Fraction of the observation when the source enters and exits the beam respectively.
    """
    from scipy.interpolate import UnivariateSpline
    time_steps = np.array(range(0, duration, dt), dtype=float)

    # For each time step record the min power so even if the source is in
    # one freq channel it's recorded
    powers_freq_min = []
    for p in powers:
        powers_freq_min.append(float(min(p) - min_z_power))

    if min(powers_freq_min) > 0.:
        enter_beam = 0.
        exit_beam = 1.
    else:
        powers_freq_min = np.array(powers_freq_min)
        logger.debug("time_steps: {}".format(time_steps))
        logger.debug("powers: {}".format(powers_freq_min))
        try:
            spline = UnivariateSpline(time_steps, powers_freq_min , s=0.)
        except:
            return None, None
        if len(spline.roots()) == 2:
            enter_beam, exit_beam = spline.roots()
            enter_beam /= duration
            exit_beam /= duration
        elif len(spline.roots()) == 1:
            if powers_freq_min[0] > powers_freq_min[-1]:
                #power declines so starts in beam then exits
                enter_beam = 0.
                exit_beam = spline.roots()[0]/duration
            else:
                enter_beam = spline.roots()[0]/duration
                exit_beam = 1.
        else:
            enter_beam = 0.
            exit_beam = 1.
    return enter_beam, exit_beam


def get_beam_power_over_time(names_ra_dec,
                             common_metadata=None,
                             dt=296, centeronly=True, verbose=False,
                             option='analytic', degrees=False,
                             start_time=0):
    """Calculates the zenith normalised power for each source over time.

    Parameters
    ----------
    names_ra_dec : `list`
        An array in the format [[source_name, RAJ, DecJ]]
    common_metadata : `list`, optional
        The list of common metadata generated from :py:meth:`vcstools.metadb_utils.get_common_obs_metadata`
    dt : `int`, optional
        The time interval of how often powers are calculated. |br| Default: 296.
    centeronly : `boolean`, optional
        Only calculates for the centre frequency. |br| Default: `True`.
    verbose : `boolean`, optional
        If `True` will not supress the output from mwa_pb. |br| Default: `False`.
    option : `str`, optional
        The primary beam model to use out of [analytic, advanced, full_EE, hyperbeam]. |br| Default: analytic.
    degrees : `boolean`, optional
        If true assumes RAJ and DecJ are in degrees. |br| Default: `False`.
    start_time : `int`, optional
        The time in seconds from the begining of the observation to start calculating at. |br| Default: 0.

    Returns
    -------
    Powers : `numpy.array`, (len(names_ra_dec), ntimes, nfreqs)
        The zenith normalised power for each source over time.
    """
    if common_metadata is None:
        common_metadata = get_common_obs_metadata(obsid)
    obsid, _, _, time, delays, centrefreq, channels = common_metadata
    names_ra_dec = np.array(names_ra_dec)
    amps = [1.0] * 16
    logger.debug("Calculating beam power for OBS ID: {0}".format(obsid))

    if option == 'hyperbeam':
        if "mwa_hyperbeam" not in sys.modules:
            logger.error("mwa_hyperbeam not installed so can not use hyperbeam to create a beam model. Exiting")
            sys.exit(1)
        beam = mwa_hyperbeam.FEEBeam(config.h5file)

    # Work out time steps to calculate over
    starttimes = np.arange(start_time, time + start_time, dt)
    stoptimes  = starttimes + dt
    stoptimes[stoptimes>time] = time
    ntimes = len(starttimes)
    midtimes = float(obsid) + 0.5 * (starttimes + stoptimes)

    # Work out frequency steps
    if centeronly:
        if centrefreq > 1e6:
            logger.warning("centrefreq is greater than 1e6, assuming input with units of Hz.")
            frequencies=np.array([centrefreq])
        else:
            frequencies=np.array([centrefreq])*1e6
        nfreqs = 1
    else:
        # in Hz
        frequencies = np.array(channels) * 1.28e6
        nfreqs = len(channels)

    # Set up np power array
    PowersX = np.zeros((len(names_ra_dec),
                        ntimes,
                        nfreqs))
    PowersY = np.zeros((len(names_ra_dec),
                        ntimes,
                        nfreqs))

    # Convert RA and Dec to desired units
    if degrees:
        RAs  = np.array(names_ra_dec[:,1], dtype=float)
        Decs = np.array(names_ra_dec[:,2], dtype=float)
    else:
        RAs, Decs = sex2deg(names_ra_dec[:,1], names_ra_dec[:,2])
    # Then check if they're valid
    if len(RAs) == 0:
        sys.stderr.write('Must supply >=1 source positions\n')
        return None
    if not len(RAs) == len(Decs):
        sys.stderr.write('Must supply equal numbers of RAs and Decs\n')
        return None

    if verbose is False:
        #Supress print statements of the primary beam model functions
        sys.stdout = open(os.devnull, 'w')
    for itime in range(ntimes):
        # this differ's from the previous ephem_utils method by 0.1 degrees
        _, Azs, Zas = mwa_alt_az_za(midtimes[itime], ra=RAs, dec=Decs, degrees=True)
        # go from altitude to zenith angle
        theta = np.radians(Zas)
        phi   = np.radians(Azs)
        for ifreq in range(nfreqs):
            #Decide on beam model
            if option == 'analytic':
                rX, rY = primary_beam.MWA_Tile_analytic(theta, phi,
                                                        freq=frequencies[ifreq], delays=delays,
                                                        zenithnorm=True, power=True)
            elif option == 'advanced':
                rX, rY = primary_beam.MWA_Tile_advanced(theta, phi,
                                                        freq=frequencies[ifreq], delays=delays,
                                                        zenithnorm=True, power=True)
            elif option == 'full_EE':
                rX, rY = primary_beam.MWA_Tile_full_EE(theta, phi,
                                                       freq=frequencies[ifreq], delays=delays,
                                                       zenithnorm=True, power=True)
            elif option == 'hyperbeam':
                jones = beam.calc_jones_array(phi, theta, int(frequencies[ifreq]), delays[0], amps, True)
                jones = jones.reshape(1, len(phi), 2, 2)
                vis = primary_beam.mwa_tile.makeUnpolInstrumentalResponse(jones, jones)
                rX, rY = (vis[:, :, 0, 0].real, vis[:, :, 1, 1].real)
        PowersX[:,itime,ifreq] = rX
        PowersY[:,itime,ifreq] = rY
    if verbose is False:
        sys.stdout = sys.__stdout__
    Powers = 0.5 * (PowersX + PowersY)
    return Powers


def source_beam_coverage(obs_list, names_ra_dec,
                         common_metadata_list=None,
                         dt_input=300, beam='analytic', min_z_power=0.3):
    """For a list of MWA observations and sources will find if the
    sources are in the beam and when they enter and exit.

    Parameters
    ----------
    obs_list : `list`
        A list of MWA Observation IDs.
    names_ra_dec : `list`
        An array in the format [[source_name, RAJ, DecJ]]
    common_metadata_list : `list` of `list`, optional
        A list of lists where each list is the common metadata generated from
        :py:meth:`vcstools.metadb_utils.get_common_obs_metadata` for each obs in the obs_list.
    dt_input : `int`, optional
        The time interval in seconds of how often powers are calculated. |br| Default: 300.
    beam : `str`, optional
        The primary beam model to use out of [analytic, advanced, full_EE, hyperbeam]. |br| Default: analytic.
    min_z_power : `float`, optional
        Zenith normalised power cut off. |br| Default: 0.3.

    Returns
    -------
    beam_coverage : `dict`
        A dictionary where the first key is the observation ID and the second is the pulsar names like so: |br|
        beam_coverage[obsid][name] = [dect_beg_norm, dect_end_norm, np.amax(source_ob_power)]
        dect_beg_norm : `float`
            Fraction of the observation when the source enters the beam.
        dect_end_norm : `float`
            Fraction of the observation when the source enters and exits the beam respectively.
    """
    beam_coverage = {}
    for on, obsid in enumerate(obs_list):
        beam_coverage[obsid] = {}
        if common_metadata_list:
            common_metadata = common_metadata_list[on]
        else:
            # No metadata supplied so make the metadata call
            common_metadata = get_common_obs_metadata(obsid)
        #common_metadata = obsid,ra_obs,dec_obs,time_obs,delays,centrefreq,channels
        obs_dur = common_metadata[3]

        if dt_input * 4 >  obs_dur:
            # If the observation time is very short then a smaller dt time is required
            # to get enough points to fit a curve
            dt = int(obs_dur / 4.)
        else:
            dt = dt_input
        logger.debug("obsid: {0}, time_obs {1} s, dt {2} s".format(obsid, obs_dur, dt))
        logger.debug("names_ra_dec: {}".format(names_ra_dec))
        powers = get_beam_power_over_time(names_ra_dec,
                                          common_metadata=common_metadata,
                                          dt=dt, centeronly=True,
                                          option=beam)
        for source_ob_power, name in zip(powers, np.array(names_ra_dec)[:,0]):
            if max(source_ob_power) > min_z_power:
                logger.debug("Running beam_enter_exit on obsid: {}".format(obsid))
                dect_beg_norm, dect_end_norm = beam_enter_exit(source_ob_power, obs_dur,
                                                                  dt=dt, min_z_power=min_z_power)
                beam_coverage[obsid][name] = [dect_beg_norm, dect_end_norm, np.amax(source_ob_power)]
    return beam_coverage



def source_beam_coverage_and_times(obsid, pulsar,
                                   p_ra=None, p_dec=None,
                                   obs_beg=None, obs_end=None,
                                   files_beg=None, files_end=None,
                                   min_z_power=0.3, dt_input=100,
                                   common_metadata=None, query=None,
                                   beam='analytic'):
    """Finds the normalised time that a pulsar is in the beam for a given obsid.
    If pulsar is not in beam, returns None, None

    Parameters
    ----------
    obsid : `int`
        The observation ID
    pulsar : `str`
        The pulsar's J name
    p_ra, p_dec : `str`, optional
        The target's right ascension and declination in sexidecimals.
        If not supplied will use the values from the ANTF.
    obs_beg, obs_end : `int`, optional
        Beginning and end GPS time of the observation.
        If not supplied will use :py:meth:`vcstools.metadb_utils.obs_max_min` to find it.
    files_beg, files_end : `int`, optional
        Beginning and end GPS time of the (fits of VCS) files.
        If not supplied will assume the full observation is available.
    min_z_power : `float`, optional
        Zenith normalised power cut off. |br| Default: 0.3.
    common_metadata : `list`, optional
        The list of common metadata generated from :py:meth:`vcstools.metadb_utils.get_common_obs_metadata`
    query : psrqpy object, optional
        A previous psrqpy query. Can be supplied to prevent performing a new query.
    beam : `str`, optional
        The primary beam model to use out of [analytic, advanced, full_EE]. |br| Default: analytic.

    Returns
    -------
    enter_files : `float`
        A float between 0 and 1 that describes the normalised time that the pulsar enters the beam
    exit_files : `float`
        A float between 0 and 1 that describes the normalised time that the pulsar exits the beam
    """
    # Perform required metadata calls
    if query is None:
        query = psrqpy.QueryATNF(psrs=pulsar, loadfromdb=data_load.ATNF_LOC).pandas
    if p_ra is None or p_dec is None:
        # Get some basic pulsar and obs info info
        query_id = list(query['PSRJ']).index(pulsar)
        p_ra = query["RAJ"][query_id]
        p_dec = query["DECJ"][query_id]
    if not common_metadata:
        common_metadata = get_common_obs_metadata(obsid)
    if obs_beg is None or obs_end is None:
        obs_beg, obs_end = obs_max_min(obsid)
    obs_dur = obs_end - obs_beg + 1
    if not files_beg:
        files_beg = obs_beg
    if not files_end:
        files_end = obs_end
    files_dur = files_end - files_beg + 1

    beam_coverage = source_beam_coverage([obsid], [[pulsar, p_ra, p_dec]],
                            common_metadata_list=[common_metadata],
                            dt_input=dt_input, beam=beam, min_z_power=min_z_power)
    if pulsar not in beam_coverage[obsid].keys():
        # Not in beam exiting
        return None, None, None, None, None, None, None, None, None
    dect_beg_norm, dect_end_norm, _ = beam_coverage[obsid][pulsar]

    # GPS times the source enters and exits beam
    dect_beg = obs_beg + obs_dur * dect_beg_norm
    dect_end = obs_beg + obs_dur * dect_end_norm

    # Normalised time the source enters/exits the beam in the files (used for Presto commands)
    files_beg_norm = (dect_beg - files_beg) / files_dur
    files_end_norm = (dect_end - files_beg) / files_dur

    if files_beg_norm > 1. or files_end_norm < 0.:
        logger.debug("source {0} is not in the beam for the files on disk".format(pulsar))
        files_beg_norm = None
        files_end_norm = None
    else:
        if files_beg_norm < 0.:
            files_beg_norm = 0.
        if files_end_norm > 1.:
            files_end_norm = 1.

    return dect_beg, dect_end, dect_beg_norm, dect_end_norm, files_beg_norm, files_end_norm, obs_beg, obs_end, obs_dur


def find_sources_in_obs(obsid_list, names_ra_dec,
                        obs_for_source=False, dt_input=300, beam='analytic',
                        min_z_power=0.3, cal_check=False, all_volt=False,
                        degrees_check=False, metadata_list=None):
    """Either creates text files for each MWA obs ID of each source within it or a text
    file for each source with each MWA obs is that the source is in.

    Parameters
    ----------
    obsid_list : `list`
        List of MWA observation IDs.
    names_ra_dec : `list`
        An array in the format [[source_name, RAJ, DecJ]]
    obs_for_source : `boolean`, optional
        If `True` creates a text file for each source with each MWA observation that the source is in.
        If `False` creates text files for each MWA obs ID of each source within it. |br| Default: `False`.
    dt_input : `int`, optional
        The time interval in seconds of how often powers are calculated. |br| Default: 300.
    beam : `str`, optional
        The primary beam model to use out of [analytic, advanced, full_EE]. |br| Default: analytic.
    min_z_power : `float`, optional
        Zenith normalised power cut off. |br| Default: 0.3.
    cal_check : `boolean`, optional
        Checks the MWA pulsar database if there is a calibration suitable for the observation ID.
    all_volt : `boolean`, optional
        Included observations with missing or incorrect voltage files. |br| Default: `False`.
    degrees_check : `boolean`, optional
        If true assumes RAJ and DecJ are in degrees. |br| Default: `False`.
    metadata_list : `list`
        List of the outputs of vcstools.metadb_utils.get_common_obs_metadata.
        If not provided, will make the metadata calls to find the data. |br| Default: `None`.

    Returns
    -------
    output_data : `dict`
        The format of output_data is dependant on obs_for_source.
        |br| If obs_for_source is `True` :
        |br|    output_data = {jname:[[obsid, duration, enter, exit, max_power],
        |br|                          [obsid, duration, enter, exit, max_power]]}
        |br| If obs_for_source is `False` :
        |br|    ouput_data = {obsid:[[jname, enter, exit, max_power],
        |br|                         [jname, enter, exit, max_power]]}
    obsid_meta : `list`
        A list of the output of get_common_obs_metadata for each obsid
    """
    import urllib.request
    #prepares metadata calls and calculates power
    powers = []
    #powers[obsid][source][time][freq]
    common_metadata_list = []
    obsid_to_remove = []

    # Loop over observations to check if there are VCS files
    for i, obsid in enumerate(obsid_list):
        # Perform the file meta data call over 10 only 10 seconds as that is suffient test
        try:
            files_meta_data = getmeta(service='data_files', params={'obs_id':obsid, 'nocache':1,
                                                                    'mintime':int(obsid) + 10,
                                                                    'maxtime':int(obsid) + 20})
        except urllib.error.HTTPError as err:
            files_meta_data = None
        if files_meta_data is None:
            logger.warning("No file metadata data found for obsid {}. Skipping".format(obsid))
            obsid_to_remove.append(obsid)
            continue

        # Check raw voltage files
        raw_available = False
        raw_deleted   = False
        for file_name in files_meta_data.keys():
            if file_name.endswith('dat'):
                deleted = files_meta_data[file_name]['deleted']
                if deleted:
                    raw_deleted   = True
                else:
                    raw_available = True

        # Check combined voltage tar files
        comb_available = False
        comb_deleted   = False
        for file_name in files_meta_data.keys():
            if file_name.endswith('tar'):
                deleted = files_meta_data[file_name]['deleted']
                if deleted:
                    comb_deleted   = True
                else:
                    comb_available = True

        if raw_available or comb_available or all_volt:
            if metadata_list:
                common_metadata, _ = metadata_list[i]
            else:
                # No metadata supplied so make the metadata call
                common_metadata = get_common_obs_metadata(obsid)
            common_metadata_list.append(common_metadata)
        elif raw_deleted and comb_deleted:
            logger.warning('Raw and combined voltage files deleted for {}'.format(obsid))
            obsid_to_remove.append(obsid)
        elif raw_deleted:
            logger.warning('Raw voltage files deleted for {}'.format(obsid))
            obsid_to_remove.append(obsid)
        elif comb_deleted:
            logger.warning('Combined voltage files deleted for {}'.format(obsid))
            obsid_to_remove.append(obsid)
        else:
            logger.warning('No raw or combined voltage files for {}'.format(obsid))
            obsid_to_remove.append(obsid)
    for otr in obsid_to_remove:
        obsid_list.remove(otr)

    # Calculate the power for all sources and obsids and find when they enter and exit the beam
    beam_coverage = source_beam_coverage(obsid_list, names_ra_dec,
                         common_metadata_list=common_metadata_list,
                         dt_input=dt_input, beam=beam, min_z_power=min_z_power)

    #chooses whether to list the source in each obs or the obs for each source
    output_data = {}
    if obs_for_source:
        for source_name in np.array(names_ra_dec)[:,0]:
            source_data = []
            for on, obsid in enumerate(obsid_list):
                if source_name in beam_coverage[obsid].keys:
                    # Source was in the beam so include it
                    _, _, _, duration, _, centre_freq, channels = common_metadata_list[on]
                    enter_beam_norm, exit_beam_norm, max_power = beam_coverage[obsid][source_name]
                    source_data.append([obsid, duration,
                                        enter_beam_norm, exit_beam_norm,
                                        max_power, centre_freq, bandwidth])
            # For each source make a dictionary key that contains a list of
            # lists of the data for each obsid
            output_data[source_name] = source_data

    else:
        #output a list of sorces for each obs
        for on, obsid in enumerate(obsid_list):
            _, _, _, duration, _, centre_freq, channels = common_metadata_list[on]
            obsid_data = []
            for source_name in np.array(names_ra_dec)[:,0]:
                if source_name in beam_coverage[obsid].keys():
                    enter_beam_norm, exit_beam_norm, max_power = beam_coverage[obsid][source_name]
                    obsid_data.append([source_name, enter_beam_norm, exit_beam_norm, max_power])
            # For each obsid make a dictionary key that contains a list of
            # lists of the data for each source/pulsar
            output_data[obsid] = obsid_data

    return output_data, common_metadata_list