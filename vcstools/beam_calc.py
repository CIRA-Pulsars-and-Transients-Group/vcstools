import os
import sys
import numpy as np
import math

#MWA scripts
from mwa_pb import primary_beam

from vcstools.pointing_utils import sex2deg
from vcstools.metadb_utils import mwa_alt_az_za, get_common_obs_metadata, getmeta

import logging
logger = logging.getLogger(__name__)


def pixel_area(ra_min, ra_max, dec_min, dec_max):
    """
    Calculate the area of a pixel on the sky from the pixel borders

    Parameters:
    -----------
    ra_min: float
        The Right Acension minimum in degrees
    ra_max: float
        The Right Acension maximum in degrees
    dec_min: float
        The Declination minimum in degrees
    dec_max: float
        The Declination maximum in degrees

    Returns:
    --------
    area: float
        Area of the pixel in square degrees
    """
    return (180./math.pi) * (ra_max - ra_min) * (math.sin(math.radians(dec_max)) - math.sin(math.radians(dec_min)))


def field_of_view(obsid,
                  beam_meta_data=None, dur=None):
    """
    Will find the field-of-view of the observation (including the drift) in square degrees.

    Parameters:
    -----------
    obsid: int
        The observation ID
    beam_meta_data: list
        OPTIONAL - the list of common metadata from vcstools.metadb_utils.get_common_obs_metadata.
        By default will download the metadata for you
    dur: int
        OPTIONAL - Duration of observation to calculate for in seconds
        By default will use the entire observation duration

    Returns:
    --------
    area: float
        The field-of-view of the observation in square degrees
    """
    if beam_meta_data is None:
        beam_meta_data = get_common_obs_metadata(obsid)

    if dur is None:
        dt = 296
    else:
        dt = 100
        # Change the dur to the inpur dur
        obsid, ra, dec, _, delays, centrefreq, channels = beam_meta_data
        beam_meta_data = [obsid, ra, dec, dur, delays, centrefreq, channels]

    # Make a pixel for each degree on the sky
    names_ra_dec = []
    for ra in range(0,360):
        for dec in range(-90,90):
            names_ra_dec.append(["sky_pos", ra+0.5, dec+0.5])

    # Get tile beam power for all pixels
    sky_powers = get_beam_power_over_time(beam_meta_data, names_ra_dec,
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


def from_power_to_gain(powers, cfreq, n, coh=True):
    from astropy.constants import c,k_B
    from math import sqrt

    obswl = c.value/cfreq
    #for coherent
    if coh:
        coeff = obswl**2*16*n/(4*np.pi*k_B.value)
    else:
        coeff = obswl**2*16*sqrt(n)/(4*np.pi*k_B.value)
    logger.debug("Wavelength {} m".format(obswl))
    logger.debug("Gain coefficient: {}".format(coeff))
    SI_to_Jy = 1e-26
    return (powers*coeff)*SI_to_Jy


def get_Trec(tab, obsfreq):
    Trec = 0.0
    for r in range(len(tab)-1):
        if tab[r][0]==obsfreq:
            Trec = tab[r][1]
        elif tab[r][0] < obsfreq < tab[r+1][0]:
            Trec = ((tab[r][1] + tab[r+1][1])/2)
    if Trec == 0.0:
        logger.debug("ERROR getting Trec")
    return Trec


def beam_enter_exit(powers, duration, dt=296, min_power=0.3):
    """
    Calculates when the source enters and exits the beam

    beam_enter_exit(min_power, powers, imax, dt):
        powers: list of powers fo the duration every dt and freq powers[times][freqs]
        dt: the time interval of how often powers are calculated
        duration: duration of the observation according to the metadata in seconds
        min_power: zenith normalised power cut off
    """
    from scipy.interpolate import UnivariateSpline
    time_steps = np.array(range(0, duration, dt), dtype=float)

    #For each time step record the min power so even if the source is in
    #one freq channel it's recorded
    powers_freq_min = []
    for p in powers:
        powers_freq_min.append(float(min(p) - min_power))

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
                #power declines so starts in beem then exits
                enter_beam = 0.
                exit_beam = spline.roots()[0]/duration
            else:
                enter_beam = spline.roots()[0]/duration
                exit_beam = 1.
        else:
            enter_beam = 0.
            exit_beam = 1.
    return enter_beam, exit_beam


def get_beam_power_over_time(beam_meta_data, names_ra_dec,
                             dt=296, centeronly=True, verbose=False,
                             option='analytic', degrees=False,
                             start_time=0):
    """
    Calulates the power (gain at coordinate/gain at zenith) for each source over time.

    get_beam_power_over_time(beam_meta_data, names_ra_dec,
                             dt=296, centeronly=True, verbose=False,
                             option = 'analytic')
    Args:
        beam_meta_data: [obsid,ra, dec, time, delays,centrefreq, channels]
                        obsid metadata obtained from meta.get_common_obs_metadata
        names_ra_dec: and array in the format [[source_name, RAJ, DecJ]]
        dt: time step in seconds for power calculations (default 296)
        centeronly: only calculates for the centre frequency (default True)
        verbose: prints extra data to (default False)
        option: primary beam model [analytic, advanced, full_EE]
        start_time: the time in seconds from the begining of the observation to
                    start calculating at
    """
    obsid, _, _, time, delays, centrefreq, channels = beam_meta_data
    names_ra_dec = np.array(names_ra_dec)
    logger.info("Calculating beam power for OBS ID: {0}".format(obsid))

    starttimes=np.arange(start_time,time+start_time,dt)
    stoptimes=starttimes+dt
    stoptimes[stoptimes>time]=time
    Ntimes=len(starttimes)
    midtimes=float(obsid)+0.5*(starttimes+stoptimes)

    if not centeronly:
        PowersX=np.zeros((len(names_ra_dec),
                             Ntimes,
                             len(channels)))
        PowersY=np.zeros((len(names_ra_dec),
                             Ntimes,
                             len(channels)))
        # in Hz
        frequencies=np.array(channels)*1.28e6
    else:
        PowersX=np.zeros((len(names_ra_dec),
                             Ntimes,1))
        PowersY=np.zeros((len(names_ra_dec),
                             Ntimes,1))
        if centrefreq > 1e6:
            logger.warning("centrefreq is greater than 1e6, assuming input with units of Hz.")
            frequencies=np.array([centrefreq])
        else:
            frequencies=np.array([centrefreq])*1e6
    if degrees:
        RAs = np.array(names_ra_dec[:,1],dtype=float)
        Decs = np.array(names_ra_dec[:,2],dtype=float)
    else:
        RAs, Decs = sex2deg(names_ra_dec[:,1],names_ra_dec[:,2])

    if len(RAs)==0:
        sys.stderr.write('Must supply >=1 source positions\n')
        return None
    if not len(RAs)==len(Decs):
        sys.stderr.write('Must supply equal numbers of RAs and Decs\n')
        return None
    if verbose is False:
        #Supress print statements of the primary beam model functions
        sys.stdout = open(os.devnull, 'w')
    for itime in range(Ntimes):
        # this differ's from the previous ephem_utils method by 0.1 degrees
        _, Azs, Zas = mwa_alt_az_za(midtimes[itime], ra=RAs, dec=Decs, degrees=True)
        # go from altitude to zenith angle
        theta = np.radians(Zas)
        phi = np.radians(Azs)
        for ifreq in range(len(frequencies)):
            #Decide on beam model
            if option == 'analytic':
                rX,rY=primary_beam.MWA_Tile_analytic(theta, phi,
                                                     freq=frequencies[ifreq], delays=delays,
                                                     zenithnorm=True,
                                                     power=True)
            elif option == 'advanced':
                rX,rY=primary_beam.MWA_Tile_advanced(theta, phi,
                                                     freq=frequencies[ifreq], delays=delays,
                                                     zenithnorm=True,
                                                     power=True)
            elif option == 'full_EE':
                rX,rY=primary_beam.MWA_Tile_full_EE(theta, phi,
                                                     freq=frequencies[ifreq], delays=delays,
                                                     zenithnorm=True,
                                                     power=True)
        PowersX[:,itime,ifreq]=rX
        PowersY[:,itime,ifreq]=rY
    if verbose is False:
        sys.stdout = sys.__stdout__
    Powers=0.5*(PowersX+PowersY)
    return Powers


def find_sources_in_obs(obsid_list, names_ra_dec,
                        obs_for_source=False, dt_input=100, beam='analytic',
                        min_power=0.3, cal_check=False, all_volt=False,
                        degrees_check=False, metadata_list=None):
    """
    Either creates text files for each MWA obs ID of each source within it or a text
    file for each source with each MWA obs is that the source is in.
    Args:
        obsid_list: list of MWA obs IDs
        names_ra_dec: [[source_name, ra, dec]]
        dt: the time step in seconds to do power calculations
        beam: beam simulation type ['analytic', 'advanced', 'full_EE']
        min_power: if above the minium power assumes it's in the beam
        cal_check: checks the MWA pulsar database if there is a calibration for the obsid
        all_volt: Use all voltages observations including some inital test data
                  with incorrect formats
        degrees_check: if false ra and dec is in hms, if true in degrees
    Output [output_data, obsid_meta]:
        output_data: The format of output_data is dependant on obs_for_source.
                     If obs_for_source is True:
                        output_data = {jname:[[obsid, duration, enter, exit, max_power],
                                              [obsid, duration, enter, exit, max_power]]}
                     If obs_for_source is False:
                        ouput_data = {obsid:[[jname, enter, exit, max_power],
                                             [jname, enter, exit, max_power]]}
        obsid_meta: a list of the output of get_common_obs_metadata for each obsid
    """
    #prepares metadata calls and calculates power
    powers = []
    #powers[obsid][source][time][freq]
    obsid_meta = []
    obsid_to_remove = []

    for i, obsid in enumerate(obsid_list):
        if metadata_list:
            beam_meta_data, _ = metadata_list[i]
        else:
            beam_meta_data = get_common_obs_metadata(obsid)
        #beam_meta_data = obsid,ra_obs,dec_obs,time_obs,delays,centrefreq,channels

        if dt_input * 4 >  beam_meta_data[3]:
            # If the observation time is very short then a smaller dt time is required
            # to get enough ower imformation
            dt = int(beam_meta_data[3] / 4.)
        else:
            dt = dt_input
        logger.debug("obsid: {0}, time_obs {1} s, dt {2} s".format(obsid, beam_meta_data[3], dt))

        #check for raw volatge files
        filedata = getmeta(service='data_files', params={'obs_id':obsid, 'nocache':1})
        if filedata is None:
            logger.warning("No file data for obsid {}. Skipping".format(obsid))
            obsid_to_remove.append(obsid)
            continue
        keys = filedata.keys()
        check = False
        for k in keys:
            if '.dat' in k: #TODO check if is still robust
                check = True
        if check or all_volt:
            powers.append(get_beam_power_over_time(beam_meta_data, names_ra_dec,
                                    dt=dt, centeronly=True, verbose=False,
                                    option=beam, degrees=degrees_check))
            obsid_meta.append(beam_meta_data)
        else:
            logger.warning('No raw voltage files for %s' % obsid)
            obsid_to_remove.append(obsid)
    for otr in obsid_to_remove:
        obsid_list.remove(otr)

    #chooses whether to list the source in each obs or the obs for each source
    output_data = {}
    if obs_for_source:
        for sn, source in enumerate(names_ra_dec):
            source_data = []
            for on, obsid in enumerate(obsid_list):
                source_ob_power = powers[on][sn]
                if max(source_ob_power) > min_power:
                    duration = obsid_meta[on][3]
                    centre_freq = obsid_meta[on][5] #MHz
                    channels = obsid_meta[on][6]
                    bandwidth = (channels[-1] - channels[0] + 1.)*1.28 #MHz
                    logger.debug("Running beam_enter_exit on obsid: {}".format(obsid))
                    enter_beam, exit_beam = beam_enter_exit(source_ob_power,duration,
                                                  dt=dt, min_power=min_power)
                    if enter_beam is not None:
                        source_data.append([obsid, duration, enter_beam, exit_beam,
                                            max(source_ob_power)[0],
                                            centre_freq, bandwidth])
            # For each source make a dictionary key that contains a list of
            # lists of the data for each obsid
            output_data[source[0]] = source_data

    else:
        #output a list of sorces for each obs
        for on, obsid in enumerate(obsid_list):
            duration = obsid_meta[on][3]
            obsid_data = []
            for sn, source in enumerate(names_ra_dec):
                source_ob_power = powers[on][sn]
                if max(source_ob_power) > min_power:
                    enter_beam, exit_beam = beam_enter_exit(source_ob_power, duration,
                                                  dt=dt, min_power=min_power)
                    obsid_data.append([source[0], enter_beam, exit_beam, max(source_ob_power)[0]])
            # For each obsid make a dictionary key that contains a list of
            # lists of the data for each source/pulsar
            output_data[obsid] = obsid_data

    return output_data, obsid_meta