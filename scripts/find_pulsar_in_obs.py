#! /usr/bin/env python3

"""
Author: Nicholas Swainston
Creation Date: /03/2016

Some of the orignal code was created by Bradley Meyers

This code is used to list the sources within the beam of observations IDs or using --obs_for_source list all the observations for each source. The sources can be input serval ways: using a list of pulsar names (--pulsar), using a complete catalogue file of pulsars (--dl_PSRCAT) or RRATs (--RRAT and --dl_RRAT), using a compatable catalogue (--in_cat with the help of --names and --coordstype) or using a RA and DEC coordinate (--coords). The observation IDs can be input (--obsid) or gathered from a directory (--FITS_dir). The default is to search all observation IDs from http://mwa-metadata01.pawsey.org.au/metadata/ that have voltages and list every known pulsar from PSRCAT in each observation ID.

Two most comon uses for this code is to search a single observation ID for pulsars like so:
> find_pulsar_in_obs.py -o 1099414416
which will output a text file called 1099414416_analytic_beam.txt
or to search all observation IDs for a single pulsar like so:
> find_pulsar_in_obs.py -p J0534+2200 --obs_for_source
which will output a text file called J0534+2200_analytic_beam.txt
"""

__author__ = 'Nicholas Swainston'
__date__ = '2016-03-21'

#TODO: (BWM) might be worth trying to only load the necessary parts of modules rather than loading an ENTIRE module. The speed up will be minimal overall, but it's also just a bit neater.
# You can also just load module/part of modules internally within functions. So if a module is only used once, just load it in the local function rather than up here (globally).

#comon
import os
import sys
import math
import argparse
import numpy as np
import csv

#astropy
from astropy.coordinates import SkyCoord
from astropy import units as u

#MWA scripts
from vcstools import data_load

import sn_flux_est as sfe
from mwa_pb import primary_beam
from mwa_metadb_utils import mwa_alt_az_za, get_common_obs_metadata,\
                             get_obs_array_phase, find_obsids_meta_pages


import logging
logger = logging.getLogger(__name__)


def yes_no(answer):
    yes = set(['Y','yes','y', 'ye', ''])
    no = set(['N','no','n'])

    while True:
        choice = input(answer).lower()
        if choice in yes:
           return True
        elif choice in no:
           return False
        else:
           logger.warning("Please respond with 'yes' or 'no'\n")


def sex2deg(ra, dec):
    """
    Convert sexagesimal coordinates to degrees.

    sex2deg( ra, dec)
    Args:
        ra: the right ascension in HH:MM:SS
        dec: the declination in DD:MM:SS
    """
    c = SkyCoord( ra, dec, frame='icrs', unit=(u.hourangle,u.deg))

    # return RA and DEC in degrees in degrees
    return [c.ra.deg, c.dec.deg]


def deg2sex(ra, dec):
    """
    Convert decimal coordingates into sexagesimal strings, i.e. hh:mm:ss.ss and dd:mm:ss.ss

    deg2sex( ra, dec)
    Args:
        ra: the right ascension in degrees
        dec: the declination in degrees
    """

    c = SkyCoord( ra, dec, frame='icrs', unit=(u.deg,u.deg))
    #coords = c.to_string('hmsdms')
    #coords = coords.replace('h',':').replace('d',':').replace('m',':').replace('s','')
    rajs = c.ra.to_string(unit=u.hour, sep=':')
    decjs = c.dec.to_string(unit=u.degree, sep=':')

    # return RA and DEC in "hh:mm:ss.ssss dd:mm:ss.ssss" form
    return rajs, decjs


def get_psrcat_ra_dec(pulsar_list=None, max_dm=1000., include_dm=False, query=None):
    """
    Uses PSRCAT to return a list of pulsar names, ras and decs. Not corrected for proper motion.
    Removes pulsars without any RA or DEC recorded
    If no pulsar_list given then returns all pulsar on the catalogue
    If include_dm is True then also ouput DM

    get_psrcat_ra_dec(pulsar_list = None)
    Args:
        pulsar_list: A space list of pulsar names eg: [J0534+2200, J0538+2817].
               (default: uses all pulsars)
    return [[Jname, RAJ, DecJ]]
    """
    import psrqpy

    #params = ['JNAME', 'RAJ', 'DECJ', 'DM']
    if query is None:
        query = psrqpy.QueryATNF(params = ['PSRJ', 'RAJ', 'DECJ', 'DM'], psrs=pulsar_list, loadfromdb=data_load.ATNF_LOC).pandas

    pulsar_ra_dec = []
    for i, _ in enumerate(query["PSRJ"]):
        # Only record if under the max_dm
        dm = query["DM"][i]
        if not math.isnan(dm):
            if float(dm) < max_dm:
                if include_dm:
                    pulsar_ra_dec.append([query["PSRJ"][i], query["RAJ"][i], query["DECJ"][i], dm])
                else:
                    pulsar_ra_dec.append([query["PSRJ"][i], query["RAJ"][i], query["DECJ"][i]])


    return pulsar_ra_dec


def grab_source_alog(source_type='Pulsar', pulsar_list=None, max_dm=1000., include_dm=False, query=None):
    """
    Will search different source catalogues and extract all the source names, RAs, Decs and, if requested, DMs

    Parameters
    ----------
    source_type: string
        The type of source you would like to get the catalogue for.
        Your choices are: ['Pulsar', 'FRB', 'rFRB', 'POI' 'RRATs', 'Fermi']
        Default: 'Pulsar'
    pulsar_list: list
        List of sources you would like to extract data for.
        If None is given then it will search for all available sources
        Default: None
    max_dm: float
        If the source_type is 'Pulsar' then you can set a maximum dm and the function will only
        return pulsars under that value.
        Default: 1000.
    include_dm: Bool
        If True the function will also return the DM if it is available at the end of the list
        Default: False

    Returns
    -------
    result: list list
        A list for each source which contains a [source_name, RA, Dec (, DM)]
        where RA and Dec are in the format HH:MM:SS
    """
    modes = ['Pulsar', 'FRB', 'rFRB', 'POI', 'RRATs', 'Fermi']
    if source_type not in modes:
        logger.error("Input source type not in known catalogues types. Please choose from: {0}".format(modes))
        return None

    #Get each source type into the format [[name, ra, dec]]
    name_ra_dec = []
    if source_type == 'Pulsar':
        name_ra_dec = get_psrcat_ra_dec(pulsar_list=pulsar_list, max_dm=max_dm, include_dm=include_dm, query=query)

    elif source_type == 'FRB':
        import json
        import urllib.request
        try:
            frb_data = json.load(urllib.request.urlopen('http://frbcat.org/products?search=&min=0&max=1000&page=1'))
        except urllib.error.HTTPError:
            logger.error('http://frbcat.org/ not available. Returning empty list')
            # putting and FRB at 90 dec which we should never be able to detect
            name_ra_dec = [['fake', "00:00:00.00", "90:00:00.00", 0.0]]
        else:
            for frb in frb_data['products']:
                name = frb['frb_name']
                logger.debug('FRB name: {}'.format(name))
                ra   = frb['rop_raj']
                dec  = frb['rop_decj']
                dm   = frb['rmp_dm'].split("&")[0]
                if include_dm:
                    name_ra_dec.append([name, ra, dec, dm])
                else:
                    name_ra_dec.append([name, ra, dec])

    elif source_type == "rFRB":
        info = get_rFRB_info(name=pulsar_list)
        if info is not None:
            for line in info:
                if include_dm:
                    name_ra_dec.append([line[0], line[1], line[2], line[3]])
                else:
                    name_ra_dec.append([line[0], line[1], line[2]])

    elif source_type == "POI":
        #POI = points of interest
        db = open(data_load.POI_CSV, "r")
        for line in db.readlines():
            if not line.startswith("#"):
                line = line.split(",")
                name_ra_dec.append([line[0], line[1], line[2]])

    elif source_type == 'RRATs':
        import urllib.request
        try:
            rrats_data = urllib.request.urlopen('http://astro.phys.wvu.edu/rratalog/rratalog.txt').read().decode()
        except urllib.error.URLError:
            logger.error('http://astro.phys.wvu.edu/rratalog/ not available. Returning empty list')
            # putting and RRAT at 90 dec which we should never be able to detect
            name_ra_dec = [['fake', "00:00:00.00", "90:00:00.00", 0.0]]
        else:
            for rrat in rrats_data.split("\n")[1:-1]:
                columns = rrat.strip().replace(" ", '\t').split('\t')
                rrat_cat_line = []
                if pulsar_list == None or (columns[0] in pulsar_list):
                    for entry in columns:
                        if entry not in ['', ' ', '\t']:
                            rrat_cat_line.append(entry.replace('--',''))
                    #Removing bad formating for the database
                    ra = rrat_cat_line[4]
                    if ra.endswith(":"):
                        ra = ra[:-1]
                    dec = rrat_cat_line[5]
                    if dec.endswith(":"):
                        dec = dec[:-1]

                    if include_dm:
                        name_ra_dec.append([rrat_cat_line[0], ra, dec, rrat_cat_line[3]])
                    else:
                        name_ra_dec.append([rrat_cat_line[0], ra, dec])

    elif source_type == 'Fermi':
        # read the fermi targets file
        try:
            fermi_loc = os.environ['FERMI_CAND_FILE']
        except:
            logger.warning("Fermi candidate file location not found. Returning nothing")
            return []
        with open(fermi_loc,"r") as fermi_file:
            csv_reader = csv.DictReader(fermi_file)
            for fermi in csv_reader:
                name = fermi['Source Name'].split()[-1]
                ra  = fermi[' RA J2000']
                dec = fermi[' Dec J2000']
                raj, decj = deg2sex(float(ra), float(dec))
                pos_u = float(fermi[' a (arcmin)'])
                if include_dm:
                    # this actually returns the position uncertainty not dm
                    name_ra_dec.append([name, raj, decj, pos_u])
                else:
                    name_ra_dec.append([name, raj, decj])

    #Remove all unwanted sources
    if pulsar_list is not None:
        rrat_cat_filtered = []
        for line in name_ra_dec:
            if line[0] in pulsar_list:
                rrat_cat_filtered.append(line)
        name_ra_dec = rrat_cat_filtered

    return name_ra_dec


def get_rFRB_info(name=None):
    """
    Gets repeating FRB info from the csv file we maintain.
    Input:
        name: a list of repeating FRB names to get info for. The default is None
              which gets all rFRBs in the catalogue.
    Output:
        [[name, ra, dec, dm, dm error]]
        A list of all the FRBs which each have a list contining name, ra, dec,
        dm and dm error
    """
    output = []
    db = open(data_load.KNOWN_RFRB_CSV, "r")
    for line in db.readlines():
        if not line.startswith("#"):
            line = line.split(",")
            #some FRBs end with a J name. We will ignore these when comparing
            #by using the first 9 characters
            FRB = line[0]
            if name is None:
                #No input FRBs so return all FRBs
                output.append(line)
            elif FRB in name:
                output.append(line)
    return output


def format_ra_dec(ra_dec_list, ra_col=0, dec_col=1):
    """
    Will format a list of lists containing RAs and Decs to uniform strings.  eg 00:00:00.00 -00:00:00.00. Will not work for numpy arrays so make sure they're list of lists
    An example input:
    format_ra_dec([[name,ra,dec]], ra_col = 1, dec_col = 2)
    """
    for i in range(len(ra_dec_list)):
        #catching errors in old psrcat RAs
        if len(ra_dec_list[i][ra_col]) >5:
            if  ra_dec_list[i][ra_col][5] == '.':
                ra_dec_list[i][ra_col] = ra_dec_list[i][ra_col][:5] + ":" +\
                        str(int(float('0'+ra_dec_list[i][ra_col][5:])*60.))
        #make sure there are two digits in the HH slot for RA or DD for Dec
        if len(ra_dec_list[i][ra_col].split(':')[0]) == 1:
            ra_dec_list[i][ra_col] = '0' + ra_dec_list[i][ra_col]
        if ra_dec_list[i][dec_col][0].isdigit():
            if len(ra_dec_list[i][dec_col].split(':')[0]) == 1:
                ra_dec_list[i][dec_col] = '0' + ra_dec_list[i][dec_col]
        else:
            if len(ra_dec_list[i][dec_col].split(':')[0]) == 2:
                ra_dec_list[i][dec_col] = ra_dec_list[i][dec_col][0] + '0' +\
                                          ra_dec_list[i][dec_col][1:]

        #make sure there is a + on positive Decs
        if not ra_dec_list[i][dec_col].startswith('-') and\
           not ra_dec_list[i][dec_col].startswith('+'):
               ra_dec_list[i][dec_col] = '+' + ra_dec_list[i][dec_col]

        #since the ra a dec are the same except the +- in the dec just loop over it
        #and change the length
        ra_dec_col = [ra_col, dec_col]
        for n in range(2):
            if len(ra_dec_list[i][ra_dec_col[n]]) == int(2 + n):
                ra_dec_list[i][ra_dec_col[n]] += ':00:00.00'
            elif len(ra_dec_list[i][ra_dec_col[n]]) == 3 + n:
                ra_dec_list[i][ra_dec_col[n]] += '00:00.00'
            elif len(ra_dec_list[i][ra_dec_col[n]]) == 4 + n:
                ra_dec_list[i][ra_dec_col[n]] += '0:00.00'
            elif len(ra_dec_list[i][ra_dec_col[n]]) == 5 + n:
                ra_dec_list[i][ra_dec_col[n]] += ':00.00'
            elif len(ra_dec_list[i][ra_dec_col[n]]) == 6 + n:
                ra_dec_list[i][ra_dec_col[n]] += '00.00'
            elif len(ra_dec_list[i][ra_dec_col[n]]) == 7 + n:
                ra_dec_list[i][ra_dec_col[n]] += '0.00'
            elif len(ra_dec_list[i][ra_dec_col[n]]) == 8 + n:
                ra_dec_list[i][ra_dec_col[n]] += '.00'
            elif len(ra_dec_list[i][ra_dec_col[n]]) == 9 + n:
                ra_dec_list[i][ra_dec_col[n]] += '00'
            elif len(ra_dec_list[i][ra_dec_col[n]]) == 10 + n:
                ra_dec_list[i][ra_dec_col[n]] += '0'

            #shorten if too long
            if len(ra_dec_list[i][ra_dec_col[n]]) > (11 + n):
                ra_dec_list[i][ra_dec_col[n]] = ra_dec_list[i][ra_dec_col[n]][:(11+n)]

    return ra_dec_list


def singles_source_search(ra, dec):
    """
    Used to creates a 30 degree box around the source to make searching for obs_ids more efficient

    singles_source_search(ra, dec)
    Args:
        ra: ra of source in degrees
        dec: dec of source in degrees
    """
    ra = float(ra)
    dec = float(dec)
    box_size = 30.
    m_o_p = False # moved over (north or south) pole

    dec_top = dec + box_size
    if dec_top > 90.:
        dec_top = 90.
        m_o_p = True

    dec_bot = dec - box_size
    if dec_top < -90.:
        dec_top = -90.
        m_o_p = True

    if m_o_p:
        obsid_list = find_obsids_meta_pages(params={'mode':'VOLTAGE_START',
                                                    'minra':0., 'maxra':360.,
                                                    'mindec':dec_bot,'maxdec':dec_top})
    else:
        ra_low = ra - 30. - box_size #30 is the how far an obs would drift in 2 hours(used as a max)
        ra_high = ra + box_size
        if ra_low < 0.:
            ra_new = 360 + ra_low
            obsid_list = find_obsids_meta_pages(params={'mode':'VOLTAGE_START',
                                                        'minra':ra_new, 'maxra':360.,
                                                        'mindec':dec_bot,'maxdec':dec_top})
            temp_obsid_list = find_obsids_meta_pages(params={'mode':'VOLTAGE_START',
                                                             'minra':0.,'maxra':ra_high,
                                                             'mindec':dec_bot,'maxdec':dec_top})
            for row in temp_obsid_list:
                obsid_list.append(row)
        elif ra_high > 360:
            ra_new = ra_high - 360
            obsid_list = find_obsids_meta_pages(params={'mode':'VOLTAGE_START',
                                                        'minra':ra_low, 'maxra':360.,
                                                        'mindec':dec_bot,'maxdec':dec_top})
            temp_obsid_list = find_obsids_meta_pages(params={'mode':'VOLTAGE_START',
                                                             'minra':0., 'maxra':ra_new,
                                                             'mindec':dec_bot,'maxdec':dec_top})
            for row in temp_obsid_list:
                obsid_list.append(row)
        else:
            obsid_list = find_obsids_meta_pages(params={'mode':'VOLTAGE_START',
                                                        'minra':ra_low, 'maxra':ra_high,
                                                        'mindec':dec_bot,'maxdec':dec_top})
    return obsid_list


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
        enter = 0.
        exit = 1.
    else:
        powers_freq_min = np.array(powers_freq_min)
        logger.debug("time_steps: {}".format(time_steps))
        logger.debug("powers: {}".format(powers_freq_min))
        try:
            spline = UnivariateSpline(time_steps, powers_freq_min , s=0.)
        except:
            return None, None
        if len(spline.roots()) == 2:
            enter, exit = spline.roots()
            enter /= duration
            exit /= duration
        elif len(spline.roots()) == 1:
            if powers_freq_min[0] > powers_freq_min[-1]:
                #power declines so starts in beem then exits
                enter = 0.
                exit = spline.roots()[0]/duration
            else:
                enter = spline.roots()[0]/duration
                exit = 1.
        else:
            enter = 0.
            exit = 1.
    return enter, exit


def cal_on_database_check(obsid):
    from mwa_pulsar_client import client
    web_address = 'https://mwa-pawsey-volt01.pawsey.org.au'
    auth = ('mwapulsar','veovys9OUTY=')
    detection_list = client.detection_list(web_address, auth)

    cal_used = False
    cal_avail = False
    for d in detection_list:
        if int(d[u'observationid']) == int(obsid):
            if d[u'calibrator'] is not None:
                cal_used = True
                #TODO add a check if there is a cal file option

    #No cal
    check_result = 'N'
    if cal_avail:
        #Cal available
        check_result = 'A'
    elif cal_used:
        #Cal used
        check_result = 'U'

    return check_result


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
            beam_meta_data, full_meta = metadata_list[i]
        else:
            beam_meta_data, full_meta = get_common_obs_metadata(obsid, return_all=True)
        #beam_meta_data = obsid,ra_obs,dec_obs,time_obs,delays,centrefreq,channels

        if dt_input * 4 >  beam_meta_data[3]:
            # If the observation time is very short then a smaller dt time is required
            # to get enough ower imformation
            dt = int(beam_meta_data[3] / 4.)
        else:
            dt = dt_input
        logger.debug("obsid: {0}, time_obs {1} s, dt {2} s".format(obsid, beam_meta_data[3], dt))

        #check for raw volatge files
        filedata = full_meta[u'files']
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
                    enter, exit = beam_enter_exit(source_ob_power,duration,
                                                  dt=dt, min_power=min_power)
                    if enter is not None:
                        source_data.append([obsid, duration, enter, exit,
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
                    enter, exit = beam_enter_exit(source_ob_power, duration,
                                                  dt=dt, min_power=min_power)
                    obsid_data.append([source[0], enter, exit, max(source_ob_power)[0]])
            # For each obsid make a dictionary key that contains a list of
            # lists of the data for each source/pulsar
            output_data[obsid] = obsid_data

    return output_data, obsid_meta

def write_output_source_files(output_data,
                              beam='analytic', min_power=0.3, cal_check=False,
                              SN_est=False, plot_est=False):
    """
    Writes an ouput file using the output of find_sources_in_obs when obs_for_source is true.
    """
    for source in output_data:
        out_name = "{0}_{1}_beam.txt".format(source, beam)
        with open(out_name,"w") as output_file:
            output_file.write('#All of the observation IDs that the {0} beam model '
                              'calculated a power of {1} or greater for the source: '
                              '{2}\n'.format(beam, min_power, source))
            output_file.write('#Column headers:\n')
            output_file.write('#Obs ID: Observation ID\n')
            output_file.write('#Dur:    The duration of the observation in seconds\n')
            output_file.write('#Enter:  The fraction of the observation when '
                                        'the source entered the beam\n')
            output_file.write('#Exit:   The fraction of the observation when '
                                        'the source exits the beam\n')
            output_file.write('#Power:  The maximum zenith normalised power of the source.\n')
            output_file.write("#OAP:    The observation's array phase where P1 is the "
                                        "phase 1 array, P2C is the phase compact array "
                                        "and P2E is the phase 2 extended array.\n")
            output_file.write("#Freq:   The centre frequency of the observation in MHz\n")
            output_file.write("#Band:   Bandwidth of the observation in MHz. If it is greater "
                                        "than 30.72 than it is a picket fence observation\n")

            if SN_est:
                output_file.write("#S/N Est: An estimate of the expected signal to noise using ANTF flux desnities\n")
                output_file.write("#S/N Err: The uncertainty of S/N Est\n")
            if cal_check:
                output_file.write('#Cal ID: Observation ID of an available '+\
                                            'calibration solution\n')
            output_file.write('#Obs ID   |Dur |Enter|Exit |Power| OAP | Freq | Band ')
            if SN_est:
                output_file.write("|S/N Est|S/N Err")

            if cal_check:
                output_file.write("|Cal ID\n")
            else:
                output_file.write('\n')
            for data in output_data[source]:
                obsid, duration, enter, leave, max_power, freq, band = data
                oap = get_obs_array_phase(obsid)
                output_file.write('{} {:4d} {:1.3f} {:1.3f} {:1.3f}  {:.3}   {:6.2f} {:6.2f}'.\
                           format(obsid, duration, enter, leave, max_power, oap, freq, band))
                if SN_est:
                    pulsar_sn, pulsar_sn_err = sfe.est_pulsar_sn(source, obsid, plot_flux=plot_est)
                    if pulsar_sn is None:
                        output_file.write('   None    None')
                    else:
                        output_file.write('{:9.2f} {:9.2f}'.format(pulsar_sn, pulsar_sn_err))

                if cal_check:
                    #checks the MWA Pulsar Database to see if the obsid has been
                    #used or has been calibrated
                    logger.info("Checking the MWA Pulsar Databse for the obsid: {0}".format(obsid))
                    cal_check_result = cal_on_database_check(obsid)
                    output_file.write("   {0}\n".format(cal_check_result))
                else:
                    output_file.write("\n")
    return


def write_output_obs_files(output_data, obsid_meta,
                           beam='analytic', min_power=0.3,
                           cal_check=False, SN_est=False, plot_est=False):
    """
    Writes an ouput file using the output of find_sources_in_obs when obs_for_source is false.
    """

    for on, obsid in enumerate(output_data):
        if SN_est:
            psr_list = [el[0] for el in output_data[obsid]]
            sn_dict = sfe.multi_psr_snfe(psr_list, obsid,\
                                         min_z_power=min_power, plot_flux=plot_est)

        oap = get_obs_array_phase(obsid)
        out_name = "{0}_{1}_beam.txt".format(obsid, beam)
        with open(out_name,"w") as output_file:
            output_file.write('#All of the sources that the {0} beam model calculated a power'
                              'of {1} or greater for observation ID: {2}\n'.format(beam,
                              min_power, obsid))
            output_file.write('#Observation data :RA(deg): {0} DEC(deg): {1} Duration(s): '
                              '{2} Array Phase: {3}\n'.format(obsid_meta[on][1],
                                  obsid_meta[on][2], obsid_meta[on][3], oap))
            if cal_check:
                #checks the MWA Pulsar Database to see if the obsid has been
                #used or has been calibrated
                logger.info("Checking the MWA Pulsar Databse for the obsid: {0}".format(obsid))
                cal_check_result = cal_on_database_check(obsid)
                output_file.write("#Calibrator Availability: {0}\n".format(cal_check_result))
            output_file.write('#Column headers:\n')
            output_file.write('#Source: Pulsar Jname\n')
            output_file.write('#Enter:  The fraction of the observation when '+\
                                        'the source entered the beam\n')
            output_file.write('#Exit:   The fraction of the observation when '+\
                                        'the source exits the beam\n')
            output_file.write('#Power:  The maximum zenith normalised power of the source.\n')
            if SN_est:
                output_file.write("#S/N Est: An estimate of the expected signal to noise using ANTF flux desnities\n")
                output_file.write("#S/N Err: The uncertainty of S/N Est\n")

            output_file.write('#Source    |Enter|Exit |Power')
            if SN_est:
                output_file.write('| S/N Est | S/N Err \n')
            else:
                output_file.write('\n')

            for data in output_data[obsid]:
                pulsar, enter, exit, max_power = data
                output_file.write('{:11} {:1.3f} {:1.3f} {:1.3f} '.format(pulsar,
                                  enter, exit, max_power))
                if SN_est:
                    beg = int(obsid) + 7
                    end = beg + int(obsid_meta[on][3])
                    pulsar_sn, pulsar_sn_err = sn_dict[pulsar]
                    if pulsar_sn is None:
                        output_file.write('   None    None\n')
                    else:
                        output_file.write('{:9.2f} {:9.2f}\n'.format(pulsar_sn, pulsar_sn_err))
                else:
                    output_file.write('\n')

    return


if __name__ == "__main__":
    # Dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING)
    beam_models = ['analytic', 'advanced', 'full_EE']
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""
    This code is used to list the sources within the beam of observations IDs or using --obs_for_source list all the observations for each source. The sources can be input serval ways: using a list of pulsar names (--pulsar), using a complete catalogue file of pulsars (--dl_PSRCAT) or RRATs (--RRAT and --dl_RRAT), using a compatable catalogue (--in_cat with the help of --names and --coordstype) or using a RA and DEC coordinate (--coords). The observation IDs can be input (--obsid) or gathered from a directory (--FITS_dir). The default is to search all observation IDs from http://mwa-metadata01.pawsey.org.au/metadata/ that have voltages and list every known pulsar from PSRCAT in each observation ID.
    """)
    parser.add_argument('--obs_for_source',action='store_true',help='Instead of listing all the sources in each observation it will list all of the observations for each source. For increased efficiency it will only search OBSIDs within the primary beam.')
    parser.add_argument('--output',type=str,help='Chooses a file for all the text files to be output to. The default is your current directory', default = './')
    parser.add_argument('-b','--beam',type=str, default = 'analytic', help='Decides the beam approximation that will be used. Options: "analytic" the analytic beam model (2012 model, fast and reasonably accurate), "advanced" the advanced beam model (2014 model, fast and slighty more accurate) or "full_EE" the full EE model (2016 model, slow but accurate). " Default: "analytic"')
    parser.add_argument('-m','--min_power',type=float,help='The minimum fraction of the zenith normalised power that a source needs to have to be recorded. Default 0.3', default=0.3)
    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO",
                                    choices=loglevels.keys(), default="INFO")
    parser.add_argument("-V", "--version", action="store_true", help="Print version and quit")

    #source options
    sourargs = parser.add_argument_group('Source options', 'The different options to control which sources are used. Default is all known pulsars.')
    sourargs.add_argument('-p','--pulsar',type=str, nargs='*',help='Searches for all known pulsars. This is the default. To search for individual pulsars list their Jnames in the format " -p J0534+2200 J0630-2834"', default = None)
    sourargs.add_argument('--max_dm',type=float, default = 250., help='The maximum DM for pulsars. All pulsars with DMs higher than the maximum will not be included in output files. Default=250.0')
    sourargs.add_argument('--source_type',type=str, default = 'Pulsar', help="An astronomical source type from ['Pulsar', 'FRB', 'rFRB', 'GC', 'RRATs', Fermi] to search for all sources in their respective web catalogue.")
    sourargs.add_argument('--in_cat',type=str,help='Location of source catalogue, must be a csv where each line is in the format "source_name, hh:mm:ss.ss, +dd:mm:ss.ss".')
    sourargs.add_argument('-c','--coords',type=str,nargs='*',help='String containing the source\'s coordinates to be searched for in the format "RA,DEC" "RA,DEC". Must be enterered as either: "hh:mm:ss.ss,+dd:mm:ss.ss" or "deg,-deg". Please only use one format.')
    #finish above later and make it more robust to incclude input as sex or deg and perhaps other coordinte systmes

    #observation options
    obargs = parser.add_argument_group('Observation ID options', 'The different options to control which observation IDs are used. Default is all observation IDs with voltages.')
    obargs.add_argument('--FITS_dir',type=str,help='Instead of searching all OBS IDs, only searchs for the obsids in the given directory. Does not check if the .fits files are within the directory. Default = /group/mwavcs/vcs')
    obargs.add_argument('-o','--obsid',type=int,nargs='*',help='Input several OBS IDs in the format " -o 1099414416 1095506112". If this option is not input all OBS IDs that have voltages will be used')
    obargs.add_argument('--all_volt',action='store_true',help='Includes observation IDs even if there are no raw voltages in the archive. Some incoherent observation ID files may be archived even though there are raw voltage files. The default is to only include files with raw voltage files.')
    obargs.add_argument('--cal_check',action='store_true',help='Check the MWA Pulsar Database to check if the obsid has every succesfully detected a pulsar and if it has a calibration solution.')
    obargs.add_argument('--sn_est',action='store_true',help='Make a expected signal to noise calculation using the flux densities from the ANTF pulsar catalogue and include them in the output file. Default: False.')
    obargs.add_argument('--plot_est',action='store_true',help='If used, will output flux estimation plots while sn_est arg is true. Default: False.')
    args=parser.parse_args()



    if args.version:
        try:
            import version
            print(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            print("Couldn't import version.py - have you installed vcstools?")
            print("ImportError: {0}".format(ie))
            sys.exit(0)

    # set up the logger for stand-alone execution
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.propagate = False

    #Parse options
    if args.in_cat and args.coords:
        logger.error("Can't use --in_cat and --coords. Please input your cooridantes "
                     "using one method. Exiting.")
        quit()
    if args.obsid and args.FITS_dir:
        logger.error("Can't use --obsid and --FITS_dir at the same time. Exiting.")
        quit()
    if args.beam not in beam_models:
        logger.error("Unknown beam model. Please use one of {0}. Exiting.".format(beam_models))
        quit()

    logger.info("Gathering sources")
    degrees_check = False
    if args.in_cat:
        names_ra_dec = []
        with open(args.in_cat,"r") as input_catalogue:
            reader = csv.reader(input_catalogue)
            for row in reader:
                names_ra_dec.append(row)
    elif args.coords:
        names_ra_dec = []
        for cn, c in enumerate(args.coords):
            names_ra_dec.append(['{0}_{1}'.format(c.split(',')[0], c.split(',')[1]),
                                c.split(',')[0], c.split(',')[1]])
            if ":" not in c:
                degrees_check = True
    else:
        names_ra_dec = grab_source_alog(source_type=args.source_type,
                                        pulsar_list=args.pulsar,
                                        max_dm=args.max_dm)

    #format ra and dec
    if not degrees_check:
        names_ra_dec = format_ra_dec(names_ra_dec, ra_col=1, dec_col=2)
    names_ra_dec = np.array(names_ra_dec)

    #Check if the user wants to use --obs for source
    if (len(names_ra_dec) ==1) and (not args.obs_for_source):
        args.obs_for_source = yes_no('You are only searching for one pulsar so it is '+\
                                     'recommened that you use --obs_for_source. Would '+\
                                     'you like to use --obs_for_source. (Y/n)')
        if args.obs_for_source:
            logger.info("Using option --obs_for_source")
        else:
            logger.info("Not using option --obs_for_source")

    #get obs IDs
    logger.info("Gathering observation IDs")
    if args.obsid:
        obsid_list = args.obsid
    elif args.FITS_dir:
        obsid_list = os.walk(args.FITS_dir).next()[1]
    elif len(names_ra_dec) == 1:
        #if there is a single pulsar simply use nearby obsids
        if degrees_check:
            ob_ra = names_ra_dec[0][1]
            ob_dec = names_ra_dec[0][2]
        else:
            ob_ra, ob_dec = sex2deg(names_ra_dec[0][1], names_ra_dec[0][2])
        obsid_list = singles_source_search(ob_ra, ob_dec)
    else:
        #use all obsids
        obsid_list = find_obsids_meta_pages({'mode':'VOLTAGE_START'})


    if args.beam == 'full_EE':
        dt = 300
    else:
        dt = 100

    logger.debug("names_ra_dec:{}".format(names_ra_dec))
    logger.info("Getting observation metadata and calculating the tile beam")
    output_data, obsid_meta = find_sources_in_obs(obsid_list, names_ra_dec,
                                obs_for_source=args.obs_for_source, dt_input=dt,
                                beam=args.beam, min_power=args.min_power,
                                cal_check=args.cal_check, all_volt=args.all_volt,
                                degrees_check=degrees_check)

    logger.info("Writing data to files")
    if args.obs_for_source:
        write_output_source_files(output_data,
                                  beam=args.beam, min_power=args.min_power,
                                  cal_check=args.cal_check,
                                  SN_est=args.sn_est, plot_est=args.plot_est)
    else:
        write_output_obs_files(output_data, obsid_meta,
                               beam=args.beam, min_power=args.min_power,
                               cal_check=args.cal_check,
                               SN_est=args.sn_est, plot_est=args.plot_est)
