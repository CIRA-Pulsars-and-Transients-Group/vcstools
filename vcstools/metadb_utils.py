from vcstools.general_utils import is_number
import logging
import os
import subprocess

logger = logging.getLogger(__name__)


def ensure_metafits(data_dir, obs_id, metafits_file):
    # TODO: To get the actual ppds file should do this with obsdownload -o <obsID> -m

    if not os.path.exists(metafits_file):
        logger.warning("{0} does not exists".format(metafits_file))
        logger.warning("Will download it from the archive. This can take a "
                      "while so please do not ctrl-C.")
        logger.warning("At the moment, even through the downloaded file is "
                       "labelled as a ppd file this is not true.")
        logger.warning("This is hopefully a temporary measure.")

        get_metafits = "wget http://ws.mwatelescope.org/metadata/fits?obs_id={0} -O {1}".format(obs_id, metafits_file)
        try:
            subprocess.call(get_metafits,shell=True)
        except:
            logger.error("Couldn't download {0}. Aborting.".\
                          format(os.basename(metafits_file)))
            sys.exit(0)
        # clean up
        #os.remove('obscrt.crt')
        #os.remove('obskey.key')
    # make a copy of the file in the product_dir if that directory exists
    # if it doesn't we might have downloaded the metafits file of a calibrator (obs_id only exists on /astro)
    # in case --work_dir was specified in process_vcs call product_dir and data_dir
    # are the same and thus we will not perform the copy
    #data_dir = data_dir.replace(comp_config['base_data_dir'], comp_config['base_product_dir']) # being pedantic
    if os.path.exists(data_dir) and not os.path.exists(metafits_file):
        logger.info("Copying {0} to {1}".format(metafits_file, data_dir))
        from shutil import copy2
        copy2("{0}".format(metafits_file), "{0}".format(data_dir))


def singles_source_search(ra, dec=None, box_size=45.):
    """
    Used to find all obsids within a box around the source to make searching through obs_ids more efficient.

    singles_source_search(ra, dec=None, box_size=45.)
    Parameters:
    ----------
    ra: float
        Right Acension of the source in degrees
    dec: float
        Declination of the source in degrees. By default will use the enitre declination range to account for grating lobes
    box_size: float
        Radius of the search box. Default: 45

    Returns:
    --------
    obsid_metadata: list
        List of of the metadata for each obsid. The metadata is in the same format as getmeta's output
    """
    ra = float(ra)
    m_o_p = False # moved over (north or south) pole

    if dec is None:
        dec_top = 90.
        dec_bot = -90.
    else:
        dec = float(dec)
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

def find_obsids_meta_pages(params=None):
    """
    Loops over pages for each page for MWA metadata calls
    """
    if params is None:
        params = {'mode':'VOLTAGE_START'}
    obsid_list = []
    temp =[]
    page = 1
    #need to ask for a page of results at a time
    while len(temp) == 200 or page == 1:
        params['page'] = page
        logger.debug("Page: {0}   params: {1}".format(page, params))
        temp = getmeta(service='find', params=params)
        if temp is not None:
            # if there are non obs in the field (which is rare) None is returned
            for row in temp:
                obsid_list.append(row[0])
        else:
            temp = []
        page += 1

    return obsid_list

def get_obs_array_phase(obsid):
    """
    For the input obsid will work out the observations array phase in the form
    of P1 for phase 1, P2C for phase 2 compact or P2E for phase to extended array
    and OTH for other.
    """
    phase_info = getmeta(service='con', params={'obs_id':obsid, 'summary':''})

    if phase_info[0] == "PHASE1":
        return "P1"
    elif phase_info[0] == "COMPACT":
        return "P2C"
    elif phase_info[0] == "LB":
        return "P2E"
    elif phase_info[0] == "OTHER":
        return "OTH"
    else:
        logger.error("Unknown phase: {0}. Exiting".format(phase_info[0]))
        exit()


def mwa_alt_az_za(obsid, ra=None, dec=None, degrees=False):
    """
    Calculate the altitude, azumith and zenith for an obsid

    Args:
        obsid  : The MWA observation id (GPS time)
        ra     : The right acension in HH:MM:SS
        dec    : The declintation in HH:MM:SS
        degrees: If true the ra and dec is given in degrees (Default:False)
    """
    from astropy.utils import iers
    iers.IERS_A_URL = 'https://datacenter.iers.org/data/9/finals2000A.all'
    logger.debug(iers.IERS_A_URL)

    from astropy.time import Time
    from astropy.coordinates import SkyCoord, AltAz, EarthLocation
    from astropy import units as u
    obstime = Time(float(obsid),format='gps')

    if ra is None or dec is None:
        #if no ra and dec given use obsid ra and dec
        ra, dec = get_common_obs_metadata(obsid)[1:3]

    if degrees:
        sky_posn = SkyCoord(ra, dec, unit=(u.deg,u.deg))
    else:
        sky_posn = SkyCoord(ra, dec, unit=(u.hourangle,u.deg))
    #earth_location = EarthLocation.of_site('Murchison Widefield Array')
    earth_location = EarthLocation.from_geodetic(lon="116:40:14.93", lat="-26:42:11.95", height=377.8)
    altaz = sky_posn.transform_to(AltAz(obstime=obstime, location=earth_location))
    Alt = altaz.alt.deg
    Az  = altaz.az.deg
    Za  = 90. - Alt
    return Alt, Az, Za


def get_common_obs_metadata(obsid, return_all=False, full_meta_data=None):
    """
    Gets needed common meta data from http://ws.mwatelescope.org/metadata/

    Parameters:
    -----------
    obsid: int
        The observation ID.
    return_all: bool
        OPTIONAL - If True will also return the full meta data dictionary. Default: False
    full_meta_data: bool
        OPTIONAL - The full meta data dictionary from getmeta. If this is not supplied will make the meta data call. Default: False.

    Returns:
    --------
    common_meta_data: list
        [obsid, ra, dec, dura, [xdelays, ydelays], centrefreq, channels]
    """
    if full_meta_data is None:
        logger.info("Obtaining metadata from http://ws.mwatelescope.org/metadata/ for OBS ID: " + str(obsid))
        full_meta_data = getmeta(service='obs', params={'obs_id':obsid})
    ra = full_meta_data[u'metadata'][u'ra_pointing'] #in sexidecimal
    dec = full_meta_data[u'metadata'][u'dec_pointing']
    dura = full_meta_data[u'stoptime'] - full_meta_data[u'starttime'] #gps time
    xdelays = full_meta_data[u'rfstreams'][u"0"][u'xdelays']
    ydelays = full_meta_data[u'rfstreams'][u"0"][u'ydelays']
    minfreq = float(min(full_meta_data[u'rfstreams'][u"0"][u'frequencies']))
    maxfreq = float(max(full_meta_data[u'rfstreams'][u"0"][u'frequencies']))
    channels = full_meta_data[u'rfstreams'][u"0"][u'frequencies']
    centrefreq = 1.28 * (minfreq + (maxfreq-minfreq)/2)

    if return_all:
        return [obsid, ra, dec, dura, [xdelays, ydelays], centrefreq, channels], full_meta_data
    else:
        return [obsid, ra, dec, dura, [xdelays, ydelays], centrefreq, channels]


def getmeta(servicetype='metadata', service='obs', params=None):
    """
    Function to call a JSON web service and return a dictionary:
    Given a JSON web service ('obs', find, or 'con') and a set of parameters as
    a Python dictionary, return a Python dictionary xcontaining the result.
    Taken verbatim from http://mwa-lfd.haystack.mit.edu/twiki/bin/view/Main/MetaDataWeb
    """
    import urllib.request
    import json

    # Append the service name to this base URL, eg 'con', 'obs', etc.
    BASEURL = 'http://ws.mwatelescope.org/'


    if params:
        # Turn the dictionary into a string with encoded 'name=value' pairs
        data = urllib.parse.urlencode(params)
    else:
        data = ''

    try:
        result = json.load(urllib.request.urlopen(BASEURL + servicetype + '/' + service + '?' + data))
    except urllib.error.HTTPError as err:
        logger.error("HTTP error from server: code=%d, response:\n %s" % (err.code, err.read()))
        return
    except urllib.error.URLError as err:
        logger.error("URL or network error: %s" % err.reason)
        return

    return result


def get_files(obsid, files_meta_data=None):
    """
    Queries the metadata to find all the file names

    Parameters:
    -----------
    obsid: str
        The ID (gps time) of the observation you are querying
    files_meta_data: dict
        The output of the getmeta function with the data_files service.
        This is an optional input that can be used if you just want to
        extract the relevant info and save a metadata call

    Output:
    -------
    files: list
        A list of all the file names
    """
    if files_meta_data is None:
        files_meta_data = getmeta(servicetype='metadata', service='data_files', params={'obs_id':str(obsid)})

    return list(files_meta_data.keys())


def calc_ta_fwhm(freq, array_phase='P2C'):
    """
    Calculates the approximate FWHM of the tied array beam in degrees.

    Parameters:
    -----------
    freq: float
        Frequency in MHz
    array_phase: string
        OPTIONAL - The different array phase (from P1, P2C, P2E) to work out the maximum baseline length. Default = 'P2C'
    Returns:
    --------
    fwhm: float
        FWHM in degrees
    """
    from scipy.constants import c
    from math import degrees

    # Work out baseline in meters
    if array_phase == 'P1':
        # True max_baseline is 2800 but due to the minimal amount of long baselines
        # the following is more realisitic
        max_baseline = 2200.
    if array_phase == 'P2C':
        # True max_baseline is 700.
        max_baseline = 360.
    elif array_phase == 'P2E':
        max_baseline = 5300.

    wavelength = c / (freq * 1e6)
    fwhm = degrees(wavelength / max_baseline)

    return fwhm


def get_channels(obsid, channels=None):
    """
    Gets the channels ids from the observation's metadata. If channels is not None assumes the
    channels have already been aquired so it doesn't do an unnecessary database call.
    """
    if channels is None:
        print("Obtaining frequency channel data from http://mwa-metadata01.pawsey.org.au/metadata/"
              "for OBS ID: {}".format(obsid))
        full_meta_data = getmeta(service='obs', params={'obs_id':obsid})
        channels = full_meta_data[u'rfstreams'][u"0"][u'frequencies']
    return channels


def obs_max_min(obsid, files_meta_data=None):
    """
    Small function to query the database and return the times of the first and last file

    Parameters:
    -----------
    obsid: str
        The ID (gps time) of the observation you are querying
    files_meta_data: dict
        The output of the getmeta function with the data_files service.
        This is an optional input that can be used if you just want to
        extract the relevant info and save a metadata call

    Output:
    -------
    obs_start_end: list
        [obs_start, obs_end]
    """
    if files_meta_data is None:
        files_meta_data = getmeta(servicetype='metadata', service='data_files', params={'obs_id':str(obsid)})

    # Make a list of gps times excluding non-numbers from list
    times = [f[11:21] for f in get_files(obsid, files_meta_data=files_meta_data) if is_number(f[11:21])]
    obs_start = int(min(times))
    obs_end = int(max(times))
    return obs_start, obs_end


def write_obs_info(obsid):
    """
    Writes obs info to a file in the current direcory
    """
    data_dict = getmeta(service='obs', params={'obs_id':str(obsid), 'nocache':1})
    filename = str(obsid)+"_info.txt"
    logger.info("Writing to file: {0}".format(filename))

    channels = data_dict["rfstreams"]["0"]["frequencies"]
    centre_freq = ( min(channels) + max(channels) ) / 2. * 1.28
    array_phase = get_obs_array_phase(obsid)
    start, stop = obs_max_min(obsid)

    f = open(filename, "w+")
    f.write("-------------------------    Obs Info    --------------------------\n")
    f.write("Obs Name:           {}\n".format(data_dict["obsname"]))
    f.write("Creator:            {}\n".format(data_dict["rfstreams"]["0"]["creator"]))
    f.write("Array phase:        {}\n".format(array_phase))
    if array_phase != 'OTH':
        f.write("~FWHM (arcminute)   {:4.2f}\n".format(calc_ta_fwhm(centre_freq,
                                                       array_phase=array_phase)*60.))
    f.write("Start time:         {}\n".format(start))
    f.write("Stop time:          {}\n".format(stop))
    f.write("Duration (s):       {}\n".format(stop-start))
    f.write("RA Pointing (deg):  {}\n".format(data_dict["metadata"]["ra_pointing"]))
    f.write("DEC Pointing (deg): {}\n".format(data_dict["metadata"]["dec_pointing"]))
    f.write("Channels:           {}\n".format(channels))
    f.write("Centrefreq (MHz):   {}\n".format(centre_freq))
    f.close()


def get_best_cal_obs(obsid):
    """
    For the input MWA observation ID find all calibration observations within 2 days
    that have the same observing channels and list them from closest in time to furthest.

    Parameters
    ----------
    obsid: int
        The MWA observation ID (gps time)

    Returns
    -------
    cal_ids: list of lists
        All calibration observations within 2 days that have the same observing channels and
        list them from closest in time to furthest
        [[obsid, mins_away, cal_target]]
    """
    from operator import itemgetter

    obs_meta = getmeta(params={'obs_id':str(obsid)})
    channels = obs_meta[u'rfstreams'][u"0"][u'frequencies']
    cenchan = channels[12]
    if channels[-1] - channels[0] == 23:
        contig = 1
    else:
        contig = 0
    two_days_secs = 2*24*60*60

    all_cals = find_obsids_meta_pages(params={'calibration':1,
                                              'mintime': obsid-two_days_secs,
                                              'maxtime': obsid+two_days_secs,
                                              'cenchan': cenchan,
                                              'contigfreq': contig,
                                              'dataquality': 126})

    cal_info = []
    for cal in all_cals:
        #get the cal metadata
        cal_meta = getmeta(params={'obs_id':str(cal), 'nocache':1})

        #check there are a factor of 24 files (no gpu boxes are down)
        gpubox_files = []
        cal_files_meta = getmeta(service='data_files', params={'obs_id':obsid})
        for f in cal_files_meta.keys():
            if 'gpubox' in f:
                gpubox_files.append(f)
        if len(gpubox_files)%24 != 0 :
            continue

        #calculate the time away from the obs and append it to the list
        cal_info.append([cal, abs(obsid-cal)/60., cal_meta['obsname']])

    #sort by time
    cal_info = sorted(cal_info, key=itemgetter(1))

    return cal_info