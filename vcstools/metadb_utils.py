from vcstools.general_utils import is_number
import logging
import os
import subprocess
from time import sleep
import numpy as np

logger = logging.getLogger(__name__)


def ensure_metafits(data_dir, obsid, metafits_file):
    """Ensure that the metafits file is in the directory, if not download it.

    Parameters
    ----------
    data_dir : `str`
        The directory to check the metafits file is in.
    obsid : `int`
        The MWA Observation ID.
    metafits_file : `str`
        The metafits file name.
    """
    # TODO: To get the actual ppds file should do this with obsdownload -o <obsID> -m

    if not os.path.exists(os.path.join(data_dir, metafits_file)):
        logger.warning("{0} does not exists".format(metafits_file))
        logger.warning("Will download it from the archive. This can take a "
                      "while so please do not ctrl-C.")
        logger.warning("At the moment, even through the downloaded file is "
                       "labelled as a ppd file this is not true.")
        logger.warning("This is hopefully a temporary measure.")

        get_metafits = "wget http://ws.mwatelescope.org/metadata/fits?obs_id={0} -O {1}".format(obsid, metafits_file)
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
    if os.path.exists(data_dir) and not os.path.exists("{}/{}".format(data_dir, metafits_file)):
        logger.info("Copying {0} to {1}".format(metafits_file, data_dir))
        from shutil import copy2
        copy2("{0}".format(metafits_file), "{0}".format(data_dir))


def singles_source_search(ra, dec=None, box_size=50., params=None):
    """Used to find all obsids within a box around the source to make searching through obs_ids more efficient.

    Parameters
    ----------
    ra : `float`
        Right Acension of the source in degrees
    dec : `float`
        Declination of the source in degrees. By default will use the enitre declination range to account for grating lobes
    box_size : `float`
        Radius of the search box. Default: 45
    params : `dict`
        The dictionary of constraints used to search for suitable observations as explained here:
        https://wiki.mwatelescope.org/display/MP/Web+Services#WebServices-Findobservations
        |br| Default: {'mode':'VOLTAGE_START'}

    Returns
    -------
    obsid_list : `list`
        List of the MWA observation IDs.
    """
    if params is None:
        # Load default params
        params={'mode':'VOLTAGE_START'}
    logger.debug("params: {}".format(params))

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
        params.update({'minra':0.,      'maxra':360.,
                       'mindec':dec_bot,'maxdec':dec_top})
    else:
        ra_low = ra - 30. - box_size #30 is the how far an obs would drift in 2 hours(used as a max)
        ra_high = ra + box_size
        if ra_low < 0.:
            ra_low = 360 + ra_low
        if ra_high > 360:
            ra_high = ra_high - 360
        params.update({'minra':ra_low,  'maxra':ra_high,
                       'mindec':dec_bot,'maxdec':dec_top})
    logger.debug("params: {}".format(params))
    obsid_list = find_obsids_meta_pages(params=params)
    return obsid_list


def find_obsids_meta_pages(params=None, pagesize=50):
    """Loops over pages for each page for MWA metadata calls

    Parameters
    ----------
    params : `dict`
        The dictionary of constraints used to search for suitable observations as explained here:
        https://wiki.mwatelescope.org/display/MP/Web+Services#WebServices-Findobservations
        |br| Default: {'mode':'VOLTAGE_START'}

    Returns
    -------
    obsid_list : `list`
        List of the MWA observation IDs.
    """
    if params is None:
        params = {'mode':'VOLTAGE_START'}
    params["pagesize"] = pagesize
    obsid_list = []
    temp = []
    page = 1
    #need to ask for a page of results at a time
    while len(temp) == pagesize or page == 1:
        params['page'] = page
        logger.debug("Page: {0}   params: {1}".format(page, params))
        temp = getmeta(service='find', params=params, retry_http_error=True)
        if temp is not None:
            # if there are non obs in the field (which is rare) None is returned
            for row in temp:
                obsid_list.append(row[0])
        else:
            temp = []
        page += 1

    return obsid_list


def get_obs_array_phase(obsid):
    """For the input obsid will work out the observations array phase in the form
    of P1 for phase 1, P2C for phase 2 compact or P2E for phase to extended array
    and OTH for other.

    Parameters
    ----------
    obsid : `int`
        The MWA Observation ID.
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
    """Calculate the altitude, azumith and zenith for an obsid

    Parameters
    ----------
    obsid : `int`
        The MWA Observation ID.
    ra : `str`
        The Right Acension in the format "HH:MM:SS.SS".
    dec : `str`
        The Declination in the format "DD:MM:SS.SS".
    degrees : `boolean`
        If `True` the ra and dec is given in degrees. |br| Default: `False`.

    Returns
    -------
    Alt : `float`
        The altitude angle in degrees.
    Az : `float`
        The azimuth angle in degrees.
    Za : `float`
        The zenith angle in degrees.
    """
    from astropy.utils import iers
    iers.IERS_A_URL = 'https://datacenter.iers.org/data/9/finals2000A.all'
    #logger.debug(iers.IERS_A_URL)

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


def get_common_obs_metadata(obsid, return_all=False, full_metadata=None):
    """Gets needed common meta data from http://ws.mwatelescope.org/metadata/

    Parameters
    ----------
    obsid : `int`
        The MWA observation ID.
    return_all : `boolean`, optional
        If `True` will also return the full meta data dictionary. |br| Default: `False`.
    full_metadata : `dict`, optional
        The dictionary of metadata generated from :py:meth:`vcstools.metadb_utils.getmeta`

    Returns
    -------
    common_metadata : `list`
        [obsid, ra, dec, dura, [xdelays, ydelays], centrefreq, channels]

        obsid : `int`
            The MWA observation ID.
        ra : `str`
            The Right Acension in the format "HH:MM:SS.SS".
        dec : `str`
            The Declination in the format "DD:MM:SS.SS".
        dura : `int`
            The duration of the observation in seconds.
        xdelays, ydelays : `list`
            The analogue delays for each antena od the tile for both the x and y polarisations.
        centrefreq : `float`
            The centre observing frequency in MHz.
        channels : `list`
            The list of observing frequency coarse channel IDs.
    """
    if full_metadata is None:
        logger.info("Obtaining metadata from http://ws.mwatelescope.org/metadata/ for OBS ID: " + str(obsid))
        full_metadata = getmeta(service='obs', params={'obs_id':obsid})
    ra = full_metadata[u'metadata'][u'ra_pointing'] #in sexidecimal
    dec = full_metadata[u'metadata'][u'dec_pointing']
    dura = full_metadata[u'stoptime'] - full_metadata[u'starttime'] #gps time
    xdelays = full_metadata[u'rfstreams'][u"0"][u'xdelays']
    ydelays = full_metadata[u'rfstreams'][u"0"][u'ydelays']
    minfreq = float(min(full_metadata[u'rfstreams'][u"0"][u'frequencies']))
    maxfreq = float(max(full_metadata[u'rfstreams'][u"0"][u'frequencies']))
    channels = full_metadata[u'rfstreams'][u"0"][u'frequencies']
    centrefreq = 1.28 * (minfreq + (maxfreq-minfreq)/2)

    if return_all:
        return [obsid, ra, dec, dura, [xdelays, ydelays], centrefreq, channels], full_metadata
    else:
        return [obsid, ra, dec, dura, [xdelays, ydelays], centrefreq, channels]


def getmeta(servicetype='metadata', service='obs', params=None, retries=3, retry_http_error=False):
    """Function to call a JSON web service to perform an MWA metadata call.
    Taken verbatim from http://mwa-lfd.haystack.mit.edu/twiki/bin/view/Main/MetaDataWeb.

    Parameters
    ----------
    servicetype : `str`
        Either the 'observation' which makes human readable html pages or 'metadata' which returns data.
        |br| Default: metadata.
    service : `str`
        The meta data service from (Defaul: obs):

        obs:
            Returns details about a single observation as explained in
            https://wiki.mwatelescope.org/display/MP/Web+Services#WebServices-Getobservation/scheduledataforanobservation

        find:
            Search the database for observations that satisfy given criteria as explained in
            https://wiki.mwatelescope.org/display/MP/Web+Services#WebServices-Findobservations

        con:
            Finds the configuration information for an observation as explained in
            https://wiki.mwatelescope.org/display/MP/Web+Services#WebServices-Gettelescopeconfigurationforagivenobservation

    params : `dict`
        A dictionary of the options to use in the metadata call which is dependent on the `service`. Examples can be found
        https://wiki.mwatelescope.org/display/MP/Web+Services#WebServices-Usingthegetmeta()function
    retries : `int`, optional
        The number of times to retry timeout errors.

    Returns
    -------
    result : `dict`
        The result for that service.
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

    # Try several times (3 by default)
    wait_time = 30
    result = None
    for x in range(0, retries):
        err = False
        try:
            result = json.load(urllib.request.urlopen(BASEURL + servicetype + '/' + service + '?' + data))
        except urllib.error.HTTPError as err:
            logger.error("HTTP error from server: code=%d, response: %s" % (err.code, err.read()))
            if retry_http_error:
                logger.error("Waiting {} seconds and trying again".format(wait_time))
                sleep(wait_time)
                pass
            else:
                raise err
                break
        except urllib.error.URLError as err:
            logger.error("URL or network error: %s" % err.reason)
            logger.error("Waiting {} seconds and trying again".format(wait_time))
            sleep(wait_time)
            pass
        else:
            break
    else:
        logger.error("Tried {} times. Exiting.".format(retries))

    return result


def get_ambient_temperature(obsid, full_metadata=None):
    """Queries the metadata to find the ambient temperature of the MWA tiles in K.
    If none recorded then assume 22.4 C.

    Parameters
    ----------
    obsid : `int`
        The MWA observation ID.
    full_metadata : `dict`, optional
        The dictionary of metadata generated from :py:meth:`vcstools.metadb_utils.getmeta`.

    Returns
    -------
    t_0 : `float`
        The ambient temperature of the MWA tiles in K.
    """
    if full_metadata is None:
        logger.info("Obtaining metadata from http://ws.mwatelescope.org/metadata/ for OBS ID: " + str(obsid))
        full_metadata = getmeta(service='obs', params={'obs_id':obsid})
    ambient_temp_dict = full_metadata[u'bftemps']
    logger.debug("ambient_temp_dict: {}".format(ambient_temp_dict))
    temperature_list_raw = []
    for tile_key in ambient_temp_dict:
        temperature_list_raw.append(ambient_temp_dict[tile_key][-1])
    # Remove Nones
    temperature_list = [i for i in temperature_list_raw if i]
    logger.debug("temperature_list: {}".format(temperature_list))
    if len(temperature_list) == 0:
        logger.warning("No ambient temperature recorded so assuming the yearly average of 22.4 C.")
        temperature_list = [22.4]
    return np.mean(temperature_list) + 273.15


def get_files(obsid, files_meta_data=None):
    """Queries the metadata to find all the file names.

    Parameters
    ----------
    obsid : `int`
        The MWA Observation ID.
    full_metadata : `dict`, optional
        The dictionary of metadata generated from :py:meth:`vcstools.metadb_utils.getmeta`

    Returns
    -------
    files : `list`
        A list of all the file names.
    """
    if files_meta_data is None:
        files_meta_data = getmeta(servicetype='metadata', service='data_files', params={'obs_id':str(obsid)})

    return list(files_meta_data.keys())


def calc_ta_fwhm(freq, array_phase='P2C'):
    """Calculates the approximate FWHM of the tied-array beam in degrees.

    Parameters
    ----------
    freq : `float`
        Frequency in MHz.
    array_phase : `str`, optional
        The different array phase (from P1, P2C, P2E) to work out the maximum baseline length. |br| Default = 'P2C'.

    Returns
    -------
    fwhm : `float`
        FWHM in degrees.
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
    """Gets the channels ids from the observation's metadata. If channels is not None assumes the
    channels have already been aquired so it doesn't do an unnecessary database call.

    Parameters
    ----------
    obsid : `int`
        The MWA Observation ID.
    channels : `list`, optional
        The list of the coarse channel frequency IDs.

    Returns
    -------
    channels : `list`
        The list of the coarse channel frequency IDs.
    """
    if channels is None:
        print("Obtaining frequency channel data from http://mwa-metadata01.pawsey.org.au/metadata/"
              "for OBS ID: {}".format(obsid))
        full_metadata = getmeta(service='obs', params={'obs_id':obsid})
        channels = full_metadata[u'rfstreams'][u"0"][u'frequencies']
    return channels


def obs_max_min(obsid, files_meta_data=None):
    """Small function to query the database and return the times of the first and last file.

    Parameters
    ----------
    obsid : `int`
        The MWA Observation ID.
    full_metadata : `dict`, optional
        The dictionary of metadata generated from :py:meth:`vcstools.metadb_utils.getmeta`

    Returns
    -------
    obs_start : `int`
        Beginning of the observation GPS time.
    obs_end : `int`
        End of the observation GPS time.
    """
    if files_meta_data is None:
        files_meta_data = getmeta(servicetype='metadata', service='data_files', params={'obs_id':str(obsid)})

    # Make a list of gps times excluding non-numbers from list
    times = [f[11:21] for f in get_files(obsid, files_meta_data=files_meta_data) if is_number(f[11:21])]
    if times:
        obs_start = int(min(times))
        obs_end   = int(max(times))
    else:
        obs_start = None
        obs_end   = None
    return obs_start, obs_end


def files_available(obsid, files_meta_data=None):
    """Query the database and return a list of all files available (remote archived and not deleted) and a list of all files

    Parameters
    ----------
    obsid : `int`
        The MWA Observation ID.
    full_metadata : `dict`, optional
        The dictionary of metadata generated from :py:meth:`vcstools.metadb_utils.getmeta`

    Returns
    -------
    available_files : `list`
        All files that have been archived and ready for download.
    all_files : `list`
        All files that are expected, once all files have been transfered/archived.
    """
    if files_meta_data is None:
        files_meta_data = getmeta(servicetype='metadata', service='data_files', params={'obs_id':str(obsid)})

    # Loop over all the files and check if they're archived and not deleted
    available_files = []
    all_files = []
    for file_name in files_meta_data.keys():
        deleted = files_meta_data[file_name]['deleted']
        remote_archived = files_meta_data[file_name]["remote_archived"]
        if remote_archived and not deleted:
            available_files.append(file_name)
        all_files.append(file_name)

    return available_files, all_files


def get_best_cal_obs(obsid):
    """For the input MWA observation ID find all calibration observations within 2 days
    that have the same observing channels and list them from closest in time to furthest.

    Parameters
    ----------
    obsid : `int`
        The MWA Observation ID.

    Returns
    -------
    cal_ids : `list`
        [[obsid, mins_away, cal_target]]

        obsid : `int`
            The MWA calibration observation ID.
        mins_away : `float`
            The mins away the calibration observation is away from the target observation.
        cal_target : `str`
            The calibration target name.
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


def combined_deleted_check(obsid, begin=None, end=None):
    """Check if the combined files are deleted (or do not exist).

    Parameters
    ----------
    obsid : `int`
        The MWA Observation ID.
    beg : `int`, optional
        Beginning of the observation GPS time to check.
    end : `int`, optional
        End of the observation GPS time to check.

    Returns
    -------
    comb_del_check : `bool`
        True if all combined files are deleted or if they do not exist.
    """
    # Work out the data_files metadata call parameters
    params = {'obs_id':obsid, 'nocache':1}
    if begin is not None:
        params['mintime'] = begin
    if end is not None:
        params['maxtime'] = end

    files_meta = getmeta(service='data_files', params=params)

    # Loop over files to search for a non deleted file
    comb_del_check = True
    for file in files_meta.keys():
        if 'combined' in file:
            deleted =         files_meta[file]['deleted']
            remote_archived = files_meta[file]["remote_archived"]
            if remote_archived and not deleted:
                # A combined file exists
                comb_del_check = False
                break

    return comb_del_check