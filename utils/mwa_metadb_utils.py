#!/usr/bin/env python3


def mwa_alt_az_za(obsid, ra=None, dec=None, degrees=False):
    """
    Calculate the altitude, azumith and zenith for an obsid

    Args:
        obsid  : The MWA observation id (GPS time)
        ra     : The right acension in HH:MM:SS
        dec    : The declintation in HH:MM:SS
        degrees: If true the ra and dec is given in degrees (Default:False)
    """
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
    earth_location = EarthLocation.of_site('Murchison Widefield Array') 
    altaz = sky_posn.transform_to(AltAz(obstime=obstime, location=earth_location)) 
    Alt = altaz.alt.deg
    Az  = altaz.az.deg 
    Za  = 90. - Alt
    return Alt, Az, Za


import logging

logger = logging.getLogger(__name__)

def get_common_obs_metadata(obs, return_all = False):
    """
    Gets needed comon meta data from http://mwa-metadata01.pawsey.org.au/metadata/
    """
    logger.info("Obtaining metadata from http://mwa-metadata01.pawsey.org.au/metadata/ for OBS ID: " + str(obs))
    #for line in txtfile:
    beam_meta_data = getmeta(service='obs', params={'obs_id':obs})
    #obn = beam_meta_data[u'obsname']
    ra = beam_meta_data[u'metadata'][u'ra_pointing'] #in sexidecimal
    dec = beam_meta_data[u'metadata'][u'dec_pointing']
    dura = beam_meta_data[u'stoptime'] - beam_meta_data[u'starttime'] #gps time
    xdelays = beam_meta_data[u'rfstreams'][u"0"][u'xdelays']
    ydelays = beam_meta_data[u'rfstreams'][u"0"][u'ydelays']
    minfreq = float(min(beam_meta_data[u'rfstreams'][u"0"][u'frequencies']))
    maxfreq = float(max(beam_meta_data[u'rfstreams'][u"0"][u'frequencies']))
    channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
    centrefreq = 1.28 * (minfreq + (maxfreq-minfreq)/2)

    if return_all:
        return [obs, ra, dec, dura, [xdelays, ydelays], centrefreq, channels], beam_meta_data
    else:
        return [obs, ra, dec, dura, [xdelays, ydelays], centrefreq, channels]


def getmeta(service='obs', params=None):
    """
    Function to call a JSON web service and return a dictionary:
    Given a JSON web service ('obs', find, or 'con') and a set of parameters as
    a Python dictionary, return a Python dictionary xcontaining the result.
    Taken verbatim from http://mwa-lfd.haystack.mit.edu/twiki/bin/view/Main/MetaDataWeb
    """
    #import urllib
    import urllib.request
    import json

    # Append the service name to this base URL, eg 'con', 'obs', etc.
    BASEURL = 'http://mwa-metadata01.pawsey.org.au/metadata/'


    if params:
        # Turn the dictionary into a string with encoded 'name=value' pairs
        data = urllib.parse.urlencode(params)
    else:
        data = ''

    if service.strip().lower() in ['obs', 'find', 'con']:
        service = service.strip().lower()
    else:
        logger.error("invalid service name: %s" % service)
        return

    try:
        result = json.load(urllib.request.urlopen(BASEURL + service + '?' + data + "&nocache"))
    except urllib.error.HTTPError as err:
        logger.error("HTTP error from server: code=%d, response:\n %s" % (err.code, err.read()))
        return
    except urllib.error.URLError as err:
        logger.error("URL or network error: %s" % err.reason)
        return

    return result


def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def obs_max_min(obs_id):
    """
    Small function to query the database and return the times of the first and last file
    """
    obsinfo = getmeta(service='obs', params={'obs_id':obs_id})
    
    # Make a list of gps times excluding non-numbers from list
    times = [f[11:21] for f in obsinfo['files'].keys() if is_number(f[11:21])]
    obs_start = int(min(times))
    obs_end = int(max(times))
    return obs_start, obs_end


