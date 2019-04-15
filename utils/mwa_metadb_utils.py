#!/usr/bin/env python

# AOCal variables
import numpy as np
import struct, os, logging
from collections import namedtuple
HEADER_FORMAT = "8s6I2d"
HEADER_SIZE = struct.calcsize(HEADER_FORMAT)
HEADER_INTRO = "MWAOCAL\0"

Header = namedtuple("header", "intro fileType structureType intervalCount antennaCount channelCount polarizationCount timeStart timeEnd")
Header.__new__.__defaults__ = (HEADER_INTRO, 0, 0, 0, 0, 0, 0, 0.0, 0.0)

"""
Collection of database related utilities that are used throughout the VCS processing pipeline
"""
class AOCal(np.ndarray):
    """
    AOCAl stored as a numpy array (with start and stop time stored as floats)

    Array is of dtype complex128 with the following dimensions:

    - calibration interval
    - antenna
    - channel
    - polarisation (order XX, XY, YX, YY)

    The following attributes are made available for convenience, however they
    are not stored explicitly, just read from the array shape.

    aocal.n_int
    aocal.n_ant
    aocal.n_chan
    aocal.n_pol
    """
    
    def __new__(cls, input_array, time_start=0.0, time_end=0.0):
        """
        See http://docs.scipy.org/doc/numpy-1.10.1/user/basics.subclassing.html
        """
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.time_start = float(time_start)
        obj.time_end = float(time_end)
        # Finally, we must return the newly created object:
        return obj
    
    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.time_start = getattr(obj, 'time_start', None)
        self.time_end = getattr(obj, 'time_end', None)

    def __getattr__(self, name):
        if name == 'n_int':
            return self.shape[0]
        elif name == 'n_ant':
            return self.shape[1]
        elif name == 'n_chan':
            return self.shape[2]
        elif name == 'n_pol':
            return self.shape[3]
        elif name == 'time_start':
            # required to avoid infinite recursion
            return object.__getattribute__(self, time_start)
        elif name == 'time_end':
            # required to avoid infinite recursion
            return object.__getattribute__(self, time_end)
        else:
            raise AttributeError, "AOCal has no Attribute %s. Dimensions can be accessed via n_int, n_ant, n_chan, n_pol" % name

    def strip_edge(self, n_chan):
        """
        return a copy of the array with edge channels removed

        useful for printing without nans but don't write out as calibration solution!
        """
        return self[:, :, n_chan:-n_chan, :]

    def tofile(self, cal_filename):
        if not (np.iscomplexobj(self) and self.itemsize == 16 and len(self.shape) == 4):
            raise TypeError, "array must have 4 dimensions and be of type complex128"
        header = Header(intervalCount=self.shape[0], antennaCount = self.shape[1], channelCount = self.shape[2], polarizationCount = self.shape[3], timeStart = self.time_start, timeEnd = self.time_end)
        with open(cal_filename, "wb") as cal_file:
            header_string = struct.pack(HEADER_FORMAT, *header)
            cal_file.write(header_string)
            logging.debug("header written")
            cal_file.seek(HEADER_SIZE, os.SEEK_SET) # skip header. os.SEEK_SET means seek relative to start of file
            np.ndarray.tofile(self, cal_file)
            logging.debug("binary file written")

    def fit(self, pols=(0, 3), mode='model', amp_order=5):
        if not (np.iscomplexobj(self) and self.itemsize == 16 and len(self.shape) == 4):
            raise TypeError, "array must have 4 dimensions and be of type complex128"
        fit_array = np.zeros(self.shape, dtype=np.complex128)
        for interval in xrange(self.shape[0]):
            for antenna in xrange(self.shape[1]):
                logging.debug("fitting antenna %d" % antenna)
                for pol in pols:
                    v = self[interval, antenna, :, pol]
                    if sum(~np.isnan(self[interval, antenna, :, pol])) > 0:
                        self[interval, antenna, :, pol] = fit_complex_gains(self[interval, antenna, :, pol])


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

def get_common_obs_metadata(obs, return_all = False):
    """
    Gets needed comon meta data from http://mwa-metadata01.pawsey.org.au/metadata/
    """
    print "Obtaining metadata from http://mwa-metadata01.pawsey.org.au/metadata/ for OBS ID: " + str(obs)
    #for line in txtfile:
    beam_meta_data = getmeta(service='obs', params={'obs_id':obs})
    #obn = beam_meta_data[u'obsname']
    ra = beam_meta_data[u'metadata'][u'ra_pointing'] #in sexidecimal
    dec = beam_meta_data[u'metadata'][u'dec_pointing']
    dura = beam_meta_data[u'stoptime'] - beam_meta_data[u'starttime'] #gps time
    xdelays = beam_meta_data[u'rfstreams'][u"0"][u'xdelays']
    minfreq = float(min(beam_meta_data[u'rfstreams'][u"0"][u'frequencies']))
    maxfreq = float(max(beam_meta_data[u'rfstreams'][u"0"][u'frequencies']))
    channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
    centrefreq = 1.28 * (minfreq + (maxfreq-minfreq)/2)

    if return_all:
        return [obs,ra,dec,dura,xdelays,centrefreq,channels],beam_meta_data
    else:
        return [obs,ra,dec,dura,xdelays,centrefreq,channels]

def getmeta(service='obs', params=None):
    """
    Function to call a JSON web service and return a dictionary:
    Given a JSON web service ('obs', find, or 'con') and a set of parameters as
    a Python dictionary, return a Python dictionary xcontaining the result.
    Taken verbatim from http://mwa-lfd.haystack.mit.edu/twiki/bin/view/Main/MetaDataWeb
    """
    import urllib
    import urllib2
    import json

    # Append the service name to this base URL, eg 'con', 'obs', etc.
    BASEURL = 'http://mwa-metadata01.pawsey.org.au/metadata/'


    if params:
        # Turn the dictionary into a string with encoded 'name=value' pairs
        data = urllib.urlencode(params)
    else:
        data = ''

    if service.strip().lower() in ['obs', 'find', 'con']:
        service = service.strip().lower()
    else:
        print "invalid service name: %s" % service
        return

    try:
        result = json.load(urllib2.urlopen(BASEURL + service + '?' + data + "&nocache"))
    except urllib2.HTTPError as error:
        print "HTTP error from server: code=%d, response:\n %s" % (error.code, error.read())
        return
    except urllib2.URLError as error:
        print "URL or network error: %s" % error.reason
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


