import os
import sys
import logging
from astropy.time import Time
from time import strptime, strftime
from astropy.utils import iers

logger = logging.getLogger(__name__)

def sfreq(freqs):

    if len(freqs) != 24:
        print("There are not 24 coarse chans defined for this obs. Got: %s" % freqs)
        return

    #freqs.sort()   # It should already be sorted, but just in case...[SET] Commenting this out because sort() is ironically putting 2-digit channels out of order
    lowchans = [f for f in freqs if f <= 128]
    highchans = [f for f in freqs if f > 128]

    highchans.reverse()
    freqs = lowchans + highchans

    return freqs


def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def gps_to_utc(gps):
    # GPS time as is done in timeconvert.py
    iers.IERS_A_URL = 'https://datacenter.iers.org/data/9/finals2000A.all'
    logger.info(iers.IERS_A_URL)
    utctime = Time(gps, format='gps', scale='utc').fits
    # remove (UTC) that some astropy versions leave on the end
    if utctime.endswith('(UTC)'):
        utctime = strptime(utctime, '%Y-%m-%dT%H:%M:%S.000(UTC)')
        utctime = strftime('%Y-%m-%dT%H:%M:%S', utctime)
    else:
        utctime = strptime(utctime, '%Y-%m-%dT%H:%M:%S.000')
        utctime = strftime('%Y-%m-%dT%H:%M:%S', utctime)
    return utctime


def mdir(path, description, gid=34858):
    """
    Simple function to create directories with the correct group permissions

    The default group ID is 'mwavcs' which is 30832 in numerical.
    Here, we try and make sure all directories created by process_vcs
    end up belonging to the user and the group 'mwavcs'.
    We also try to give rwx (read, write, execute) permissions and
    set the sticky bit for both user and group.
    """
    try:
        # TODO: this doesn't carry over permissions correctly to "combined" for some reason...
        os.makedirs(path)
        # we leave the uid unchanged but change gid to mwavcs
        os.chown(path, -1, gid)
        os.chmod(path, 0o771)
        os.system("chmod -R g+s {0}".format(path))
    except:
        if (os.path.exists(path)):
            logger.info("{0} Directory Already Exists\n".format(description))
        else:
            logger.error("{0} Could not make new directories\n".format(description))
            sys.exit(0)


