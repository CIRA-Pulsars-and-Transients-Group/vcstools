import os
import sys
import logging
from astropy.time import Time
from time import strptime, strftime
from astropy.utils import iers
import numpy as np

logger = logging.getLogger(__name__)

def sfreq(freqs):
    """Sore te coarse frequency channel IDs into the strange MWA format that reverses the order of channels above 128.

    Parameters
    ----------
    freqs : `list`
        List of coarse frequency channel IDs.

    Returns
    -------
    freqs : `list`
        List of coarse frequency channel IDs in the new order.
    """
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
    """Simple is check to see if a string can be converted to an interger.
    
    Parameters
    ----------
    s : `str`
        String to check.

    Returns
    -------
    result : `boolean`
        Boolean if it can be converted to an int.
    """
    try:
        int(s)
        return True
    except ValueError:
        return False


def gps_to_utc(gps):
    # GPS time as is done in timeconvert.py
    iers.IERS_A_URL = 'https://datacenter.iers.org/data/9/finals2000A.all'
    #logger.debug(iers.IERS_A_URL)
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
    """Simple function to create directories with the correct group permissions (771).
    
    Parameters
    ----------
    path : `str`
        The path of the directory we want to create.
    description : `str`
        The description of the directory to be printed to logger.
    gid : `int`, optional
        The group ID to apply to the directory. |br| Default: 34858 which the mwavcs.
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



def create_link(data_dir, target_dir, product_dir, link):
    """Creates a symbolic link product_dir/link that points to data_dir/target_dir.

    Parameters
    ----------
    data_dir : `str`
        The absolute path to the base directory of the true location of the files.
        For our uses this is often a scratch partition like /astro on Galaxy
    target_dir : `str`
        The folder you would like to be linked to
    product_dir : `str`
        The absolute path of the link you would like to create
    link : `str`
        The name of the link you would like to create. Often the same as target_dir
    """
    data_dir = os.path.abspath(data_dir)
    product_dir = os.path.abspath(product_dir)
    if data_dir == product_dir:
        # base directories are the same so doing nothing
        return

    # add product_dir and data_dir to link and target_dir respectively
    link = link.replace(product_dir, '') # just in case...
    link = link.replace('/', '')
    link = os.path.join(product_dir, link)
    target_dir = target_dir.replace(data_dir,'')
    if target_dir.startswith("/"):
        target_dir = target_dir[1:]
    target_dir = os.path.join(data_dir, target_dir)

    # check if link exists and whether it already points to where we'd like it to
    if os.path.exists(link):
        if os.path.islink(link):
            if os.readlink(link) == target_dir:
                return
            else:
                logger.warning("The link {0} already exists but points at {1} while you "
                               "asked it to point at {2}. Deleting the link and creating"
                               "a new one".format(link, os.readlink(link), target_dir))
                os.unlink(link)
                os.symlink(target_dir, link)
        else:
            logger.error("{0} is an existing directory and cannot be turned into a link. Aborting...".format(link))
            sys.exit(0)
    else:
        logger.info("Trying to link {0} against {1}".format(link, target_dir))
        os.symlink(target_dir, link)

def setup_logger(logger, log_level=logging.INFO):
    """Setup the logger with the format we prefer and apply it to all import vcstools modules.

    Parameters
    ----------
    logger : logger object
        The logger object to modify.
    log_level : logging class
        The logging level to apply. |br| Default: logging.INFO
    
    Returns
    -------
    logger : logger object
        The modified logger object.
    """
    # set up the logger for stand-alone execution
    formatter = logging.Formatter('%(asctime)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    # Set up local logger
    logger.setLevel(log_level)
    logger.addHandler(ch)
    logger.propagate = False
    # Loop over imported vcstools modules and set up their loggers
    for imported_module in sys.modules.keys():
        if imported_module.startswith('vcstools'):
            logging.getLogger(imported_module).setLevel(log_level)
            logging.getLogger(imported_module).addHandler(ch)
            logging.getLogger(imported_module).propagate = False
    return logger

def split_remove_remainder(array, nchunks):
    """Split an array into nchunks and remove any remaining elements.

    Parameters
    ----------
    array : np.array, (N)
        A single dimension array.
    nchunks :`int`
        The number of sub arrays to split array into.

    Returns
    -------
    array_chunks : `list`
        A list containing `nhcunks` arrays of equal size.
    """
    # Work out size of the sub arrays
    size = len(array) // nchunks
    array_chunks = np.split(array, np.arange(size,len(array),size))
    if len(array) % nchunks != 0:
        # Remove remainder
        logger.debug("Removing {} size array from end of array_chunks.".format(len(array_chunks[-1])))
        array_chunks = array_chunks[:-1]
    return array_chunks