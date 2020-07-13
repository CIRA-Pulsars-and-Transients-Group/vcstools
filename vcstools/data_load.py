"""
Loads all the data required by vcstools from the data directory.
"""

import os

datadir = os.path.join(os.path.dirname(__file__), 'data')

# Hard code the path of the MWA tile receiver temperature file
TRCVR_FILE = os.path.join(datadir, 'MWA_Trcvr_tile_56.csv')

# Hard code the paths of pulsar search files
KNOWN_RFRB_CSV = os.path.join(datadir, 'known_repeating_FRBs.csv')
POI_CSV = os.path.join(datadir, 'points_of_interest.csv')

# Hard code the path of the ATNF psrcat database file
ATNF_LOC = os.path.join(datadir, 'psrcat.db')
# Check if the file exists, if not download the latest zersion
if not os.path.exists(ATNF_LOC):
    # Importing download functions here to avoid unnessiary imports when the file is available
    import urllib.request
    import gzip
    import shutil
    import tarfile
    print("The ANTF psrcat database file does not exist. Downloading it from www.atnf.csiro.au")
    # Download the file
    psrcat_zip_dir = urllib.request.urlretrieve('https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz')[0]
    # Unzip it
    with gzip.open(psrcat_zip_dir,  'rb') as f_in:
        with open('psrcat_pkg.tar', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    # Untar the file we require
    psrcat_tar = tarfile.open(psrcat_zip_dir)
    # Do some python magic to no download the file within it's subdirectory from https://stackoverflow.com/questions/8405843/python-extract-using-tarfile-but-ignoring-directories
    member = psrcat_tar.getmember('psrcat_tar/psrcat.db')
    member.name = os.path.basename(member.name)
    psrcat_tar.extract(member, path=datadir)
    print("Download complete")


