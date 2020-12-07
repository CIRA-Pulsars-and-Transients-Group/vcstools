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