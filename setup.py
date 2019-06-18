#! /usr/bin/env python3
"""
Setup for vcstools
"""
import os
import sys
from setuptools import setup

def read(fname):
    """Read a file"""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# Since we mostly run this on supercomputers it probably isn't correct to 
# pip install all these modules
"""
reqs = ['python>=3.6.3',
        'argparse>=1.4.0',
        'numpy>=1.13.3',
        'matplotlib>=2.1.0',
        'astropy>=2.0.2',
        'distribute>=0.7.3',
        'setuptools>=38.2.1',
        'psycopg2>=2.7.3.2',
        'requests>=2.18.4']
"""

setup(name="vcstools",
      #version=get_version(),
      description="Scripts used to process the Murchison Widefield Array's Voltage Capture System data",
      url="https://github.com/ICRAR/vcstools",
      long_description=read('README.md'),
      #packages=['vcstools'],
      #install_requires=reqs,
      scripts=['scripts/checks.py', 'scripts/calibrate_vcs.py', 'scripts/create_psrfits.sh',
               'scripts/file_maxmin.py', 'scripts/find_pulsar_in_obs.py', 'scripts/obs_query.py',
               'scripts/process_vcs.py', 'scripts/recombine.py', 'scripts/rename_corr_output.py',
               'scripts/reorder_chans.py', 'scripts/rts2ao.py', 'scripts/untar.sh',
               'scripts/cleanup.py',
               'database/submit_to_database.py', 'database/database_vcs.py',
               'utils/zapchan.py', 'utils/calc_ephem.py', 'utils/check_disk_usage.sh',
               'utils/check_quota.sh', 'utils/mdir.py', 'utils/mwa_metadb_utils.py',
               'utils/job_submit.py', 'utils/plotFlatTileBeam.py', 'utils/plotPolarTileBeam.py',
               'utils/plotTiedArrayBeam.py', 'utils/get_MWAobsInfo.py', 'utils/aocal.py',
               'utils/config.py'],
      #data_files=[('AegeanTools', [os.path.join(data_dir, 'MOC.fits')]) ],
      setup_requires=['pytest-runner'],
      tests_require=['pytest']#, 'nose']
)
