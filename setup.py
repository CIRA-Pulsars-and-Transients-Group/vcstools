#! /usr/bin/env python3
"""
Setup for vcstools
"""
import os
import sys
from setuptools import setup
from subprocess import check_output

def read(fname):
    """Read a file"""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

#The following two functions were taken from the repo: https://github.com/pyfidelity/setuptools-git-version/blob/master/setuptools_git_version.py
def format_version(version, fmt='{tag}.{commitcount}_{gitsha}'):
    parts = version.split('-')
    if len(parts) == 1:
        return parts[0]
    assert len(parts) in (3, 4)
    dirty = len(parts) == 4
    tag, count, sha = parts[:3]
    if count == '0' and not dirty:
        return tag
    return fmt.format(tag=tag, commitcount=count, gitsha=sha.lstrip('g'))

def get_git_version():
    git_version = check_output('git describe --tags --long --dirty --always'.split()).decode('utf-8').strip()
    return format_version(version=git_version)

vcstools_version = get_git_version()

# Since we mostly run this on supercomputers it probably isn't correct to 
# pip install all these modules
reqs = ['python>=3.6.3',
        'argparse>=1.4.0',
        'numpy>=1.13.3',
        'matplotlib>=2.1.0',
        'astropy>=2.0.2']

#make a temporary version file to be installed then delete it
with open('version.py', 'a') as the_file:
    the_file.write('__version__ = "{}"\n'.format(vcstools_version))

setup(name="vcstools",
      version=vcstools_version,
      description="Scripts used to process the Murchison Widefield Array's Voltage Capture System data",
      url="https://github.com/ICRAR/vcstools",
      long_description=read('README.md'),
      #packages=['vcstools'],
      #install_requires=reqs,
      scripts=['scripts/checks.py', 'scripts/calibrate_vcs.py', 'scripts/create_psrfits.sh',
               'scripts/file_maxmin.py', 'scripts/find_pulsar_in_obs.py', 'scripts/obs_query.py',
               'scripts/process_vcs.py', 'scripts/recombine.py', 'scripts/rename_corr_output.py',
               'scripts/reorder_chans.py', 'scripts/rts2ao.py', 'scripts/untar.sh',
               'scripts/cleanup.py', 'scripts/create_ics_psrfits.py',
               'scripts/splice.sh', 'scripts/auto_plot.bash',
               'database/submit_to_database.py', 'database/database_vcs.py',
               'database/known_repeating_FRBs.csv',
               'utils/zapchan.py', 'utils/calc_ephem.py', 'utils/check_disk_usage.sh',
               'utils/check_quota.sh', 'utils/mdir.py', 'utils/mwa_metadb_utils.py',
               'utils/job_submit.py', 'utils/plotFlatTileBeam.py', 'utils/plotPolarTileBeam.py',
               'utils/plotTiedArrayBeam.py', 'utils/aocal.py',
               'utils/config.py', 'utils/sn_flux_est.py',
               'version.py'],
      #data_files=[('AegeanTools', [os.path.join(data_dir, 'MOC.fits')]) ],
      setup_requires=['pytest-runner'],
      tests_require=['pytest']#, 'nose']
)

os.remove('version.py')
