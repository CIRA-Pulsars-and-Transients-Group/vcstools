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
def format_version(version, fmt='{tag}.{commitcount}'):
    parts = version.split('-')
    if len(parts) == 1:
        return parts[0]
    assert len(parts) in (3, 4)
    dirty = len(parts) == 4
    tag, count, sha = parts[:3]
    if count == '0' and not dirty:
        return tag
    return fmt.format(tag=tag, commitcount=count)

def get_git_version():
    git_version = check_output('git describe --tags --long --dirty --always'.split()).decode('utf-8').strip()
    return format_version(version=git_version)


reqs = ['astropy>=3.2.1',
        'argparse>=1.4.0',
        'numpy>=1.13.3',
        'matplotlib>=2.1.0',
        'psrqpy>=1.0.5',
        'mpi4py',
        #mwa software
        'mwa-voltage',
        'mwa_pb',
        'mwa-hyperbeam',
        'pulsar_spectra',
       ]

# Download the ANTF_pulsar_database_file file if it doesn't exist
datadir = os.path.join(os.path.dirname(__file__), 'vcstools', 'data')


#vcstools_version = get_git_version()
vcstools_version = "2.7.1"
#make a temporary version file to be installed then delete it
with open('version.py', 'a') as the_file:
    the_file.write(f'__version__ = "{vcstools_version}"\n')

setup(name="mwa_vcstools",
      version=vcstools_version,
      description="Scripts used to process the Murchison Widefield Array's Voltage Capture System data",
      url="https://github.com/CIRA-Pulsars-and-Transients-Group/vcstools.git",
      #long_description=read('README.md'),
      packages=['vcstools'],
      package_data={'vcstools':['data/*.csv', 'data/*.db']},
      python_requires='>=3.7',
      install_requires=reqs,
      scripts=[# bash
               'scripts/untar.sh', 'scripts/create_psrfits.sh', 'scripts/splice.sh',
               'scripts/check_disk_usage.sh', 'scripts/auto_plot.bash',
               # python
               'scripts/checks.py', 'scripts/calibrate_vcs.py', 'scripts/submit_to_database.py',
               'scripts/find_pulsar_in_obs.py', 'scripts/RVM_fit.py', 'scripts/mwa_metadb_utils.py',
               'scripts/process_vcs.py', 'scripts/recombine.py', 'scripts/rename_corr_output.py',
               'scripts/reorder_chans.py', 'scripts/rts2ao.py', 'scripts/splice_wrapper.py',
               'scripts/cleanup.py', 'scripts/create_ics_psrfits.py', 'scripts/rm_synthesis.py',
               'scripts/zapchan.py', 'scripts/sn_flux_est.py', 'scripts/prof_estimate.py',
               'scripts/pabeam.py',
               # plotting scripts
               'scripts/plotting/plotPolarTileBeam.py', 'scripts/plotting/plotFlatTileBeam.py',
               'scripts/plotting/plotTiedArrayBeam.py', 'scripts/plotting/plotSkyMap.py',
               'scripts/plotting/plot_BPcal_128T.py', 'scripts/plotting/calc_ephem.py',
               # temporary automatically generated version file
               'version.py'],
      setup_requires=['pytest-runner'],
      tests_require=['pytest']
)

# remove files
os.remove('version.py')
if os.path.isfile('record.txt'):
    os.remove('record.txt')
