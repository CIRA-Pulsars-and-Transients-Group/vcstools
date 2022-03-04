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

def download_ANTF_pulsar_database_file(datadir, version="v1.65"):
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
        psrcat_zip_dir = urllib.request.urlretrieve('https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.{}.tar.gz'.format(version))[0]
        # Unzip it
        with gzip.open(psrcat_zip_dir,  'rb') as f_in:
            with open('psrcat_pkg.{}.tar'.format(version), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        # Untar the file we require
        psrcat_tar = tarfile.open(psrcat_zip_dir)
        # Do some python magic to no download the file within it's subdirectory from
        # https://stackoverflow.com/questions/8405843/python-extract-using-tarfile-but-ignoring-directories
        member = psrcat_tar.getmember('psrcat_tar/psrcat.db')
        member.name = os.path.basename(member.name)
        psrcat_tar.extract(member, path=datadir)
        print("Download complete")
        os.remove('psrcat_pkg.{}.tar'.format(version))

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
download_ANTF_pulsar_database_file(datadir)


#vcstools_version = get_git_version()
vcstools_version = "2.7"
#make a temporary version file to be installed then delete it
with open('version.py', 'a') as the_file:
    the_file.write('__version__ = "{}"\n'.format(vcstools_version))

setup(name="mwa_vcstools",
      version=vcstools_version,
      description="Scripts used to process the Murchison Widefield Array's Voltage Capture System data",
      url="https://github.com/CIRA-Pulsars-and-Transients-Group/vcstools.git",
      #long_description=read('README.md'),
      packages=['vcstools'],
      package_data={'vcstools':['data/*.csv', 'data/*.db']},
      python_requires='>=3.6',
      install_requires=reqs,
      scripts=[# bash
               'scripts/untar.sh', 'scripts/create_psrfits.sh', 'scripts/splice.sh',
               'scripts/check_disk_usage.sh', 'scripts/check_quota.sh', 'scripts/auto_plot.bash',
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
