#! /usr/bin/env python

import os
import sys
import requests
import json
import socket

from vcstools.pointing_utils import sex2deg
from vcstools.metadb_utils import getmeta

import logging
logger = logging.getLogger(__name__)

auth=('admin', 'testing123')
base_url = 'http://localhost:8000'
urls = {'search_parameters': base_url+'/search_parameters_create/',
        'ml_parameters': base_url + '/ml_parameters_create/',
        'ml_candidates': base_url + '/ml_candidates_create/',
        'candidates': base_url + '/candidates_create/',
        'ratings': base_url + '/ratings_create/',
        'beams': base_url + '/beams_create/',
        'pulsar': base_url + '/pulsar_create/',
        'supercomputers': base_url + '/supercomputers_create/',
        'users': base_url + '/users_create/',
        'observation_setting': base_url + '/observation_setting_create/',
        'followups': base_url + '/followups_create/',
        'calibrator': base_url + '/calibrator_create/',
        'detection': base_url + '/detection_create/',
        'detection_file':base_url + '/detection_file_create/'}


def upload_wrapper(url, data, files=None):
    """
    Wrapper for the upload to the MWA pulsar database

    Parameters
    ----------
    url: string
        Website url for the relivent table
    data: dict
        Data dict in the format for the table
    files: dict
        Files dict in the format for files in the table. Default None
    """
    # set up a session with the given auth
    session = requests.session()
    session.auth = auth

    #logger.debug(data)
    print(data)
    if files is None:
        r = session.post(url, data=data)
    else:
        r = session.post(url, data=data, files=files)
    if r.status_code > 201:
        logger.error("Error uploading. Status code: {}".format(r.status_code))
        sys.exit(r.status_code)
    session.close()


def find_supercomputer_id():
    """
    Work out which supercomputer you are on and return its ID
    """
    hostname = socket.gethostname()
    if hostname.startswith('john') or hostname.startswith('farnarkle'):
        return 1
    elif hostname.startswith('mwa') or hostname.startswith('garrawarla'):
        return 2
    elif hostname.startswith('x86')  or hostname.startswith('arm'):
        return 3
    else:
        logger.error('Unknown computer {}. Exiting'.format(hostname))
        sys.exit(1)

def upload_beam(name_obs_pos, begin, end,
                vcstools_version=None, mwa_search_version=None,
                mwa_search_command=None, supercomputer_id=None):
    """
    Upload the a pointing to the beams table of the MWA pulsar database

    Parameters
    ----------
    name_obs_pos: string
        String used to sort data in search pipeline in the format name_obsid_raj_decj
    begin: int
        Begin time of the beam in GPS format
    end: int
        End time of the beam in GPS format
    vcstools_version: string
        The vcstools version. By default it will use vcstools.version.__version__
    mwa_search_version: string
        The mwa_search version. By default it will use mwa_search.version.__version__
    mwa_search_command: string
        The mwa_search_pipeline.nf command
    supercompter_id: int
        The ID of the supercomputer. By default will use the supercomputer you're currently on
    """

    # Get info from input name
    name, obsid, raj, decj = name_obs_pos.split("_")
    rad, decd = sex2deg(raj, decj)

    # Get versions
    if vcstools_version is None:
        import vcstools.version
        vcstools_version = vcstools.version.__version__
    if mwa_search_version is None:
        import mwa_search.version
        mwa_search_version = mwa_search.version.__version__

    if supercomputer_id is None:
        supercomputer_id = find_supercomputer_id()

    # Set up data to upload
    data = {'ra_degrees': rad,
            'dec_degrees': decd,
            'begin_gps_seconds': begin,
            'end_gps_seconds': end,
            'username': auth[0],
            'supercomputer_id': supercomputer_id,
            'vcstools_version': vcstools_version,
            'mwa_search_version': mwa_search_version,
            'mwa_search_command': mwa_search_command,
            'observation_id': obsid}

    upload_wrapper(urls['beams'], data)


def upload_obsid(obsid):
    """
    Upload an observation to the observation_settings table of the MWA pulsar database

    Parameters
    ----------
    obsid: int or list of ints
        MWA observation ID or list of IDs
    """
    # Make sure it's a list
    obsids = [obsid] if not isinstance(obsid, list) else obsid

    # Set up metadb keys to upload
    base_keys = ['starttime', 'stoptime', 'obsname', 'creator',
                 'modtime', 'mode', 'unpowered_tile_name', 'mode_params',
                 'dec_phase_center', 'ra_phase_center', 'projectid',
                 'int_time', 'freq_res']#, 'checked_user']
    metadata_keys = ['elevation_pointing','azimuth_pointing',
                     'ra_pointing', 'dec_pointing', 'sun_elevation',
                     'sun_pointing_distance', 'jupiter_pointing_distance',
                     'moon_pointing_distance', 'sky_temp', 'calibration',
                     'gridpoint_name', 'gridpoint_number',
                     'local_sidereal_time_deg', 'calibrators']
    rfstreams_keys = ['azimuth', 'elevation', 'ra', 'dec',
                     'frequencies', 'frequency_type', 'walsh_mode',
                     'gain_control_type', 'gain_control_value', 'dipole_exclusion']
    all_keys = base_keys + metadata_keys + rfstreams_keys


    for obs in obsids:
        data_dict = {'observation_id':int(obs)}
        data = getmeta(params={'obsid':obs})
        for keys, db in zip((base_keys, metadata_keys, rfstreams_keys),
                            (data, data['metadata'],data['rfstreams']['0'])):
            for key in keys:
                if key not in db.keys():
                    data_dict[key] = ""
                elif db[key] is not None:
                    # no 'checked_user' keys anywhere
                    data_dict[key] = db[key]
        upload_wrapper(urls['observation_setting'], data_dict)


def upload_cand(beam_id, search_type='Blind', search_params_id=1, png_file=None, pfd_file=None):
    """
    Upload a candidate to the candidates table of the MWA pulsar database

    Parameters
    ----------
    beam_id: int
        The beam ID of the pointing for the candidate
    search_type: string
        A label of the type of search that produced this candidate. Default Blind
    search_params_id: int
        The search params ID. This links to the search_paramaters table which
        lists the search parameters such as DM range
    png_file: string
        Path to the png file of the candidate
    pfd_file: string
        Path to the pfd file of the candidate
    """
    # Open files
    files = {'png': None, 'pfd': None}
    if png_file is not None:
        png = open(png_file, 'rb')
        files['png'] = png
    if pfd_file is not None:
        pfd = open(pfd_file, 'rb')
        files['pfd'] = pfd

    # Set up data to upload
    data = {'search_beam_id': beam_id,
            'username': auth[0],
            'search_type': search_type,
            'search_params_id': search_params_id}

    upload_wrapper(urls['candidates'], data, files=files)

    # Close files
    if png_file is not None:
        png.close()
    if pfd_file is not None:
        pfd.close()
