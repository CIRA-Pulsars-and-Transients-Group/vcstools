#! /usr/bin/env python

import os
import sys
import requests

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

def upload_beam(name_obs_pos, begin, end,
                vcstools_version=None, mwa_search_version=None, mwa_search_command=None):
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
    """
    # set up a session with the given auth
    session = requests.session()
    session.auth = auth

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

    # Set up data to upload
    data = {'ra_degrees': rad,
            'dec_degrees': decd,
            'begin_gps_seconds': begin,
            'end_gps_seconds': end,
            'supercomputer': ,
            'vcstools_version': vcstools_version,
            'mwa_search_version': mwa_search_version,
            'mwa_search_command': mwa_search_command,
            'observation': obsid}

    r = session.post(urls['beams'], data=data)
