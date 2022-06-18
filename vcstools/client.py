#! /usr/bin/env python

import os
import sys
import requests
from requests.exceptions import HTTPError
import json
import socket
from pathlib import Path
from urllib.parse import urljoin
import concurrent.futures
import time
import numpy as np

from vcstools.pointing_utils import sex2deg
from vcstools.metadb_utils import getmeta

from presto.prepfold import pfd

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

TABLE_TO_PATH = {
    'search_parameters': './search_parameters_create/',
    'ml_parameters': './ml_parameters_create/',
    'ml_candidates': './ml_candidates_create/',
    'candidates': './candidates_create/',
    'ratings': './ratings_create/',
    'beams': './beams_create/',
    'pulsar': './pulsar_create/',
    'supercomputers': './supercomputers_create/',
    'users': './users_create/',
    'observation_setting': './observation_setting_create/',
    'followups': './followups_create/',
    'calibrator': './calibrator_create/',
    'detection': './detection_create/',
    'detection_file': './detection_file_create/'
}
MAX_THREADS = 8
DEFAULT_TARGET = 'https://apps.datacentral.org.au/smart/'
TOKEN_ENV_VAR_NAME = "SMART_TOKEN"
BASE_URL_ENV_VAR_NAME = "SMART_BASE_URL"


class TokenAuth(requests.auth.AuthBase):
    def __init__(self, token):
        self.token = token

    def __call__(self, r):
        r.headers['Authorization'] = "Token {}".format(self.token)
        return r


class RowUploader:
    def __init__(self, *, session, table, directory, base_url):
        """
        parameters:
        -----------
        session : `requests.session`
            A connection session
        table : str
            The name of the table that we are uploading to
        directory : str
            The base directory that we prepend to all filenames
        """
        self.session = session
        self.table = table
        self.directory = Path(directory)
        self.url = urljoin(base_url, TABLE_TO_PATH[table])

    def __call__(self, row):
        """
        parameters:
        -----------
        row : dict
            The data that is being uploaded
        """
        data = {}
        for key, val in row.items():
            if isinstance(val, str):
                val = val.strip()
            data[key.lower().strip()] = val

        if self.table == 'candidates':
            # fidget the data/files to include the files for upload
            files = {}
            if 'pfd_path' in data:
                pfd =  open(self.directory.joinpath(data['pfd_path']), 'rb')
                files['pfd'] = pfd
                del data['pfd_path']
            if 'png_path' in data:
                png = open(self.directory.joinpath(data['png_path']), 'rb')
                files['png'] = png
                del data['png_path']
            r = self.session.post(self.url, data=data, files=files)
            if 'pfd' in files:
                pfd.close()
            if 'png' in files:
                png.close()
        elif self.table == 'calibrator':
            # fidget the data/files to include the files for upload
            with open(self.directory.joinpath(data['cal_file']), 'rb') as cal:
                files = {'cal_file': cal}
                del data['cal_file']
                r = self.session.post(self.url, data=data, files=files)
        elif self.table == 'detection_file':
            # fidget the data/files to include the files for upload
            with open(self.directory.joinpath(data['path']), 'rb') as df:
                files = {'path': df}
                del data['path']
                r = self.session.post(self.url, data=data, files=files)
        else:
            r = self.session.post(self.url, data=data)

        print("Response: %s", r.text) # Because I cannot get logger to output anything
        logger.debug("Response: %s", r.text)
        try:
            r.raise_for_status()
        except HTTPError:
            logger.warning(r.text)
            # wait a few second and try again
            time.sleep(2)
            r = self.session.post(self.url, data=data)
            r.raise_for_status()

        return "{0} : {1}".format(data, r.status_code)


def upload_wrapper(
        data_list,
        table,
        directory='.',
        base_url=DEFAULT_TARGET,
        max_threads=MAX_THREADS,
    ):
    """
    Wrapper for the upload to the MWA pulsar database

    Parameters
    ----------
    data_list : `list` of `dict`
        A list of each data dict to upload in the format for the table
    table : `str``
        Destination table name
    files : `dict`
        Files dict in the format for files in the table. Default None
    """
    token = os.environ.get(TOKEN_ENV_VAR_NAME)
    if token is None:
        logger.error("Token not found, set {}".format(TOKEN_ENV_VAR_NAME))

    # set up a session with the given auth
    session = requests.session()
    session.auth = TokenAuth(token)
    row_uploader = RowUploader(
        session=session, table=table, directory=directory, base_url=base_url,
    )

    logger.info("Sending rows to {} using {} threads".format(table, max_threads))
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=max_threads
    ) as executor:
        for result in executor.map(row_uploader, data_list):
            logger.info(result)
    logger.info("Completed")
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


def upload_beam(
        pointing_list,
        obsid,
        begin,
        end,
        vcstools_version=None,
        mwa_search_version=None,
        mwa_search_command=None,
        supercomputer_id=None
    ):
    """
    Upload the a pointing to the beams table of the MWA pulsar database

    Parameters
    ----------
    pointing_list : `list`
        A list of pointings in the format "HH:MM:SS.SS_+DD:MM:SS.SS".
    begin : `int`
        Begin time of the beam in GPS format
    end : `int`
        End time of the beam in GPS format
    vcstools_version : `str`, optional
        The vcstools version. By default it will use vcstools.version.__version__
    mwa_search_version : `str`, optional
        The mwa_search version. By default it will use mwa_search.version.__version__
    mwa_search_command : `str`, optional
        The mwa_search_pipeline.nf command
    supercompter_id : `int`, optional
        The ID of the supercomputer. By default will use the supercomputer you're currently on
    """
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
    data_list = []
    for point in pointing_list:
        raj, decj = point.split("_")
        rad, decd = sex2deg(raj, decj)
        data = {'ra_degrees': rad,
                'dec_degrees': decd,
                'begin_gps_seconds': begin,
                'end_gps_seconds': end,
                'supercomputer_id': supercomputer_id,
                'vcstools_version': vcstools_version,
                'mwa_search_version': mwa_search_version,
                'mwa_search_command': mwa_search_command,
                'observation_id': obsid}
        data_list.append(data)
        print("Adding beam (RA {}, Dec {})".format(data['ra_degrees'], data['dec_degrees']))

    upload_wrapper(data_list, 'beams')


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

    data_list = []
    for obs in obsids:
        data_dict = {'observation_id':str(obs)}
        data = getmeta(params={'obsid':obs})
        for keys, db in zip((base_keys, metadata_keys, rfstreams_keys),
                            (data, data['metadata'],data['rfstreams']['0'])):
            for key in keys:
                if key not in db.keys():
                    data_dict[key] = ""
                elif key == 'sky_temp' and db[key] is None:
                    # If no skytemp in the getmeta make an approximation
                    data_dict[key] = 240.
                elif db[key] is not None:
                    # no 'checked_user' keys anywhere
                    data_dict[key] = db[key]
        data_list.append(data_dict)

    upload_wrapper(data_list, 'observation_setting')


def upload_cand(
        pfd_file_list,
        obsid,
        search_type='Blind',
        search_params_id=1,
    ):
    """
    Upload a candidate to the candidates table of the MWA pulsar database

    Parameters
    ----------
    pfd_file_list : `list`
        List of pfd file locations to upload.
    obsid: int
        The observation ID of the candidate/beam
    search_type: `string`
        A label of the type of search that produced this candidate. Default 'Blind'.
    search_params_id: `int`
        The search params ID. This links to the search_paramaters table which
        lists the search parameters such as DM range. Default: 1.
    """
    # Get data from pfd
    data_list = []
    for pfd_file in pfd_file_list:
        PFDObject = pfd(pfd_file)
        raj = PFDObject.rastr.decode("utf-8")
        decj = PFDObject.decstr.decode("utf-8")
        rad, decd = sex2deg(raj, decj)
        data = {
            'rad': rad,
            'decd': decd,
            'obsid': obsid,
            'pfd_path': pfd_file,
            'period':PFDObject.bary_p1,
            'dm':PFDObject.bestdm,
            'sigma':PFDObject.calc_sigma(),
            'search_type': search_type,
            'search_params_id': search_params_id,
        }

        # Check for png file in same location
        if os.path.isfile(pfd_file+".png"):
            data['png_path'] = pfd_file+".png"
        data_list.append(data)
        print("Adding candidate (RA {}, Dec {})".format(raj, decj))

    upload_wrapper(data_list, 'candidates')


def upload_supercomputer(supercomputers, id_list=None):
    """Upload new supercomputers to the supercomputer table of the MWA pulsar database.

    Parameters
    ----------
    supercomputers : `list`
        List of supercomputer strings.
    """
    data_list = []
    for si, sc in enumerate(supercomputers):
        data = {'Name': sc}
        if id_list is None:
            data['id'] = si
        data_list.append(data)

    upload_wrapper(data_list, 'supercomputers')


def upload_defaults():
    """Upload the defaults to the MWA pulsar database for ml_parameters, supercomputers.
    """
    # ml_parameters
    data_list = [
                 {'ID':1,
                  'Name':'LOTAAS_periodic_classifier',
                  'ML_version':'v1.0',
                  'Training_set_version':'v1.0',
                  'Command':'LOTAASClassifier.jar',
                  'Commnents':"LOTAAS's original version for periodic canidates"},
                 {'ID':2,
                  'Name':'LOTAAS_single_pulse_classifier',
                  'ML_version':'v1.0',
                  'Training_set_version':'v1.0',
                  'Command':'single_pulse_searcher.py',
                  'Commnents':"LOTAAS's original version for single pulse canidates"},
                ]
    upload_wrapper(data_list, 'ml_parameters')

    # search_parameters
    data_list = [{'ID':1,
                  'DM_min':1,
                  'DM_max':250,
                  'DM_min_step':0.01,
                  'Acceleration_search_max':0.0,
                  'Accelsearch_sigma_cutoff':10.}]
    upload_wrapper(data_list, 'search_parameters')

    # supercomputers
    upload_supercomputer(['Ozstar', 'Garrawarla (Pawsey)', 'SHAO'])

def upload_pulsars():
    """Upload pulsar parameters using the ATNF pulsar catalogue.
    """
    data_list = []
    import psrqpy
    query = psrqpy.QueryATNF().pandas
    for qid in range(len(query["PSRJ"])):
        RA, Dec = sex2deg(query["RAJ"][qid], query["DECJ"][qid])
        dm = query["DM"][qid]
        # Turn dm and period nans into Nones
        if np.isnan(dm):
            dm = None
        period = query["P0"][qid]
        if np.isnan(period):
            period = None
        data_list.append({"id":qid,
                          "name":query["PSRJ"][qid],
                          "ra":RA,
                          "dec":Dec,
                          "period":period,
                          "dm":dm,
                          "new":False,
                         })
    upload_wrapper(data_list, 'pulsar')
