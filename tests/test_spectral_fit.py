#! /usr/bin/env python
"""
Tests the spectral_fit.py script
"""
import os
import numpy as np
from numpy.testing import assert_almost_equal
import psrqpy
import pandas as pd

from vcstools import data_load
from vcstools.spectral_fit import find_best_spectral_fit
from vcstools.sn_flux_utils import flux_from_atnf

import logging
logger = logging.getLogger(__name__)

'''
def test_fit_spectral_model():
    """Tests the fit_spectral_model funtion with examples from table Table Ci in Jankowski et al. 2018
    """
    pulsars = ['J0820-1350', 'J1644-4559']
    query = psrqpy.QueryATNF(loadfromdb=data_load.ATNF_LOC).pandas
    for pulsar in pulsars:
        freq_all, flux_all, flux_err_all, spind, spind_err = flux_from_atnf(pulsar, query=query)
        # convert to GHz
        #freq_all = np.array(freq_all) / 1e9
        print(f"\nFitting {pulsar}")
        print(freq_all, flux_all, flux_err_all)

        find_best_spectral_fit(pulsar, freq_all, flux_all, flux_err_all)
            
    #if a != expected_a or a_err != expected_a_err:
    #    raise AssertionError()
    #assert_almost_equal(a,     expected_a,    decimal=1)
    #sassert_almost_equal(a_err, expected_a_err, decimal=1)
'''

def test_fit_spectral_model_with_chris_data():
    """Tests the fit_spectral_model funtion with examples from table Table Ci in Jankowski et al. 2018
    """
    pulsars = ['J0953+0755', 'J1645-0317']
    for pulsar in pulsars:
        df = pd.read_csv("test_files/{}.csv".format(pulsar))
        print(f"\nFitting {pulsar}")
        freq_all = np.array(df['FREQ'].tolist())*1e6
        flux_all = df['S'].tolist()
        flux_err_all = df['SERR'].tolist()
        print(freq_all, flux_all, flux_err_all)
        models, fit_results = find_best_spectral_fit(pulsar, freq_all, flux_all, flux_err_all, plot=True)
        # convert to GHz
        #freq_all = np.array(freq_all) / 1e9
        print(models)

if __name__ == "__main__":
    """
    Tests the relevant functions in spectral_fit.py
    """
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
