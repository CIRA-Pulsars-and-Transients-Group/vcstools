#! /usr/bin/env python
"""
Tests the radiometer_equation.py script
"""
import os
from numpy.testing import assert_almost_equal
from vcstools.metadb_utils import get_common_obs_metadata
import psrqpy

from vcstools import data_load
from vcstools.radiometer_equation import find_t_sys_gain, find_pulsar_w50,\
                                         est_pulsar_flux, est_pulsar_sn

import logging
logger = logging.getLogger(__name__)


#Global obsid information to speed things up
md_dict={}
obsid_list = [1222697776, 1226062160, 1225713560, 1117643248]
for obs in obsid_list:
    md_dict[str(obs)] = get_common_obs_metadata(obs, return_all=True)

query = psrqpy.QueryATNF(loadfromdb=data_load.ATNF_LOC).pandas

def test_find_t_sys_gain():
    """
    Tests the find_t_sys_gain function
    """
    print("find_t_sys_gain")
    obsid=1222697776
    pulsar= "J2330-2005"
    # Full obs
    obs_beg = 1222697783
    obs_end = 1222702696

    md = md_dict[str(obsid)][0]
    full_md = md_dict[str(obsid)][1]

    test_cases = []
    test_cases.append((pulsar, obsid, obs_beg, obs_end, md, full_md,
                295.7376549611132, 5.914753099222264, 0.083530955951217, 0.01882702843010659))


    for pulsar, obsid, obs_beg, obs_end, metadata, full_metadata,\
        exp_t_sys, exp_t_sys_err, exp_gain, exp_gain_err in test_cases:
        t_sys, t_sys_err, gain, gain_err = find_t_sys_gain(pulsar, obsid,
                                                           obs_beg=obs_beg, obs_end=obs_end,
                                                           query=query,
                                                           common_metadata=metadata,
                                                           full_metadata=full_metadata)
        assert_almost_equal(t_sys,     exp_t_sys,     decimal=4)
        assert_almost_equal(t_sys_err, exp_t_sys_err, decimal=4)
        assert_almost_equal(gain,      exp_gain,      decimal=4)
        assert_almost_equal(gain_err,  exp_gain_err,  decimal=4)

def test_find_pulsar_w50():
    """
    Tests the find_pulsar_w50 function
    """
    print("find_pulsar_w50")
    test_cases=[]
    test_cases.append(("J1614-2230", 2.3581508027839326e-06, 1.9651256689866104e-06))

    for psr, exp_w50, exp_w50_err in test_cases:
        w50, w50_err = find_pulsar_w50(psr, query=query)
        assert_almost_equal(w50,     exp_w50,     decimal=4)
        assert_almost_equal(w50_err, exp_w50_err, decimal=4)


def test_est_pulsar_flux():
    """
    Tests the est_pulsar_flux function
    """
    print("est_pulsar_flux")

    test_cases = []
    # below was true with ATNF catalogue version 1.63
    #test_cases.append(("J2330-2005", 1226062160, 0.13693009478998122, 0.02425180108093204))
    test_cases.append(("J2330-2005", 1226062160, 0.3253646594291489, 0.10930747950834567))

    for psr, obsid, exp_flux, exp_flux_err in test_cases:
        flux, flux_err = est_pulsar_flux(psr, obsid, query=query)
        print(flux, flux_err)
        assert_almost_equal(flux,     exp_flux,     decimal=4)
        assert_almost_equal(flux_err, exp_flux_err, decimal=4)


def test_est_pulsar_sn():
    """
    I have tested this function various pulsar/observation combinations and compared to the estimated S/N from bestrof files. Here are the results:

    1225462936, -b 1225462943 -e 1225463543 psr J0152-1637: The function estimated a S/N of 100.92 +/- 22.48.
    I compared this to what was estimated from the profile analyse_pulse_prof() and attained a S/N of 75.26 +/- 1.70

    1226062160, -b 1226062167 -e 1226062767 psr J2330-2005: Function estimated a S/N of 412.26 +/- 103.80.
    I compared this to the profile estimation and attained S/N of 51.51 +/- 1.16

    1225713560, -b 1225713567 -e 1225714167, psr J2241-5236: Function estimated a S/N of 119.82 +/- 41.61.
    I Compared this to the profile estimation and attained S/N of 17.36 +/- 2.41

    """
    print("est_pulsar_sn")
    test_cases=[]
    #Has 3 fluxes on database
    obsid=1225713560
    md = md_dict[str(obsid)][0]
    full_md = md_dict[str(obsid)][1]
    # Below is true for ATNF version 1.62
    #test_cases.append(("J2241-5236", obsid, 1225713567, 1225714167, None, None, md, full_md,
    #                   175.5310369703948, 63.791707034461005, 0.2845788591631629, 0.03456583343266214))
    #                  psr, obsid, beg, end, ra, dec, md, full_md
    test_cases.append(("J2241-5236", obsid, 1225713567, 1225714167, None, None, md, full_md,
    #                  exp_sn, exp_sn_err, exp_flux, exp_flux_err
                       164.1305464077595, 67.58112355899925, 0.2788556217863279, 0.0628397803003623))
    #Has 6 fluxes on database
    obsid=1226062160
    md = md_dict[str(obsid)][0]
    full_md = md_dict[str(obsid)][1]
    test_cases.append(("J2330-2005", obsid, 1226062167, 1226062767, "23:00:00", "-20:00:00", md, full_md,
                       209.2713035503201, 89.60097245078296, 0.3253646594291489, 0.10930747950834567))
    #adding a test that gave a none type error for alpha_bound, c_bound
    obsid=1117643248
    md = md_dict[str(obsid)][0]
    full_md = md_dict[str(obsid)][1]
    test_cases.append(("J1623-2631", obsid, 1117643268, 1117645615, None, None, md, full_md,
                       29.467042642173908, float("nan"), 0.3666299299976159, float("nan")))

    for psr, obsid, beg, end, ra, dec, md, full_md, exp_sn, exp_sn_err, exp_flux, exp_flux_err in test_cases:
        sn, sn_err, flux, flux_err = est_pulsar_sn(psr, obsid, dect_beg=beg, dect_end=end, p_ra=ra, p_dec=dec, common_metadata=md, full_metadata=full_md, query=query, plot_flux=True)
        print(sn, sn_err, flux, flux_err)
        assert_almost_equal(sn,       exp_sn,       decimal=2)
        assert_almost_equal(sn_err,   exp_sn_err,   decimal=2)
        assert_almost_equal(flux,     exp_flux,     decimal=2)
        assert_almost_equal(flux_err, exp_flux_err, decimal=2)

if __name__ == "__main__":
    """
    Tests the relevant functions in radiometer_equation.py
    Uses psrcat version 1.59. Values may change for different versions
    """
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
