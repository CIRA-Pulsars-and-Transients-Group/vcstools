#! /usr/bin/env python3
"""
Tests the sn_flux_est.py script
"""
import os
from numpy.testing import assert_almost_equal
from vcstools.metadb_utils import get_common_obs_metadata
import psrqpy

from vcstools import data_load

import vcstools.sn_flux_est as snfe

import logging
logger = logging.getLogger(__name__)


#Global obsid information to speed things up
md_dict={}
obsid_list = [1223042480, 1222697776, 1226062160, 1225713560, 1117643248]
for obs in obsid_list:
    md_dict[str(obs)] = get_common_obs_metadata(obs, return_all=True)

def test_pulsar_beam_coverage():
    """
    Tests the pulsar_beam_coverage funtion
    """
    print("test_pulsar_beam_coverage")
    obsid = 1223042480
    pulsar = "J2234+2114"
    test_cases = []
    test_cases.append((1223042487, 1223042587, 0.0, 1.0))

    test_none_cases = []
    test_none_cases.append((-1e10, -1e4, 0.0, 0.072731816627943369))

    md = md_dict[str(obsid)][0]
    full_md = md_dict[str(obsid)][1]

    for beg, end, exp_enter, exp_exit in test_cases:
        out_enter, out_exit = snfe.pulsar_beam_coverage(obsid, pulsar, beg, end, metadata=md, full_meta=full_md)
        assert_almost_equal(exp_enter, out_enter, decimal=4)
        assert_almost_equal(exp_exit, out_exit, decimal=4)

    #seperately test the none tests because assert_almost_equal doesn't handle Nones
    for beg, end, exp_enter, exp_exit in test_none_cases:
        out_enter, out_exit = snfe.pulsar_beam_coverage(obsid, pulsar, beg, end, metadata=md, full_meta=full_md)
        if out_enter is not None or out_exit is not None:
            raise AssertionError()


def test_find_times():
    """
    Tests the find_times function
    """
    print("find_times")
    obsid = 1223042480
    pulsar = "J2234+2114"
    test_cases = []
    test_cases.append((1223042487, 1223042587, 1223042487, 1223042587, 100.0))

    md = md_dict[str(obsid)][0]
    full_md = md_dict[str(obsid)][1]

    for beg, end, exp_beg, exp_end, exp_int in test_cases:
        out_beg, out_end, out_int = snfe.find_times(obsid, pulsar, beg, end, metadata=md, full_meta=full_md)
        assert_almost_equal(exp_beg, out_beg, decimal=4)
        assert_almost_equal(exp_end, out_end, decimal=4)
        assert_almost_equal(exp_int, out_int, decimal=4)


def test_find_t_sys_gain():
    """
    Tests the find_t_sys_gain function
    """
    print("find_t_sys_gain")
    obsid=1222697776
    pulsar= "J2330-2005"
    beg = 1222697776 + 7
    end = beg + 3599

    md = md_dict[str(obsid)][0]
    full_md = md_dict[str(obsid)][1]

    q = psrqpy.QueryATNF(psrs=pulsar, loadfromdb=data_load.ATNF_LOC)
    test_cases = []
    test_cases.append((pulsar, obsid, beg, end, q, md, full_md,\
                295.7376549611132, 5.914753099222264, 0.31480323776600816, 0.05136230390439281))


    for pulsar, obsid, beg, end query, metadata, full_meta,\
        exp_t_sys, exp_t_sys_err, exp_gain, exp_gain_err in test_cases:
        t_sys, t_sys_err, gain, gain_err = snfe.find_t_sys_gain(pulsar, obsid, beg, end,\
                query=query, obs_metadata=metadata, full_meta=full_meta)
        assert_almost_equal(exp_t_sys, t_sys, decimal=4)
        assert_almost_equal(exp_t_sys_err, t_sys_err, decimal=4)
        assert_almost_equal(exp_gain, gain, decimal=4)
        assert_almost_equal(exp_gain_err, gain_err, decimal=4)

def test_find_pulsar_w50():
    """
    Tests the find_pulsar_w50 function
    """
    print("find_pulsar_w50")
    test_cases=[]
    test_cases.append(("J1614-2230", 2.3581508027839326e-06, 1.9651256689866104e-06))

    for psr, exp_w50, exp_w50_err in test_cases:
        w50, w50_err = snfe.find_pulsar_w50(psr)
        assert_almost_equal(exp_w50, w50, decimal=4)
        assert_almost_equal(exp_w50_err, w50_err, decimal=4)


def test_est_pulsar_flux():
    """
    Tests the est_pulsar_flux function
    """
    print("est_pulsar_flux")

    test_cases = []
    test_cases.append(("J2330-2005", 1226062160, 0.13693009478998122, 0.02425180108093204))

    for psr, obsid, exp_flux, exp_flux_err in test_cases:
        flux, flux_err = snfe.est_pulsar_flux(psr, obsid)
        assert_almost_equal(exp_flux, flux, decimal=4)
        assert_almost_equal(exp_flux_err, flux_err, decimal=4)


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
    test_cases.append(("J2241-5236", obsid, None, None, None, None, md, full_md,\
                    280.9035270638417, 102.95827699010243))
    #Has 6 fluxes on database
    obsid=1226062160
    md = md_dict[str(obsid)][0]
    full_md = md_dict[str(obsid)][1]
    test_cases.append(("J2330-2005", obsid, None, None, "23:00:00", "-20:00:00", md, full_md,\
                    142.0052581169168, 45.315968765282484))
    #adding a test that gave a none type error for alpha_bound, c_bound
    obsid=1117643248
    md = md_dict[str(obsid)][0]
    full_md = md_dict[str(obsid)][1]
    test_cases.append(("J1623-2631", obsid, 1117643268, 1117645615, None, None, md, full_md,\
                    49.80414834502304, 12.304114528369126))

    for psr, obsid, beg, end, ra, dec, md, full_md, exp_sn, exp_sn_err in test_cases:
        sn, sn_err = snfe.est_pulsar_sn(psr, obsid, beg=beg, end=end, p_ra=ra, p_dec=dec, obs_metadata=md, full_meta=full_md)
        assert_almost_equal(exp_sn, sn, decimal=2)
        assert_almost_equal(exp_sn_err, sn_err, decimal=2)


if __name__ == "__main__":
    """
    Tests the relevant functions in sn_flux_est.py
    Uses psrcat version 1.59. Values may change for different versions
    """

    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
