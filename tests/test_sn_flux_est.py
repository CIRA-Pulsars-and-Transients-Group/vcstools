#! /usr/bin/env python3
"""
Tests the sn_flux_est.py script
"""
import os
from numpy.testing import assert_almost_equal
import mwa_metadb_utils
import psrqpy
import sn_flux_est as snfe

import logging
logger = logging.getLogger(__name__)

try:
    ATNF_LOC = os.environ['PSRCAT_FILE']
except KeyError:
    print("ATNF database could not be found on disk.")
    ATNF_LOC = None


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
    test_none_cases.append((None, None, 0.0, 0.072731816627943369))
    test_none_cases.append((-1e10, -1e4, 0.0, 0.072731816627943369))

    for beg, end, exp_enter, exp_exit in test_cases:
        out_enter, out_exit = snfe.pulsar_beam_coverage(obsid, pulsar, beg, end)
        assert_almost_equal(exp_enter, out_enter, decimal=6)
        assert_almost_equal(exp_exit, out_exit, decimal=6)

    #seperately test the none tests because assert_almost_equal doesn't handle Nones
    for beg, end, exp_enter, exp_exit in test_none_cases:
        out_enter, out_exit = snfe.pulsar_beam_coverage(obsid, pulsar, beg, end)
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
    test_cases.append((None, None, 1223042487.0, 1223042844.4768934, 357.4768934249878))
    test_cases.append((1223042487, 1223042587, 1223042487, 1223042587, 100.0))

    for beg, end, exp_beg, exp_end, exp_int in test_cases:
        out_beg, out_end, out_int = snfe.find_times(obsid, pulsar, beg=beg, end=end)
        assert_almost_equal(exp_beg, out_beg, decimal=6)
        assert_almost_equal(exp_end, out_end, decimal=6)
        assert_almost_equal(exp_int, out_int, decimal=6)


def test_find_t_sys_gain():
    """
    Tests the find_t_sys_gain function
    """
    print("find_t_sys_gain")
    obsid_1=1223042480
    pulsar_1= "J2234+2114"
    obsid_2=1222697776
    pulsar_2= "J2330-2005"
    md_1 = mwa_metadb_utils.get_common_obs_metadata(obsid_1) #this is to speed up the tests
    md_2 = mwa_metadb_utils.get_common_obs_metadata(obsid_2)
    q_1 = psrqpy.QueryATNF(psrs=[pulsar_1], loadfromdb=ATNF_LOC)
    q_2 = psrqpy.QueryATNF(psrs=[pulsar_2], loadfromdb=ATNF_LOC)
    test_cases = []
    test_cases.append((pulsar_1, obsid_1, 1223042887, 600, q_1, md_1,\
                325.17339020027805, 6.5034678040055613, 0.12547449408456718, 0.095633857616224477))
    test_cases.append((pulsar_2, obsid_2, None, None, q_2, md_2,\
                295.73944858721245, 5.9147889717442492, 0.29127750035795497, 0.049367868667752078))


    for pulsar, obsid, beg, t_int, query, metadata,\
        exp_t_sys, exp_t_sys_err, exp_gain, exp_gain_err in test_cases:
        t_sys, t_sys_err, gain, gain_err = snfe.find_t_sys_gain(pulsar, obsid, beg=beg,
                t_int=t_int, query=query, obs_metadata=metadata,
                trcvr='database/MWA_Trcvr_tile_56.csv')
        assert_almost_equal(exp_t_sys, t_sys, decimal=6)
        assert_almost_equal(exp_t_sys_err, t_sys_err, decimal=6)
        assert_almost_equal(exp_gain, gain, decimal=6)
        assert_almost_equal(exp_gain_err, gain_err, decimal=6)

def test_find_pulsar_w50():
    """
    Tests the find_pulsar_w50 function
    """
    print("find_pulsar_w50")
    test_cases=[]
    test_cases.append(("J1614-2230", 2.3581508027839326e-06, 1.9651256689866104e-06))

    for psr, exp_w50, exp_w50_err in test_cases:
        w50, w50_err = snfe.find_pulsar_w50(psr)
        assert_almost_equal(exp_w50, w50, decimal=6)
        assert_almost_equal(exp_w50_err, w50_err, decimal=6)


def test_est_pulsar_flux():
    """
    Tests the est_pulsar_flux function
    """
    print("est_pulsar_flux")

    test_cases = []
    test_cases.append(("J2330-2005", 1226062160, 0.15964747164099558, 0.02807867696250783))

    for psr, obsid, exp_flux, exp_flux_err in test_cases:
        flux, flux_err = snfe.est_pulsar_flux(psr, obsid)
        assert_almost_equal(exp_flux, flux, decimal=6)
        assert_almost_equal(exp_flux_err, flux_err, decimal=6)


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
    test_cases.append(("J2241-5236", 1225713560, None, None, None, None, None, None,\
                    317.75361997000766, 111.37467565933801))
    #Has 6 fluxes on database
    test_cases.append(("J2330-2005", 1226062160, None, None, "23:00:00", "-20:00:00", None, None,\
                    293.1886829824123, 93.36097783427014))
    #adding a test that gave a none type error for alpha_bound, c_bound
    test_cases.append(("J1623-2631", 1117643248, 1117643268, 1117645615, None, None, 0.000, 0.842,\
                    49.81658085212985, 12.307207767010034))


    for psr, obsid, beg, end, ra, dec, enter, exit, exp_sn, exp_sn_err in test_cases:
        sn, sn_err = snfe.est_pulsar_sn(psr, obsid, beg=beg, end=end, p_ra=ra, p_dec=dec,\
                    o_enter=enter, o_exit=exit, trcvr='database/MWA_Trcvr_tile_56.csv')
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
