#! /usr/bin/env python3
"""Compares the results of the FEE beam model using mwa_pb and mwa_hyperbeam."""
import numpy as np
from numpy.testing import assert_almost_equal
import time

from mwa_pb import primary_beam
from mwa_pb import config
import mwa_hyperbeam

def test_hyperbeam_vs_pb():
    """Compares the results of the FEE beam model using mwa_pb and mwa_hyperbeam."""
    beam = mwa_hyperbeam.FEEBeam(config.h5file)

    # Set up fake data.
    n = 1000000
    az = np.linspace(0, 0.9 * np.pi, n)
    za = np.linspace(0.1, 0.9 * np.pi / 2, n)
    freq = 167000000
    delays = [0] * 16
    amps = [1.0] * 16
    test_decimals = 4

    # Jones --------------------------------------------------
    # mwa_pb method
    start_time = time.perf_counter()
    pb_jones = primary_beam.MWA_Tile_full_EE(za, az,
                                                freq=freq, delays=np.array(delays),
                                                zenithnorm=True,
                                                power=True,
                                                jones=True,
                                                interp=False)
    print("mwa_pb benchmark: {} s".format(time.perf_counter()-start_time))

    # hyperbeam method
    #print(freq, delays, amps)
    start_time = time.perf_counter()
    hb_jones = beam.calc_jones_array(az, za, freq, delays, amps, True)
    print("hyperbeam benchmark: {} s".format(time.perf_counter()-start_time))
    hb_jones = hb_jones.reshape(n, 2, 2)

    # Compare Jones
    assert_almost_equal(pb_jones, hb_jones, decimal=test_decimals)

    # Power ---------------------------------------------------
    # mwa_pb method
    pb_xx, pb_yy = primary_beam.MWA_Tile_full_EE(za, az,
                                               freq=freq, delays=np.array(delays),
                                               zenithnorm=True,
                                               power=True)

    # hyperbeam method
    hb_jones = hb_jones.reshape(1, n, 2, 2)
    vis = primary_beam.mwa_tile.makeUnpolInstrumentalResponse(hb_jones, hb_jones)
    hb_xx, hb_yy = (vis[:, :, 0, 0].real, vis[:, :, 1, 1].real)

    # Compare power
    assert_almost_equal(pb_xx, hb_xx[0], decimal=test_decimals)
    assert_almost_equal(pb_yy, hb_yy[0], decimal=test_decimals)


if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
