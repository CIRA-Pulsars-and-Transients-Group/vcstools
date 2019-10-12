#! /usr/bin/env python3
"""
Tests the find_pulsar_in_obs.py script
"""
import checks
from numpy.testing import assert_approx_equal, assert_almost_equal


def test_get_files_and_sizes():
    """Test the get_files_and_sizes function"""
    #obsid, mode, file_loc, suffix, number
    tests = [[1118463408, 'raw',    'tests/test_files/checks_get_files_1118463408_raw.txt', '.dat', 253440000],
             [1221832280, 'tar_ics','tests/test_files/checks_get_files_1221832280_tar_ics.txt','.tar', 7865368576],
             [1221832280, 'ics',    'tests/test_files/checks_get_files_1221832280_ics.txt', '_ics.dat', 30720000]]
    for test in tests:
        obsid, mode, file_loc, suffix, number = test
        # Read in file list that was too long to store in this function
        with open(file_loc, 'r') as file_read:
            file_list = []
            for line in file_read:
                file_list.append(line[:-1])
        expected_ans = (file_list, suffix, number)
        ans = checks.get_files_and_sizes(obsid, mode)
        if ans != expected_ans:
            raise AssertionError()


if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
