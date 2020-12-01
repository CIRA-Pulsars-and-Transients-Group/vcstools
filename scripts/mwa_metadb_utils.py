#!/usr/bin/env python3
import logging
import argparse

from vcstools.metadb_utils import getmeta, write_obs_info, get_obs_array_phase, obs_max_min, calc_ta_fwhm, get_best_cal_obs

logger = logging.getLogger(__name__)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Returns information on an input Obs ID""")
    parser.add_argument("obsid", type=int, help="Input Observation ID")
    parser.add_argument("-w", "--write", action="store_true", help="OPTIONAL - Use to write results to file.")
    parser.add_argument("-c", "--cal_best", action="store_true", help="If this option is used it will list "
                        "calibration observations within 2 days that have the same observing channels and "
                        "list them from closest in time to furthest.")
    args = parser.parse_args()

    if args.write:
        write_obs_info(args.obsid)
    else:
        data_dict = getmeta(params={"obsid":args.obsid, 'nocache':1})
        channels = data_dict["rfstreams"]["0"]["frequencies"]
        centre_freq = ( min(channels) + max(channels) ) / 2. * 1.28
        array_phase = get_obs_array_phase(args.obsid)
        start, stop = obs_max_min(args.obsid, meta=data_dict)

        print("-------------------------    Obs Info    --------------------------")
        print("Obs Name:           {}".format(data_dict["obsname"]))
        print("Creator:            {}".format(data_dict["rfstreams"]["0"]["creator"]))
        print("Array phase:        {}".format(array_phase))
        if array_phase != 'OTH':
            print("~FWHM (arcminute)   {:4.2f}".format(calc_ta_fwhm(centre_freq,
                                                       array_phase=array_phase)*60.))
        print("Start time:         {}".format(start))
        print("Stop time:          {}".format(stop))
        print("Duration (s):       {}".format(stop-start))
        print("RA Pointing (deg):  {}".format(data_dict["metadata"]["ra_pointing"]))
        print("DEC Pointing (deg): {}".format(data_dict["metadata"]["dec_pointing"]))
        print("Channels:           {}".format(data_dict["rfstreams"]["0"]["frequencies"]))
        print("Centrefreq (MHz):   {}".format(centre_freq))

    if args.cal_best:
        all_cals = get_best_cal_obs(args.obsid)
        print()
        print("{:14}|{:8}|{}".format("Calibration ID", "Hrs away", "Cal Target"))
        for cal in all_cals:
            print("{:14}|{:8.2f}|{}".format(cal[0], cal[1]/60., cal[2]))