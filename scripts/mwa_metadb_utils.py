#!/usr/bin/env python3
import logging
import argparse

from vcstools.metadb_utils import write_obs_info, get_obs_array_phase, obs_max_min, calc_ta_fwhm,\
                                  get_best_cal_obs, get_common_obs_metadata, files_available
from vcstools.beam_calc import field_of_view

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
        beam_common_data, data_dict = get_common_obs_metadata(args.obsid, return_all=True)
        channels = data_dict["rfstreams"]["0"]["frequencies"]
        centre_freq = ( min(channels) + max(channels) ) / 2. * 1.28
        array_phase = get_obs_array_phase(args.obsid)
        start, stop = obs_max_min(args.obsid)
        available_files, all_files = files_available(args.obsid)
        fov = field_of_view(args.obsid, beam_meta_data=beam_common_data)

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
        print("Files available:    {:4.1f}%  ({}/{})".format(len(available_files)/len(all_files)*100,
                                                         len(available_files), len(all_files)))
        print("RA Pointing (deg):  {}".format(data_dict["metadata"]["ra_pointing"]))
        print("DEC Pointing (deg): {}".format(data_dict["metadata"]["dec_pointing"]))
        print("Channels:           {}".format(data_dict["rfstreams"]["0"]["frequencies"]))
        print("Centrefreq (MHz):   {}".format(centre_freq))
        print("FoV (square deg):   {:5.1f}".format(fov))

    if args.cal_best:
        all_cals = get_best_cal_obs(args.obsid)
        print()
        print("{:14}|{:8}|{}".format("Calibration ID", "Hrs away", "Cal Target"))
        for cal in all_cals:
            print("{:14}|{:8.2f}|{}".format(cal[0], cal[1]/60., cal[2]))