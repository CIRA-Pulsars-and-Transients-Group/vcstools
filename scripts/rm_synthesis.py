#! /usr/bin/env python3
import logging
import argparse

from vcstools.rm_synth_utils import rm_synth_pipe
from vcstools.general_utils import setup_logger

logger = logging.getLogger(__name__)


def rm_synth_main(kwargs):
    rm_dict, _ = rm_synth_pipe(kwargs)
    if rm_dict:
        for i in rm_dict.keys():
            logger.info("For phase range: {0} - {1}".format(*rm_dict[i]["phase_range"]))
            logger.info("RM: {0:7.3f} +/- {1:6.3f}".format(rm_dict[i]["rm"], rm_dict[i]["rm_e"]))


if __name__ == '__main__':

    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR = logging.ERROR)

    parser = argparse.ArgumentParser(description="A script that performs the RM synthesis technique",\
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group("Required Inputs:")
    required.add_argument("-a", "--archive", required=True, type=str, help="The name of the archvie file to work with.")

    fitting = parser.add_argument_group("Fitting Options:")
    fitting.add_argument("--phase_ranges", type=float, nargs="+", help="The phase range(s) to fit the RM to. If unsupplied, will find the on-pulse and fit that range.\
                         Supports multiple ranges. eg. 0.1 0.15 0.55 0.62 will fit from 0.1 to 0.15 and from 0.55 or 0.62.")
    fitting.add_argument("--phi_steps", type=int, default=10000, help="The number of rm steps to use for synthesis.")
    fitting.add_argument("--phi_range", type=float, default=(-300, 300), nargs="+", help="The range of RMs so synthsize. Giving a smaller window will speed up operations.")
    fitting.add_argument("--force_single", action="store_true", help="use this tag to force using only a single phase range (if phase_ranges is unsupplied)")

    output = parser.add_argument_group("Output Options:")
    output.add_argument("--label", type=str, help="A label for the output.")
    output.add_argument("--write", action="store_true", help="Use this tag to write the results to a labelled file")
    output.add_argument("--plot_rm", action="store_true", help="Use this tag to plot the resulting RM fit.")
    output.add_argument("--plot_gfit", action="store_true", help="Use this tag to plot the gaussian fit (if phase ranges not supplied).")
    output.add_argument("--keep_QUV", action="store_true", help="Use this tag to keep the QUVflux.out file from rmfit.")

    gfit = parser.add_argument_group("Gaussian Fit Options")
    gfit.add_argument("--cliptype", type=str, default="regular", help="Verbosity of clipping for gaussian fitting. Options - regular, noisy, verbose")

    optional = parser.add_argument_group("additional Inputs:")
    optional.add_argument("-d", "--work_dir", type=str, default="./", help="The directory to work in")
    optional.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level. Default: INFO", choices=loglevels.keys())

    args=parser.parse_args()

    # set up the logger for stand-alone execution
    logger = setup_logger(logger, log_level=loglevels[args.loglvl])

    kwargs = vars(args)
    rm_synth_main(kwargs)