#!/usr/bin/env python3
import logging
import argparse
import sys

from vcstools.sn_flux_utils import est_pulsar_sn, ATNF_spectral_data_plot
from vcstools.general_utils import setup_logger

logger = logging.getLogger(__name__)


def snfe_main(kwargs):
    if kwargs["mode"] == "SNFE":
        raj = None
        decj = None
        if kwargs["pointing"]:
            raj = kwargs["pointing"].split("_")[0]
            decj = kwargs["pointing"].split("_")[1]
        SN, SN_e, flux, flux_e = est_pulsar_sn(kwargs["pulsar"], kwargs["obsid"], kwargs["beg"], kwargs["end"],
                                plot_flux=kwargs["plot_est"], p_ra=raj, p_dec=decj)
        logger.info("{0} Flux: {1} +/- {2} Jy".format(kwargs["pulsar"], flux, flux_e))
        logger.info("{0} S/N: {1} +/- {2}".format(kwargs["pulsar"], SN, SN_e))
    elif kwargs["mode"] == "ATNF":
        ATNF_spectral_data_plot(kwargs["pulsar"])
    else:
        logger.error("Valid mode not selected. Please refer to documentation for options")
        sys.exit(1)


if __name__ == "__main__":

    loglevels = dict(DEBUG=logging.DEBUG,\
                    INFO=logging.INFO,\
                    WARNING=logging.WARNING,\
                    ERROR=logging.ERROR)

    parser = argparse.ArgumentParser(description="""A utility file for estimating the S/N of a pulsar in an obsid""")

    parser.add_argument("-o", "--obsid", type=int, required=True, help="The Observation ID (e.g. 1221399680)")
    parser.add_argument("-p", "--pulsar", type=str, required=True, help="The pulsar's name (e.g. J2241-5236).\
                        Takes a single argument for flux estimation or multiple for ATNF plotting")
    parser.add_argument("-b", "--beg", type=int, required=True, help="The beginning time of observation.\
                        If None, will use beginning given by a metadata call. Default: None")
    parser.add_argument("-e", "--end", type=int, required=True, help="The end time of observation.\
                        If None, will use the end given by a metadata call. Default: None")
    parser.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbostity level. Default: INFO")
    parser.add_argument("--pointing", type=str, default=None, help="The pointing of the target in the format '12:34:56_98:76:54'.\
                        If None, will obtain from a call to ATNF. Default: None")
    parser.add_argument("--min_z_power", type=float, default=0.3, help="The minimum zenith normalised power used to determine if the pulsar\
                        is in the beam or not")
    parser.add_argument("--plot_est", action="store_true", help="Use this tag to create a plot of flux estimation.")
    parser.add_argument("--mode", type=str, help="""MODES: 'SNFE' = Estimate S/N and flux for a single pulsar in an obsid\n
                                                'ATNF' = Plot the spectral energy distribution for any number of pulsars using data from ATNF.""")
    args = parser.parse_args()

    # set up the logger for stand-alone execution
    logger = setup_logger(logger, log_level=loglevels[args.loglvl])

    kwargs=vars(args)
    snfe_main(kwargs)