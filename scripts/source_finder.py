#!/usr/bin/env python

"""
Code to find sources in MWA VCS observations.
"""

__author__ = 'C. P. Lee'
__date__ = '2023-10-11'
__version__ = '0.1'

import sys
import logging
import argparse

import numpy as np
import psrqpy


class InvalidSourceError(Exception):
    """Raise when a source is not valid for any reason."""
    pass


class UnknownPulsarError(Exception):
    """Raise when a pulsar is not found in a catalogue."""
    pass

class FileReadError(Exception):
    """Raise when a file cannot be read."""
    pass


def parse_source_list(sources):
    pointings = []
    query_flag = False
    for source in sources:
        if source.startswith(('J', 'B')):
            query_flag = True
            break
    if query_flag:
        logging.debug('Querying the pulsar catalogue...')
        query = psrqpy.QueryATNF(params=['PSRJ', 'PSRB', 'RAJ', 'DECJ'])
        logging.info(f'Using ATNF pulsar catalogue version {query.get_version}')
    for source in sources:
        if source.startswith(('J', 'B')):
            logging.debug(f'Treating "{source}" as a pulsar')
            raj, decj = get_pulsar_coords(source, query)
        elif '_' in source:
            logging.debug(f'Treating "{source}" as coordinates')
            raj, decj = interpret_coords(source)
        else:
            raise InvalidSourceError(f'Source not recognised: {source}')
        pointing = dict(Name=source, RAJ=raj, DECJ=decj)
        pointings.append(pointing)
    return pointings


def read_sources_file(sources_file):
    try:
        sources = np.loadtxt(sources_file, dtype=str, unpack=True)
    except FileReadError:
        print('The file provided cannot be interpreted.')
    return list(sources)


def get_pulsar_coords(pulsar, query):
    if pulsar.startswith('B'):
        try:
            pid = list(query['PSRB']).index(pulsar)
        except UnknownPulsarError:
            print(f'Pulsar not found in catalogue: {pulsar}')
        pulsar = query['PSRJ'][pid]
    if pulsar not in list(query['PSRJ']):
        raise UnknownPulsarError(f'Pulsar not found in catalogue: {pulsar}')
    psrs = query.get_pulsars()
    raj = psrs[pulsar].RAJ
    decj = psrs[pulsar].DECJ
    return raj, decj


def get_all_pulsars():
    logging.debug('Querying the pulsar catalogue...')
    query = psrqpy.QueryATNF(params=['PSRJ', 'RAJ', 'DECJ'])
    logging.info(f'Using ATNF pulsar catalogue version {query.get_version}')
    psrjs = list(query.table['PSRJ'])
    decjs = list(query.table['DECJ'])
    rajs = list(query.table['RAJ'])
    pointings = []
    for psrj, raj, decj in zip(psrjs, rajs, decjs):
        raj = format_sexigesimal(raj)
        decj = format_sexigesimal(decj, add_sign=True)
        pointing = dict(Name=psrj, RAJ=raj, DECJ=decj)
        pointings.append(pointing)
    return pointings


def format_sexigesimal(coord, add_sign=False):
    if coord.startswith(('-', '–')):
        sign = '-'
        coord = coord[1:]
    elif coord.startswith('+'):
        sign = '+'
        coord = coord[1:]
    else:
        sign = '+'
    parts = coord.split(':')
    if len(parts) == 3:
        deg = int(parts[0])
        minute = int(parts[1])
        second = float(parts[2])
    elif len(parts) == 2:
        deg = int(parts[0])
        minute = int(float(parts[1]))
        second = 0.
    elif len(parts) == 1:
        deg = int(float(parts[0]))
        minute = 0
        second = 0.
    else:
        raise ValueError(f'Cannot interpret coordinate: {coord}.')
    if add_sign:
        formatted_coord = f'{sign}{int(deg):02d}:{int(minute):02d}:{second:05.2f}'
    else:
        formatted_coord = f'{int(deg):02d}:{int(minute):02d}:{second:05.2f}'
    return formatted_coord


def interpret_coords(coords):
    # Split up the coordinates
    raj = coords.split('_')[0]
    decj = coords.split('_')[1]
    # Check how the coordinates are formatted
    if ':' in raj or ':' in decj:
        if ':' not in raj or ':' not in decj:
            raise InvalidSourceError(f'Coordinates not valid: {coords}')
        decimal_flag = False
        if decj[0].isdigit():
            decj = f'+{decj}'
    elif is_float(raj) and is_float(decj):
        decimal_flag = True
    else:
        raise InvalidSourceError(f'Coordinates not valid: {coords}')
    # Convert coordinates to sexigesimal if required
    if decimal_flag:
        raj, decj = decimal_to_sexigesimal(raj, decj)
    # Check that the coordinates are valid
    if not is_sexigesimal(raj, 'RA') or not is_sexigesimal(decj, 'DEC'):
        raise InvalidSourceError(f'Coordinates not valid: {coords}')
    return raj, decj


def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def is_int(string):
    if string.isnumeric():
        return True
    else:
        return False


def is_sexigesimal(coord, mode):
    logging.debug(f'Checking if {coord} is valid sexigesimal')
    if mode == 'RA':
        hours, minutes, seconds = coord.split(':')
        logging.debug(f'{hours=}')
        logging.debug(f'{minutes=}')
        logging.debug(f'{seconds=}')
        if not is_int(hours):
            return False
        if int(hours) < 0 or int(hours) > 24:
            return False
    elif mode == 'DEC':
        degrees, minutes, seconds = coord.split(':')
        logging.debug(f'{degrees=}')
        logging.debug(f'{minutes=}')
        logging.debug(f'{seconds=}')
        if degrees.startswith(('-', '–', '+')):
            degrees = degrees[1:]
        if not is_int(degrees):
            return False
        if int(degrees) < 0 or int(degrees) > 90:
            return False
    if not is_int(minutes):
        return False
    if int(minutes) < 0 or int(minutes) > 60:
        return False
    if not is_float(seconds):
        return False
    if float(seconds) < 0. or float(seconds) > 60.:
        return False
    return True


def decimal_to_sexigesimal(ra_decimal, dec_decimal):
    # Convert right ascension
    ra_decimal = float(ra_decimal)
    if ra_decimal < 0 or ra_decimal >= 360:
        raise InvalidSourceError(f'Right ascension out of bounds: {ra_decimal}')
    ra_hrs_round, ra_hrs_remainder = divmod(ra_decimal, 15)
    ra_minutes, ra_seconds = divmod(ra_hrs_remainder*240, 60)
    ra_sexigesimal = f'{int(ra_hrs_round):02d}:{int(ra_minutes):02d}:{ra_seconds:05.2f}'
    # Convert declination
    dec_decimal = float(dec_decimal)
    if dec_decimal < 0:
        sign = '-'
    else:
        sign = '+'
    if abs(dec_decimal) > 90:
        raise InvalidSourceError(f'Declination is out of bounds: {dec_deci_deg}')
    dec_deg_round, dec_deg_remainder = divmod(abs(dec_decimal), 1)
    dec_minutes, dec_seconds = divmod(dec_deg_remainder*3600, 60)
    dec_sexigesimal = f'{sign}{int(dec_deg_round):02d}:{int(dec_minutes):02d}:{dec_seconds:05.2f}'
    return ra_sexigesimal, dec_sexigesimal


def get_input_args():
    parser = argparse.ArgumentParser(
        usage='%(prog)s [options]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Code to find sources in MWA VCS observations.',
        add_help=False
    )
    loglevels = dict(
        DEBUG=logging.DEBUG,
        INFO=logging.INFO,
        WARNING=logging.WARNING
    )
    # Optional arguments
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-h', '--help', action='help',
        help='Show this help information and exit.')
    optional.add_argument('-V', '--version', action='version',
        version='%(prog)s {}'.format(__version__), help='Print version and exit.')
    optional.add_argument('-L', '--loglvl', type=str, choices=loglevels,
        default='INFO', help='Logger verbosity level.')
    # Source arguments
    source_args = parser.add_argument_group('Source arguments',
        'Options to specify the source(s) to find. The default is all pulsars.')
    source_args.add_argument('-s', '--sources', type=str, nargs='*', default = None,
        help='A list of sources to find. Sources can be specified as pulsar names ' + \
        'or equatorial coordinates separated by an underscore. Coordinates can be ' + \
        'in either decimal or sexigesimal format. For example, the following ' + \
        'arguments are all valid: "B1451-68" "J1456-6843" "14:55:59.92_-68:43:39.50".')
    source_args.add_argument('-f', '--sources_file', type=str, default = None,
        help='A file containing a list of sources to find. Each source should be ' + \
        'listed on a new line. The source format is the same as the -s option.')
    # Parse arguments
    args = parser.parse_args()
    # Set the logging level
    logging.basicConfig(level=loglevels[args.loglvl])
    return args


def main():
    args = get_input_args()

    if args.sources or args.sources_file:
        logging.info('Parsing the provided source list...')
        sources = args.sources
        if args.sources_file:
            sources_from_file = read_sources_file(args.sources_file)
            if sources:
                sources += sources_from_file
            else:
                sources = sources_from_file
        pointings = parse_source_list(sources)
        for pointing in pointings:
            logging.info(f"Source: {pointing['Name']:30} RAJ: {pointing['RAJ']:14} DECJ: {pointing['DECJ']:15}")
    else:
        pointings = get_all_pulsars()
        logging.debug(f'{len(pointings)} pulsars parsed from the catalogue')
    

if __name__ == "__main__":
    main()