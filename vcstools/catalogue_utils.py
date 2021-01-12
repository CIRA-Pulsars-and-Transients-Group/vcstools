import os
import math
import csv

from vcstools import data_load
from vcstools.pointing_utils import deg2sex

import logging
logger = logging.getLogger(__name__)


def get_psrcat_ra_dec(pulsar_list=None, max_dm=250., include_dm=False, query=None):
    """
    Uses PSRCAT to return a list of pulsar names, ras and decs. Not corrected for proper motion.
    Removes pulsars without any RA or DEC recorded
    If no pulsar_list given then returns all pulsar on the catalogue
    If include_dm is True then also ouput DM

    get_psrcat_ra_dec(pulsar_list = None)
    Args:
        pulsar_list: A space list of pulsar names eg: [J0534+2200, J0538+2817].
               (default: uses all pulsars)
    return [[Jname, RAJ, DecJ]]
    """
    import psrqpy

    #params = ['JNAME', 'RAJ', 'DECJ', 'DM']
    if query is None:
        query = psrqpy.QueryATNF(params = ['PSRJ', 'RAJ', 'DECJ', 'DM'], psrs=pulsar_list, loadfromdb=data_load.ATNF_LOC).pandas

    pulsar_ra_dec = []
    for i, _ in enumerate(query["PSRJ"]):
        # Only record if under the max_dm
        dm = query["DM"][i]
        if not math.isnan(dm):
            if float(dm) < max_dm:
                if include_dm:
                    pulsar_ra_dec.append([query["PSRJ"][i], query["RAJ"][i], query["DECJ"][i], dm])
                else:
                    pulsar_ra_dec.append([query["PSRJ"][i], query["RAJ"][i], query["DECJ"][i]])


    return pulsar_ra_dec


def grab_source_alog(source_type='Pulsar', pulsar_list=None, max_dm=1000., include_dm=False, query=None):
    """
    Will search different source catalogues and extract all the source names, RAs, Decs and, if requested, DMs

    Parameters
    ----------
    source_type: string
        The type of source you would like to get the catalogue for.
        Your choices are: ['Pulsar', 'FRB', 'rFRB', 'POI' 'RRATs', 'Fermi']
        Default: 'Pulsar'
    pulsar_list: list
        List of sources you would like to extract data for.
        If None is given then it will search for all available sources
        Default: None
    max_dm: float
        If the source_type is 'Pulsar' then you can set a maximum dm and the function will only
        return pulsars under that value.
        Default: 1000.
    include_dm: Bool
        If True the function will also return the DM if it is available at the end of the list
        Default: False

    Returns
    -------
    result: list list
        A list for each source which contains a [source_name, RA, Dec (, DM)]
        where RA and Dec are in the format HH:MM:SS
    """
    modes = ['Pulsar', 'FRB', 'rFRB', 'POI', 'RRATs', 'Fermi']
    if source_type not in modes:
        logger.error("Input source type not in known catalogues types. Please choose from: {0}".format(modes))
        return None

    #Get each source type into the format [[name, ra, dec]]
    name_ra_dec = []
    if source_type == 'Pulsar':
        name_ra_dec = get_psrcat_ra_dec(pulsar_list=pulsar_list, max_dm=max_dm, include_dm=include_dm, query=query)

    elif source_type == 'FRB':
        import urllib.request
        from io import StringIO
        try:
            frb_csv = urllib.request.urlopen("https://www.wis-tns.org/search?&include_frb=1&objtype[]=130&num_page=500&format=csv").read().decode('utf-8')
        except urllib.error.HTTPError:
            logger.error('frbcat (https://www.wis-tns.org) not available. Returning empty list')
            # putting and FRB at 90 dec which we should never be able to detect
            name_ra_dec = [['fake', "00:00:00.00", "90:00:00.00", 0.0]]
        else:
            f = StringIO(frb_csv)
            reader = csv.reader(f, delimiter=',')
            for frb in reader:
                name = frb[1]
                logger.debug('FRB name: {}'.format(name))
                ra   = frb[2]
                dec  = frb[3]
                dm   = frb[6]
                if include_dm:
                    name_ra_dec.append([name, ra, dec, dm])
                else:
                    name_ra_dec.append([name, ra, dec])

    elif source_type == "rFRB":
        info = get_rFRB_info(name=pulsar_list)
        if info is not None:
            for line in info:
                if include_dm:
                    name_ra_dec.append([line[0], line[1], line[2], line[3]])
                else:
                    name_ra_dec.append([line[0], line[1], line[2]])

    elif source_type == "POI":
        #POI = points of interest
        db = open(data_load.POI_CSV, "r")
        for line in db.readlines():
            if not line.startswith("#"):
                line = line.split(",")
                name_ra_dec.append([line[0], line[1], line[2]])

    elif source_type == 'RRATs':
        import urllib.request
        try:
            rrats_data = urllib.request.urlopen('http://astro.phys.wvu.edu/rratalog/rratalog.txt').read().decode()
        except urllib.error.URLError:
            logger.error('http://astro.phys.wvu.edu/rratalog/ not available. Returning empty list')
            # putting and RRAT at 90 dec which we should never be able to detect
            name_ra_dec = [['fake', "00:00:00.00", "90:00:00.00", 0.0]]
        else:
            for rrat in rrats_data.split("\n")[1:-1]:
                columns = rrat.strip().replace(" ", '\t').split('\t')
                rrat_cat_line = []
                if pulsar_list == None or (columns[0] in pulsar_list):
                    for entry in columns:
                        if entry not in ['', ' ', '\t']:
                            rrat_cat_line.append(entry.replace('--',''))
                    #Removing bad formating for the database
                    ra = rrat_cat_line[4]
                    if ra.endswith(":"):
                        ra = ra[:-1]
                    dec = rrat_cat_line[5]
                    if dec.endswith(":"):
                        dec = dec[:-1]

                    if include_dm:
                        name_ra_dec.append([rrat_cat_line[0], ra, dec, rrat_cat_line[3]])
                    else:
                        name_ra_dec.append([rrat_cat_line[0], ra, dec])

    elif source_type == 'Fermi':
        # read the fermi targets file
        try:
            fermi_loc = os.environ['FERMI_CAND_FILE']
        except:
            logger.warning("Fermi candidate file location not found. Returning nothing")
            return []
        with open(fermi_loc,"r") as fermi_file:
            csv_reader = csv.DictReader(fermi_file)
            for fermi in csv_reader:
                name = fermi['Source Name'].split()[-1]
                ra  = fermi[' RA J2000']
                dec = fermi[' Dec J2000']
                raj, decj = deg2sex(float(ra), float(dec))
                pos_u = float(fermi[' a (arcmin)'])
                if include_dm:
                    # this actually returns the position uncertainty not dm
                    name_ra_dec.append([name, raj, decj, pos_u])
                else:
                    name_ra_dec.append([name, raj, decj])

    #Remove all unwanted sources
    if pulsar_list is not None:
        rrat_cat_filtered = []
        for line in name_ra_dec:
            if line[0] in pulsar_list:
                rrat_cat_filtered.append(line)
        name_ra_dec = rrat_cat_filtered

    return name_ra_dec


def get_rFRB_info(name=None):
    """
    Gets repeating FRB info from the csv file we maintain.
    Input:
        name: a list of repeating FRB names to get info for. The default is None
              which gets all rFRBs in the catalogue.
    Output:
        [[name, ra, dec, dm, dm error]]
        A list of all the FRBs which each have a list contining name, ra, dec,
        dm and dm error
    """
    output = []
    db = open(data_load.KNOWN_RFRB_CSV, "r")
    for line in db.readlines():
        if not line.startswith("#"):
            line = line.split(",")
            #some FRBs end with a J name. We will ignore these when comparing
            #by using the first 9 characters
            FRB = line[0]
            if name is None:
                #No input FRBs so return all FRBs
                output.append(line)
            elif FRB in name:
                output.append(line)
    return output