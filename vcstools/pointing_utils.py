from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import units as u
from astropy.time import Time
import numpy as np

def sex2deg(ra, dec):
    """
    Convert sexagesimal coordinates to degrees.

    sex2deg( ra, dec)
    Args:
        ra: the right ascension in HH:MM:SS
        dec: the declination in DD:MM:SS
    """
    c = SkyCoord( ra, dec, frame='icrs', unit=(u.hourangle,u.deg))

    # return RA and DEC in degrees in degrees
    return [c.ra.deg, c.dec.deg]


def deg2sex(ra, dec):
    """
    Convert decimal coordingates into sexagesimal strings, i.e. hh:mm:ss.ss and dd:mm:ss.ss

    deg2sex( ra, dec)
    Args:
        ra: the right ascension in degrees
        dec: the declination in degrees
    """

    c = SkyCoord( ra, dec, frame='icrs', unit=(u.deg,u.deg))
    #coords = c.to_string('hmsdms')
    #coords = coords.replace('h',':').replace('d',':').replace('m',':').replace('s','')
    rajs = c.ra.to_string(unit=u.hour, sep=':')
    decjs = c.dec.to_string(unit=u.degree, sep=':')

    # return RA and DEC in "hh:mm:ss.ssss dd:mm:ss.ssss" form
    return rajs, decjs


def format_ra_dec(ra_dec_list, ra_col=0, dec_col=1):
    """
    Will format a list of lists containing RAs and Decs to uniform strings.  eg 00:00:00.00 -00:00:00.00. Will not work for numpy arrays so make sure they're list of lists
    An example input:
    format_ra_dec([[name,ra,dec]], ra_col = 1, dec_col = 2)
    """
    for i in range(len(ra_dec_list)):
        #catching errors in old psrcat RAs
        if len(ra_dec_list[i][ra_col]) >5:
            if  ra_dec_list[i][ra_col][5] == '.':
                ra_dec_list[i][ra_col] = ra_dec_list[i][ra_col][:5] + ":" +\
                        str(int(float('0'+ra_dec_list[i][ra_col][5:])*60.))
        #make sure there are two digits in the HH slot for RA or DD for Dec
        if len(ra_dec_list[i][ra_col].split(':')[0]) == 1:
            ra_dec_list[i][ra_col] = '0' + ra_dec_list[i][ra_col]
        if ra_dec_list[i][dec_col][0].isdigit():
            if len(ra_dec_list[i][dec_col].split(':')[0]) == 1:
                ra_dec_list[i][dec_col] = '0' + ra_dec_list[i][dec_col]
        else:
            if len(ra_dec_list[i][dec_col].split(':')[0]) == 2:
                ra_dec_list[i][dec_col] = ra_dec_list[i][dec_col][0] + '0' +\
                                          ra_dec_list[i][dec_col][1:]

        #make sure there is a + on positive Decs
        if not ra_dec_list[i][dec_col].startswith('-') and\
           not ra_dec_list[i][dec_col].startswith('+'):
               ra_dec_list[i][dec_col] = '+' + ra_dec_list[i][dec_col]

        #since the ra a dec are the same except the +- in the dec just loop over it
        #and change the length
        ra_dec_col = [ra_col, dec_col]
        for n in range(2):
            if len(ra_dec_list[i][ra_dec_col[n]]) == int(2 + n):
                ra_dec_list[i][ra_dec_col[n]] += ':00:00.00'
            elif len(ra_dec_list[i][ra_dec_col[n]]) == 3 + n:
                ra_dec_list[i][ra_dec_col[n]] += '00:00.00'
            elif len(ra_dec_list[i][ra_dec_col[n]]) == 4 + n:
                ra_dec_list[i][ra_dec_col[n]] += '0:00.00'
            elif len(ra_dec_list[i][ra_dec_col[n]]) == 5 + n:
                ra_dec_list[i][ra_dec_col[n]] += ':00.00'
            elif len(ra_dec_list[i][ra_dec_col[n]]) == 6 + n:
                ra_dec_list[i][ra_dec_col[n]] += '00.00'
            elif len(ra_dec_list[i][ra_dec_col[n]]) == 7 + n:
                ra_dec_list[i][ra_dec_col[n]] += '0.00'
            elif len(ra_dec_list[i][ra_dec_col[n]]) == 8 + n:
                ra_dec_list[i][ra_dec_col[n]] += '.00'
            elif len(ra_dec_list[i][ra_dec_col[n]]) == 9 + n:
                ra_dec_list[i][ra_dec_col[n]] += '00'
            elif len(ra_dec_list[i][ra_dec_col[n]]) == 10 + n:
                ra_dec_list[i][ra_dec_col[n]] += '0'

            #shorten if too long
            if len(ra_dec_list[i][ra_dec_col[n]]) > (11 + n):
                ra_dec_list[i][ra_dec_col[n]] = ra_dec_list[i][ra_dec_col[n]][:(11+n)]

    return ra_dec_list


def getTargetAZZA(ra, dec, time, lat=-26.7033, lon=116.671, height=377.827):
    """
    Function to get the target position in alt/az at a given EarthLocation and Time.

    Default lat,lon,height is the centre of  MWA.

    Input:
      ra - target right ascension in astropy-readable format
      dec - target declination in astropy-readable format
      time - time of observation in UTC (i.e. a string on form: yyyy-mm-dd hh:mm:ss.ssss)
      lat - observatory latitude in degrees
      lon - observatory longitude in degrees

    Returns:
      a list containing four elements in the following order:
        list[0] = target azimuth in radians
        list[1] = target zenith angle in radians
        list[2] = target azimuth in degrees
        list[3] = target zenith angle in degrees
    """

    # Create Earth location for MWA
    location = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=height*u.m)

    # Create sky coordinates for target
    coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

    # Create a time object for desired observing time
    obstime = Time(time)

    # Convert from RA/Dec to Alt/Az
    altaz = coord.transform_to(AltAz(obstime=obstime, location=location))

    az = altaz.az.rad
    azdeg = altaz.az.deg

    za = np.pi / 2 - altaz.alt.rad
    zadeg = 90 - altaz.alt.deg

    return az, za, azdeg, zadeg