from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import units as u
from astropy.time import Time
import numpy as np

def sex2deg(ra, dec):
    """Convert sexagesimal coordinates to degrees.

    Parameters
    ----------
    ra : `str`
        The Right Acension in the format "HH:MM:SS.SS".
    dec : `str`
        The Declination in the format "DD:MM:SS.SS".

    Returns
    -------
    ra : `float`
        The Right Acension in degrees.
    dec : `float`
        The Declination in degrees.
    """
    c = SkyCoord( ra, dec, frame='icrs', unit=(u.hourangle,u.deg))

    # return RA and DEC in degrees in degrees
    return [c.ra.deg, c.dec.deg]


def deg2sex(ra, dec):
    """Convert decimal coordingates into sexagesimal strings.

    Parameters
    ----------
    ra : `float`
        The Right Acension in degrees.
    dec : `float`
        The Declination in degrees.

    Returns
    -------
    ra : `str`
        The Right Acension in the format "HH:MM:SS.SS".
    dec : `str`
        The Declination in the format "DD:MM:SS.SS".
    """
    c = SkyCoord( ra, dec, frame='icrs', unit=(u.deg,u.deg))
    #coords = c.to_string('hmsdms')
    #coords = coords.replace('h',':').replace('d',':').replace('m',':').replace('s','')
    rajs = c.ra.to_string(unit=u.hour, sep=':')
    decjs = c.dec.to_string(unit=u.degree, sep=':')

    # return RA and DEC in "hh:mm:ss.ssss dd:mm:ss.ssss" form
    return rajs, decjs


def format_ra_dec(ra_dec_list, ra_col=0, dec_col=1):
    """Will format a list of lists containing RAs and Decs to uniform strings. eg 00:00:00.00 -00:00:00.00.
    Will not work for numpy arrays so make sure they're list of lists
    An example input:
    format_ra_dec([[name,ra,dec]], ra_col = 1, dec_col = 2)

    Parameters
    ----------
    ra_dec_list : `list` of `lists`
        A list of lists where each row is a different source with an RA and declination.
    ra_col : `int`, optional
        The coloumn ID that contains the RA. |br| Default: 0.
    dec_col : `int`, optional
        The coloumn ID that contains the Declination. |br| Default: 1.

    Returns
    -------
    ra_dec_list : `list` of `lists`
        The formated ra_dec_list.
    
    Examples
    --------
    >>> format_ra_dec([["J2351+8533", "23:51:03", "+85:33:20.6"]], ra_col = 1, dec_col = 2)
    [["J2351+8533", "23:51:03.00", "+85:33:20.60"]]
    """
    for i in range(len(ra_dec_list)):
        #catching errors in old psrcat RAs and Decs
        if len(ra_dec_list[i][ra_col]) >5:
            if  ra_dec_list[i][ra_col][5] == '.':
                ra_dec_list[i][ra_col] = ra_dec_list[i][ra_col][:5] + ":" +\
                        str(int(float('0'+ra_dec_list[i][ra_col][5:])*60.))
        if ra_dec_list[i][dec_col].count(":") == 1 and ra_dec_list[i][dec_col].count(".") == 1:
            # Only arcminutes so split it into arcseconds
            ra_dec_list[i][dec_col] = ra_dec_list[i][dec_col].split(".")[0] + ":" +\
                        str(float("0." + str(ra_dec_list[i][dec_col].split(".")[-1]))*60 )[:5]
        if len(ra_dec_list[i][dec_col].split(":")[-1]) == 1:
            # Some final digits do not have leading zeros
            ra_dec_list[i][dec_col] = ra_dec_list[i][dec_col].rsplit(':', 1)[0] + ":0" +\
                        ra_dec_list[i][dec_col].split(":")[-1]
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


def getTargetAZZA(ra, dec, time, lat=-26.7033, lon=116.671, height=377.827, units=(u.hourangle, u.deg)):
    """Function to get the target position in alt/az at a given EarthLocation and Time. The default lat,lon,height is the centre of MWA.

    Parameters
    ----------
    ra : `str`
        The target right ascension in astropy-readable format.
    dec : `str`
        The target declination in astropy-readable format.
    time : `str`
        The time of observation in UTC (i.e. a string on form: yyyy-mm-dd hh:mm:ss.ssss)
    lat : `float`, optional
        The observatory latitude in degrees. |br| Default: -26.7033.
    lon : `float`, optional
        The observatory longitude in degrees. |br| Default: 116.671.
    height : `float`, optional
        The observatory height from sea level in meters. |br| Default: 377.827.
    units : `tuple`, optional
        The astropy units of the ra and dec. |br| Default: (u.hourangle, u.deg)

    Returns
    -------
    az : `float`
        Target azimuth in radians.
    za : `float`
        Target zenith angle in radians.
    azdeg : `float`
        Target azimuth in degrees.
    zadeg : `float`
        Target zenith angle in degrees.
    """
    # Create Earth location for MWA
    location = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=height*u.m)

    # Create sky coordinates for target
    coord = SkyCoord(ra, dec, unit=units)

    # Convert from RA/Dec to Alt/Az
    altaz = coord.transform_to(AltAz(obstime=Time(time), location=location))

    az = altaz.az.rad
    azdeg = altaz.az.deg

    za = np.pi / 2 - altaz.alt.rad
    zadeg = 90. - altaz.alt.deg

    return az, za, azdeg, zadeg


def getTargetRADec(az, za, time, lat=-26.7033, lon=116.671, height=377.827, units=(u.deg, u.deg)):
    """Function to get the target position in ra dec at a given EarthLocation and Time. The default lat,lon,height is the centre of MWA.

    Parameters
    ----------
    az : `str`
        The target azimuth in astropy-readable format.
    za : `str`
        The target zenith angle in astropy-readable format.
    time : `str`
        The time of observation in UTC (i.e. a string on form: yyyy-mm-dd hh:mm:ss.ssss)
    lat : `float`, optional
        The observatory latitude in degrees. |br| Default: -26.7033.
    lon : `float`, optional
        The observatory longitude in degrees. |br| Default: 116.671.
    height : `float`, optional
        The observatory height from sea level in meters. |br| Default: 377.827.
    units : `tuple`, optional
        The astropy units of the ra and dec. |br| Default: (u.deg, u.deg)

    Returns
    -------
    ra : `float`
        Target right ascension in radians.
    dec : `float`
        Target declination in radians.
    radeg : `float`
        Target right ascension in degrees.
    decdeg : `float`
        Target declination in degrees.
    """
    # Create Earth location for MWA
    location = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=height*u.m)

    # Create sky coordinates for target
    coord = SkyCoord(az=az, alt=(90.*u.deg - za*units[0]), unit=units, frame='altaz', obstime=Time(time), location=location)
    radec = coord.icrs

    return radec.ra.rad, radec.dec.rad, radec.ra.deg, radec.dec.deg
