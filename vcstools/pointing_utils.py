from astropy.coordinates import SkyCoord
from astropy import units as u

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
