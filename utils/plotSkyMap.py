#!/usr/bin/env python 

import os
import sys
import argparse
import urllib
import urllib2
import json
import subprocess
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astropy.time import Time
from mwapy.pb import primary_beam
import ephem
from mwapy import ephem_utils,metadata
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.tri as tri
import matplotlib.cm as cm

from mwa_metadb_utils import getmeta


def get_beam_power(obsid_data,
                   sources,
                   dt=296,
                   centeronly=True,
                   verbose=False,
                   min_power=0.6):
    """
    obsid_data = [obsid,ra, dec, time, delays,centrefreq, channels]
    sources=[names,coord1,coord2] #astropy table coloumn names

    Calulates the power (gain at coordinate/gain at zenith) for each source and if it is above
    the min_power then it outputs it to the text file.

    """
    #print "Calculating beam power"
    obsid,ra, dec, time, delays,centrefreq, channels = obsid_data
    
    #starttimes=np.arange(0,time,dt)
    starttimes=np.arange(0,time,time)
    stoptimes=starttimes+dt
    stoptimes[stoptimes>time]=time
    Ntimes=len(starttimes)
    midtimes=float(obsid)+0.5*(starttimes+stoptimes)

    mwa = ephem_utils.Obs[ephem_utils.obscode['MWA']]
    # determine the LST
    observer = ephem.Observer()
    # make sure no refraction is included
    observer.pressure = 0
    observer.long = mwa.long / ephem_utils.DEG_IN_RADIAN
    observer.lat = mwa.lat / ephem_utils.DEG_IN_RADIAN
    observer.elevation = mwa.elev

    if not centeronly:
        PowersX=np.zeros((len(sources),
                             Ntimes,
                             len(channels)))
        PowersY=np.zeros((len(sources),
                             Ntimes,
                             len(channels)))
        # in Hz
        frequencies=np.array(channels)*1.28e6
    else:
        PowersX=np.zeros((len(sources),
                             Ntimes,1))
        PowersY=np.zeros((len(sources),
                             Ntimes,1))
        frequencies=[centrefreq]

    RAs=np.array([x[0] for x in sources])
    Decs=np.array([x[1] for x in sources])
    if len(RAs)==0:
        sys.stderr.write('Must supply >=1 source positions\n')
        return None
    if not len(RAs)==len(Decs):
        sys.stderr.write('Must supply equal numbers of RAs and Decs\n')
        return None
    for itime in xrange(Ntimes):
        obstime = Time(midtimes[itime],format='gps',scale='utc')
        observer.date = obstime.datetime.strftime('%Y/%m/%d %H:%M:%S')
        LST_hours = observer.sidereal_time() * ephem_utils.HRS_IN_RADIAN

        HAs = -RAs + LST_hours * 15
        Azs, Alts = ephem_utils.eq2horz(HAs, Decs, mwa.lat)
        # go from altitude to zenith angle
        theta=np.radians(90-Alts)
        phi=np.radians(Azs)
        
        for ifreq in xrange(len(frequencies)):
            rX,rY=primary_beam.MWA_Tile_analytic(theta, phi,
                                                 freq=frequencies[ifreq], delays=delays,
                                                 zenithnorm=True,
                                                 power=True)
            PowersX[:,itime,ifreq]=rX
            PowersY[:,itime,ifreq]=rY

    #Power [#sources, #times, #frequencies]
    Powers=0.5*(PowersX+PowersY)
                 
    return Powers

def plotFigure():


    fig, ax = plt.subplots(projection="mollweide")


    Path = mpath.Path
    #MWA could reach
    newra_s = -np.pi
    newra_e = np.pi
    newdec_s = -0.5*np.pi
    newdec_e = (30)/180.0*np.pi
    path_data = [
                (Path.MOVETO, (newra_s, newdec_s)),
                (Path.LINETO, (newra_s, newdec_e)),
                (Path.LINETO, (newra_e, newdec_e)),
                (Path.LINETO, (newra_e, newdec_s)),
                (Path.CLOSEPOLY, (newra_e, newdec_s)),
                ]   
    codes, verts = zip(*path_data)
    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(path, facecolor='g', alpha=0.07)
    #ax.add_patch(patch)

    #only MWA could reach
    newra_s = -np.pi
    newra_e = np.pi
    newdec_s = -0.5*np.pi
    newdec_e = (-50)/180.0*np.pi
    path_data = [
                (Path.MOVETO, (newra_s, newdec_s)),
                (Path.LINETO, (newra_s, newdec_e)),
                (Path.LINETO, (newra_e, newdec_e)),
                (Path.LINETO, (newra_e, newdec_s)),
                (Path.CLOSEPOLY, (newra_e, newdec_s)),
                ]   
    codes, verts = zip(*path_data)
    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(path, linewidth='0.0', facecolor='g', alpha=0.13)
    #ax.add_patch(patch)


    psrcat = Table.read("psrcat.csv", delimiter=";")
    psrcat_coords = SkyCoord(psrcat["RAJ"], psrcat["DECJ"], unit=(u.hourangle, u.deg))

    rrats = Table.read("RRATs.csv", delimiter=",")
    rrat_obsids = np.genfromtxt("RRAT_obsids.txt", dtype=str)
    rrat_coords = SkyCoord(rrats["RAJ"], rrats["DECJ"], unit=(u.hourangle, u.deg))

    #sys.exit(0)


    levels = np.arange(0.25, 1., 0.05)
    #colors=['0.25', '0.5', '0.5', '0.5', '0.5', '0.5', '0.5', '0.5', '0.5', '0.5', '0.5', '0.5', '0.5', '0.5', '0.5']
    cmap = plt.get_cmap("gray_r")
    linewidths=[0.7, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]
    #linewidths=[0.7, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]



    for line in rrat_obsids:
        #words=line.split()
        words=line.split(',')
        ob=words[0]
        if ob[0]!=('1'):
            continue
        print "Obtaining metadata from http://mwa-metadata01.pawsey.org.au/metadata/ for OBS ID: "
        beam_meta_data = getmeta(service='obs', params={'obs_id':ob})
            
        ra = beam_meta_data[u'metadata'][u'ra_pointing']
        dec = beam_meta_data[u'metadata'][u'dec_pointing']
        time = beam_meta_data[u'stoptime'] - beam_meta_data[u'starttime'] #gps time
        delays = beam_meta_data[u'rfstreams'][u'0'][u'xdelays']

        minfreq = float(min(beam_meta_data[u'rfstreams'][u"0"][u'frequencies']))
        maxfreq = float(max(beam_meta_data[u'rfstreams'][u"0"][u'frequencies']))
        centrefreq = 1.28e6 * (minfreq + (maxfreq-minfreq)/2) #in MHz
        channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']

        cord = [ob, ra, dec, time, delays, centrefreq, channels]
        #print ra
        #print dec
        #Dec=range(-90,90,1)
        RA=[] ; Dec=[] ; z=[] ; x=[] ; y=[]
        
        for i in range(-87,88,3):
            for j in range(0,361,3):
                Dec.append(i)
                RA.append(j)


        powout=get_beam_power(cord, zip(RA,Dec), dt=600)
        for i in range(len(RA)):
            temppower=powout[i,0,0]
            for t in range(0,(time+361)/720):
                if i%121 >= t:
                    power_ra = powout[i-t,0,0] #ra+t*15./3600 3deg
                else : 
                    power_ra = powout[i+121-t,0,0]
                if power_ra > temppower:
                    temppower = power_ra
                #print temppower, power_ra
            z.append(temppower)
            if RA[i] > 180:
                x.append(-RA[i]/180.*np.pi+2*np.pi)
            else:
                x.append(-RA[i]/180.*np.pi)

            y.append(Dec[i]/180.*np.pi)
        plt.tricontour(x, y, z, levels=levels,
                       cmap=cmap,
                       linewidths=linewidths, zorder=1.1)

    plt.colorbar(fraction=0.02, pad=0.03, label="Zenith normalised power")

    #p1=ax.scatter(ra_PCAT_N, dec_PCAT_N,0.19, marker='o', color ='gray', label="Known pulsars beyond MWA Dec limitation")
    #p1=ax.scatter(ra_PCAT, dec_PCAT,0.2, marker='o', color ='blue', label="Known pulsars MWA could reach")
    maskGood = psrcat_coords.dec.wrap_at(180*u.deg).deg < 30
    maskBad = psrcat_coords.dec.wrap_at(180*u.deg).deg >= 30

    ax.scatter(-psrcat_coords.ra.wrap_at(180*u.deg).rad[maskGood], psrcat_coords.dec.wrap_at(180*u.deg).rad[maskGood], 0.01, marker='x', color ='C0', zorder=1.6)
    ax.scatter(-psrcat_coords.ra.wrap_at(180*u.deg).rad[maskBad], psrcat_coords.dec.wrap_at(180*u.deg).rad[maskBad], 0.01, marker='x', color ='0.5', zorder=1.6)
    #newra_s = -np.pi
    #newra_e = np.pi
    #newdec_s = -0.5*np.pi
    #newdec_e = (-29.9)/180.0*np.pi
    ax.scatter(-rrat_coords.ra.wrap_at(180*u.deg).rad, rrat_coords.dec.wrap_at(180*u.deg).rad, 10, marker="x", color="red", zorder=1.7)
    #p1 = ax.plot([newra_s,newra_e], [newdec_e,newdec_e],color ='cyan', label="LOFAR Dec limitation")
    #p2=ax.scatter(ra_BOTH, dec_BOTH,9, marker='*', color ='tomato',zorder=52, label="Both Tara and us detected")
    #p2=ax.scatter(ra_MWA, dec_MWA,9, marker='*', color ='crimson',zorder=51, label="MWA VCS pulsars")#We detected while Tara didn't
    #p1=ax.scatter(ra_MENO, dec_MENO,7, marker='v', color ='darkcyan',zorder=50, label="Tara detected, we fold but not got it")
    #p1=ax.scatter(ra_TARA, dec_TARA,7, marker='v', color ='green',zorder=49, label="Tara detected, we haven't fold yet")
    """
    p1=ax.plot(ra_TARA, dec_TARA,ms=2, marker='v', linestyle='none', color ='green',zorder=49, label="Tara detected, we haven't fold yet")
    #p1=ax.plot(ra_PCAT, dec_PCAT,ms=0.5, marker='o', linestyle='none', color ='blue', label="2536 Known pulsars")
    p2=ax.plot(ra_MWA, dec_MWA,ms=3, marker='*', linestyle='none', color ='crimson',zorder=51, label="We detected while Tara didn't")
    p2=ax.plot(ra_BOTH, dec_BOTH,ms=3, marker='*', linestyle='none', color ='gold', zorder=52,label="Both Tara and us detected")
    p1=ax.plot(ra_MENO, dec_MENO,ms=2, marker='v', linestyle='none', color ='cyan',zorder=50, label="Tara detected, we fold but not got it")
    """
    #print 'PSR both',len(ra_BOTH)

    xtick_labels = ['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h']
    ax.set_xticklabels(xtick_labels)
    ax.set_xlabel("Right Ascension")
    ax.set_ylabel("Declination")
    ax.grid(True)


    #handles, labels = ax.get_legend_handles_labels()
    plt.legend(bbox_to_anchor=(0.65, 1.02,0.4,0.2), loc=3,numpoints=1,
               ncol=1, mode="expand", borderaxespad=0., fontsize=6)

    #plt.savefig('all_contour_drift_TARA.png', figsize=(5,3), dpi=600)
    #plt.savefig('all_archive_drift_TARA_a.png', figsize=(5,3), dpi=600)
    #plt.savefig('all_archive_drift_TARA_b.png', figsize=(5,3), dpi=600)
    #plt.savefig('all_contour_archive.png', figsize=(5,3), dpi=600)
    #plt.savefig('all_contour_testtt.png', figsize=(5,3), dpi=600)
    #plt.savefig('rrats.eps', format="eps", bbox_inches="tight")
    plt.savefig('rrats.png', format="png", dpi=300, bbox_inches="tight")
    #plt.savefig('all_contour_combine_limit.png', figsize=(5,3), dpi=600)
    #plt.savefig('all_archive_drift_TARA_a.ps', figsize=(5,3))
    #plt.savefig('all_archive_drift_TARA_b.ps', figsize=(5,3))
    #plt.savefig('all_contour_archive.png.ps', figsize=(5,3))
    #plt.savefig('all_contour_combine_limit.ps', figsize=(5,3))
    #plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--obsfile", type=str, help="File containing a list of observation IDs, one ID per line")
    
    parser.add_argument("-t", "--targetfile", type=str, help="File containing a list of target sources, with positions in columns labelled 'RAJ' and 'DECJ'")


    parser.add_argument("--show_psrcat", action="store_true", help="Whether to show background pulsars on map")
    parser.add_argument("--show_mwa_sky", action="store_true", help="Whether to split the sky based on MWA visibility")

    args = parser.parse_args()
