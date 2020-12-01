#!/usr/bin/env python
#
# script that reads in a TLE file, estimates the power output by the beamformer, and plots the result.
#
# Requires:
#   - PyEphem
#   - Matplotlib
#

import matplotlib
matplotlib.use('Agg')
import sys, getopt, string, re
from pylab import *

#----------------------------------------------------------------------------------------------------------------------#
# check for command-line options

def help_message():
    print('''
Options: --file=name -- Optional. bandpass file to plot.
Options: --phases    -- Optional. plot phases rather than amplitudes.
Options: --chan      -- Optional. use channel number as x axis.
Options: --raw       -- Optional. use channel number as x axis, with no ordering.
                                  i.e., if above pfb channel 128.
Options: --outname   -- Optional. save plot to a file instead of showing it.
    ''')
    sys.exit(0)

#----------------------------------------------------------------------------------------------------------------------#
# initialise command-line variables.

bp_file = 'BandpassCalibration_node001.dat';
sel_offset = 1
plot_chan = 0
plot_raw  = 0
save_plot = 0

#----------------------------------------------------------------------------------------------------------------------#
# read command-line arguments

try:
    options, xarguments = getopt.getopt(sys.argv[1:], '', ['file=','phases','chan','raw','outname='])
except getopt.error:
    help_message()

for a in options[:]:
    if a[0] == '--file' and a[1] != '':
        bp_file = a[1]
        options.remove(a)
        break
    elif a[0] == '--file' and a[1] == '':
        print('--file expects an argument')
        sys.exit(0)
for a in options[:]:
    if a[0] == '--phases':
        sel_offset = 2
        options.remove(a)
        break
for a in options[:]:
    if a[0] == '--chan':
        plot_chan = 1
        options.remove(a)
        break
for a in options[:]:
    if a[0] == '--raw':
        plot_raw  = 1
        plot_chan = 1
        options.remove(a)
        break
for a in options[:]:
    if a[0] == '--outname' and a[1] != '':
        save_plot = 1
        outname = a[1]
        options.remove(a)
        break
    elif a[0] == '--outname' and a[1] == '':
        print('--outname expects an argument')
        sys.exit(0)
        break

#----------------------------------------------------------------------------------------------------------------------#
# MAIN

# open the tle file
try:
    fid = open(bp_file, 'r')
except IOError:
    print('Can\'t open file \'' + bp_file + '\' for reading.')
    sys.exit(0)

# read the file into an array
lines = fid.readlines()

# check that something was read in
if len(lines)==0:
    print("Erorr reading bandpass file: no lines read from the file.")
    sys.exit(0)
elif len(lines)%8!=1:
    print("Erorr reading bandpass file: number of lines should be 1 plus a multiple of 8.")
    sys.exit(0)

# close the file
fid.close()

# initialise the body list
PX_lsq = []
PY_lsq = []
QX_lsq = []
QY_lsq = []
PX_fit = []
PY_fit = []
QX_fit = []
QY_fit = []

tmp1 = str.split( lines[0], ',' )
tmp2 = []
for val in tmp1:
    tmp2.append( float(val) )
N_ch = len(tmp2)

freq = zeros(N_ch)
for k in range(0,N_ch):
    freq[k] = tmp2[k]

if plot_raw:
    ch = range(0,N_ch)
else:
    ch = zeros(N_ch)
    for k in range(0,N_ch):
        ch[k] = freq[k]/0.04

freq_idx = argsort(freq)

chan_sel = sel_offset + 2*array(freq_idx)

for lineIndex in range(1, len(lines), 8):  # get 0, 8, 16, ...

    tmp = str.split( lines[lineIndex+0], ',' ); PX_lsq.append([]); PX_lsq[-1]=zeros(N_ch*2+1)
    for k in range(0,len(tmp)): PX_lsq[-1][k] = float(tmp[k])
    tmp = str.split( lines[lineIndex+1], ',' ); PX_fit.append([]); PX_fit[-1]=zeros(N_ch*2+1)
    for k in range(0,len(tmp)): PX_fit[-1][k] = float(tmp[k])
    tmp = str.split( lines[lineIndex+2], ',' ); PY_lsq.append([]); PY_lsq[-1]=zeros(N_ch*2+1)
    for k in range(0,len(tmp)): PY_lsq[-1][k] = float(tmp[k])
    tmp = str.split( lines[lineIndex+3], ',' ); PY_fit.append([]); PY_fit[-1]=zeros(N_ch*2+1)
    for k in range(0,len(tmp)): PY_fit[-1][k] = float(tmp[k])
    tmp = str.split( lines[lineIndex+4], ',' ); QX_lsq.append([]); QX_lsq[-1]=zeros(N_ch*2+1)
    for k in range(0,len(tmp)): QX_lsq[-1][k] = float(tmp[k])
    tmp = str.split( lines[lineIndex+5], ',' ); QX_fit.append([]); QX_fit[-1]=zeros(N_ch*2+1)
    for k in range(0,len(tmp)): QX_fit[-1][k] = float(tmp[k])
    tmp = str.split( lines[lineIndex+6], ',' ); QY_lsq.append([]); QY_lsq[-1]=zeros(N_ch*2+1)
    for k in range(0,len(tmp)): QY_lsq[-1][k] = float(tmp[k])
    tmp = str.split( lines[lineIndex+7], ',' ); QY_fit.append([]); QY_fit[-1]=zeros(N_ch*2+1)
    for k in range(0,len(tmp)): QY_fit[-1][k] = float(tmp[k])

N_ant = len(PX_lsq)

print("plotting Jones matrices for %d frequency channels and %d antennas" % (N_ch, N_ant))

# Step through the channels and print spectra

fig1 = figure(figsize=(12,8))

if sel_offset == 1:

    figure(1); fig1.canvas.set_window_title('Amp. of the Jones matrices fits') 

elif sel_offset == 2:

    figure(1); fig1.canvas.set_window_title('Phase of the Jones matrices fits') 

    # unwrap 2pi phase jumps
    for ant in range(0, N_ant):
        last_PX_lsq_phase = PX_lsq[ant][2]
        last_PY_lsq_phase = PY_lsq[ant][2]
        last_QX_lsq_phase = QX_lsq[ant][2]
        last_QY_lsq_phase = QY_lsq[ant][2]
        last_PX_fit_phase = PX_fit[ant][2]
        last_PY_fit_phase = PY_fit[ant][2]
        last_QX_fit_phase = QX_fit[ant][2]
        last_QY_fit_phase = QY_fit[ant][2]
        PX_lsq[ant][2] *= 180/pi
        PY_lsq[ant][2] *= 180/pi
        QX_lsq[ant][2] *= 180/pi
        QY_lsq[ant][2] *= 180/pi
        PX_fit[ant][2] *= 180/pi
        PY_fit[ant][2] *= 180/pi
        QX_fit[ant][2] *= 180/pi
        QY_fit[ant][2] *= 180/pi

        for offset in range(4,2*N_ch+1,2):

            #if PX_lsq[ant][offset]-last_PX_lsq_phase > pi: PX_lsq[ant][offset] -= 2*pi
            #if PY_lsq[ant][offset]-last_PY_lsq_phase > pi: PY_lsq[ant][offset] -= 2*pi
            #if QX_lsq[ant][offset]-last_QX_lsq_phase > pi: QX_lsq[ant][offset] -= 2*pi
            #if QY_lsq[ant][offset]-last_QY_lsq_phase > pi: QY_lsq[ant][offset] -= 2*pi
            #if PX_fit[ant][offset]-last_PX_fit_phase > pi: PX_fit[ant][offset] -= 2*pi
            #if PY_fit[ant][offset]-last_PY_fit_phase > pi: PY_fit[ant][offset] -= 2*pi
            #if QX_fit[ant][offset]-last_QX_fit_phase > pi: QX_fit[ant][offset] -= 2*pi
            #if QY_fit[ant][offset]-last_QY_fit_phase > pi: QY_fit[ant][offset] -= 2*pi

            #if PX_lsq[ant][offset]-last_PX_lsq_phase < -pi: PX_lsq[ant][offset] += 2*pi
            #if PY_lsq[ant][offset]-last_PY_lsq_phase < -pi: PY_lsq[ant][offset] += 2*pi
            #if QX_lsq[ant][offset]-last_QX_lsq_phase < -pi: QX_lsq[ant][offset] += 2*pi
            #if QY_lsq[ant][offset]-last_QY_lsq_phase < -pi: QY_lsq[ant][offset] += 2*pi
            #if PX_fit[ant][offset]-last_PX_fit_phase < -pi: PX_fit[ant][offset] += 2*pi
            #if PY_fit[ant][offset]-last_PY_fit_phase < -pi: PY_fit[ant][offset] += 2*pi
            #if QX_fit[ant][offset]-last_QX_fit_phase < -pi: QX_fit[ant][offset] += 2*pi
            #if QY_fit[ant][offset]-last_QY_fit_phase < -pi: QY_fit[ant][offset] += 2*pi

            last_PX_lsq_phase = PX_lsq[ant][offset]
            last_PY_lsq_phase = PY_lsq[ant][offset]
            last_QX_lsq_phase = QX_lsq[ant][offset]
            last_QY_lsq_phase = QY_lsq[ant][offset]
            last_PX_fit_phase = PX_fit[ant][offset]
            last_PY_fit_phase = PY_fit[ant][offset]
            last_QX_fit_phase = QX_fit[ant][offset]
            last_QY_fit_phase = QY_fit[ant][offset]

            PX_lsq[ant][offset] *= 180/pi
            PY_lsq[ant][offset] *= 180/pi
            QX_lsq[ant][offset] *= 180/pi
            QY_lsq[ant][offset] *= 180/pi
            PX_fit[ant][offset] *= 180/pi
            PY_fit[ant][offset] *= 180/pi
            QX_fit[ant][offset] *= 180/pi
            QY_fit[ant][offset] *= 180/pi

# --------------------------------------------------------------------------------- #

if plot_chan==1: freq=ch

band_start = freq[0]
quart_bw = ( freq[-1] - band_start ) / 4.0

rc('font', size=8)
legendfont = matplotlib.font_manager.FontProperties(size=8)

# --------------------------------------------------------------------------------- #

if sel_offset == 1: thresh = 0.35
#if sel_offset == 1: thresh = 0.75
if sel_offset == 2: thresh = 30

#fig, axs = plt.subplots(2,2)
#ax1 = axs[0,0]
#ax2 = axs[0,1]
#ax3 = axs[1,0]
#ax4 = axs[1,1]
ax1 = plt.subplot(2, 2, 1)
ax2 = plt.subplot(2, 2, 2)
ax3 = plt.subplot(2, 2, 3)
ax4 = plt.subplot(2, 2, 4)


for ant in range(0, N_ant):

    tile = PX_lsq[ant][0]

    if ant==0:
      print(ant, tile)

    val = max(fabs(PX_lsq[ant][chan_sel]))
    if fabs( val - 1.0 ) < thresh:
        ax1.plot(freq[freq_idx],PX_lsq[ant][chan_sel], '-b.')
    else:
        ax1.plot(freq[freq_idx],PX_lsq[ant][chan_sel], '-c.', label='ID=%3d, max=%f (flag %3d?)'%(tile,val,tile-1))
        print('Possible PP flag: ID=%3d, max=%f (flag %3d?)'%(tile,val,tile-1))

    val = max(fabs(PY_lsq[ant][chan_sel]))
    if val < thresh or sel_offset == 2:
        ax2.plot(freq[freq_idx],PY_lsq[ant][chan_sel], '-b.')
    else:
        ax2.plot(freq[freq_idx],PY_lsq[ant][chan_sel], '-c.', label='ID=%3d, max=%f (flag %3d?)'%(tile,val,tile-1))
        if sel_offset == 1:
            print('Possible PQ flag: ID=%3d, max=%f (flag %3d?)'%(tile,val,tile-1))

    val = max(fabs(QX_lsq[ant][chan_sel]))
    if val < thresh or sel_offset == 2:
        ax3.plot(freq[freq_idx],QX_lsq[ant][chan_sel], '-b.')
    else:
        ax3.plot(freq[freq_idx],QX_lsq[ant][chan_sel], '-c.', label='ID=%3d, max=%f (flag %3d?)'%(tile,val,tile-1))
        if sel_offset == 1:
            print('Possible QP flag: ID=%3d, max=%f (flag %3d?)'%(tile,val,tile-1))

    val = max(fabs(QY_lsq[ant][chan_sel]))
    if fabs( val - 1.0 ) < thresh:
        ax4.plot(freq[freq_idx],QY_lsq[ant][chan_sel], '-b.')
    else:
        ax4.plot(freq[freq_idx],QY_lsq[ant][chan_sel], '-c.', label='ID=%3d, max=%f (flag %3d?)'%(tile,val,tile-1))
        print('Possible QQ flag: ID=%3d, max=%f (flag %3d?)'%(tile,val,tile-1))

# --------------------------------------------------------------------------------- #

for ant in range(0, N_ant):

    tile = PX_lsq[ant][0]

    val = max(fabs(PX_lsq[ant][chan_sel]))
    if fabs( val - 1.0 ) < thresh:
        ax1.plot(freq[freq_idx],PX_fit[ant][chan_sel], '--r')
    else:
        ax1.plot(freq[freq_idx],PX_fit[ant][chan_sel], '-k')

    val = max(fabs(PY_lsq[ant][chan_sel]))
    if val < thresh or sel_offset == 2:
        ax2.plot(freq[freq_idx],PY_fit[ant][chan_sel], '--r')
    else:
        ax2.plot(freq[freq_idx],PY_fit[ant][chan_sel], '-k')

    val = max(fabs(QX_lsq[ant][chan_sel]))
    if val < thresh or sel_offset == 2:
        ax3.plot(freq[freq_idx],QX_fit[ant][chan_sel], '--r')
    else:
        ax3.plot(freq[freq_idx],QX_fit[ant][chan_sel], '-k')

    val = max(fabs(QY_lsq[ant][chan_sel]))
    if fabs( val - 1.0 ) < thresh:
        ax4.plot(freq[freq_idx],QY_fit[ant][chan_sel], '--r')
    else:
        ax4.plot(freq[freq_idx],QY_fit[ant][chan_sel], '-k')

# --------------------------------------------------------------------------------- #

#legend(loc=[0,0.85],prop=legendfont)

ax1.legend(prop=legendfont)
ax1.set_title( "Jones matrix element P <- X" )
ax2.set_title( "Jones matrix element P <- Y" )
ax3.set_title( "Jones matrix element Q <- X" )
ax4.legend(prop=legendfont)
ax4.set_title( "Jones matrix element Q <- Y" )


for ax in [ax1, ax2, ax3, ax4]:
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(5))
    if plot_chan==1: ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
    if plot_chan==0: ax.set_xlabel( "MHz" )
    if plot_chan==1: ax.set_xlabel( "dumb chan # (doesn't skip flagged channels)" )
    if sel_offset == 1: ax.set_ylabel( "gain relative to band average" )
    if sel_offset == 2: ax.set_ylabel( "degrees" )
    #if sel_offset == 2: ylim(-180,180)
    ax.grid(True)

if save_plot:
    savefig(outname)
else:
    show()