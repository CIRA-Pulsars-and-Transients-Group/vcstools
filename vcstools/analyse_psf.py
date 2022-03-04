"""
Scripts to analyise the output PSFs
"""
import numpy as np
from itertools import chain
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from astropy.time import Time
from astropy.constants import c
from astropy.io import fits

from vcstools.pointing_utils import getTargetRADec, deg2sex

import logging
logger = logging.getLogger(__name__)


def read_vcsbeam_psf(psf_file):
    """Read in the PSF data file output from vcsbeam's mwa_tied_array_beam_psf

    Parameters
    ----------
    psf_file : `str`
        The PSF data file output from vcsbeam's mwa_tied_array_beam_psf.

    Returns
    -------
    ra : `numpy.array`, (Ny, Nx)
        The RA (degrees) in the meshgrid format.
    dec : `numpy.array`, (Ny, Nx)
        The declination (degrees) in the meshgrid format.
    power : `numpy.array`, (Ny, Nx)
        The the data values of the fits file in the meshgrid format.
    """
    input_array = np.loadtxt(psf_file, dtype=float)
    ra    = input_array[:,0]
    dec   = input_array[:,1]
    power = input_array[:,2]
    ra.shape = dec.shape = power.shape = (int(np.sqrt(input_array.shape[0])),
                                          int(np.sqrt(input_array.shape[0])))
    return ra, dec, power


def plot_vcsbeam_psf(psf_file, output_name="vcsbeam_psf.png", normalise=False, vmin=None):
    """Plot the PSF data file output from vcsbeam's mwa_tied_array_beam_psf

    Parameters
    ----------
    psf_file : `str`
        The PSF data file output from vcsbeam's mwa_tied_array_beam_psf.
    output_name : `str`, optional
        Output plot name. |br| Default: vcsbeam_psf.png.
    normalise : `boolean`, optional
        Normalise the PSF. |br| Default: `False`.
    vmin : `float`, optional
        Minimum value of plot. |br| Default: `None`.
    """
    ra, dec, power = read_vcsbeam_psf(psf_file)

    if normalise:
        power = power / np.amax(power)

    ax = plt.subplot(1, 1, 1,)
    im = ax.pcolormesh(ra, dec, power, norm=colors.LogNorm(), cmap='plasma', vmin=vmin)

    plt.xlabel(r"Right Ascension (hours)")
    plt.ylabel(r"Declination ($^{\circ}$)")
    plt.colorbar(im,#spacing='uniform', shrink = 0.65, #ticks=[2., 10., 20., 30., 40., 50.],
                 label="Normalised array factor power")
    plt.savefig(output_name)


def plot_track_beam_response(response_file, output_name="vcsbeam_response.png", time_max=-1):
    """Plot the data file output from vcsbeam's mwa_track_primary_beam_response

    Parameters
    ----------
    response_file : `str`
        The PSF data file output from vcsbeam's mwa_track_primary_beam_response.
    output_name : `str`, optional
        Output plot name. |br| Default: vcsbeam_psf.png.
    time_max : `int`, optional
        Maximum number of seconds to process. |br| Default: -1.
    """
    input_array = np.loadtxt(response_file, dtype=float)
    sec          = input_array[:,0]
    freq         = input_array[:,1]
    stokes_I     = input_array[:,5]
    array_factor = input_array[:,9]
    response = stokes_I * array_factor**2

    # Split into freq chunks
    nsec = max(sec)
    nfreq = input_array.shape[0] // nsec
    sec_map   = np.array(np.array_split(sec,      nfreq))[:,:time_max]
    freq_map  = np.array(np.array_split(freq,     nfreq))[:,:time_max]
    power_map = np.array(np.array_split(response, nfreq))[:,:time_max]

    # Sum over frequency
    power_time = np.average(power_map, axis=0)
    sec_time = sec_map[0]

    # Sum over time
    power_freq = np.average(power_map, axis=1)
    freq_freq = freq_map[:,0]

    size = 3
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(2*size,size))

    # Plot rack over time
    ax1.plot(sec_time, power_time)
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Array Factor")

    # Plot rack over time
    ax2.plot(freq_freq, power_freq)
    ax2.set_xlabel("Frequency (MHz)")
    ax2.set_ylabel("Array Factor")
    plt.tight_layout()
    plt.savefig(output_name)


def plot_pabeam(dat_file, output_name="pabeam_psf.png"):
    """Makes a polar plot over the sky response with the data file output from pabeam.py

    Parameters
    ----------
    dat_file : `str`
        The data file output from pabeam.py.
    output_name : `str`, optional
        Output plot name. |br| Default: pabeam_psf.png.
    """
    with open(dat_file) as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]
        nza = int(lines[7].split(" ")[-1])
        naz = int(lines[8].split(" ")[-1])
    input_array = np.loadtxt(dat_file, dtype=float)
    za    = input_array[:,0]
    az    = input_array[:,1]
    power = input_array[:,2]
    za.shape = az.shape = power.shape = (naz, nza)

    ax = plt.subplot(1, 1, 1, projection='polar')
    im = plt.pcolormesh(np.radians(az), za, power)
    plt.xlabel(r"Azimuth ($^{\circ}$)")
    plt.ylabel(r"Zenith ($^{\circ}$)")
    plt.colorbar(im, label="Normalised array factor power")
    plt.savefig(output_name, dpi=500)


def read_pabeam_ra_dec(dat_file):
    """Read in the PSF data file output from pabeam.py when using 
    the --ra_dec_projection option.

    Parameters
    ----------
    dat_file : `str`
        The PSF data file output from pabeam.py

    Returns
    -------
    ra : `numpy.array`, (Nx, Ny)
        The RA (degrees) in the meshgrid format
    dec : `numpy.array`, (Nx, Ny)
        The declination (degrees) in the meshgrid format
    power : `numpy.array`, (Nx, Ny)
        The the data values of the fits file in the meshgrid format
    """
    with open(dat_file) as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]
        date  = lines[3].split("'")[-2]
        time = Time(date, format='iso', scale='utc')
        nra   = int(lines[7].split(" ")[-1])
        ndec  = int(lines[8].split(" ")[-1])
    input_array = np.loadtxt(dat_file, dtype=float)
    za    = input_array[:,0]
    az    = input_array[:,1]
    power = input_array[:,2]

    # Convert to RA and Dec
    _, _, ra, dec = getTargetRADec(az, za, time)
    ra.shape = dec.shape = power.shape = (nra, ndec)
    return ra, dec, power

def plot_pabeam_ra_dec(dat_file, output_name="pabeam_psf.png",
                       normalise=False, vmin=None):
    """Plots the PSF data file output pabeam.py when using 
    the --ra_dec_projection option.

    Parameters
    ----------
    dat_file : `str`
        The PSF data file output pabeam.py when using 
        the --ra_dec_projection option
    output_name : `str`, optional
        Output plot name. |br| Default: pabeam_psf.png.
    normalise : `boolean`, optional
        Normalise the PSF. |br| Default: `False`.
    vmin : `float`, optional
        Minimum value of plot. |br| Default: `None`.
    """
    ra, dec, power = read_pabeam_ra_dec(dat_file)
    ax = plt.subplot(1, 1, 1,)

    if normalise:
        power = power / np.amax(power)

    if fft_abs:
        # Take the fft of the PSF and plot the amplitude
        fits_fft = np.fft.fftshift(np.fft.fft2(power))
        im = plt.imshow(np.abs(fits_fft), cmap='plasma', norm=colors.LogNorm())
    elif fft_angle:
        # Take the fft of the PSF and plot the phase
        fits_fft = np.fft.fftshift(np.fft.fft2(power))
        im = plt.imshow(np.angle(fits_fft), cmap='htp',
                        vmin=-np.pi, vmax=np.pi)
    else:
        # Create plot
        im = plt.pcolormesh(ra, dec, power, cmap='plasma',
                            vmin=vmin, norm=colors.LogNorm(), vmax=0.3)
        plt.xlabel(r"Right Acension ($^{\circ}$)")
        plt.ylabel(r"Declination ($^{\circ}$)")

    plt.colorbar(im)
    plt.savefig(output_name, dpi=500)


def read_psf_fits(fits_file, centre_size=None):
    """Read in an imaging PSF fits file (from WSCLEAN for example).

    Parameters
    ----------
    fits_file : `str`
        The PSF fits file gained from imaging (WSCLEAN).
    centre_size : `int`, optional
        The radius in pixels to plot from the centre. |br| Default: `None`.

    Returns
    -------
    rav : `numpy.array`, (Nx, Ny)
        The RA (degrees) in the meshgrid format.
    decv : `numpy.array`, (Nx, Ny)
        The declination (degrees) in the meshgrid format.
    fits_data : `numpy.array`, (Nx, Ny)
        The the data values of the fits file in the meshgrid format.
    """
    hdul = fits.open(fits_file)

    # RA header info
    ra_pix_res   = hdul[0].header['CDELT1']
    ra_ref_pixel = int(hdul[0].header['CRPIX1'])
    ra_ref_pos   = hdul[0].header['CRVAL1']
    if ra_ref_pos < 0.:
        ra_ref_pos = 360. + ra_ref_pos

    # Dec header info
    dec_pix_res   = hdul[0].header['CDELT2']
    dec_ref_pixel = int(hdul[0].header['CRPIX2'])
    dec_ref_pos   = hdul[0].header['CRVAL2']

    # Other header info
    centre_freq = hdul[0].header['CRVAL3'] #Hz
    if len(hdul[0].data.shape) == 4:
        _, _, ra_npixels, dec_npixels = hdul[0].data.shape
        fits_data = hdul[0].data[0][0]
    else:
        ra_npixels, dec_npixels = hdul[0].data.shape
        fits_data = hdul[0].data
    raj, decj = deg2sex(float(ra_ref_pos), float(dec_ref_pos))
    pointing = "{}_{}".format(raj,decj)

    # Calculate RA Dec range
    ra_start = ra_ref_pos + ra_ref_pixel * ra_pix_res
    ra_stop  = ra_start - ra_npixels * ra_pix_res
    ra_range = np.linspace(ra_stop, ra_start, ra_npixels)

    dec_start = dec_ref_pos - dec_ref_pixel * dec_pix_res
    dec_stop  = dec_start + dec_npixels * dec_pix_res
    dec_range = np.linspace(dec_start, dec_stop, dec_npixels)

    # Plot data
    rav, decv = np.meshgrid(ra_range, dec_range)

    if centre_size:
        # Grab the centre pixels
        rav = rav[ra_ref_pixel-centre_size:ra_ref_pixel+centre_size,
                  dec_ref_pixel-centre_size:dec_ref_pixel+centre_size]
        decv = decv[ra_ref_pixel-centre_size:ra_ref_pixel+centre_size,
                    dec_ref_pixel-centre_size:dec_ref_pixel+centre_size]
        fits_data = fits_data[ra_ref_pixel-centre_size:ra_ref_pixel+centre_size,
                              dec_ref_pixel-centre_size:dec_ref_pixel+centre_size]

    return rav, decv, fits_data

def plot_imaging_psf(fits_file, output_name="imaging_psf.png",
                     normalise=False, vmin=None, centre_size=None,
                     fft_abs=False, fft_angle=False, centre=False):
    """Plots the imaging PSF (from WSCLEAN for example).

    Parameters
    ----------
    fits_file : `str`
        The PSF fits file gained from imaging (WSCLEAN).
    output_name : `str`, optional
        Output plot name. |br| Default: imaging_psf.png.
    normalise : `boolean`, optional
        Normalise the PSF. |br| Default: `False`.
    vmin : `float`, optional
        Minimum value of plot. |br| Default: `None`.
    centre_size : `int`, optional
        The radius in pixels to plot from the centre. |br| Default: `None`.
    fft_abs : `boolean`, optional
        Take the fft of the PSF and plot the amplitude. |br| Default: `False`.
    fft_angle : `boolean`, optional
        Take the fft of the PSF and plot the phase. |br| Default: `False`.
    """
    rav, decv, fits_data = read_psf_fits(fits_file, centre_size=centre_size)
    fig, ax = plt.subplots()
    if normalise:
        fits_data = fits_data / np.amax(fits_data)
    if centre:
        xcentre = rav.shape[0] // 2
        ycentre = rav.shape[1] // 2
        rav  = (rav  - rav[xcentre][ycentre]) * 3600
        decv = (decv - decv[xcentre][ycentre]) * 3600
        circle1 = plt.Circle((0, 0), 30, color='r', fill=False)
        ax.add_patch(circle1)

    if fft_abs:
        # Take the fft of the PSF and plot the amplitude
        fits_fft = np.fft.fftshift(np.fft.fft2(fits_data))
        im = plt.imshow(np.abs(fits_fft), cmap='plasma', norm=colors.LogNorm())
    elif fft_angle:
        # Take the fft of the PSF and plot the phase
        fits_fft = np.fft.fftshift(np.fft.fft2(fits_data))
        im = plt.imshow(np.angle(fits_fft), cmap='htp',
                        vmin=-np.pi, vmax=np.pi)
    else:
        # Create plot
        im = plt.pcolormesh(rav, decv, fits_data, cmap='plasma',
                            vmin=vmin, norm=colors.LogNorm())
        #image_data = fits.getdata(fits_file, ext=0)
        #plt.imshow(image_data, cmap='plasma', extent=[np.amin(rav), np.amax(rav),
        #                                              np.amin(decv), np.amax(decv)])
        #plt.colorbar()#, label="Normalised array factor power")
        plt.xlabel(r"Right Acension ($^{\circ}$)")
        plt.ylabel(r"Declination ($^{\circ}$)")
    
    plt.colorbar(im)
    plt.savefig(output_name, dpi=500)

    #plt.clf()
    # Replot with squared power
    #im = plt.pcolormesh(rav, decv, fits_data**2, cmap='plasma',
    #                    vmin=vmin, norm=colors.LogNorm())
    #plt.colorbar(im)
    #plt.savefig("squared_{}".format(output_name), dpi=500)


def plot_psf_comparison(imaging_psf,
                        pabeam_psf,
                        c_pabeam_psf):
    """Compare the PSFs FWHM along RA and declination from imaging, pabeam.py and vcsbeam's mwa_tied_array_beam_psf.
    Outputs a plot called comparing_psf_cuts.png.

    Parameters
    ----------
    imaging_psf : `str`
        The PSF fits file gained from imaging (WSCLEAN).
    pabeam_psf : `str`
        The PSF data file output from pabeam.py.
    c_pabeam_psf : `str`
        The PSF data file output from vcsbeam's mwa_tied_array_beam_psf.
    """
    # read image psf
    rai, deci, ipower = read_psf_fits(imaging_psf)
    rai_middle, deci_middle = rai.shape
    rai_middle  = rai_middle // 2
    deci_middle = deci_middle // 2
    # read pabeam psf
    rap, decp, ppower = read_pabeam_ra_dec(pabeam_psf)
    rap_middle, decp_middle = rap.shape
    rap_middle  = rap_middle // 2
    decp_middle = decp_middle // 2
    # read Sammy's C pabeam psf
    rac, decc, cpower = read_vcsbeam_psf(c_pabeam_psf)
    rac = rac * 15
    rac_middle, decc_middle = rac.shape
    rac_middle  = rac_middle // 2
    decc_middle = decc_middle // 2
    # normalise
    ppower = ppower / np.amax(ppower)
    cpower = cpower / np.amax(cpower)
    ipower = ipower / np.amax(ipower)

    size = 3
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(2*size,size))
    pixel_radius = 30

    # remove offsets
    rap_diff = rap[rap_middle][rap_middle] - rai[rai_middle][rai_middle]
    rac_diff = rac[rac_middle][rac_middle] - rai[rai_middle][rai_middle]

    # plot RA cuts -----------------------------------------------------------------
    image_ra = list(rai[rai_middle][rai_middle-pixel_radius:rai_middle+pixel_radius])
    image_ra.reverse()
    spline = UnivariateSpline(image_ra,
                              ipower[rai_middle][rai_middle-pixel_radius:rai_middle+pixel_radius] - 0.5, s=0)
    enter_beam, exit_beam = spline.roots()
    fwhm = exit_beam - enter_beam
    ax1.plot(rai[rai_middle][rai_middle-pixel_radius:rai_middle+pixel_radius],
             ipower[rai_middle][rai_middle-pixel_radius:rai_middle+pixel_radius],
             label=r"Image FWHM: {:.2f}$^\prime$".format(fwhm*60))

    spline = UnivariateSpline(rap[rap_middle][rap_middle-pixel_radius:rap_middle+pixel_radius],
                              ppower[rap_middle][rap_middle-pixel_radius:rap_middle+pixel_radius] - 0.5, s=0)
    enter_beam, exit_beam = spline.roots()
    fwhm = exit_beam - enter_beam
    ax1.plot(rap[rap_middle][rap_middle-pixel_radius:rap_middle+pixel_radius] - rap_diff,
             ppower[rap_middle][rap_middle-pixel_radius:rap_middle+pixel_radius],
             label=r"Phased Array FWHM: {:.2f}$^\prime$".format(fwhm*60))

    psfc_slice = []
    rac_slice = []
    for i in range(rac_middle-pixel_radius, rac_middle+pixel_radius):
        psfc_slice.append(cpower[i][decc_middle])
        rac_slice.append(rac[i][rac_middle])
    psfc_slice = np.array(psfc_slice)
    rac_slice = np.array(rac_slice)
    #print(rac_slice)
    spline = UnivariateSpline(rac_slice, psfc_slice - 0.5, s=0)
    enter_beam, exit_beam = spline.roots()
    fwhm = exit_beam - enter_beam
    ax1.plot(rac_slice - rac_diff, psfc_slice,
             label=r"C Phased Array FWHM: {:.2f}$^\prime$".format(fwhm*60))

    ax1.set_xlabel(r"Right Acension ($^{\circ}$)")
    ax1.legend(loc="upper right", fontsize=6)

    # plot dec cuts -----------------------------------------------------------------
    psfi_slice = []
    deci_slice = []
    for i in range(deci_middle-pixel_radius, deci_middle+pixel_radius):
        psfi_slice.append(ipower[i][deci_middle])
        deci_slice.append(deci[i][deci_middle])
    psfi_slice = np.array(psfi_slice)
    deci_slice = np.array(deci_slice)
    spline = UnivariateSpline(deci_slice, psfi_slice - 0.5, s=0)
    enter_beam, exit_beam = spline.roots()
    fwhm = exit_beam - enter_beam
    deci_diff = deci_slice[pixel_radius] - decc[decc_middle][decc_middle]
    ax2.plot(deci_slice - deci_diff, psfi_slice,
             label=r"Image FWHM: {:.2f}$^\prime$".format(fwhm*60))

    psfp_slice = []
    decp_slice = []
    for i in range(decp_middle-pixel_radius, decp_middle+pixel_radius):
        psfp_slice.append(ppower[i][decp_middle])
        decp_slice.append(decp[i][decp_middle])
    psfp_slice = np.array(psfp_slice)
    decp_slice = np.array(decp_slice)
    spline = UnivariateSpline(decp_slice, psfp_slice - 0.5, s=0)
    enter_beam, exit_beam = spline.roots()
    fwhm = exit_beam - enter_beam
    ax2.plot(decp_slice, psfp_slice,
             label=r"Phased Array FWHM: {:.2f}$^\prime$".format(fwhm*60))

    spline = UnivariateSpline(decc[decc_middle][decc_middle-pixel_radius:decc_middle+pixel_radius],
                              cpower[decc_middle][decc_middle-pixel_radius:decc_middle+pixel_radius] - 0.5, s=0)
    enter_beam, exit_beam = spline.roots()
    fwhm = exit_beam - enter_beam
    ax2.plot(decc[decc_middle][decc_middle-pixel_radius:decc_middle+pixel_radius],
             cpower[decc_middle][decc_middle-pixel_radius:decc_middle+pixel_radius],
             label=r"C Phased Array FWHM: {:.2f}$^\prime$".format(fwhm*60))

    ax2.legend(loc="upper right", fontsize=6)
    ax2.set_xlabel(r"Declination ($^{\circ}$)")
    plt.savefig("comparing_psf_cuts.png", dpi=500, bbox_inches="tight")