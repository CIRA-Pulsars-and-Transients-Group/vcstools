#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "slalib.h"
#include "slamac.h"
#include "psrfits.h"
#include "beam_common.h"
#include "beam_psrfits.h"
#include "mycomplex.h"

void printf_psrfits( struct psrfits *pf ) {
    fprintf(stdout, "\nPSRFITS:\n");
    fprintf(stdout, "Basename of output file     [%s]\n", pf->basefilename);
    fprintf(stdout, "Filename of output file     [%s]\n", pf->filename);
    fprintf(stdout, "CFITSIO file pointer        [%p]\n", pf->fptr);

    fprintf(stdout, "\nPSRFITS HDRINFO:\n");
    fprintf(stdout, "Obs mode                    [%s]\n", pf->hdr.obs_mode);
    fprintf(stdout, "Telescope                   [%s]\n", pf->hdr.telescope);
    fprintf(stdout, "Observer                    [%s]\n", pf->hdr.observer);
    fprintf(stdout, "Source                      [%s]\n", pf->hdr.source);
    fprintf(stdout, "Front End                   [%s]\n", pf->hdr.frontend);
    fprintf(stdout, "Back End                    [%s]\n", pf->hdr.backend);
    fprintf(stdout, "Project ID                  [%s]\n", pf->hdr.project_id);
    fprintf(stdout, "Date Obs                    [%s]\n", pf->hdr.date_obs);
    fprintf(stdout, "RA (string)                 [%s]\n", pf->hdr.ra_str);
    fprintf(stdout, "Dec (string)                [%s]\n", pf->hdr.dec_str);
    fprintf(stdout, "Pol recorded (LIN or CIRC)  [%s]\n", pf->hdr.poln_type);
    fprintf(stdout, "Order of pols               [%s]\n", pf->hdr.poln_order);
    fprintf(stdout, "Track mode                  [%s]\n", pf->hdr.track_mode);
    fprintf(stdout, "Cal mode                    [%s]\n", pf->hdr.cal_mode);
    fprintf(stdout, "Feed mode                   [%s]\n", pf->hdr.feed_mode);
    fprintf(stdout, "Start MJD                   [%Lf]\n", pf->hdr.MJD_epoch);
    fprintf(stdout, "Sample Time (s)             [%lf]\n", pf->hdr.dt);
    fprintf(stdout, "Centre Frequency (MHz)      [%lf]\n", pf->hdr.fctr);
    fprintf(stdout, "Orig freq spacing (MHz)     [%lf]\n", pf->hdr.orig_df);
    fprintf(stdout, "Freq spacing (MHz)          [%lf]\n", pf->hdr.df);
    fprintf(stdout, "Bandwidth (MHz)             [%lf]\n", pf->hdr.BW);
    fprintf(stdout, "RA (2000) (deg)             [%lf]\n", pf->hdr.ra2000);
    fprintf(stdout, "Dec (2000) (deg)            [%lf]\n", pf->hdr.dec2000);
    fprintf(stdout, "Azimuth (deg)               [%lf]\n", pf->hdr.azimuth);
    fprintf(stdout, "Zenith Angle (deg)          [%lf]\n", pf->hdr.zenith_ang);
    fprintf(stdout, "Beam FWHM (deg)             [%lf]\n", pf->hdr.beam_FWHM);

    fprintf(stdout, "Length of scan in this file [%lf]\n", pf->hdr.scanlen);
    fprintf(stdout, "Seconds past 00h LST        [%lf]\n", pf->hdr.start_lst);
    fprintf(stdout, "Seconds past 00h UTC        [%lf]\n", pf->hdr.start_sec);

    fprintf(stdout, "Start MJD (whole day)       [%d]\n", pf->hdr.start_day);
    fprintf(stdout, "Scan Number                 [%d]\n", pf->hdr.scan_number);
    fprintf(stdout, "Number of bits per sample   [%d]\n", pf->hdr.nbits);

    fprintf(stdout, "Number of Channels          [%d]\n", pf->hdr.nchan);
    fprintf(stdout, "Number of polarisations     [%d]\n", pf->hdr.npol);
    fprintf(stdout, "Number of spectra per row   [%d]\n", pf->hdr.nsblk);

    fprintf(stdout, "Summed Polarisations? [1/0] [%d]\n", pf->hdr.summed_polns);
    fprintf(stdout, "Receiver Polarisation       [%d]\n", pf->hdr.rcvr_polns);
    fprintf(stdout, "Offset Subint               [%d]\n", pf->hdr.offset_subint);
    fprintf(stdout, "Dwnsmpl fact in time        [%d]\n", pf->hdr.ds_time_fact);
    fprintf(stdout, "Dwnsmpl fact in freq        [%d]\n", pf->hdr.ds_freq_fact);
    fprintf(stdout, "Only Stokes I?              [%d]\n", pf->hdr.onlyI);

    fprintf(stdout, "\nPSRFITS SUBINT:\n");
    fprintf(stdout, "Length of subint (sec)      [%lf]\n", pf->sub.tsubint);
    fprintf(stdout, "Offset (sec)                [%lf]\n", pf->sub.offs);
    fprintf(stdout, "LST (sec)                   [%lf]\n", pf->sub.lst);
    fprintf(stdout, "RA (J2000) (deg)            [%lf]\n", pf->sub.ra);
    fprintf(stdout, "Dec (J2000) (deg)           [%lf]\n", pf->sub.dec);
    fprintf(stdout, "Gal. long. (deg)            [%lf]\n", pf->sub.glon);
    fprintf(stdout, "Gal. lat. (deg)             [%lf]\n", pf->sub.glat);
    fprintf(stdout, "Feed angle (deg)            [%lf]\n", pf->sub.feed_ang);
    fprintf(stdout, "Pos angle of feed (deg)     [%lf]\n", pf->sub.pos_ang);
    fprintf(stdout, "Parallactic angle           [%lf]\n", pf->sub.par_ang);
    fprintf(stdout, "Telescope azimuth           [%lf]\n", pf->sub.tel_az);
    fprintf(stdout, "Telescope zenith angle      [%lf]\n", pf->sub.tel_zen);
    fprintf(stdout, "Bytes per row of raw data   [%d]\n", pf->sub.bytes_per_subint);
    fprintf(stdout, "FITS data typecode          [%d]\n", pf->sub.FITS_typecode);

}

void populate_psrfits_header(
        struct psrfits *pf,
        char           *metafits,
        char           *obsid,
        char           *time_utc,
        unsigned int    sample_rate,
        long int        frequency,
        int             nchan,
        long int        chan_width,
        int             outpol,
        char           *rec_channel,
        struct delays  *delay_vals )
{
    int is_coherent = 0;
    if (outpol == 4)
        is_coherent = 1;
    else if (outpol != 1)
    {
        fprintf( stderr, "warning: populate_psrfits_header: "
                         "unusual number of output pols = %d\n", outpol );
    }


    fitsfile *fptr = NULL;
    int status     = 0;

    fits_open_file(&fptr, metafits, READONLY, &status);
    fits_read_key(fptr, TSTRING, "PROJECT", pf->hdr.project_id, NULL, &status);
    fits_close_file(fptr, &status);

    // Now set values for our hdrinfo structure
    strcpy(pf->hdr.obs_mode,  "SEARCH");
    strcpy(pf->hdr.observer,  "MWA User");
    strcpy(pf->hdr.telescope, "MWA");
    strncpy(pf->hdr.source, obsid, 23);
    pf->hdr.scanlen = 1.0; // in sec

    strcpy(pf->hdr.frontend, "MWA-RECVR");
    snprintf(pf->hdr.backend, 24*sizeof(char), "vcstools %s", VERSION_BEAMFORMER );

    // Now let us finally get the time right
    strcpy(pf->hdr.date_obs,   time_utc);
    strcpy(pf->hdr.poln_type,  "LIN");
    strcpy(pf->hdr.track_mode, "TRACK");
    strcpy(pf->hdr.cal_mode,   "OFF");
    strcpy(pf->hdr.feed_mode,  "FA");

    pf->hdr.dt   = 1.0/sample_rate;                            // (sec)
    pf->hdr.fctr = (frequency + (nchan/2.0)*chan_width)/1.0e6; // (MHz)
    pf->hdr.BW   = (nchan*chan_width)/1.0e6;                   // (MHz)

    // npols + nbits and whether pols are added
    pf->filenum       = 0;       // This is the crucial one to set to initialize things
    pf->rows_per_file = 200;     // I assume this is a max subint issue

    pf->hdr.npol         = outpol;
    pf->hdr.nchan        = nchan;
    pf->hdr.onlyI        = 0;

    pf->hdr.scan_number   = 1;
    pf->hdr.rcvr_polns    = 2;
    pf->hdr.offset_subint = 0;

    if (is_coherent)
        pf->hdr.summed_polns = 0;
    else
        pf->hdr.summed_polns = 1;

    pf->hdr.df         = chan_width/1.0e6; // (MHz)
    pf->hdr.orig_nchan = pf->hdr.nchan;
    pf->hdr.orig_df    = pf->hdr.df;
    pf->hdr.nbits      = 8;
    pf->hdr.nsblk      = sample_rate;  // block is always 1 second of data

    pf->hdr.ds_freq_fact = 1;
    pf->hdr.ds_time_fact = 1;

    // some things that we are unlikely to change
    pf->hdr.fd_hand  = 1;
    pf->hdr.fd_sang  = 45.0;
    pf->hdr.fd_xyph  = 0.0;
    pf->hdr.be_phase = 0;
    pf->hdr.chan_dm  = 0.0;

    // Now set values for our subint structure
    pf->tot_rows     = 0;
    pf->sub.tsubint  = roundf(pf->hdr.nsblk * pf->hdr.dt);
    pf->sub.offs     = roundf(pf->tot_rows * pf->sub.tsubint) + 0.5*pf->sub.tsubint;

    pf->sub.feed_ang = 0.0;
    pf->sub.pos_ang  = 0.0;
    pf->sub.par_ang  = 0.0;

    // Specify psrfits data type
    pf->sub.FITS_typecode = TBYTE;

    pf->sub.bytes_per_subint = (pf->hdr.nbits * pf->hdr.nchan *
                                pf->hdr.npol  * pf->hdr.nsblk) / 8;

    // Create and initialize the subint arrays
    pf->sub.dat_freqs   = (float *)malloc(sizeof(float) * pf->hdr.nchan);
    pf->sub.dat_weights = (float *)malloc(sizeof(float) * pf->hdr.nchan);

    double dtmp = pf->hdr.fctr - 0.5 * pf->hdr.BW + 0.5 * pf->hdr.df;
    int i;
    for (i = 0 ; i < pf->hdr.nchan ; i++) {
        pf->sub.dat_freqs[i] = dtmp + i * pf->hdr.df;
        pf->sub.dat_weights[i] = 1.0;
    }

    // the following is definitely needed for 8 bit numbers
    pf->sub.dat_offsets = (float *)malloc(sizeof(float) * pf->hdr.nchan * pf->hdr.npol);
    pf->sub.dat_scales  = (float *)malloc(sizeof(float) * pf->hdr.nchan * pf->hdr.npol);
    for (i = 0 ; i < pf->hdr.nchan * pf->hdr.npol ; i++) {
        pf->sub.dat_offsets[i] = 0.0;
        pf->sub.dat_scales[i]  = 1.0;
    }

    pf->sub.data    = (unsigned char *)malloc(pf->sub.bytes_per_subint);
    pf->sub.rawdata = pf->sub.data;

    int ch = atoi(rec_channel);
    if (!is_coherent)
    {
        sprintf(pf->basefilename, "%s_%s_ch%03d_incoh",
                pf->hdr.project_id, pf->hdr.source, ch);
    }
    else
    {
        sprintf(pf->basefilename, "%s_%s_ch%03d",
                pf->hdr.project_id, pf->hdr.source, ch);
    }

    // Update values that depend on get_delays()
    if (delay_vals != NULL) {

        pf->hdr.ra2000  = delay_vals->mean_ra  * DR2D;
        pf->hdr.dec2000 = delay_vals->mean_dec * DR2D;

        dec2hms(pf->hdr.ra_str,  pf->hdr.ra2000/15.0, 0);
        dec2hms(pf->hdr.dec_str, pf->hdr.dec2000,     1);

        pf->hdr.azimuth    = delay_vals->az*DR2D;
        pf->hdr.zenith_ang = 90.0 - (delay_vals->el*DR2D);

        pf->hdr.beam_FWHM = 0.25;
        pf->hdr.start_lst = delay_vals->lmst * 60.0 * 60.0;        // Local Apparent Sidereal Time in seconds
        pf->hdr.start_sec = roundf(delay_vals->fracmjd*86400.0);   // this will always be a whole second
        pf->hdr.start_day = delay_vals->intmjd;
        pf->hdr.MJD_epoch = delay_vals->intmjd + delay_vals->fracmjd;

        // Now set values for our subint structure
        pf->sub.lst      = pf->hdr.start_lst;
        pf->sub.ra       = pf->hdr.ra2000;
        pf->sub.dec      = pf->hdr.dec2000;
        slaEqgal(pf->hdr.ra2000*DD2R, pf->hdr.dec2000*DD2R,
                 &pf->sub.glon, &pf->sub.glat);
        pf->sub.glon    *= DR2D;
        pf->sub.glat    *= DR2D;
        pf->sub.tel_az   = pf->hdr.azimuth;
        pf->sub.tel_zen  = pf->hdr.zenith_ang;
    }
}


void correct_psrfits_stt( struct psrfits *pf )
{
    /* now we have to correct the STT_SMJD/STT_OFFS as they will have been broken by the write_psrfits*/
    int    itmp    = 0;
    int    itmp2   = 0;
    double dtmp    = 0;
    int    status  = 0;

    //fits_open_file(&(pf.fptr), pf.filename, READWRITE, &status);

    fits_read_key(pf->fptr, TDOUBLE, "STT_OFFS", &dtmp,  NULL, &status);
    fits_read_key(pf->fptr, TINT,    "STT_SMJD", &itmp,  NULL, &status);
    fits_read_key(pf->fptr, TINT,    "STT_IMJD", &itmp2, NULL, &status);

    if (dtmp > 0.5) {
        itmp = itmp+1;
        if (itmp == 86400) {
            itmp = 0;
            itmp2++;
        }
    }
    dtmp = 0.0;

    fits_update_key(pf->fptr, TINT, "STT_SMJD", &itmp, NULL, &status);
    fits_update_key(pf->fptr, TINT, "STT_IMJD", &itmp2, NULL, &status);
    fits_update_key(pf->fptr, TDOUBLE, "STT_OFFS", &dtmp, NULL, &status);

}


void psrfits_write_second( struct psrfits *pf, float *data_buffer, int nchan,
        int outpol )
{
    int8_t *out_buffer_8 = (int8_t *)malloc( outpol*nchan*pf->hdr.nsblk * sizeof(int8_t) );

    flatten_bandpass( pf->hdr.nsblk, nchan, outpol, data_buffer);

    float_to_unit8(data_buffer, pf->hdr.nsblk * nchan * outpol, out_buffer_8);

    memcpy( pf->sub.data, out_buffer_8, pf->sub.bytes_per_subint );

    if (psrfits_write_subint(pf) != 0)
    {
        fprintf(stderr, "error: Write subint failed. File exists?\n");
        exit(EXIT_FAILURE);
    }

    pf->sub.offs = roundf(pf->tot_rows * pf->sub.tsubint) + 0.5*pf->sub.tsubint;
    pf->sub.lst += pf->sub.tsubint;

    free( out_buffer_8 );
}


