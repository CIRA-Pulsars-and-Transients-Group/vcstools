/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "vdifio.h"
#include "psrfits.h"
#include "slamac.h"
#include "beam_common.h"
#include "beam_vdif.h"
#include "mwa_header.h"
#include "vdifio.h"
#include "ascii_header.h"
#include "filter.h"
#include "mycomplex.h"

#ifndef HAVE_CUDA
#include <omp.h>
#endif

void vdif_write_second( struct vdifinfo *vf, vdif_header *vhdr,
                        float *data_buffer_vdif, float *gain )
{

    // Set level occupancy
    float rmean, imean;
    ComplexFloat cmean, stddev;

    get_mean_complex(
            (ComplexFloat *)data_buffer_vdif,
            vf->sizeof_buffer/2.0,
            &rmean, &imean, &cmean );

    stddev = get_std_dev_complex(
            (ComplexFloat *)data_buffer_vdif,
            vf->sizeof_buffer/2.0 );

    if (fabsf(rmean) > 0.001)
    {
        fprintf( stderr, "warning: vdif_write_second: significantly "
                         "non-zero mean (%f), adjusting data\n", rmean );
        unsigned int i;
        for (i = 0; i < vf->sizeof_buffer/2; i++)
        {
            data_buffer_vdif[2*i+0] -= CRealf(cmean);
            data_buffer_vdif[2*i+1] -= CImagf(cmean);
        }
    }

    vf->b_scales[0] = CRealf(stddev);
    vf->b_scales[1] = CRealf(stddev);

    vf->got_scales = 1;
    set_level_occupancy(
            (ComplexFloat *)data_buffer_vdif,
            vf->sizeof_buffer/2.0, gain);

    // Normalise
    normalise_complex(
            (ComplexFloat *)data_buffer_vdif,
            vf->sizeof_buffer/2.0,
            1.0/(*gain) );

    float *data_buffer_ptr = data_buffer_vdif;
    size_t offset_out_vdif = 0;

    int8_t *out_buffer_8_vdif = (int8_t *)malloc(vf->block_size);

    while  (offset_out_vdif < vf->block_size) {

        // Add the current header
        memcpy( (out_buffer_8_vdif + offset_out_vdif), vhdr, VDIF_HEADER_SIZE );

        // Offset into the output array
        offset_out_vdif += VDIF_HEADER_SIZE;

        // Convert from float to int8
        float2int8_trunc( data_buffer_ptr, vf->sizeof_beam, -126.0, 127.0,
                          (out_buffer_8_vdif + offset_out_vdif) );
        to_offset_binary( (out_buffer_8_vdif + offset_out_vdif),
                          vf->sizeof_beam );

        offset_out_vdif += vf->frame_length - VDIF_HEADER_SIZE; // increment output offset
        data_buffer_ptr += vf->sizeof_beam;
        nextVDIFHeader( vhdr, vf->frame_rate );
    }

    // Write a full second's worth of samples
    vdif_write_data( vf, out_buffer_8_vdif );

    free( out_buffer_8_vdif );
}

void vdif_write_data( struct vdifinfo *vf, int8_t *output )
{
    // form the filename
    // there is a standard naming convention
    char  filename[1030];
    sprintf( filename, "%s.vdif", vf->basefilename );

    //fprintf(stderr,"Attempting to open VDIF file for writing: %s\n",filename);
    FILE *fs = fopen( filename, "a" );
    fwrite( output, vf->block_size, 1, fs );
    fclose( fs );

    // write a CPSR2 test header for DSPSR
    char ascii_header[MWA_HEADER_SIZE] = MWA_HEADER_INIT;
    //ascii_header_set( ascii_header, "UTC_START", "%s", vf->date_obs  );
    ascii_header_set( ascii_header, "DATAFILE",   "%s", filename      );
    ascii_header_set( ascii_header, "INSTRUMENT", "%s", "VDIF"        );
    ascii_header_set( ascii_header, "TELESCOPE",  "%s", vf->telescope );
    ascii_header_set( ascii_header, "MODE",       "%s", vf->obs_mode  );
    ascii_header_set( ascii_header, "FREQ",       "%f", vf->fctr      );

    ascii_header_set( ascii_header, "BW",         "%f", vf->BW        );
    ascii_header_set( ascii_header, "RA",         "%s", vf->ra_str    );
    ascii_header_set( ascii_header, "DEC",        "%s", vf->dec_str   );
    ascii_header_set( ascii_header, "SOURCE",     "%s", vf->source    );

    sprintf( filename, "%s.hdr", vf->basefilename );
    fs = fopen( filename,"w" );
    fwrite( ascii_header, MWA_HEADER_SIZE, 1, fs );
    fclose( fs );

}


void populate_vdif_header(
        struct vdifinfo *vf,
        vdif_header     *vhdr,
        char            *metafits,
        char            *obsid,
        char            *time_utc,
        int              sample_rate,
        long int         frequency,
        int              nchan, 
        long int         chan_width,
        char            *rec_channel,
        struct delays   *delay_vals,
        int              npointing )
{
    for ( int p=0; p<npointing; p++ )
    {
        // First how big is a DataFrame
        vf[p].bits              = 8;   // this is because it is all the downstream apps support (dspsr/diFX)
        vf[p].iscomplex         = 1;   // (it is complex data)
        vf[p].nchan             = 2;   // I am hardcoding this to 2 channels per thread - one per pol
        vf[p].samples_per_frame = 128; // also hardcoding to 128 time-samples per frame
        vf[p].sample_rate       = sample_rate*128;  // = 1280000 (also hardcoding this to the raw channel rate)
        vf[p].BW                = 1.28;

        vf[p].frame_length  = (vf[p].nchan * (vf[p].iscomplex+1) * vf[p].samples_per_frame) +
                            VDIF_HEADER_SIZE;                                         // = 544
        vf[p].threadid      = 0;
        sprintf( vf[p].stationid, "mw" );

        vf[p].frame_rate = sample_rate;                                                 // = 10000
        vf[p].block_size = vf[p].frame_length * vf[p].frame_rate;                           // = 5440000

        // A single frame (128 samples). Remember vf.nchan is kludged to npol
        vf[p].sizeof_beam = vf[p].samples_per_frame * vf[p].nchan * (vf[p].iscomplex+1);      // = 512

        // One full second (1.28 million 2 bit samples)
        vf[p].sizeof_buffer = vf[p].frame_rate * vf[p].sizeof_beam;                         // = 5120000

        createVDIFHeader( vhdr, vf[p].frame_length, vf[p].threadid, vf[p].bits, vf[p].nchan,
                                vf[p].iscomplex, vf[p].stationid);

        // Now we have to add the time
        uint64_t start_day = delay_vals->intmjd;
        uint64_t start_sec = roundf( delay_vals->fracmjd * 86400.0 );
        uint64_t mjdsec    = (start_day * 86400) + start_sec; // Note the VDIFEpoch is strange - from the standard

        setVDIFEpoch( vhdr, start_day );
        setVDIFMJDSec( vhdr, mjdsec );
        setVDIFFrameNumber( vhdr, 0 );

        // Get the project ID directly from the metafits file
        fitsfile *fptr = NULL;
        int status     = 0;

        fits_open_file(&fptr, metafits, READONLY, &status);
        fits_read_key(fptr, TSTRING, "PROJECT", vf[p].exp_name, NULL, &status);
        fits_close_file(fptr, &status);

        strncpy( vf[p].scan_name, obsid, 17 );

        vf[p].b_scales   = (float *)malloc( sizeof(float) * vf[p].nchan );
        vf[p].b_offsets  = (float *)malloc( sizeof(float) * vf[p].nchan );
        vf[p].got_scales = 1;

        strncpy( vf[p].telescope, "MWA", 24);
        strncpy( vf[p].obs_mode,  "PSR", 8);

        // Determine the RA and Dec strings
        double ra2000  = delay_vals[p].mean_ra  * DR2D;
        double dec2000 = delay_vals[p].mean_dec * DR2D;

        dec2hms(vf[p].ra_str,  ra2000/15.0, 0); // 0 = no '+' sign
        dec2hms(vf[p].dec_str, dec2000,     1); // 1 = with '+' sign

        strncpy( vf[p].date_obs, time_utc, 24);

        vf[p].MJD_epoch = delay_vals->intmjd + delay_vals->fracmjd;
        vf[p].fctr      = (frequency + (nchan/2.0)*chan_width)/1.0e6; // (MHz)
        strncpy( vf[p].source, "unset", 24 );

        // The output file basename
        int ch = atoi(rec_channel);
        sprintf( vf[p].basefilename, "%s_%s_%s_%s_ch%03d",
                 vf[p].exp_name, vf[p].scan_name, vf[p].ra_str, vf[p].dec_str, ch);
    }
}


ComplexFloat get_std_dev_complex(ComplexFloat *input, int nsamples)
{
    // assume zero mean
    float rtotal = 0;
    float itotal = 0;
    float isigma = 0;
    float rsigma = 0;
    int i;

    for (i=0;i<nsamples;i++){
         rtotal = rtotal+(CRealf(input[i])*CRealf(input[i]));
         itotal = itotal+(CImagf(input[i])*CImagf(input[i]));

     }
    rsigma = sqrtf((1.0/(nsamples-1))*rtotal);
    isigma = sqrtf((1.0/(nsamples-1))*itotal);

    return CMakef( rsigma, isigma );
}

void set_level_occupancy(ComplexFloat *input, int nsamples, float *new_gain)
{
    //float percentage = 0.0;
    //float occupancy = 17.0;
    float limit = 0.00001;
    float step = 0.001;
    int i = 0;
    float gain = *new_gain;

    float percentage_clipped = 100;
    while (percentage_clipped > 0 && percentage_clipped > limit) {
        int count = 0;
        int clipped = 0;
        for (i = 0; i < nsamples; i++) {
            if (isnan(CRealf(input[i])) || isnan(CImagf(input[i])))
            {
                fprintf( stderr, "error: set_level_occupancy: input[%d] = "
                                 "NaN\n", i );
                exit(EXIT_FAILURE);
            }
            if (gain*CRealf(input[i]) >= 0 && gain*CRealf(input[i]) < 64)
            {
                count++;
            }
            if (fabs(gain*CRealf(input[i])) > 127)
            {
                clipped++;
            }
        }
        percentage_clipped = ((float) clipped/nsamples) * 100;
        if (percentage_clipped < limit) {
            gain = gain + step;
        }
        else {
            gain = gain - step;
        }
        //percentage = ((float)count/nsamples)*100.0;
        //fprintf(stdout,"Gain set to %f (linear)\n",gain);
        //fprintf(stdout,"percentage of samples in the first 64 (+ve) levels - %f percent \n",percentage);
        //fprintf(stdout,"percentage clipped %f percent\n",percentage_clipped);
    }
    *new_gain = gain;
}


void get_mean_complex( ComplexFloat *input, int nsamples, float *rmean,
                       float *imean, ComplexFloat *cmean)
{
    int i;

    float rtotal = 0;
    float itotal = 0 ;

    ComplexFloat ctotal = CMakef( 0.0, 0.0 );

    for (i = 0; i < nsamples; i++)
    {
//if (isnan(CRealf(input[i])) || isnan(CImagf(input[i]))) { fprintf(stderr, "\ninput[%d] = %e + %e*I\n\n", i, CRealf(input[i]), CImagf(input[i])); exit(1); }
        rtotal += CRealf( input[i] );
        itotal += CImagf( input[i] );
        ctotal  = CAddf( ctotal, input[i] );
    }

    *rmean = rtotal / nsamples;
    *imean = itotal / nsamples;
    *cmean = CSclf( ctotal, 1.0 / (float)nsamples );
}

void normalise_complex(ComplexFloat *input, int nsamples, float scale)
{
    int i=0;

    for (i=0;i<nsamples;i++){
        input[i] = CSclf( input[i], scale );
    }
}


void to_offset_binary(int8_t *i, int n)
{
    int j;
    for (j = 0; j < n; j++) {
        i[j] = i[j] ^ 0x80;
    }
}

#ifndef HAVE_CUDA

void invert_pfb_ifft( ComplexDouble ***detected_beam, int file_no,
                      int nsamples, int nchan, int npol,
                      float *data_buffer_vdif )
/* "Invert the PFB" by simply applying an inverse FFT.
 * This function expects "detected_beam" to be structured as follows:
 *
 *   detected_beam[2*nsamples][nchan][npol]
 *
 * Although detected_samples potentially contains 2 seconds' worth of data,
 * this function only FFTs one second. The appropriate second is worked out
 * using file_no: if it is even, the first half of detected_beam is used,
 * if odd, the second half.
 *
 * The output of the inverse FFT is packed back into data_buffer_vdif, a 1D
 * array whose ordering is as follows:
 *
 *   time, pol, complexity
 *
 * This ordering is suited for immediate output to the VDIF format.
 */
{
    // Allocate FFTW arrays
    int arr_size = nsamples * nchan * npol;
    fftwf_complex *in  = (fftwf_complex *)fftwf_malloc( arr_size * sizeof(fftwf_complex) );

    // Create a plan for doing column-wise 1D transforms
    int rank     = 1;
    int n[]      = { nchan };
    int howmany  = nsamples * npol;
    int idist    = nchan;
    int odist    = nchan;
    int istride  = 1;
    int ostride  = 1;
    int *inembed = n, *onembed = n;
    fftwf_plan p = fftwf_plan_many_dft( rank, n, howmany,
                                        in, inembed, istride, idist,
                                        in, onembed, ostride, odist,
                                        FFTW_BACKWARD, FFTW_ESTIMATE );

    // Populate the FFTW arrays such that the middle channel of detected_beam
    // is placed nearest the DC term.

    int s;    // sample index

#pragma omp parallel for
    for (s = 0; s < nsamples; s ++)
    {
        int ds, ch, pol;
        int ii;   // "in" index
        int chi;  // corrected channel index for "in" array

        // Calculate the proper sample index for this second
        ds = (file_no % 2)*nsamples + s;

        for (ch  = 0; ch  < nchan; ch++ )
        for (pol = 0; pol < npol;  pol++)
        {
            // Swap the two halves of the array
            chi = (ch < nchan/2 ? ch + (nchan/2) : ch - (nchan/2));

            // Calculate the "in" index
            ii = nchan * npol * s +
                 nchan * pol +
                 chi;

            // Copy across the data (but set DC bin to 0)
            in[ii] = (chi == 0 ? 0.0 : detected_beam[ds][ch][pol]);
        }
    }

/*
    fprintf( stderr, "  First column to be iFFT'd (inside invert_pfb_ifft()): [\n" );
    for (s = 0; s < nchan; s++)
       fprintf( stderr, " %f + %f*I\n", creal(in[s]), cimag(in[s]) );
    fprintf( stderr, "]\n" );
*/

    // Execute the FFT
    fftwf_execute( p );

    // Pack result into the output array

#pragma omp parallel for
    for (s = 0; s < nsamples; s ++)
    {
        int ch, pol;
        int ii, oi; // "in" index & "out" index

        for (ch  = 0; ch  < nchan; ch++ )
        for (pol = 0; pol < npol;  pol++)
        {
            // Calculate the "in" index
            ii = nchan * npol * s +
                 nchan * pol +
                 ch;

            // Calculate the "out" index ("ch" here turns into a subdivision
            // of time)
            oi = 2 * npol * nchan * s +
                 2 * npol * ch +
                 2 * pol;

            // Copy data across, dividing by nchan to account for the lack of
            // normalisation in the FFTW library.
            data_buffer_vdif[oi]   = crealf(in[ii]) / (double)nchan;
            data_buffer_vdif[oi+1] = cimagf(in[ii]) / (double)nchan;
        }
    }

    // Clean up
    fftwf_free( in );
    fftwf_destroy_plan( p );
}

void invert_pfb_ord( ComplexDouble ***detected_beam, int file_no,
                      int nsamples, int nchan, int npol,
                      ComplexDouble **fils, int fil_size,
                      float *data_buffer_uvdif )
/* "Invert the PFB" by applying a resynthesis filter.
 * This function expects "detected_beam" to be structured as follows:
 *
 *   detected_beam[2*nsamples][nchan][npol]
 *
 * Although detected_samples potentially contains 2 seconds' worth of data,
 * this function only inverts 1 second. The appropriate second is worked out
 * using file_no: if it is even, the first half of detected_beam is used,
 * if odd, the second half.
 *
 * The output of the inversion is packed back into data_buffer_vdif, a 1D
 * array whose ordering is as follows:
 *
 *   time, pol, complexity
 *
 * This ordering is suited for immediate output to the VDIF format.
 *
 * Finally, fils points to a 2D array of filter coefficients, each row of
 * which has been "rotated" with phase ramps of different amounts. It is
 * assumed that fils has size:
 *
 *   fils[nchan][fil_size]
 */
{
    // Set the output buffer to zeros
    int s;
#pragma omp parallel for
    for (s = 0; s < npol*nchan*nsamples*2; s++)
    {
        data_buffer_uvdif[s] = 0.0;
    }

    // Loop over (output) sample -- embarassingly parallel
#pragma omp parallel for
    for (s = 0; s < nchan*nsamples; s++)
    {
        //fprintf( stderr, "  Thread num: %d, s = %d\n", omp_get_thread_num(), s );
        int U        = nchan;        // upsampling factor = number of channels
        int i0;                      // The index of the first input sample to
                                     // be included in the output sum
        int f0;                      // The index of the first filter coeffi-
                                     // cient to be included in the output sum
        int N        = nsamples * U; // The total number of output samples
        int ch, f, i, pol, oi;       // Various loop counters
        ComplexDouble part;

        for (pol = 0; pol < npol; pol++)
        {
            // Calculate the output index for data_buffer_uvdif
            oi = 2*npol*s + 2*pol;

            // First take care of the corner case = the very first second
            if (file_no == 0 && s < fil_size - 1)
            {
                //data_buffer_uvdif[oi  ] = 0.0; // "real"
                //data_buffer_uvdif[oi+1] = 0.0; // "imag"
                continue;
            }

            // Calculate the first input idx to be included in this out sample
            if (file_no % 2 == 0)
                i0 = ((s + 2*N - fil_size + U) / U) % (2*nsamples);
            else // file_no % 2 == 1
                i0 = (s + 1*N - fil_size + U) / U;

            // Calculate the first filter coefficient index
            f0 = (U - (s % U) - 1) % U;

            // Loop over channels and filter coefficients to calculate output
            for (ch = 0; ch < nchan; ch++)
            //for (ch = 3; ch < 4; ch++)
            {
                i = i0;
                for (f = f0; f < fil_size; f += U)
                {
                    part = CMuld( fils[ch][(fil_size-1) - f], detected_beam[i][ch][pol] );
                    data_buffer_uvdif[oi  ] += CReald(part);
                    data_buffer_uvdif[oi+1] += CImagd(part);

                    // Update input index simultaneously with filter coeff
                    i++;
                    if (i == 2*nsamples)  i = 0; // (i.e. loop back around to
                                                 //  the other second)
                } // Loop over relevant filter coefficients
            } // Loop over channels

            // Normalise the result
            data_buffer_uvdif[oi  ] /= nchan;
            data_buffer_uvdif[oi+1] /= nchan;

        } // Loop over X/Y pol
    } // Loop over samples
}

#endif