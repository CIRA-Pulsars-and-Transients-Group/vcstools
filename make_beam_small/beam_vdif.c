#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "vdifio.h"
#include "psrfits.h"
#include "slamac.h"
#include "mwac_utils.h"
#include "beam_common.h"
#include "beam_vdif.h"
#include "mwa_header.h"
#include "vdifio.h"
#include "ascii_header.h"


void vdif_write_second( struct vdifinfo *vf, vdif_header *vhdr,
        float *data_buffer_vdif, float *gain )
{
    if (vf->got_scales == 0) {

        float rmean, imean;
        complex float cmean, stddev;

        get_mean_complex(
                (complex float *)data_buffer_vdif,
                vf->sizeof_buffer/2.0,
                &rmean, &imean, &cmean );

        stddev = get_std_dev_complex(
                (complex float *)data_buffer_vdif,
                vf->sizeof_buffer/2.0 );

        if (fabsf(rmean) > 0.001) {
            fprintf( stderr, "error: vdif_write_second: significantly "
                             "non-zero mean (%f)", rmean );
            exit(EXIT_FAILURE);
        }

        vf->b_scales[0] = crealf(stddev);
        vf->b_scales[1] = crealf(stddev);

        vf->got_scales = 1; // TODO: find out if this is ever meant to be reset to 0
        set_level_occupancy(
                (complex float *)data_buffer_vdif,
                vf->sizeof_buffer/2.0, gain);

    }

    normalise_complex(
            (complex float *)data_buffer_vdif,
            vf->sizeof_buffer/2.0,
            1.0/(*gain) );

    float *data_buffer_ptr = data_buffer_vdif;
    size_t offset_out_vdif = 0;

    int8_t *out_buffer_8_vdif = (int8_t *)malloc(vf->block_size);

    while  (offset_out_vdif < vf->block_size) {

        memcpy( (out_buffer_8_vdif + offset_out_vdif), vhdr, VDIF_HEADER_SIZE ); // add the current header
        offset_out_vdif += VDIF_HEADER_SIZE; // offset into the output array

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
        struct delays   *delay_vals )
{
    // First how big is a DataFrame
    vf->bits              = 8;   // this is because it is all the downstream apps support (dspsr/diFX)
    vf->iscomplex         = 1;   // (it is complex data)
    vf->nchan             = 2;   // I am hardcoding this to 2 channels per thread - one per pol
    vf->samples_per_frame = 128; // also hardcoding to 128 time-samples per frame
    vf->sample_rate       = sample_rate*128;  // also hardcoding this to the raw channel rate
    vf->BW                = 1.28;

    vf->frame_length  = (vf->nchan * (vf->iscomplex+1) * (vf->bits) * vf->samples_per_frame) +
                        VDIF_HEADER_SIZE;
    vf->threadid      = 0;
    sprintf( vf->stationid, "mw" );

    vf->frame_rate = sample_rate;
    vf->block_size = vf->frame_length * vf->frame_rate;

    // A single frame (128 samples). Remember vf.nchan is kludged to npol
    vf->sizeof_beam = vf->samples_per_frame * vf->nchan * (vf->iscomplex+1);

    // One full second (1.28 million 2 bit samples)
    vf->sizeof_buffer = vf->frame_rate * vf->sizeof_beam;

    createVDIFHeader( vhdr, vf->frame_length, vf->threadid, vf->bits, vf->nchan,
                            vf->iscomplex, vf->stationid);

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
    fits_read_key(fptr, TSTRING, "PROJECT", vf->exp_name, NULL, &status);
    fits_close_file(fptr, &status);

    strncpy( vf->scan_name, obsid, 17 );

    vf->b_scales   = (float *)malloc( sizeof(float) * vf->nchan );
    vf->b_offsets  = (float *)malloc( sizeof(float) * vf->nchan );
    vf->got_scales = 0;

    strncpy( vf->telescope, "MWA", 24);
    strncpy( vf->obs_mode,  "PSR", 8);

    // Determine the RA and Dec strings
    double ra2000  = delay_vals->mean_ra  * DR2D;
    double dec2000 = delay_vals->mean_dec * DR2D;

    dec2hms(vf->ra_str,  ra2000/15.0, 0); // 0 = no '+' sign
    dec2hms(vf->dec_str, dec2000,     1); // 1 = with '+' sign

    strncpy( vf->date_obs, time_utc, 24);

    vf->MJD_epoch = delay_vals->intmjd + delay_vals->fracmjd;
    vf->fctr      = (frequency + (nchan/2.0)*chan_width)/1.0e6; // (MHz)
    strncpy( vf->source, "unset", 24 );

    // The output file basename
    int ch = atoi(rec_channel);
    sprintf( vf->basefilename, "%s_%s_ch%03d",
             vf->exp_name, vf->scan_name, ch);
}


complex float get_std_dev_complex(complex float *input, int nsamples)
{
    // assume zero mean
    float rtotal = 0;
    float itotal = 0;
    float isigma = 0;
    float rsigma = 0;
    int i;

    for (i=0;i<nsamples;i++){
         rtotal = rtotal+(crealf(input[i])*crealf(input[i]));
         itotal = itotal+(cimagf(input[i])*cimagf(input[i]));

     }
    rsigma = sqrtf((1.0/(nsamples-1))*rtotal);
    isigma = sqrtf((1.0/(nsamples-1))*itotal);

    return rsigma+I*isigma;
}

void set_level_occupancy(complex float *input, int nsamples, float *new_gain)
{
    float percentage = 0.0;
    //float occupancy = 17.0;
    float limit = 0.00001;
    float step = 0.001;
    int i = 0;
    float gain = *new_gain;

    float percentage_clipped = 100;
    while (percentage_clipped > 0 && percentage_clipped > limit) {
        int count = 0;
        int clipped = 0;
        for (i=0;i<nsamples;i++) {
            if (gain*creal(input[i]) >= 0 && gain*creal(input[i]) < 64) {
                count++;
            }
            if (fabs(gain*creal(input[i])) > 127) {
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
        percentage = ((float)count/nsamples)*100.0;
        fprintf(stdout,"Gain set to %f (linear)\n",gain);
        fprintf(stdout,"percentage of samples in the first 64 (+ve) levels - %f percent \n",percentage);
        fprintf(stdout,"percentage clipped %f percent\n",percentage_clipped);
    }
    *new_gain = gain;
}


void get_mean_complex(complex float *input, int nsamples, float *rmean,float *imean, complex float *cmean)
{

    int i=0;
    float rtotal = 0;
    float itotal = 0 ;
    complex float ctotal = 0 + I*0.0;
    for (i=0;i<nsamples;i++){
        rtotal = rtotal+crealf(input[i]);
        itotal = itotal+cimagf(input[i]);
        ctotal = ctotal + input[i];
    }
    *rmean=rtotal/nsamples;
    *imean=itotal/nsamples;
    *cmean=ctotal/nsamples;

}

void normalise_complex(complex float *input, int nsamples, float scale)
{
    int i=0;

    for (i=0;i<nsamples;i++){
        input[i]=input[i]*scale;
    }
}


void to_offset_binary(int8_t *i, int n)
{
    int j;
    for (j = 0; j < n; j++) {
        i[j] = i[j] ^ 0x80;
    }
}

void invert_pfb_ifft( complex float *array, int nchan,
                      fftwf_plan p, fftwf_complex *in, fftwf_complex *out )
/* "Invert the PFB" by simply applying an inverse FFT.
 * This function expects "array" to contain at least nchan elements. It also
 * expects that the arrays in the plan also contain at least nchan elements.
 * The output of the inverse FFT is packed back into "array".
 * "in" and "out" must be the same arrays as those used to make the plan "p"
 */
{
    // Populate the FFTW arrays
    int ch;
    for (ch = 0; ch < nchan; ch++)
    {
        if (ch < nchan/2)
        {
            // these are the negative frequencies
            // pack them in the second half of the array
            // skip over the edges
            in[(nchan/2) + ch] = array[ch];
        }
        else if (ch > nchan/2)
        {
            // positive frequencies -- shift them to the first half
            in[ch-(nchan/2)] = array[ch];
        }
        else // if (ch == nchan/2)
        {
            // Nyquist bin - give it a zero mean
            in[0] = 0.0;
        }
    }

    // Make it so!
    fftwf_execute( p );

    // Pack result into the output array
    for (ch = 0; ch < nchan; ch++)
        array[ch] = out[ch];

}

