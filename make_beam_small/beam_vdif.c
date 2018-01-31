#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <omp.h>
#include "vdifio.h"
#include "psrfits.h"
#include "slamac.h"
#include "mwac_utils.h"
#include "beam_common.h"
#include "beam_vdif.h"
#include "mwa_header.h"
#include "vdifio.h"
#include "ascii_header.h"
#include "filter.h"


void vdif_write_second( struct vdifinfo *vf, vdif_header *vhdr,
        float *data_buffer_vdif, float *gain )
{

    // Set level occupancy
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
        printf( "warning: vdif_write_second: significantly "
                "non-zero mean (%f), adjusting data\n", rmean );
        unsigned int i;
        for (i = 0; i < vf->sizeof_buffer/2; i++)
        {
            data_buffer_vdif[2*i+0] -= creal(cmean);
            data_buffer_vdif[2*i+1] -= cimag(cmean);
        }
    }

    vf->b_scales[0] = crealf(stddev);
    vf->b_scales[1] = crealf(stddev);

    vf->got_scales = 1;
    set_level_occupancy(
            (complex float *)data_buffer_vdif,
            vf->sizeof_buffer/2.0, gain);

    // Normalise
    normalise_complex(
            (complex float *)data_buffer_vdif,
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
    vf->sample_rate       = sample_rate*128;  // = 1280000 (also hardcoding this to the raw channel rate)
    vf->BW                = 1.28;

    vf->frame_length  = (vf->nchan * (vf->iscomplex+1) * vf->samples_per_frame) +
                        VDIF_HEADER_SIZE;                                         // = 544
    vf->threadid      = 0;
    sprintf( vf->stationid, "mw" );

    vf->frame_rate = sample_rate;                                                 // = 10000
    vf->block_size = vf->frame_length * vf->frame_rate;                           // = 5440000

    // A single frame (128 samples). Remember vf.nchan is kludged to npol
    vf->sizeof_beam = vf->samples_per_frame * vf->nchan * (vf->iscomplex+1);      // = 512

    // One full second (1.28 million 2 bit samples)
    vf->sizeof_buffer = vf->frame_rate * vf->sizeof_beam;                         // = 5120000

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
    vf->got_scales = 1;

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
        //percentage = ((float)count/nsamples)*100.0;
        //fprintf(stdout,"Gain set to %f (linear)\n",gain);
        //fprintf(stdout,"percentage of samples in the first 64 (+ve) levels - %f percent \n",percentage);
        //fprintf(stdout,"percentage clipped %f percent\n",percentage_clipped);
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


void invert_pfb_ifft( complex float ***detected_beam, int file_no,
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

            // Calculate the "out" index ("ch" here turns into a subvidivion
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


void invert_pfb_ord( complex float ***detected_beam, int file_no,
                      int nsamples, int nchan, int npol,
                      filter fils[], float *data_buffer_uvdif )
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
 * Finally, fils points to an array of filter coefficients, each row of which
 * has been "rotated" with phase ramps of different amounts.
 */
{
    // Set the output buffer to zeros
    int s;
#pragma omp parallel for
    for (s = 0; s < 2*npol*nchan*nsamples; s++)
    {
        data_buffer_uvdif[s] = 0.0;
    }

    // Loop over (output) sample -- embarassingly parallel
#pragma omp parallel for
    for (s = 0; s < nchan*nsamples; s++)
    {
        //fprintf( stderr, "  Thread num: %d, s = %d\n", omp_get_thread_num(), s );
        int U        = nchan;        // upsampling factor = number of channels
        int fil_size = fils[0].size; // All filters should have the same size
        int i0;                      // The index of the first input sample to
                                     // be included in the output sum
        int f0;                      // The index of the first filter coeffi-
                                     // cient to be included in the output sum
        int N        = nsamples * U; // The total number of output samples
        int ch, f, i, pol, oi;       // Various loop counters
        complex float part;

        for (pol = 0; pol < npol; pol++)
        {
            // Calculate the output index for data_buffer_uvdif
            oi = 2 * npol * s +
                 2 * pol;

            // First take care of the corner case = the very first second
            if (file_no == 0 && s < fil_size - 1)
            {
                data_buffer_uvdif[oi] = 0.0;
                continue;
            }

            // Calculate the first input idx to be included in this out sample
            if (file_no % 2 == 0)
                i0 = ((s + 2*N - fil_size + U - 1) / U) % (2*nsamples);
            else // file_no % 2 == 1
                i0 = (s + 1*N - fil_size + U - 1) / U;

            // Calculate the first filter coefficient index
            f0 = (U - (s % U)) % U;

            // Loop over channels and filter coefficients to calculate output
            for (ch = 0; ch < nchan; ch++)
            {
                i = i0;
                for (f = f0; f < fil_size; f += U)
                {
                    part = fils[ch].coeffs[f] * detected_beam[i][ch][pol];
                    data_buffer_uvdif[oi  ] += creal(part);
                    data_buffer_uvdif[oi+1] += cimag(part);

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

/*********************
 * Testing functions *
 *********************/

int test_invert_pfb_ifft()
{
    int test_success = 1;

    // Set up arrays
    int s, c, p, i; // counters
    int nsamples = 5;
    int nchan    = 4;
    int npol     = 2;

    int noutsamples = nsamples * nchan * npol * 2; // The "2" is for complexity

    complex float ***in;
    in = (complex float ***)malloc( 2*nsamples * sizeof(complex float **) ); // The "2" here is for nfiles
    for (s = 0; s < 2*nsamples; s++)
    {
        in[s] = (complex float **)malloc( nchan * sizeof(complex float *) );
        for (c = 0; c < nchan; c++)
        {
            in[s][c] = (complex float *)malloc( npol * sizeof(complex float) );
        }
    }

    float *out = (float *)malloc( noutsamples * sizeof(float) );
    complex float in1D[] = { 0.06449101+0.47745452*I,  0.25715001-0.64475143*I,
                             0.51664289+1.76467373*I, -1.44484082+1.03794392*I,
                            -0.41289831+0.27941595*I,  0.3068476 +1.37876269*I,
                            -1.9354946 +0.32219191*I, -0.42653615-1.69119755*I,
                             2.2170282 +0.80472379*I,  0.12864692-0.47128854*I,
                            -1.16028624-0.37376863*I, -0.38224916-0.93524965*I,
                            -1.22904786-0.76877662*I,  0.78255331-0.61905262*I,
                            -0.39498729+1.78955492*I,  1.06901631-1.11650729*I,
                            -0.0222534 +0.52556743*I, -0.36377478-1.1323636 *I,
                            -0.06430987+0.21966019*I,  2.03954773-0.79814271*I,
                            -0.66810883+0.0962631 *I, -1.19963108-0.72895494*I,
                            -0.01451912+0.07647457*I, -0.24382818+0.95664832*I,
                             0.12938902-0.30502573*I,  1.1381252 -0.21292474*I,
                            -1.19407573-0.77315504*I,  1.17001473+0.30662279*I,
                             0.56896738-1.64906108*I,  0.75396181+1.4724647 *I,
                             0.53814834-1.16005429*I, -0.53033866-1.21012213*I,
                             0.93440265+0.15479393*I,  0.14508885-0.08175225*I,
                             2.34790655+0.41494762*I,  2.38678981+0.5477534 *I,
                             2.01970984+2.30349442*I, -0.23616837+1.6253156 *I,
                             1.37385358+0.96571704*I, -1.09512013+2.26146954*I,
                             0.22483308+0.2637857 *I, -0.010462  -0.25085408*I,
                            -0.14401447+0.32246797*I, -2.98975638-0.72696144*I,
                             0.60517084+0.00792186*I,  1.694062  +1.23106275*I,
                             0.67741499-0.18467748*I,  1.50002306-0.49210533*I,
                             0.70718233-0.1475476 *I,  0.21912304+0.49289078*I,
                            -1.38471223+0.47236468*I, -0.17955162-1.17805713*I,
                             0.51350124-1.77279338*I,  0.55928752+1.67283547*I,
                             1.26361809+0.25470412*I, -1.95402836-0.26665733*I,
                            -0.68856163+0.73175239*I,  1.45928063+0.005512  *I,
                            -1.5173021 -1.72574429*I, -0.40862674+0.83482037*I,
                             1.57055389+0.51548405*I, -0.31160221+0.80255972*I,
                             0.79664571-0.23661208*I, -0.21066871-1.22762975*I,
                             1.65369268-1.46882212*I,  0.12355657+0.32145898*I,
                            -1.30112467+0.7759431 *I, -0.37798103-0.52303737*I,
                            -1.65956614+1.46966904*I,  2.05323508-0.98037911*I,
                             0.78154538-0.455914  *I, -0.689024  +0.92592841*I,
                            -1.69215217+0.33113928*I, -0.81932695+1.16758996*I,
                             1.46291286-2.12344018*I,  1.09994985-1.48012812*I,
                             0.66778226+0.46633416*I,  1.30965981+1.03930631*I,
                             0.65389738-0.29542409*I,  1.02400447+0.37062682*I };

    float answer[] = {
        -0.2807999477632719, 0.21766175477734398, -0.13314671627421729, -0.008254728810845727,
        -0.43679659049406105, -0.5448613376151576, 0.8707146363974354, -0.2531999225033581,
        0.3130454544049577, 0.021065506591868982, -0.07330244099319, 0.14796270152238572,
        0.40455108385237526, 0.30613407624594474, -0.6642654791300282, 0.1134919497918181,
        0.4908564905519233, -0.15045360028479235, 0.15563044145654487, -0.626084131743714,
        -0.4382667741214965, -0.3289049656857378, 0.43162563185324815, 0.26380990478809074,
        0.6176576095197633, 0.552815494953696, -0.770154369074445, 0.24169582249586719,
        -0.6702473259501902, -0.07345692898316583, 0.18289829576465194, 0.1205784044597561,
        0.4133798892328422, -0.35123472201809364, -0.5278920243598099, 0.08098911671684547,
        -0.07799187244419897, 0.4694387708809175, -0.25437360684907573, 0.21488495081193662,
        -0.4245065876277181, 0.6140184361870462, 0.19383760861849197, -0.03285756890669147,
        0.08911857083907493, -0.73222248504987, 0.5884280225903936, -0.26301649862209064,
        0.6093822388480371, -0.05283192079945933, 0.19814763402090865, -0.34667962657334295,
        -0.1622341384218563, 0.0842288140753566, 0.5284048614445283, 0.09119015366913336,
        -0.5446877266503222, -0.09968094456870844, 0.08633605642588385, -0.4778509147592732,
        0.09753962622414143, 0.06828405129281118, -0.8128885518913207, 0.7333403876634828,
        0.86657032833051, 0.15519876924511317, 0.17210533317835713, 1.5475698906460242,
        -0.39097707565514134, 0.5217267562424636, -0.6639659436520617, -0.790611543652574,
        -0.39936900343713755, -0.07780180281795758, 0.8377495849591154, -0.39582268015858124,
        -0.07622424923823112, -0.5991237226696192, -0.3458889744854108, -0.3611356668348691,
        -0.6938463235936019, -0.17850745682795408, 0.9498139739946144, 0.18671982091661202,
        0.062818569580091, -0.8107700219539551, 0.27949930847891546, -0.05049019887142139,
        0.8062628646273196, 0.31040030776670413, -0.6472285528213756, -0.18275889210432347,
        -0.1752351106138086, 0.6788771710152051, -0.5820847296521542, 0.04652927005913285,
        0.18668843705522598, -0.2081784876207473, -0.22030990039413167, -0.09165381200038308,
        0.24094139613451562, -0.06278176535369435, 0.3564978909112116, -0.18513062217636156,
        0.16690272596740224, 0.1344046896676525, 0.4770605191073424, -0.7947428791857158,
        -0.5945325591571439, 0.13655556330678914, -0.6132485096244222, 1.0715273133624605,
        0.09052306390714529, 0.3930211912718532, 0.2620707415326092, 0.0226035054235785,
        -0.03518668499549801, -0.6499149408556733, 0.11490889464080073, -0.10363763803837592,
        -0.4348038767051533, -0.02714499455251776, 0.523206203092407, 0.23513852148970785,
        0.37946749779350597, 0.2840387441363378, -0.9001858392658169, -0.15410438887491043,
        0.3498170563565887, -0.4176001290938339, -0.0738387644691415, 0.3538045819616782,
        -0.20229908201731714, 0.24182113136415367, -0.061685343761384404, -1.0529820289082108,
        0.477029283700344, -0.316810933347468, -0.755944307946986, 0.381029936341159,
        -0.6245472580396155, 0.49258993107714827, 0.8914684161775119, 0.31814751060537355,
        -0.3528823198000375, 0.004650281050137414, 0.7503616354000141, 0.4690668221692462,
        1.084967563660022, 0.397034381612424, 0.00022430848263746417, -0.18799737439440709,
        -0.49319376721137786, 0.1609193568887628, -0.4164705066776212, -0.2358997404721381,
        -0.2388914766486066, -0.5626040195513242, -0.33411543720503034, -0.04516970730270106 };

    // Set up input data
    for (s = 0; s < 2*nsamples; s++)
    for (p = 0; p < npol;       p++)
    for (c = 0; c < nchan;      c++)
    {
        i = nchan*npol*s + nchan*p + c;
        in[s][c][p] = in1D[i];
    }

    // Run it through invert_pfb_ifft() (for each "file")
    int f;
    for (f = 0; f < 2; f++)
    {
        //printf( "     Answer,  calculated\n" );
        invert_pfb_ifft( in, f, nsamples, nchan, npol, out );

        for (s = 0; s < noutsamples; s++)
        {
            i = noutsamples * f + s;
            //printf( "     %.12e,  %.12e\n", answer[i], out[s] );
            if (fabs( answer[i] - out[s] ) > 1.0e-6) // Seems to only succeed to single precision...
            {
                test_success = 0;
                break;
            }
        }
    }

    // Free memory
    for (s = 0; s < 2*nsamples; s++)
    {
        for (c = 0; c < nchan; c++)
        {
            free( in[s][c] );
        }
        free( in[s] );
    }
    free( in );
    free( out );

    return test_success;
}
