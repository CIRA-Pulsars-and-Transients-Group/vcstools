#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "beam_vdif.h"
#include "filter.h"
#include "tests.h"

int main()
{
    // TEST: invert_pfb_ifft()
    if (test_invert_pfb_ifft())
        printf( "# invert_pfb_ifft(): test successful\n" );
    else
        printf( "# invert_pfb_ifft(): test FAILED\n" );

    // TEST: invert_pfb_ord()
    if (test_invert_pfb_ord())
        printf( "# invert_pfb_ord(): test successful\n" );
    else
        printf( "# invert_pfb_ord(): test FAILED\n" );

}

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
    complex float in1D[]   = INVERT_PFB_IFFT_IN;
    float         answer[] = INVERT_PFB_IFFT_OUT;

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


int test_invert_pfb_ord()
{
    int test_success = 1;

    int s, c, p; // Some usefule loop counters

    // Set up the input array
    int nsamples = 46;
    int nchan    = 4;
    int npol     = 1;

    complex float input_ch0[] = INVERT_PFB_ORD_IN_CH0;
    complex float input_ch1[] = INVERT_PFB_ORD_IN_CH1;
    complex float input_ch2[] = INVERT_PFB_ORD_IN_CH2;
    complex float input_ch3[] = INVERT_PFB_ORD_IN_CH3;

    complex float ***detected_beam;
    detected_beam = (complex float ***)malloc( 2*nsamples * sizeof(complex float **) );
    for (s = 0; s < 2*nsamples; s++)
    {
        detected_beam[s] = (complex float **)malloc( nchan * sizeof(complex float *) );
        for (c = 0; c < nchan; c++)
        {
            detected_beam[s][c] = (complex float *)malloc( npol * sizeof(complex float) );
        }
    }

    for (s = 0; s < 2*nsamples; s++)
    for (p = 0; p < npol      ; p++)
    {
        detected_beam[s][0][p] = input_ch2[s];
        detected_beam[s][1][p] = input_ch3[s];
        detected_beam[s][2][p] = input_ch0[s];
        detected_beam[s][3][p] = input_ch1[s];
    }

    // Set up the filters
    filter fil;
    fil.ntaps  = 8;
    fil.size   = nchan * fil.ntaps;
    fil.coeffs = (complex double *)malloc( fil.size * sizeof(complex double) );

    double filter_coeffs[] = INVERT_PFB_ORD_FILTER_COEFFS;

    for (s = 0; s < fil.size; s++)
    {
        fil.coeffs[s] = filter_coeffs[s];
    }

    filter fil_ramps[nchan];
    for (c = 0; c < nchan; c++)
        create_filter( &(fil_ramps[c]), fil.size, fil.ntaps );
    apply_mult_phase_ramps( &fil, nchan, fil_ramps );

    // Set up the output array
    float data_buffer_uvdif[2*nsamples*npol*nchan*2];

    // Set up the answer array
    float answer[] = INVERT_PFB_ORD_OUT;

    // (I forgot to record the normalised data in tests.h, so I have to
    // normalise the answers here, now)
    //for (s = 0; s < nsamples*nchan*npol*2*2; s++)
    //    answer[s] /= nchan;

    // Test the invert_pfb_ord() function!
    int nfiles = 2;
    int f;
    for (f = 0; f < nfiles; f++)
    {
        invert_pfb_ord( detected_beam, f, nsamples, nchan, npol,
                        fil_ramps, data_buffer_uvdif );
    }

    // A mini test--does fir_filter_1D() work as expected?
    complex double ch0[2*nsamples*nchan];
    for (s = 0; s < 2*nsamples*nchan; s++)
    {
        if (s % nchan == 0)
            ch0[s] = (complex double)input_ch2[s/nchan];
        else
            ch0[s] = 0.0 + 0.0*I;
    }
    complex double inv_channel[2*nsamples*nchan];
    fir_filter_1D( &(fil_ramps[0]), ch0, 2*nsamples*nchan, inv_channel );
    for (s = 0; s < 2*nsamples*nchan; s++)
        inv_channel[s] /= nchan;

    // Compare the results
    printf( "#         ch0                   data_buffer_uvdif               inv_channel\n" );
    printf( "#---------------------        ---------------------        ---------------------\n" );
    for (s = 0; s < 4*nsamples*npol*nchan; s++)
    {
        if (fabs( answer[s] - data_buffer_uvdif[s] ) > 1e-6)
        {
            test_success = 0;
            //break;
        }
        //printf( "%f  %f\n", answer[s], data_buffer_uvdif[s] );
        if (s%2 == 0)
            printf( "%10f  %10f  -->  %10f  %10f   ==? %10f  %10f\n",
                    creal(ch0[s/2]), cimag(ch0[s/2]),
                    data_buffer_uvdif[s], data_buffer_uvdif[s+1],
                    creal(inv_channel[s/2]), cimag(inv_channel[s/2]) );
    }

    // Free memory
    destroy_filter( &fil );

    for (c = 0; c < nchan; c++)
        destroy_filter( &(fil_ramps[c]) );

    for (s = 0; s < 2*nsamples; s++)
    {
        for (c = 0; c < nchan; c++)
            free( detected_beam[s][c] );
        free( detected_beam[s] );
    }
    free( detected_beam );

    return test_success;
}
