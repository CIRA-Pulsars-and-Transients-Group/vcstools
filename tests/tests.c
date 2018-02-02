#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>
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

    int f, s, c, p; // Some useful loop counters

    // Set up the input array
    int nsamples = 46;
    int nchan    = 4;
    int npol     = 2;

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

    // Pack the input values into the detected_beam array. Artificially make
    // a second polarisation that is 180 deg out of phase with the original
    // polarisation.
    for (s = 0; s < 2*nsamples; s++)
    {
        // First polarisation
        p = 0;
        detected_beam[s][0][p] = input_ch2[s];
        detected_beam[s][1][p] = input_ch3[s];
        detected_beam[s][2][p] = input_ch0[s];
        detected_beam[s][3][p] = input_ch1[s];
        // Second polarisation
        p = 1;
        detected_beam[s][0][p] = -input_ch2[s];
        detected_beam[s][1][p] = -input_ch3[s];
        detected_beam[s][2][p] = -input_ch0[s];
        detected_beam[s][3][p] = -input_ch1[s];
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
    int nfiles = 2;
    float data_buffer_uvdif[nfiles][2*nsamples*npol*nchan];

    // Set up the answer array
    float answer[] = INVERT_PFB_ORD_OUT;

    // Test the invert_pfb_ord() function!
    for (f = 0; f < nfiles; f++)
    {
        invert_pfb_ord( detected_beam, f, nsamples, nchan, npol,
                        fil_ramps, data_buffer_uvdif[f] );
    }

    // Compare the results
    int ai, di; // indexes for answer and data_buffer_uvdif respectively
    int x;      // loop counter for complexity
    float ans, calc;
    for (f = 0; f < nfiles; f++)
    for (s = 0; s < nsamples*nchan; s++)
    for (p = 0; p < npol; p++)
    for (x = 0; x < 2; x++)
    {
        ai = 2*nsamples*nchan*f + 2*s + x;
        di = 2*npol*s + 2*p + x;

        ans = (p == 0 ? answer[ai] : -answer[ai]);
        calc = data_buffer_uvdif[f][di];

        if (fabs( ans - calc ) > 1e-5)
        {
            test_success = 0;
        }
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
