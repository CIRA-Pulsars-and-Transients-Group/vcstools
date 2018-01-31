#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "beam_vdif.h"
#include "tests.h"

int main()
{
    // TEST: invert_pfb_ifft()
    if (test_invert_pfb_ifft())
        printf( "invert_pfb_ifft(): test successful\n" );
    else
        printf( "invert_pfb_ifft(): test FAILED\n" );
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

    return test_success;
}
