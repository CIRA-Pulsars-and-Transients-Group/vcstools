/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>

extern "C" {
#include "mycomplex.h"
#include "ipfb.h"
}

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    /* Wrapper function for GPU/CUDA error handling. Every CUDA call goes through
       this function. It will return a message giving your the error string,
       file name and line of the error. Aborts on error. */

    if (code != 0)
    {
        fprintf(stderr, "GPUAssert:: %s - %s (%d)\n", cudaGetErrorString(code), file, line);
        if (abort)
        {
            exit(code);
        }
    }
}

// define a macro for accessing gpuAssert
#define gpuErrchk(ans) {gpuAssert((ans), __FILE__, __LINE__, true);}

__global__ void filter_kernel( float *in_real, float *in_imag,
                               float *ft_real, float *ft_imag,
                               int ntaps, int npol, float *out )
/* This kernel computes the synthesis filter:

              1              K-1
   xhat[n] = --- SUM f[n-mM] SUM X_k[m] e^(2πjkn/K)
              K   m          k=0

   The sum over m is nominally over all integers, but in practice only
   involves a few terms because of the finiteness of the filter, f. To be
   precise, there are precisely ntaps non-zero values.

   X_k[m] represents the complex-valued inputs, in_real and in_imag.
   Every possible value of f[n]*e^(2πjkn/K) is provided in ft_real and
   ft_imag.

   K is the number of channels, and because this is a critically sampled
   PFB, M = K. We will also use P to mean the number of taps in the synthesis
   filter, and N = KP to mean the size of the filter.

   The polarisations are computed completely independently.

   And, of course, xhat[n] is represented by the out array.

   Here, we don't perform the 1/K scaling at the very end -- this is left to
   the caller to do, if desired.
 */
{
    // First, set a generic variable for this thread
    int idx = blockDim.x*blockIdx.x + threadIdx.x;

    // The polarisation for this thread is
    int pol = idx % npol;

    // and the time index of the output (i.e. the index for xhat) is
    int n = idx / npol;

    // Other constants we'll need are:
    int K = blockDim.x / npol;  // Total number of channels (should be 128)
    int M = K;
    int P = ntaps;
    int F = P*K;

    // Because we must have 0 <= n-mM < F, the smallest allowed value of m
    // is:
    int m0 = (n - F)/M + 1;

    // out[] must also include both polarisations, so its indices are simply
    int o_real = 2*idx;
    int o_imag = o_real + 1;

    // Initialise the output sample to zero
    out[o_real] = 0.0;
    out[o_imag] = 0.0;

    // Perform the double sum
    int m, k, f, tw, ft, i;
    for (m = m0; m < m0 + P; m0++)
    {
        // With m now known, we can get the index for the filter
        f = n - m*M;

        for (k = 0; k < K; k++)
        {
            // The index for the twiddle factor is
            tw = (k*n) % K;

            // The index into the ft (= filter/twiddle) array is
            ft = F*tw + f;

            // The "in" index
            // (See cu_invert_pfb_ord() for how the in[] arrays were packed)
            i = npol*K*m + npol*k + pol;

            // Complex multiplication
            out[o_real] += in_real[i] * ft_real[ft] -
                           in_imag[i] * ft_imag[ft];
            out[o_imag] += in_real[i] * ft_imag[ft] +
                           in_imag[i] * ft_real[ft];
        }
    }

    __syncthreads();
}

extern "C"
void cu_invert_pfb_ord( ComplexDouble ***detected_beam, int file_no,
                        int nsamples, int nchan, int npol,
                        struct gpu_ipfb_arrays *g, float *data_buffer_uvdif )
/* "Invert the PFB" by applying a resynthesis filter, using GPU
 * acceleration.
 *
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
 * It is assumed that the inverse filter coefficients have already been loaded
 * to the GPU.
 */
{
    // Setup input values:
    // The starting sample index is "ntaps" places from the end of the second
    // half of detected_beam if the file number is even, and "ntaps" places
    // from the end of the first half of detected_beam if the file number is
    // odd.
    int start_s = (file_no % 2 == 0 ? 2*nsamples - g->ntaps : nsamples - g->ntaps);

    int s_in, s, ch, pol, i;
    for (s_in = 0; s_in < nsamples + g->ntaps; s_in++)
    {
        s = (start_s + s_in) % (2*nsamples);
        for (ch = 0; ch < nchan; ch++)
        {
            for (pol = 0; pol < npol; pol++)
            {
                // Calculate the index for in_real and in_imag;
                i = npol*nchan*s_in + npol*ch + pol;

                // Copy the data across - taking care of the file_no = 0 case
                if (file_no == 0 && s_in < g->ntaps)
                {
                    g->in_real[i] = 0.0;
                    g->in_imag[i] = 0.0;
                }
                else
                {
                    g->in_real[i] = CReald( detected_beam[s][ch][pol] );
                    g->in_imag[i] = CImagd( detected_beam[s][ch][pol] );
                }
            }
        }
    }

    // Copy the data to the device
    gpuErrchk(cudaMemcpy( g->d_in_real, g->in_real, g->in_size, cudaMemcpyHostToDevice ));
    gpuErrchk(cudaMemcpy( g->d_in_imag, g->in_imag, g->in_size, cudaMemcpyHostToDevice ));

    // Call the kernel
    filter_kernel<<<nsamples, nchan*npol>>>( g->d_in_real, g->d_in_imag,
                                             g->d_ft_real, g->d_ft_imag,
                                             g->ntaps, npol, g->d_out );
    cudaDeviceSynchronize();

    // Copy the result back into host memory
    gpuErrchk(cudaMemcpy( data_buffer_uvdif, g->d_out, g->out_size, cudaMemcpyDeviceToHost ));
}


void cu_load_filter( double *coeffs, ComplexDouble *twiddles, struct gpu_ipfb_arrays *g,
        int nchan )
/* This function loads the inverse filter coefficients and the twiddle factors
   into GPU memory. If they were loaded separately (as floats), then the
   multiplication of the filter coefficients and the twiddle factors will be
   less precise than if a single array containing every combination of floats
   and twiddle factors is calculated in doubles, and then demoted to floats.
   Hence, this pre-calculation is done in this function before cudaMemcpy is
   called.

   The result is 2x 1D arrays loaded onto the GPU (one for real, one for imag)
   where the ith element is equal to

   result[i] = f[n] * exp(2πjk/K),
   n = i % N  (N is the filter size, "fil_size")
   k = i / N
   and K is the number of channels (nchan).

   This function assumes that the size of resulting array has already been
   calculated (in bytes) and set in
     g->ft_size
   and that the number of elements in twiddles is
     nchan
*/
{
    int ch, f, i;
    int fil_size = g->ft_size / nchan / sizeof(float);

    // Setup filter values:
    ComplexDouble ft; // pre-calculated filter coeffs times twiddle factor
    ComplexDouble cf; // temp variable for complex version of filter coeffs
    for (f = 0; f < fil_size; f++)
    {
        cf = CMaked( coeffs[f], 0.0 );
        for (ch = 0; ch < nchan; ch++)
        {
            i = fil_size*ch + f;
            ft = CMuld( twiddles[ch], cf );
            g->ft_real[i] = CReald( ft );
            g->ft_imag[i] = CImagd( ft );
        }
    }

    gpuErrchk(cudaMemcpy( g->d_ft_real, g->ft_real, g->ft_size, cudaMemcpyHostToDevice ));
    gpuErrchk(cudaMemcpy( g->d_ft_imag, g->ft_imag, g->ft_size, cudaMemcpyHostToDevice ));
}


void malloc_ipfb( struct gpu_ipfb_arrays *g, int ntaps, int nsamples,
        int nchan, int npol, int fil_size )
{
    // Flatten the input array (detected_array) for GPU.
    // We only need one second's worth, plus 12 time samples tacked onto the
    // beginning (from the previous second)

    g->ntaps     = ntaps;
    g->in_size   = ((nsamples + ntaps) * nchan * npol) * sizeof(float);
    g->ft_size   = fil_size * nchan * sizeof(float);
    g->out_size  = nsamples * nchan * npol * 2 * sizeof(float);

    // Allocate memory on the device
    gpuErrchk(cudaMalloc( (void **)&g->d_in_real, g->in_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_in_imag, g->in_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_ft_real, g->ft_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_ft_imag, g->ft_size ));

    gpuErrchk(cudaMalloc( (void **)&g->d_out, g->out_size ));

    // Allocate memory for host copies of the same
    g->in_real = (float *)malloc( g->in_size );
    g->in_imag = (float *)malloc( g->in_size );
    g->ft_real = (float *)malloc( g->ft_size );
    g->ft_imag = (float *)malloc( g->ft_size );

}


void free_ipfb( struct gpu_ipfb_arrays *g )
{
    // Free memory on host and device
    free( g->in_real );
    free( g->in_imag );
    free( g->ft_real );
    free( g->ft_imag );
    cudaFree( g->d_in_real );
    cudaFree( g->d_in_imag );
    cudaFree( g->d_ft_real );
    cudaFree( g->d_ft_imag );
    cudaFree( g->d_out );
}
