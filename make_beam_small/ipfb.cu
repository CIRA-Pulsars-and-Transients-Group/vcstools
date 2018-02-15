#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>

extern "C" {
#include "mycomplex.h"
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

__global__ void filter_kernel( float   *in_real, float   *in_imag,
                               float *fils_real, float *fils_imag,
                               int ntaps, int npol, float *out )
{
    int idx = blockDim.x * blockIdx.x + threadIdx.x;

    // Calculate the number of channels
    int nchan = blockDim.x / npol;

    // Calculate the first "out" index
    int o_real = 2*idx;
    int o_imag = o_real + 1;

    // Calculate the first "in" index for this thread
    int i0 = ((idx + blockDim.x - npol) / blockDim.x) * blockDim.x +
             (threadIdx.x % npol);

    // Calculate the "fils" first column index
    //int f0 = nchan - ((idx/npol + nchan - 1) % nchan) - 1;
    int f0 = (idx/npol) % nchan;

    // Initialise the output sample to zero
    out[o_real] = 0.0;
    out[o_imag] = 0.0;

    // Multiply the filter with the in data and sum
    int tap, ch;
    int i, f;
    for (tap = 0; tap < ntaps; tap++)
    {
        for (ch = 0; ch < nchan; ch++)
        {
            // The "in" index
            i = i0 + blockDim.x*tap + npol*ch;

            // The "fils" index
            f = f0 + nchan*(ntaps-1-tap) + nchan*ntaps*ch;

            // Complex multiplication
            out[o_real] += in_real[i] * fils_real[f] -
                           in_imag[i] * fils_imag[f];
            out[o_imag] += in_real[i] * fils_imag[f] +
                           in_imag[i] * fils_real[f];
        }
    }

    // Normalise
    out[o_real] /= (float)nchan;
    out[o_imag] /= (float)nchan;

    __syncthreads();
}

extern "C"
void cu_invert_pfb_ord( ComplexDouble ***detected_beam, int file_no,
                        int nsamples, int nchan, int npol,
                        ComplexDouble **fils, int fil_size,
                        float *data_buffer_uvdif )
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
 * Finally, fils points to a 2D array of filter coefficients, each row of
 * which has been "rotated" with phase ramps of different amounts. It is
 * assumed that fils has size:
 *
 *   fils[nchan][fil_size]
 */
{
    // Flatten the input array (detected_array) for GPU.
    // We only need one second's worth, plus 12 time samples tacked onto the
    // beginning (from the previous second)

    int ntaps = 12;
    int   in_size = ((nsamples + ntaps) * nchan * npol) * sizeof(float);
    int fils_size = nchan * fil_size * sizeof(float);
    int  out_size = nsamples * nchan * npol * 2 * sizeof(float);

    float     *in_real,     *in_imag;
    float   *d_in_real,   *d_in_imag;
    float   *fils_real,   *fils_imag;
    float *d_fils_real, *d_fils_imag;

    float *d_out; // This should match the format/structure of
                  // data_buffer_uvdif exactly

    // Allocate memory on the device
    gpuErrchk(cudaMalloc( (void **)&d_in_real, in_size ));
    gpuErrchk(cudaMalloc( (void **)&d_in_imag, in_size ));
    gpuErrchk(cudaMalloc( (void **)&d_fils_real, fils_size ));
    gpuErrchk(cudaMalloc( (void **)&d_fils_imag, fils_size ));

    gpuErrchk(cudaMalloc( (void **)&d_out, out_size ));

    // Allocate memory for host copies of the same
    in_real = (float *)malloc( in_size );
    in_imag = (float *)malloc( in_size );
    fils_real = (float *)malloc( fils_size );
    fils_imag = (float *)malloc( fils_size );

    // Setup input values:
    // The starting sample index is "ntaps" places from the end of the second
    // half of detected_beam if the file number is even, and "ntaps" places
    // from the end of the first half of detected_beam if the file number is
    // odd.
    int start_s = (file_no % 2 == 0 ? 2*nsamples-ntaps : nsamples-ntaps);

    int s_in, s, ch, pol, i, f;
    for (s_in = 0; s_in < nsamples + ntaps; s_in++)
    {
        s = (start_s + s_in) % (2*nsamples);
        for (ch = 0; ch < nchan; ch++)
        {
            for (pol = 0; pol < npol; pol++)
            {
                // Calculate the index for in_real and in_imag;
                i = npol*nchan*s_in + npol*ch + pol;

                // Copy the data across - taking care of the file_no = 0 case
                if (file_no == 0 && s_in < ntaps)
                {
                    in_real[i] = 0.0;
                    in_imag[i] = 0.0;
                }
                else
                {
                    in_real[i] = CReald( detected_beam[s][ch][pol] );
                    in_imag[i] = CImagd( detected_beam[s][ch][pol] );
                }
            }
        }
    }

    // Setup filter values:
    for (ch = 0; ch < nchan; ch++)
    for (f = 0; f < fil_size; f++)
    {
        i = fil_size*ch + f;
        fils_real[i] = CReald( fils[ch][f] );
        fils_imag[i] = CImagd( fils[ch][f] );
    }

    // Copy the data to the device
    gpuErrchk(cudaMemcpy( d_in_real, in_real, in_size, cudaMemcpyHostToDevice ));
    gpuErrchk(cudaMemcpy( d_in_imag, in_imag, in_size, cudaMemcpyHostToDevice ));
    gpuErrchk(cudaMemcpy( d_fils_real, fils_real, fils_size, cudaMemcpyHostToDevice ));
    gpuErrchk(cudaMemcpy( d_fils_imag, fils_imag, fils_size, cudaMemcpyHostToDevice ));

    // Call the kernel
    filter_kernel<<<nsamples, nchan*npol>>>( d_in_real, d_in_imag,
                                             d_fils_real, d_fils_imag,
                                             ntaps, npol, d_out );
    cudaDeviceSynchronize();

    // Copy the result back into host memory
    gpuErrchk(cudaMemcpy( data_buffer_uvdif, d_out, out_size, cudaMemcpyDeviceToHost ));

    // Free memory on host and device
    free( in_real );
    free( in_imag );
    free( fils_real );
    free( fils_imag );
    cudaFree( d_in_real );
    cudaFree( d_in_imag );
    cudaFree( d_fils_real );
    cudaFree( d_fils_imag );
    cudaFree( d_out );
}

