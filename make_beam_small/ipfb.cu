#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include <cuComplex.h>

__global__ void filter_kernel( float   *in_real, float   *in_imag,
                               float *fils_real, float *fils_imag,
                               int ntaps, float *out )
{
    // Called the kernel with
    //   (10000, 128, 2) blocks  = (nsamples, nchan, npol)
    //   (    1,   1, 1) threads

    // Calculate the first "in" index for this block
    int i0 = blockDim.z * blockDim.y * blockIdx.x +
             blockDim.z * blockIdx.y +
             blockIdx.z;

    // Calculate the "out" index
    int o_real = 2 * i;
    int o_imag = o_real + 1;

    // Calculate the "fils" first column index
    int f0 = blockDim.y - ((blockIdx.x-1) % blockDim.y) - 1;

    // Initialise the output sample to zero
    out[o_real] = 0.0;
    out[o_imag] = 0.0;

    // Multiply the filter with the in data and sum
    int tap, ch;
    int i, f;
    float tmp;
    for (tap = 0; tap < ntaps; tap++)
    {
        for (ch = 0; ch < blockDim.y; ch++)
        {
            // The "in" index
            i = blockDim.z * blockDim.y * (blockIdx.x + tap) +
                blockDim.z * ch +
                blockIdx.z;

            // The "fils" index
            f = blockDim.y * ntaps * ch +
                blockDim.y * tap +
                f0;

            // Complex multiplication
            out[o_real] += in_real[i] * fils_real[i] -
                           in_imag[i] * fils_imag[i];
            out[o_imag] += in_real[i] * fils_imag[i] +
                           in_imag[i] * fils_real[i];
        }
    }

    // Normalise
    out[o_real] /= (float)blockDim.y;
    out[o_imag] /= (float)blockDim.y;

    __syncthreads();
}

void cu_invert_pfb_ord( ComplexFloat ***detected_beam, int file_no,
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
    cudaMalloc( (void **)&d_in_real, in_size );
    cudaMalloc( (void **)&d_in_imag, in_size );
    cudaMalloc( (void **)&d_fils_real, fils_size );
    cudaMalloc( (void **)&d_fils_imag, fils_size );

    cudaMalloc( (void **)&d_out, out_size );

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

    int s_in;
#pragma omp parallel for
    for (s_in = 0; s_in < nsamples + ntaps; s_in++)
    {
        int s, ch, pol, i;
        s = (start_s + s_in) % (2*nsamples);
        for (ch = 0; ch < nchan; ch++)
        {
            for (pol = 0; pol < npol; pol++)
            {
                // Calculate the index for in_real and in_imag;
                i = npol*nchan*s_in + npol*ch + pol;

                // Copy the data across
                in_real[i] = creal( detected_beam[s][ch][pol] );
                in_imag[i] = cimag( detected_beam[s][ch][pol] );
            }
        }
    }

    // Setup filter values:
    int ch, f, i;
#pragma omp parallel for
    for (ch = 0; ch < nchan; ch++)
    for (f = 0; f < fil_size; f++)
    {
        i = fil_size*ch + f;
        fils_real[i] = creal( fils[ch][f] );
        fils_imag[i] = cimag( fils[ch][f] );
    }

    // Copy the data to the device
    cudaMemcpy( d_in_real, in_real, in_size, cudaMemcpyHostToDevice );
    cudaMemcpy( d_in_imag, in_imag, in_size, cudaMemcpyHostToDevice );
    cudaMemcpy( d_fils_real, fils_real, fils_size, cudaMemcpyHostToDevice );
    cudaMemcpy( d_fils_imag, fils_imag, fils_size, cudaMemcpyHostToDevice );

    // Call the kernel with
    //   (10000, 128, 2) blocks  = (nsamples, nchan, npol)
    //   (    1,   1, 1) threads
    dim3 blocks( nsamples, nchan, npol );
    dim3 threads( 1, 1, 1 );

    filter_kernel<<<blocks, threads>>>( d_in_real, d_in_imag,
                                        d_fils_real, d_fils_imag,
                                        ntaps, d_out );

    // Copy the result back into host memory
    cudaMemcpy( data_buffer_uvdif, d_out, out_size, cudaMemcpyDeviceToHost );

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

