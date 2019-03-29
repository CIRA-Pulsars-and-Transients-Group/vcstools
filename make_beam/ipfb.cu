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

__global__ void filter_kernel( float   *in_real, float   *in_imag,
                               float *fils_real, float *fils_imag,
                               int ntaps, int npol, float *out )
{
    int s = blockIdx.x;
    int nsamples = gridDim.x;
    //blockDim.x = nchan * npol 
    //threadIdx.x = npol*ch + pol
 
    // Calculate the number of channels
    int nchan = blockDim.x / npol;
    int p = threadIdx.y;
    int idx = nchan * npol * nsamples * p + nchan * npol * s + threadIdx.x;
 
    //idx = npol*nchan*(nsamples+ntaps)*p +
    //      npol*nchan*s_in +
    //      npol*ch +         pol;

    // Calculate the first "out" index
    int o_real = 2*idx;
    int o_imag = o_real + 1;

    // Calculate the first "in" index for this thread
    int i0 = ((idx + nchan * npol - npol) / (nchan * npol)) * nchan * npol 
             + (threadIdx.x % npol) + p * ntaps * nchan * npol;
    //f = fil_size*ch + tap
    
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
            i = i0 + nchan * npol*tap + npol*ch;

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
void cu_invert_pfb_ord( ComplexDouble ****detected_beam, int file_no,
                        int npointing, int nsamples, int nchan, int npol,
                        struct gpu_ipfb_arrays **g, float *data_buffer_uvdif )
/* "Invert the PFB" by applying a resynthesis filter, using GPU
 * acceleration.
 *
 * This function expects "detected_beam" to be structured as follows:
 *
 *   detected_beam[3*nsamples][nchan][npol]
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
    int start_s;
    
    if (file_no % 3 == 0)      start_s = 3*nsamples - (*g)->ntaps;
    else if (file_no % 3 == 1) start_s = nsamples - (*g)->ntaps;
    else                       start_s = 2*nsamples - (*g)->ntaps;

    int p, s_in, s, ch, pol, i;
    for (p = 0; p < npointing; p++)
    for (s_in = 0; s_in < nsamples + (*g)->ntaps; s_in++)
    {
        s = (start_s + s_in) % (3*nsamples);
        for (ch = 0; ch < nchan; ch++)
        {
            for (pol = 0; pol < npol; pol++)
            {
                // Calculate the index for in_real and in_imag;
                i = npol*nchan*(nsamples + (*g)->ntaps)*p +
                    npol*nchan*s_in + 
                    npol*ch + 
                    pol;
                // Copy the data across - taking care of the file_no = 0 case
                // The s_in%(npol*nchan*nsamples) does this for each pointing
                if (file_no == 0 && (s_in%(npol*nchan*nsamples)) < (*g)->ntaps)
                {
                    (*g)->in_real[i] = 0.0;
                    (*g)->in_imag[i] = 0.0;
                }
                else
                {
                    (*g)->in_real[i] = CReald( detected_beam[p][s][ch][pol] );
                    (*g)->in_imag[i] = CImagd( detected_beam[p][s][ch][pol] );
                }
            }
        }
    }
    
    // Copy the data to the device
    gpuErrchk(cudaMemcpy( (*g)->d_in_real, (*g)->in_real, (*g)->in_size, cudaMemcpyHostToDevice ));
    gpuErrchk(cudaMemcpy( (*g)->d_in_imag, (*g)->in_imag, (*g)->in_size, cudaMemcpyHostToDevice ));
    
    // Call the kernel
    dim3 n_cpol_p(nchan*npol, npointing);
    filter_kernel<<<nsamples, n_cpol_p>>>( (*g)->d_in_real, (*g)->d_in_imag,
                                             (*g)->d_fils_real, (*g)->d_fils_imag,
                                             (*g)->ntaps, npol, (*g)->d_out );
    cudaDeviceSynchronize();

    // Copy the result back into host memory
    gpuErrchk(cudaMemcpy( data_buffer_uvdif, (*g)->d_out, (*g)->out_size, cudaMemcpyDeviceToHost ));
}


void cu_load_filter( ComplexDouble **fils, struct gpu_ipfb_arrays **g,
        int nchan )
/* This function loads the inverse filter coefficients into GPU memory.
   It assumes that the filter size has already been set in
     (*g)->fils_size
*/
{
    int ch, f, i;
    int fil_size = (*g)->fils_size / nchan / sizeof(float);

    // Setup filter values:
    for (ch = 0; ch < nchan; ch++)
    for (f = 0; f < fil_size; f++)
    {
        i = fil_size*ch + f;
        (*g)->fils_real[i] = CReald( fils[ch][f] );
        (*g)->fils_imag[i] = CImagd( fils[ch][f] );
    }

    gpuErrchk(cudaMemcpy( (*g)->d_fils_real, (*g)->fils_real, (*g)->fils_size, cudaMemcpyHostToDevice ));
    gpuErrchk(cudaMemcpy( (*g)->d_fils_imag, (*g)->fils_imag, (*g)->fils_size, cudaMemcpyHostToDevice ));
}


void malloc_ipfb( struct gpu_ipfb_arrays **g, int ntaps, int nsamples,
        int nchan, int npol, int fil_size, int npointing )
{
    // Flatten the input array (detected_array) for GPU.
    // We only need one second's worth, plus 12 time samples tacked onto the
    // beginning (from the previous second)

    (*g)->ntaps     = ntaps;
    (*g)->in_size   = npointing * ((nsamples + ntaps) * nchan * npol) * sizeof(float);
    // fils_size = nchan * nchan * ntaps = 128 * 128 * 12
    (*g)->fils_size = nchan * fil_size * sizeof(float);
    (*g)->out_size  = npointing * nsamples * nchan * npol * 2 * sizeof(float);

    // Allocate memory on the device
    gpuErrchk(cudaMalloc( (void **)&(*g)->d_in_real,   (*g)->in_size ));
    gpuErrchk(cudaMalloc( (void **)&(*g)->d_in_imag,   (*g)->in_size ));
    gpuErrchk(cudaMalloc( (void **)&(*g)->d_fils_real, (*g)->fils_size ));
    gpuErrchk(cudaMalloc( (void **)&(*g)->d_fils_imag, (*g)->fils_size ));

    gpuErrchk(cudaMalloc( (void **)&(*g)->d_out, (*g)->out_size ));

    // Allocate memory for host copies of the same
    (*g)->in_real   = (float *)malloc( (*g)->in_size );
    (*g)->in_imag   = (float *)malloc( (*g)->in_size );
    (*g)->fils_real = (float *)malloc( (*g)->fils_size );
    (*g)->fils_imag = (float *)malloc( (*g)->fils_size );

}


void free_ipfb( struct gpu_ipfb_arrays **g )
{
    // Free memory on host and device
    free( (*g)->in_real );
    free( (*g)->in_imag );
    free( (*g)->fils_real );
    free( (*g)->fils_imag );
    cudaFree( (*g)->d_in_real );
    cudaFree( (*g)->d_in_imag );
    cudaFree( (*g)->d_fils_real );
    cudaFree( (*g)->d_fils_imag );
    cudaFree( (*g)->d_out );
}
