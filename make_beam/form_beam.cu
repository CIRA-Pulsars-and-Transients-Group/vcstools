#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <cuda_runtime.h>

extern "C" {
#include "beam_common.h"
#include "form_beam.h"
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


// define constants to be used in the kernel
#define NSTATION  128
#define NPOL      2
#define NSTOKES   4


__global__ void beamform_kernel( uint8_t *data,
                                 ComplexDouble *W,
                                 ComplexDouble *J,
                                 double invw,
                                 ComplexDouble *Bd,
                                 float *C,
                                 float *I )
/* Layout for input arrays:
 *   data [nsamples] [nchan] [NPFB] [NREC] [NINC] -- see docs
 *   W    [NSTATION] [nchan] [NPOL]               -- weights array
 *   J    [NSTATION] [nchan] [NPOL] [NPOL]        -- jones matrix
 * Layout for output arrays:
 *   Bd   [nsamples] [nchan]   [NPOL]             -- detected beam
 *   C    [nsamples] [NSTOKES] [nchan]            -- coherent full stokes
 *   I    [nsamples] [nchan]                      -- incoherent
 */
{
    // Translate GPU block/thread numbers into meaningful names
    int sample = blockIdx.x;
    int nchan  = blockDim.x;
    int ch     = threadIdx.x;

    // Calculate the indices for the input arrays
    int Di[NSTATION][NPOL];
    int Wi[NSTATION][NPOL];
    int Ji[NSTATION][NPOL][NPOL];

    int ant, pol, pol2, st;
    int pfb, rec, inc;
    for (ant = 0; ant < NSTATION; ant++)
    {
        pfb = ant / 32;
        inc = (ant / 8) % 4;
        for (pol = 0; pol < NPOL; pol++)
        {
            rec = (2*ant+pol) % 16;

            Di[ant][pol] = sample * (NINC*NREC*NPFB*nchan) +
                           ch     * (NINC*NREC*NPFB)       +
                           pfb    * (NINC*NREC)            +
                           rec    * (NINC)                 +
                           inc;

            Wi[ant][pol] = ant * (NPOL*nchan) +
                           ch  * (NPOL)       +
                           pol;

            for (pol2 = 0; pol2 < NPOL; pol2++)
            {
                Ji[ant][pol][pol2] = ant  * (NPOL*NPOL*nchan) +
                                     ch   * (NPOL*NPOL)       +
                                     pol  * (NPOL)            +
                                     pol2;
            }
        }
    }

    // Calculate the indices for the output arrays
    int Bdi[NPOL];
    int Ci[NSTOKES];
    int Ii;

    for (pol = 0; pol < NPOL; pol++)
        Bdi[pol] = sample * (NPOL*nchan) +
                   ch     * (NPOL)       +
                   pol;

    for (st = 0; st < NSTOKES; st++)
        Ci[st] = sample * (nchan*NSTOKES) +
                 st     * (nchan)         +
                 ch;

    Ii = sample*nchan + ch;

    // Calculate the beam and the noise floor
    ComplexDouble B[NPOL];
    ComplexDouble D[NPOL];
    ComplexDouble WD[NPOL];
    ComplexDouble N[NPOL][NPOL];


    /* Fix from Maceij regarding NaNs in output when running on Athena, 13 April 2018.
       Apparently the different compilers and architectures are treating what were 
       unintialised variables very differently */
    for(pol = 0; pol < NPOL; pol++)
    {
        B[pol]  = CMaked( 0.0, 0.0 );
        D[pol]  = CMaked( 0.0, 0.0 );
        WD[pol] = CMaked( 0.0, 0.0 );
    }


    I[Ii] = 0.0;
    for (pol = 0; pol < NPOL; pol++)
    {
        // Initialise beams and noise floor to zero
        Bd[Bdi[pol]] = CMaked( 0.0, 0.0 );
        for (pol2 = 0; pol2 < NPOL; pol2++)
            N[pol][pol2] = CMaked( 0.0, 0.0 );

        for (ant = 0; ant < NSTATION; ant++)
        {
            // Calculate the coherent beam (B = J*W*D)
            B[pol]  = CMaked( 0.0, 0.0 );
            D[pol]  = UCMPLX4_TO_CMPLX_FLT(data[Di[ant][pol]]);
            WD[pol] = CMuld( W[Wi[ant][pol]], D[pol] );

            // (... and along the way, calculate the incoherent beam...)
            I[Ii] += DETECT(D[pol]);

            for (pol2 = 0; pol2 < NPOL; pol2++)
            {
                B[pol] = CAddd( B[pol], CMuld( J[Ji[ant][pol][pol2]],
                                               WD[pol2] ) );
            }

            // Detect the coherent beam
            Bd[Bdi[pol]] = CAddd( Bd[Bdi[pol]], B[pol] );

            // Calculate the noise floor (N = B*B')
            for (pol2 = 0; pol2 < NPOL; pol2++)
            {
                N[pol][pol2] = CAddd( N[pol][pol2],
                                      CMuld( B[pol], CConjd( B[pol2] ) ) );
            }
        }
    }

    // Form the stokes parameters for the coherent beam
    float bnXX = DETECT(Bd[Bdi[0]]) - CReald(N[0][0]);
    float bnYY = DETECT(Bd[Bdi[1]]) - CReald(N[1][1]);
    ComplexDouble bnXY = CSubd( CMuld( Bd[Bdi[0]], CConjd( Bd[Bdi[1]] ) ),
                                N[0][1] );

    // Stokes I, Q, U, V:
    C[Ci[0]] = invw*(bnXX + bnYY);
    C[Ci[1]] = invw*(bnXX - bnYY);
    C[Ci[2]] =  2.0*invw*CReald( bnXY );
    C[Ci[3]] = -2.0*invw*CImagd( bnXY );

    __syncthreads();
}

void cu_form_beam( uint8_t *data, struct make_beam_opts *opts,
                   ComplexDouble ***complex_weights_array,
                   ComplexDouble ****invJi, int file_no, int nstation, int nchan,
                   int npol, int outpol_coh, int outpol_incoh, double invw,
                   ComplexDouble ***detected_beam, float *coh, float *incoh )
/* The CPU version of the beamforming operations, using OpenMP for
 * parallelisation.
 *
 * Inputs:
 *   data    = array of 4bit+4bit complex numbers. For data order, refer to the
 *             documentation.
 *   opts    = passed option parameters, containing meta information about the
 *             obs and the data
 *   W       = complex weights array. [nstation][nchan][npol]
 *   J       = inverse Jones matrix. [nstation][nchan][npol][npol]
 *   file_no = number of file we are processing, starting at 0.
 *   nstation     = 128
 *   nchan        = 128
 *   npol         = 2 (X,Y)
 *   outpol_coh   = 4 (I,Q,U,V)
 *   outpol_incoh = 1 (I)
 *   invw         = the reciprocal of the sum of the antenna weights
 *
 * Outputs:
 *   detected_beam = result of beamforming operation, summed over antennas
 *                   [2*nsamples][nchan][npol]
 *   coh           = result in Stokes parameters (minus noise floor)
 *                   [nsamples][nstokes][nchan]
 *   incoh         = result (just Stokes I)
 *                   [nsamples][nchan]
 *
 * Assumes "coh" and "incoh" contain only zeros.
 */
{
    // Calculate array sizes for host and device
    size_t coh_size   = opts->sample_rate * outpol_coh   * nchan * sizeof(float);
    size_t incoh_size = opts->sample_rate * outpol_incoh * nchan * sizeof(float);
    size_t data_size  = opts->sample_rate * nstation * nchan * npol * sizeof(uint8_t);
    size_t Bd_size    = opts->sample_rate * nchan * npol * sizeof(ComplexDouble);
    size_t W_size     = nstation * nchan * npol          * sizeof(ComplexDouble);
    size_t J_size     = nstation * nchan * npol * npol   * sizeof(ComplexDouble);

    // Arrays to be passed to the GPU kernel
    // (We don't need to allocate host memory for data, coh, or incoh -- we
    // assume this is allocated before these pointers were passed into this
    // function)
    ComplexDouble *W, *d_W;
    ComplexDouble *J, *d_J;
    ComplexDouble *Bd, *d_Bd;
    uint8_t *d_data;
    float   *d_coh;
    float   *d_incoh;

    // Allocate host memory
    W  = (ComplexDouble *)malloc( W_size );
    J  = (ComplexDouble *)malloc( J_size );
    Bd = (ComplexDouble *)malloc( Bd_size );


    // Allocate device memory
    gpuErrchk(cudaMalloc( (void **)&d_W,     W_size ));
    gpuErrchk(cudaMalloc( (void **)&d_J,     J_size ));
    gpuErrchk(cudaMalloc( (void **)&d_Bd,    Bd_size ));
    gpuErrchk(cudaMalloc( (void **)&d_data,  data_size ));
    gpuErrchk(cudaMalloc( (void **)&d_coh,   coh_size ));
    gpuErrchk(cudaMalloc( (void **)&d_incoh, incoh_size ));

    // Setup input values (= populate W and J)
    int s, ant, ch, pol, pol2;
    int Wi, Ji;
    for (ant = 0; ant < nstation; ant++)
    for (ch  = 0; ch  < nchan   ; ch++ )
    for (pol = 0; pol < npol    ; pol++)
    {
        Wi = ant * (npol*nchan) +
             ch  * (npol) +
             pol;
        W[Wi] = complex_weights_array[ant][ch][pol];

        for (pol2 = 0; pol2 < npol; pol2++)
        {
            Ji = Wi*npol + pol2;
            J[Ji] = invJi[ant][ch][pol][pol2];
        }
    }

    // Copy the data to the device
    gpuErrchk(cudaMemcpy( d_data, data, data_size, cudaMemcpyHostToDevice ));
    gpuErrchk(cudaMemcpy( d_W,    W,    W_size,    cudaMemcpyHostToDevice ));
    gpuErrchk(cudaMemcpy( d_J,    J,    J_size,    cudaMemcpyHostToDevice ));

    // Call the kernel
    beamform_kernel<<<opts->sample_rate, nchan>>>(
            d_data, d_W, d_J, invw, d_Bd, d_coh, d_incoh );
    cudaDeviceSynchronize();

    // Copy the results back into host memory
    gpuErrchk(cudaMemcpy( coh,   d_coh,   coh_size,   cudaMemcpyDeviceToHost ));
    gpuErrchk(cudaMemcpy( incoh, d_incoh, incoh_size, cudaMemcpyDeviceToHost ));
    gpuErrchk(cudaMemcpy( Bd,    d_Bd,    Bd_size,    cudaMemcpyDeviceToHost ));

    // Copy the data back from Bd back into the detected_beam array
    // Make sure we put it back into the correct half of the array, depending
    // on whether this is an even or odd second.
    int offset, i;
    if (file_no % 2 == 0)
        offset = 0;
    else
        offset = opts->sample_rate;

    for (s   = 0; s   < opts->sample_rate; s++  )
    for (ch  = 0; ch  < nchan            ; ch++ )
    for (pol = 0; pol < npol             ; pol++)
    {
        i = s  * (npol*nchan) +
            ch * (npol)       +
            pol;

        detected_beam[s+offset][ch][pol] = Bd[i];
    }

    // Free memory on host and device
    free( W );
    free( J );
    free( Bd );
    cudaFree( d_W );
    cudaFree( d_J );
    cudaFree( d_Bd );
    cudaFree( d_data );
    cudaFree( d_coh );
    cudaFree( d_incoh );
}

