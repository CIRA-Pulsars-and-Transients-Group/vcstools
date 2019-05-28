/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <cuda_runtime.h>

extern "C" {
#include "beam_common.h"
#include "form_beam.h"
#include "mycomplex.h"
}


#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)

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
// maximum number of pointings (currently)
#define NPOINTING 4



__global__ void beamform_kernel( uint8_t *data,
                                 ComplexDouble *W,
                                 ComplexDouble *J,
                                 double invw,
                                 ComplexDouble *Bd,
                                 float *C,
                                 float *I,
                                 int p)
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
    // Translate GPU block/thread numbers into meaning->l names
    int s    = blockIdx.x;  /* The (s)ample number */
    int ns   = gridDim.x;   /* The (n)umber of (s)amples (=10000)*/
    int c    = blockIdx.y;  /* The (c)hannel number */
    int nc   = gridDim.y;   /* The (n)umber of (c)hannels (=128) */
    
    int ant  = threadIdx.x; /* The (ant)enna number */
    int nant = blockDim.x;  /* The (n_umber of (ant)ennas */

    // GPU profiling
    clock_t start, stop;
    double setup_t, detect_t, sum_t, stokes_t;
    if ((p == 0) && (ant == 0) && (c == 0) && (s == 0)) start = clock();
    
    // Calculate the beam and the noise floor
    __shared__ double Ia[NSTATION];
    __shared__ ComplexDouble Bx[NSTATION], By[NSTATION];
    ComplexDouble Dx, Dy;
    ComplexDouble WDx, WDy;

    __shared__ ComplexDouble Nxx[NSTATION], Nxy[NSTATION],
                             Nyy[NSTATION];//Nyx[NSTATION]


    /* Fix from Maceij regarding NaNs in output when running on Athena, 13 April 2018.
       Apparently the different compilers and architectures are treating what were 
       unintialised variables very differently */
    Bx[ant]  = CMaked( 0.0, 0.0 );
    By[ant]  = CMaked( 0.0, 0.0 );

    Dx  = CMaked( 0.0, 0.0 );
    Dy  = CMaked( 0.0, 0.0 );

    WDx = CMaked( 0.0, 0.0 );
    WDy = CMaked( 0.0, 0.0 );

    Nxx[ant] = CMaked( 0.0, 0.0 );
    Nxy[ant] = CMaked( 0.0, 0.0 );
    //Nyx[ant] = CMaked( 0.0, 0.0 );
    Nyy[ant] = CMaked( 0.0, 0.0 );

    if (p == 0) Ia[ant] = 0.0;

    // Calculate beamform products for each antenna, and then add them together
    // Calculate the coherent beam (B = J*W*D)
    if ((p == 0) && (ant == 0) && (c == 0) && (s == 0))
    {
        stop = clock();
        setup_t = (double)(stop - start) / CLOCKS_PER_SEC * NPOINTING * NANT;
        start = clock();
    }

    Dx  = UCMPLX4_TO_CMPLX_FLT(data[D_IDX(s,c,ant,0,nc)]);
    Dy  = UCMPLX4_TO_CMPLX_FLT(data[D_IDX(s,c,ant,1,nc)]);

    if (p == 0) Ia[ant] = DETECT(Dx) + DETECT(Dy);

    WDx = CMuld( W[W_IDX(p,ant,c,0,nc)], Dx );
    WDy = CMuld( W[W_IDX(p,ant,c,1,nc)], Dy );

    Bx[ant] = CAddd( CMuld( J[J_IDX(p,ant,c,0,0,nc)], WDx ),
                               CMuld( J[J_IDX(p,ant,c,1,0,nc)], WDy ) );
    By[ant] = CAddd( CMuld( J[J_IDX(p,ant,c,0,1,nc)], WDx ),
                               CMuld( J[J_IDX(p,ant,c,1,1,nc)], WDy ) );

    Nxx[ant] = CMuld( Bx[ant], CConjd(Bx[ant]) );
    Nxy[ant] = CMuld( Bx[ant], CConjd(By[ant]) );
    //Nyx[ant] = CMuld( By[ant], CConjd(Bx[ant]) );
    Nyy[ant] = CMuld( By[ant], CConjd(By[ant]) );

    // Detect the coherent beam
    // A summation over an array is faster on a GPU if you add half on array 
    // to its other half as than can be done in parallel. Then this is repeated
    // with half of the previous array until the array is down to 1.
    __syncthreads();

    if ((p == 0) && (ant == 0) && (c == 0) && (s == 0))
    {
        stop = clock();
        detect_t = (double)(stop - start) / CLOCKS_PER_SEC * NPOINTING * NANT;
        start = clock();
    }
    for ( int h_ant = nant / 2; h_ant > 0; h_ant = h_ant / 2 )
    {
        if (ant < h_ant)
        {
            if (p == 0) Ia[ant] += Ia[ant+h_ant];
            Bx[ant]  = CAddd( Bx[ant],  Bx[ant  + h_ant] );
            By[ant]  = CAddd( By[ant],  By[ant  + h_ant] );
            Nxx[ant] = CAddd( Nxx[ant], Nxx[ant + h_ant] );
            Nxy[ant] = CAddd( Nxy[ant], Nxy[ant + h_ant] );
            //Nyx[ant]=CAddd( Nyx[ant], Nyx[ant + h_ant] );
            Nyy[ant] = CAddd( Nyy[ant], Nyy[ant + h_ant] );
        }
        // below makes no difference so removed
        // else return;
        __syncthreads();
    }

    if ((p == 0) && (ant == 0) && (c == 0) && (s == 0))
    {
        stop = clock();
        sum_t = (double)(stop - start) / CLOCKS_PER_SEC * NPOINTING * NANT;
        start = clock();

    }
    
    // Form the stokes parameters for the coherent beam
    // Only doing it for ant 0 so that it only prints once
    if (ant == 0)
    {
        float bnXX = DETECT(Bx[0]) - CReald(Nxx[0]);
        float bnYY = DETECT(By[0]) - CReald(Nyy[0]);
        ComplexDouble bnXY = CSubd(
                                 CMuld( Bx[0], CConjd( By[0] ) ),
                                 Nxy[0] );

        // The incoherent beam
        I[I_IDX(s,c,nc)] = Ia[0];

        // Stokes I, Q, U, V:
        C[C_IDX(p,s,0,c,ns,nc)] = invw*(bnXX + bnYY);
        C[C_IDX(p,s,1,c,ns,nc)] = invw*(bnXX - bnYY);
        C[C_IDX(p,s,2,c,ns,nc)] =  2.0*invw*CReald( bnXY );
        C[C_IDX(p,s,3,c,ns,nc)] = -2.0*invw*CImagd( bnXY );

        // The beamformed products
        Bd[B_IDX(p,s,c,0,ns,nc)] = Bx[0];
        Bd[B_IDX(p,s,c,1,ns,nc)] = By[0];
    }
    if ((p == 0) && (ant == 0) && (c == 0) && (s == 0))
    {
        stop = clock();
        stokes_t = (double)(stop - start) / CLOCKS_PER_SEC * NPOINTING * NANT;
        printf("Time:  setup: % f detect: %f    sum: %f     stokes: %f\n", setup_t, detect_t, sum_t, stokes_t);
    }
    
}

__global__ void flatten_bandpass_I_kernel(float *I,
                                     int nstep)/* uint8_t *Iout ) */
{
    // For just doing stokes I
    // One block
    // 128 threads each thread will do one channel
    // (we have already summed over all ant)

    // For doing the C array (I,Q,U,V)
    // ... figure it out later.

    // Translate GPU block/thread numbers into meaningful names
    int chan = threadIdx.x; /* The (c)hannel number */
    int nchan = blockDim.x; /* The total number of channels */
    float band;

    int new_var = 32; /* magic number */
    int i;

    float *data_ptr;

    // initialise the band 'array'
    band = 0.0;

    // accumulate abs(data) over all time samples and save into band
    data_ptr = I + I_IDX(0, chan, nchan);
    for (i=0;i<nstep;i++) { // time steps
        band += fabsf(*data_ptr);
        data_ptr = I + I_IDX(i,chan,nchan);
    }

    // now normalise the incoherent beam
    data_ptr = I + I_IDX(0, chan, nchan);
    for (i=0;i<nstep;i++) { // time steps
        *data_ptr = (*data_ptr)/( (band/nstep)/new_var );
        data_ptr = I + I_IDX(i,chan,nchan);
    }

}


__global__ void flatten_bandpass_C_kernel(float *C,
                                          int nstep)/* uint8_t *Iout ) */
{
    // For just doing stokes I
    // One block
    // 128 threads each thread will do one channel
    // (we have already summed over all ant)

    // For doing the C array (I,Q,U,V)
    // ... figure it out later.

    // Translate GPU block/thread numbers into meaningful names
    int chan   = threadIdx.x; /* The (c)hannel number */
    int nchan  = blockDim.x;  /* The (n)umber of (c)hannels */
    int stokes = threadIdx.y; /* The (stokes) number */
    
    int p      = blockIdx.x;  /* The (p)ointing number */

    float band;

    int new_var = 32; /* magic number */
    int i;

    float *data_ptr;

    // initialise the band 'array'
    band = 0.0;

    // accumulate abs(data) over all time samples and save into band
    //data_ptr = C + C_IDX(0,stokes,chan,nchan);
    for (i=0;i<nstep;i++) { // time steps
        data_ptr = C + C_IDX(p,i,stokes,chan,nstep,nchan);
        band += fabsf(*data_ptr);
    }

    // now normalise the coherent beam
    //data_ptr = C + C_IDX(0,stokes,chan,nchan);
    for (i=0;i<nstep;i++) { // time steps
        data_ptr = C + C_IDX(p,i,stokes,chan,nstep,nchan);
        *data_ptr = (*data_ptr)/( (band/nstep)/new_var );
    }

}


void cu_form_beam( uint8_t *data, struct make_beam_opts *opts,
                   ComplexDouble ****complex_weights_array,
                   ComplexDouble *****invJi, int file_no, 
                   int npointing, int nstation, int nchan,
                   int npol, int outpol_coh, double invw,
                   struct gpu_formbeam_arrays *g,
                   ComplexDouble ****detected_beam, float *coh, float *incoh,
                   cudaStream_t *streams )
/* The CPU version of the beamforming operations, using OpenMP for
 * parallelisation.
 *
 * Inputs:
 *   data    = array of 4bit+4bit complex numbers. For data order, refer to the
 *             documentation.
 *   opts    = passed option parameters, containing meta information about the
 *             obs and the data
 *   W       = complex weights array. [npointing][nstation][nchan][npol]
 *   J       = inverse Jones matrix. [npointing][nstation][nchan][npol][npol]
 *   file_no = number of file we are processing, starting at 0.
 *   nstation     = 128
 *   nchan        = 128
 *   npol         = 2 (X,Y)
 *   outpol_coh   = 4 (I,Q,U,V)
 *   invw         = the reciprocal of the sum of the antenna weights
 *   g            = struct containing pointers to various arrays on
 *                  both host and device
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
    // Setup input values (= populate W and J)
    int p, ant, ch, pol, pol2;
    int Wi, Ji;
    for (p   = 0; p   < npointing; p++  )
    for (ant = 0; ant < nstation ; ant++)
    for (ch  = 0; ch  < nchan    ; ch++ )
    for (pol = 0; pol < npol     ; pol++)
    {
        Wi = p   * (npol*nchan*nstation) +
             ant * (npol*nchan) +
             ch  * (npol) +
             pol;
        g->W[Wi] = complex_weights_array[p][ant][ch][pol];

        for (pol2 = 0; pol2 < npol; pol2++)
        {
            Ji = Wi*npol + pol2;
            g->J[Ji] = invJi[p][ant][ch][pol][pol2];
        }
    }
    // Copy the data to the device
    gpuErrchk(cudaMemcpyAsync( g->d_data, data, g->data_size, cudaMemcpyHostToDevice ));
    gpuErrchk(cudaMemcpyAsync( g->d_W,    g->W, g->W_size,    cudaMemcpyHostToDevice ));
    gpuErrchk(cudaMemcpyAsync( g->d_J,    g->J, g->J_size,    cudaMemcpyHostToDevice ));

    // Call the kernels
    // sammples_chan(index=blockIdx.x  size=gridDim.x,
    //               index=blockIdx.y  size=gridDim.y)
    // stat_point   (index=threadIdx.x size=blockDim.x,
    //               index=threadIdx.y size=blockDim.y)
    dim3 samples_chan(opts->sample_rate, nchan);
    dim3 stat(NSTATION);
    // Send off a parrellel cuda stream for each pointing
    for ( p = 0; p < npointing; p++ )
    {    
        beamform_kernel<<<samples_chan, stat, 0, streams[p]>>>( g->d_data,
                            g->d_W, g->d_J, invw, 
                            g->d_Bd, g->d_coh, g->d_incoh, p );
            
        gpuErrchk( cudaPeekAtLastError() );
        //cudaDeviceSynchronize();
        // sync not required between kernel queues since each stream acts like a FIFO queue
        // so all instances of the above kernel will complete before we move to the next
        // we are using the "default" stream since we don't specify any stream id

        // 1 block per pointing direction, hence the 1 for now
        // TODO check if these actually work, can't see them return values.
        // The incoh kernal also takes 40 second for some reason so commenting out
        if ( p == 0 )
            flatten_bandpass_I_kernel<<<1, nchan, 0, streams[p]>>>(g->d_incoh,
                                                     opts->sample_rate);
        //cudaDeviceSynchronize();

        // now do the same for the coherent beam
        dim3 chan_stokes(nchan, outpol_coh);
        // This doesn't seem to change anything some commenting out
        flatten_bandpass_C_kernel<<<npointing, chan_stokes, 0, streams[p]>>>(g->d_coh, 
                                                               opts->sample_rate);
        //cudaDeviceSynchronize(); // Memcpy acts as a synchronize step so don't sync here
    }
    gpuErrchk( cudaDeviceSynchronize() );
    
    // Copy the results back into host memory
    gpuErrchk(cudaMemcpyAsync( g->Bd, g->d_Bd,    g->Bd_size,    cudaMemcpyDeviceToHost ));
    gpuErrchk(cudaMemcpyAsync( incoh, g->d_incoh, g->incoh_size, cudaMemcpyDeviceToHost ));
    gpuErrchk(cudaMemcpyAsync( coh,   g->d_coh,   g->coh_size,   cudaMemcpyDeviceToHost ));
    
    // Copy the data back from Bd back into the detected_beam array
    // Make sure we put it back into the correct half of the array, depending
    // on whether this is an even or odd second.
    int offset, i;
    offset = file_no % 3 * opts->sample_rate;
    
    for ( int p   = 0; p   < npointing        ; p++  )
    for ( int s   = 0; s   < opts->sample_rate; s++  )
    for ( int ch  = 0; ch  < nchan            ; ch++ )
    for ( int pol = 0; pol < npol             ; pol++)
    {
        i = p  * (npol*nchan*opts->sample_rate) +
            s  * (npol*nchan)                   +
            ch * (npol)                         +
            pol;

        detected_beam[p][s+offset][ch][pol] = g->Bd[i];
    }
}

void malloc_formbeam( struct gpu_formbeam_arrays *g, unsigned int sample_rate,
        int nstation, int nchan, int npol, int outpol_coh, int outpol_incoh, int npointing)
{
    // Calculate array sizes for host and device
    g->coh_size   = npointing * sample_rate * outpol_coh   * nchan * sizeof(float);
    g->incoh_size = sample_rate * outpol_incoh * nchan * sizeof(float);
    g->data_size  = sample_rate * nstation * nchan * npol * sizeof(uint8_t);
    g->Bd_size    = npointing * sample_rate * nchan * npol * sizeof(ComplexDouble);
    g->W_size     = npointing * nstation * nchan * npol * sizeof(ComplexDouble);
    g->J_size     = npointing * nstation * nchan * npol * npol * sizeof(ComplexDouble);

    // Allocate host memory
    //g->W  = (ComplexDouble *)malloc( g->W_size );
    //g->J  = (ComplexDouble *)malloc( g->J_size );
    //g->Bd = (ComplexDouble *)malloc( g->Bd_size );
    cudaMallocHost( &g->W, g->W_size );
    cudaCheckErrors("cudaMallocHost W fail");
    cudaMallocHost( &g->J, g->J_size );
    cudaCheckErrors("cudaMallocHost J fail");
    cudaMallocHost( &g->Bd, g->Bd_size );
    cudaCheckErrors("cudaMallocHost Bd fail");

    // Allocate device memory
    gpuErrchk(cudaMalloc( (void **)&g->d_W,     g->W_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_J,     g->J_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_Bd,    g->Bd_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_data,  g->data_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_coh,   g->coh_size ));
    gpuErrchk(cudaMalloc( (void **)&g->d_incoh, g->incoh_size ));

    printf("%d GB GPU memory allocated\n", (g->W_size + g->J_size + 
                                            g->Bd_size + g->data_size +
                                            g->coh_size + g->incoh_size)
                                            /1000000000 );
}

void free_formbeam( struct gpu_formbeam_arrays *g )
{
    // Free memory on host and device
    cudaFreeHost( g->W );
    cudaFreeHost( g->J );
    cudaFreeHost( g->Bd );
    cudaFree( g->d_W );
    cudaFree( g->d_J );
    cudaFree( g->d_Bd );
    cudaFree( g->d_data );
    cudaFree( g->d_coh );
    cudaFree( g->d_incoh );
}

float *create_pinned_data_buffer_psrfits( size_t size )
{
    float *ptr;
    cudaMallocHost( &ptr, size * sizeof(float) );
    //cudaError_t status = cudaHostRegister((void**)&ptr, size * sizeof(float),
    //                                      cudaHostRegisterPortable );
    cudaCheckErrors("cudaMallocHost data_buffer_psrfits fail");
    return ptr;
}

float *create_pinned_data_buffer_vdif( size_t size )
{
    float *ptr;
    cudaMallocHost( &ptr, size * sizeof(float) );
    //cudaError_t status = cudaHostRegister((void**)&ptr, size * sizeof(float),
    //                                      cudaHostRegisterPortable );
    cudaCheckErrors("cudaMallocHost data_buffer_vdif fail");
    return ptr;
}

