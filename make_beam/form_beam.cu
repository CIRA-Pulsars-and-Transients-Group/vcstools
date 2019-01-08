/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

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
    int s   = blockIdx.x;  /* The (s)ample number */
    int ns  = gridDim.x;   /* The (n)umber of (s)amples (=10000)*/
    int c   = blockIdx.y;  /* The (c)hannel number */
    int p   = blockIdx.z;  /* The (p)ointing number */
    int nc  = gridDim.y;   /* The (n)umber of (c)hannels (=128) */
    int ant = threadIdx.x; /* The (ant)enna number */
    
    // Calculate the beam and the noise floor
    __shared__ double Ia[NSTATION];
    __shared__ ComplexDouble Bx[NSTATION], By[NSTATION];
    ComplexDouble Dx, Dy;
    ComplexDouble WDx, WDy;

    __shared__ ComplexDouble Nxx[NSTATION], Nxy[NSTATION],
                             Nyx[NSTATION], Nyy[NSTATION];


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
    Nyx[ant] = CMaked( 0.0, 0.0 );
    Nyy[ant] = CMaked( 0.0, 0.0 );

    Ia[ant] = 0.0;

    // Calculate beamform products for each antenna, and then add them together
    // Calculate the coherent beam (B = J*W*D)
    Dx  = UCMPLX4_TO_CMPLX_FLT(data[D_IDX(s,c,ant,0,nc)]);
    Dy  = UCMPLX4_TO_CMPLX_FLT(data[D_IDX(s,c,ant,1,nc)]);

    Ia[ant] = DETECT(Dx) + DETECT(Dy);

    WDx = CMuld( W[W_IDX(p,ant,c,0,nc)], Dx );
    WDy = CMuld( W[W_IDX(p,ant,c,1,nc)], Dy );

    Bx[ant] = CAddd( CMuld( J[J_IDX(p,ant,c,0,0,nc)], WDx ),
                     CMuld( J[J_IDX(p,ant,c,1,0,nc)], WDy ) );
    By[ant] = CAddd( CMuld( J[J_IDX(p,ant,c,0,1,nc)], WDx ),
                     CMuld( J[J_IDX(p,ant,c,1,1,nc)], WDy ) );

    Nxx[ant] = CMuld( Bx[ant], CConjd(Bx[ant]) );
    Nxy[ant] = CMuld( Bx[ant], CConjd(By[ant]) );
    Nyx[ant] = CMuld( By[ant], CConjd(Bx[ant]) );
    Nyy[ant] = CMuld( By[ant], CConjd(By[ant]) );

    // Detect the coherent beam
    __syncthreads();
    if (ant < 64)
    {
        Ia[ant] += Ia[ant+64];
        Bx[ant] = CAddd( Bx[ant], Bx[ant+64] );
        By[ant] = CAddd( By[ant], By[ant+64] );
        Nxx[ant] = CAddd( Nxx[ant], Nxx[ant+64] );
        Nxy[ant] = CAddd( Nxy[ant], Nxy[ant+64] );
        Nyx[ant] = CAddd( Nyx[ant], Nyx[ant+64] );
        Nyy[ant] = CAddd( Nyy[ant], Nyy[ant+64] );
    }
    __syncthreads();
    if (ant < 32)
    {
        Ia[ant] += Ia[ant+32];
        Bx[ant] = CAddd( Bx[ant], Bx[ant+32] );
        By[ant] = CAddd( By[ant], By[ant+32] );
        Nxx[ant] = CAddd( Nxx[ant], Nxx[ant+32] );
        Nxy[ant] = CAddd( Nxy[ant], Nxy[ant+32] );
        Nyx[ant] = CAddd( Nyx[ant], Nyx[ant+32] );
        Nyy[ant] = CAddd( Nyy[ant], Nyy[ant+32] );
    }
    if (ant < 16)
    {
        Ia[ant] += Ia[ant+16];
        Bx[ant] = CAddd( Bx[ant], Bx[ant+16] );
        By[ant] = CAddd( By[ant], By[ant+16] );
        Nxx[ant] = CAddd( Nxx[ant], Nxx[ant+16] );
        Nxy[ant] = CAddd( Nxy[ant], Nxy[ant+16] );
        Nyx[ant] = CAddd( Nyx[ant], Nyx[ant+16] );
        Nyy[ant] = CAddd( Nyy[ant], Nyy[ant+16] );
    }
    if (ant < 8)
    {
        Ia[ant] += Ia[ant+8];
        Bx[ant] = CAddd( Bx[ant], Bx[ant+8] );
        By[ant] = CAddd( By[ant], By[ant+8] );
        Nxx[ant] = CAddd( Nxx[ant], Nxx[ant+8] );
        Nxy[ant] = CAddd( Nxy[ant], Nxy[ant+8] );
        Nyx[ant] = CAddd( Nyx[ant], Nyx[ant+8] );
        Nyy[ant] = CAddd( Nyy[ant], Nyy[ant+8] );
    }
    if (ant < 4)
    {
        Ia[ant] += Ia[ant+4];
        Bx[ant] = CAddd( Bx[ant], Bx[ant+4] );
        By[ant] = CAddd( By[ant], By[ant+4] );
        Nxx[ant] = CAddd( Nxx[ant], Nxx[ant+4] );
        Nxy[ant] = CAddd( Nxy[ant], Nxy[ant+4] );
        Nyx[ant] = CAddd( Nyx[ant], Nyx[ant+4] );
        Nyy[ant] = CAddd( Nyy[ant], Nyy[ant+4] );
    }
    if (ant < 2)
    {
        Ia[ant] += Ia[ant+2];
        Bx[ant] = CAddd( Bx[ant], Bx[ant+2] );
        By[ant] = CAddd( By[ant], By[ant+2] );
        Nxx[ant] = CAddd( Nxx[ant], Nxx[ant+2] );
        Nxy[ant] = CAddd( Nxy[ant], Nxy[ant+2] );
        Nyx[ant] = CAddd( Nyx[ant], Nyx[ant+2] );
        Nyy[ant] = CAddd( Nyy[ant], Nyy[ant+2] );
    }
    if (ant < 1)
    {
        Ia[ant] += Ia[ant+1];
        Bx[ant] = CAddd( Bx[ant], Bx[ant+1] );
        By[ant] = CAddd( By[ant], By[ant+1] );
        Nxx[ant] = CAddd( Nxx[ant], Nxx[ant+1] );
        Nxy[ant] = CAddd( Nxy[ant], Nxy[ant+1] );
        Nyx[ant] = CAddd( Nyx[ant], Nyx[ant+1] );
        Nyy[ant] = CAddd( Nyy[ant], Nyy[ant+1] );
    }
    __syncthreads();

    // Form the stokes parameters for the coherent beam
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
        if ( C_IDX(p,s,0,c,ns,nc) > 40960000)
            printf("C_IDX: %d\n",C_IDX(p,s,0,c,ns,nc));
        C[C_IDX(p,s,0,c,ns,nc)] = invw*(bnXX + bnYY);
        C[C_IDX(p,s,1,c,ns,nc)] = invw*(bnXX - bnYY);
        C[C_IDX(p,s,2,c,ns,nc)] =  2.0*invw*CReald( bnXY );
        C[C_IDX(p,s,3,c,ns,nc)] = -2.0*invw*CImagd( bnXY );

        // The beamformed products
        if (B_IDX(p,s,c,0,ns,nc) > 10240000)
            printf("B_IDX: %d\n", B_IDX(p,s,c,0,ns,nc));
        Bd[B_IDX(p,s,c,0,ns,nc)] = Bx[0];
        Bd[B_IDX(p,s,c,1,ns,nc)] = By[0];
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
    int nchan  = blockDim.x; /* The total number of channels */
    int p      = blockIdx.x;
    int stokes = threadIdx.y;
    //int nstokes = blockDim.y;
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
                   struct gpu_formbeam_arrays **g,
                   ComplexDouble ****detected_beam, float *coh, float *incoh )
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
    int p, s, ant, ch, pol, pol2;
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
        (*g)->W[Wi] = complex_weights_array[p][ant][ch][pol];

        for (pol2 = 0; pol2 < npol; pol2++)
        {
            Ji = Wi*npol + pol2;
            (*g)->J[Ji] = invJi[p][ant][ch][pol][pol2];
        }
    }

    // events for timing
    cudaEvent_t startEvent, stopEvent; 
    cudaEventCreate(&startEvent);
    cudaEventCreate(&stopEvent);
    float time;

    // Copy the data to the device
    cudaEventRecord(startEvent, 0);
    gpuErrchk(cudaMemcpy( (*g)->d_data, data,    (*g)->data_size, cudaMemcpyHostToDevice ));
    cudaEventRecord(stopEvent, 0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&time, startEvent, stopEvent);
    printf(" d_data memmcpy time : %7.3f ms  transfer: %7.3f GB/s \n",time, (*g)->data_size * 1e-6 / time);

    cudaEventRecord(startEvent, 0);
    gpuErrchk(cudaMemcpy( (*g)->d_W,    (*g)->W, (*g)->W_size,    cudaMemcpyHostToDevice ));
    cudaEventRecord(stopEvent, 0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&time, startEvent, stopEvent);
    printf(" d_w    memmcpy time : %7.3f ms  transfer: %7.3f GB/s \n",time, (*g)->W_size * 1e-6 / time);

    cudaEventRecord(startEvent, 0);
    gpuErrchk(cudaMemcpy( (*g)->d_J,    (*g)->J, (*g)->J_size,    cudaMemcpyHostToDevice ));
    cudaEventRecord(stopEvent, 0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&time, startEvent, stopEvent);
    printf(" d_j    memmcpy time : %7.3f ms  transfer: %7.3f GB/s \n",time, (*g)->J_size * 1e-6 / time);

    // Call the kernels
    dim3 samples_chan(opts->sample_rate, nchan, npointing);
    #define NOW  ((double)clock()/(double)CLOCKS_PER_SEC)
    double begintime = NOW;
    beamform_kernel<<<samples_chan, NSTATION>>>(
            (*g)->d_data, (*g)->d_W, (*g)->d_J, invw, (*g)->d_Bd, (*g)->d_coh, (*g)->d_incoh );
    //cudaDeviceSynchronize();
    // sync not required between kernel queues since each stream acts like a FIFO queue
    // so all instances of the above kernel will complete before we move to the next
    // we are using the "default" stream since we don't specify any stream id

    // 1 block per pointing direction, hence the 1 for now
    // TODO check if these actually work, can't see them return values.
    // The incoh kernal also takes 40 second for some reason so commenting out
    //flatten_bandpass_I_kernel<<<1, nchan>>>((*g)->d_incoh, opts->sample_rate);
    //cudaDeviceSynchronize();

    // now do the same for the coherent beam
    dim3 chan_stokes(nchan, outpol_coh);
    //flatten_bandpass_C_kernel<<<npointing, chan_stokes>>>((*g)->d_coh, opts->sample_rate);
    //cudaDeviceSynchronize(); // Memcpy acts as a synchronize step so don't sync here
    // Copy the results back into host memory
    cudaEventRecord(startEvent, 0);
    gpuErrchk(cudaMemcpy( (*g)->Bd, (*g)->d_Bd,    (*g)->Bd_size,    cudaMemcpyDeviceToHost ));
    cudaEventRecord(stopEvent, 0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&time, startEvent, stopEvent);
    printf(" d_Bd   memmcpy time : %7.3f ms  transfer: %7.3f GB/s \n",time, (*g)->Bd_size * 1e-6 / time);

    cudaEventRecord(startEvent, 0);
    gpuErrchk(cudaMemcpy( incoh,    (*g)->d_incoh, (*g)->incoh_size, cudaMemcpyDeviceToHost ));
    cudaEventRecord(stopEvent, 0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&time, startEvent, stopEvent);
    printf(" incoh  memmcpy time : %7.3f ms  transfer: %7.3f GB/s \n",time, (*g)->incoh_size * 1e-6 / time);

    cudaEventRecord(startEvent, 0);
    gpuErrchk(cudaMemcpy( coh,      (*g)->d_coh,   (*g)->coh_size,   cudaMemcpyDeviceToHost ));
    cudaEventRecord(stopEvent, 0);
    cudaEventSynchronize(stopEvent);
    cudaEventElapsedTime(&time, startEvent, stopEvent);
    printf(" coh    memmcpy time : %7.3f ms  transfer: %7.3f GB/s \n",time, (*g)->coh_size * 1e-6 / time);
    // Copy the data back from Bd back into the detected_beam array
    // Make sure we put it back into the correct half of the array, depending
    // on whether this is an even or odd second.
    int offset, i;
    offset = file_no % 3 * opts->sample_rate;
    
    for (p   = 0; p   < npointing        ; p++  )
    for (s   = 0; s   < opts->sample_rate; s++  )
    for (ch  = 0; ch  < nchan            ; ch++ )
    for (pol = 0; pol < npol             ; pol++)
    {
        i = p  * (npol*nchan*opts->sample_rate) +
            s  * (npol*nchan)                   +
            ch * (npol)                         +
            pol;

        detected_beam[p][s+offset][ch][pol] = (*g)->Bd[i];
    }
}

void malloc_formbeam( struct gpu_formbeam_arrays **g, unsigned int sample_rate,
        int nstation, int nchan, int npol, int outpol_coh, int outpol_incoh, int npointing)
{
    // Calculate array sizes for host and device
    (*g)->coh_size   = npointing * sample_rate * outpol_coh   * nchan * sizeof(float);
    (*g)->incoh_size = sample_rate * outpol_incoh * nchan * sizeof(float);
    (*g)->data_size  = npointing * sample_rate * nstation * nchan * npol * sizeof(uint8_t);
    (*g)->Bd_size    = npointing * sample_rate * nchan * npol * sizeof(ComplexDouble);
    (*g)->W_size     = npointing * nstation * nchan * npol * sizeof(ComplexDouble);
    (*g)->J_size     = npointing * nstation * nchan * npol * npol * sizeof(ComplexDouble);

    // Allocate host memory
    (*g)->W  = (ComplexDouble *)malloc( (*g)->W_size );
    (*g)->J  = (ComplexDouble *)malloc( (*g)->J_size );
    (*g)->Bd = (ComplexDouble *)malloc( (*g)->Bd_size );


    // Allocate device memory
    gpuErrchk(cudaMalloc( (void **)&(*g)->d_W,     (*g)->W_size ));
    gpuErrchk(cudaMalloc( (void **)&(*g)->d_J,     (*g)->J_size ));
    gpuErrchk(cudaMalloc( (void **)&(*g)->d_Bd,    (*g)->Bd_size ));
    gpuErrchk(cudaMalloc( (void **)&(*g)->d_data,  (*g)->data_size ));
    gpuErrchk(cudaMalloc( (void **)&(*g)->d_coh,   (*g)->coh_size ));
    gpuErrchk(cudaMalloc( (void **)&(*g)->d_incoh, (*g)->incoh_size ));

    printf("%d B GPU memory allocated\n", (*g)->W_size + (*g)->J_size + (*g)->Bd_size + 
                                      (*g)->data_size + (*g)->coh_size + (*g)->incoh_size );
}

void free_formbeam( struct gpu_formbeam_arrays **g )
{
    // Free memory on host and device
    free( (*g)->W );
    free( (*g)->J );
    free( (*g)->Bd );
    cudaFree( (*g)->d_W );
    cudaFree( (*g)->d_J );
    cudaFree( (*g)->d_Bd );
    cudaFree( (*g)->d_data );
    cudaFree( (*g)->d_coh );
    cudaFree( (*g)->d_incoh );
}
