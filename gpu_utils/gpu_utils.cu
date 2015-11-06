#include <cuda_runtime.h>
#include <cufft.h>
#include <helper_cuda.h>
#include <helper_functions.h>

#include "gpu_utils.h"

static __device__ __host__ inline Complex ComplexAdd(Complex, Complex);
static __device__ __host__ inline Complex ComplexScale(Complex, float);
static __device__ __host__ inline Complex ComplexMul(Complex, Complex);
static __global__ void ComplexPointwiseMulAndScale(Complex *, const Complex *, int, float);
static __global__ void UpSample(Complex *, Complex *, int, int);
static __global__ void ElementWiseMultiply(Complex *, Complex *, int );
static __global__ void ShiftandReduce(Complex *signal, int offset, int shift, int nelements); 
void Filter1D( Complex  * input, Complex * output, Complex * filter, int nchan_in, int ntap, int nsamples )
{
    int mem_size = sizeof(Complex) * nsamples;

    // Allocate device memory for signal
    Complex *d_signal;
    checkCudaErrors(cudaMalloc((void **)&d_signal, mem_size));
    // Copy host memory to device
    checkCudaErrors(cudaMemcpy(d_signal, input, mem_size,cudaMemcpyHostToDevice));

    // Allocate device memory for filter kernel
    Complex *d_filter_kernel;
    checkCudaErrors(cudaMalloc((void **)&d_filter_kernel, mem_size));

    // Copy host memory to device
    checkCudaErrors(cudaMemcpy(d_filter_kernel, filter, mem_size,cudaMemcpyHostToDevice));

    // CUFFT plan
    cufftHandle plan;
    checkCudaErrors(cufftPlan1d(&plan, nsamples, CUFFT_C2C, 1));

    // Transform signal and kernel
    printf("Transforming signal cufftExecC2C\n");
    checkCudaErrors(cufftExecC2C(plan, (cufftComplex *)d_signal, (cufftComplex *)d_signal, CUFFT_FORWARD));
    checkCudaErrors(cufftExecC2C(plan, (cufftComplex *)d_filter_kernel, (cufftComplex *)d_filter_kernel, CUFFT_FORWARD));

    // Multiply the coefficients together and normalize the result
    printf("Launching ComplexPointwiseMulAndScale<<< >>>\n");
    ComplexPointwiseMulAndScale<<<32, 256>>>(d_signal, d_filter_kernel, nsamples, 1.0f / nsamples);

    // Check if kernel execution generated and error
    getLastCudaError("Kernel execution failed [ ComplexPointwiseMulAndScale ]");

    // Transform signal back
    printf("Transforming signal back cufftExecC2C\n");
    checkCudaErrors(cufftExecC2C(plan, (cufftComplex *)d_signal, (cufftComplex *)d_signal, CUFFT_INVERSE));

    // Copy device memory to host
    Complex *h_convolved_signal = output;
    checkCudaErrors(cudaMemcpy(h_convolved_signal, d_signal, mem_size,cudaMemcpyDeviceToHost));

    // Clean up

    checkCudaErrors(cudaFree(d_signal));
    checkCudaErrors(cudaFree(d_filter_kernel));

}
void cuda_invert_pfb(Complex  * input, Complex * output, Complex * h_filter, int nchan_in, int ntap, int nsamples) {
// A more specific version of the filter that tries to invert the PFB a bit more efficiently
// input:
// nsamples time samples of an nchan filter bank packed in the following order:
//
// [time][ch][complexity]
//
// Each pol is processed separately 
// all this routine has to do is to apply the filter

// filter: padded to nsamples in length
// 
// Algorithm   
// we first: Upsample the data by a factor of nchan by adding nchan-1 zeros to the streal
// Then we low pass filter each channel
// I then rotate the phase of each sample - this mimics rotating the phase of the synthesis filter
// which I will try doing if this works...


    int mem_size = sizeof(Complex) * nsamples * nchan_in * nchan_in;
    int nelements = nsamples * nchan_in;
    int verbose = 0;
    static int planned = 0;    
    
    if (verbose)
        fprintf(stdout,"Allocate memory for buffer\n");

    Complex *d_buffer;
    checkCudaErrors(cudaMalloc((void **)&d_buffer, nelements*sizeof(Complex)));
    if (verbose)
        fprintf(stdout,"Moving host memeory to device\n");

    checkCudaErrors(cudaMemcpy(d_buffer, input, nelements*sizeof(Complex),cudaMemcpyHostToDevice));
    if (verbose)
        fprintf(stdout,"Allocate device memory for signal matrix\n");
    Complex *d_signal;
    checkCudaErrors(cudaMalloc((void **)&d_signal, mem_size));
    getLastCudaError("Kernel execution failed [ cudaMalloc ]");
    // default to zero
    if (verbose)
        fprintf(stdout,"Initialising signal memory\n");

    checkCudaErrors(cudaMemset((void *) d_signal,0,mem_size));

    getLastCudaError("Kernel execution failed [ cudaMemset ]");
    // Allocate device memory for filter kernel
    if (verbose)
        fprintf(stdout,"Allocating %d complex elements to filter\n",nelements);

    Complex *d_filter_kernel;
    checkCudaErrors(cudaMalloc((void **)&d_filter_kernel, nelements*sizeof(Complex)));
  
    //Copy the filter over
    checkCudaErrors(cudaMemcpy(d_filter_kernel,h_filter,nelements*sizeof(Complex),cudaMemcpyHostToDevice));

    getLastCudaError("Kernel execution failed [ cudaMalloc ]");
    // 
    // Copy th input filter to the device (assume it has already
    // been padded and converted to a complex type.

    // Kernel to copy host signal to upsampled device matrix
    // and simultaneously perform the transpose
    
    if (verbose)
        fprintf(stdout,"Upsampling\n");

    UpSample<<<nchan_in,256>>>((Complex *) d_signal, (Complex *) d_buffer, nchan_in,nsamples);
    getLastCudaError("Kernel execution failed [ UpSample ]");
    // CUFFT plan
    // Batch FFT to do all the FFTs of the rows of the upsampled device matrix
    static cufftHandle plan;
    if (!planned) {
        checkCudaErrors(cufftPlan1d(&plan, nelements, CUFFT_C2C, nchan_in));
        planned = 1;
    }

       // Transform signal and kernel
    printf("cufftComplex size %lu\n",sizeof(cufftComplex));
    checkCudaErrors(cufftExecC2C(plan, (cufftComplex *)d_signal, (cufftComplex *)d_signal, CUFFT_FORWARD));

    // FFT padded filter
    
    checkCudaErrors(cufftExecC2C(plan, (cufftComplex *)d_filter_kernel, (cufftComplex *)d_filter_kernel, CUFFT_FORWARD));
    // multiply each row of the matrix by the filter vector
    ElementWiseMultiply<<<nchan_in,256>>>((Complex *) d_signal, (Complex *) d_filter_kernel, (int) nelements);
    getLastCudaError("Kernel execution failed [ElementWiseMultiply]");
    // Fourier Transform Back 
   
    // printf("Transforming signal back cufftExecC2C\n");
    checkCudaErrors(cufftExecC2C(plan, (cufftComplex *)d_signal, (cufftComplex *)d_signal, CUFFT_INVERSE));
 // 
//    Complex *h_buffer;
//    h_buffer = (Complex *) malloc(nelements*sizeof(Complex));
 //   checkCudaErrors(cudaMemcpy(h_buffer, d_signal, nelements*sizeof(Complex),cudaMemcpyDeviceToHost));
//    for (int i=0;i<nelements;i++){
//        fprintf(stdout,"After: %f %f\n",h_buffer[i].x, h_buffer[i].y);
//    }
//    free(h_buffer);

    
 
    // Permute each row by nsamp entries
    // and Reduce down the columns
//    for (int c = 0; c < nchan_in-1; c++){
 //       int offset = c*nelements;
 //       int shift = c*nsamples;
 //       ShiftandReduce<<<nelements,1>>>((Complex *) d_signal, (int) offset, (int) shift, (int) nelements );
 //       getLastCudaError("Kernel execution failed [ShiftandReduce]");
 //   } 
    // Copy last row to host

    int offset = (63)*nelements;

    // Copy device memory to host
    Complex *h_convolved_signal = output;
    Complex *d_outptr = &d_signal[offset];
    checkCudaErrors(cudaMemcpy(h_convolved_signal, d_outptr, nelements*sizeof(Complex),cudaMemcpyDeviceToHost));

    // Clean up

    checkCudaErrors(cudaFree(d_signal));
    checkCudaErrors(cudaFree(d_buffer));
    checkCudaErrors(cudaFree(d_filter_kernel));


}
static __device__ __host__ inline Complex ComplexAdd(Complex a, Complex b)
{
    Complex c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    return c;
}

// Complex scale
static __device__ __host__ inline Complex ComplexScale(Complex a, float s)
{
    Complex c;
    c.x = s * a.x;
    c.y = s * a.y;
    return c;
}

// Complex multiplication
static __device__ __host__ inline Complex ComplexMul(Complex a, Complex b)
{
    Complex c;
    c.x = a.x * b.x - a.y * b.y;
    c.y = a.x * b.y + a.y * b.x;
    return c;
}

// Complex pointwise multiplication
static __global__ void ComplexPointwiseMulAndScale(Complex *a, const Complex *b, int size, float scale)
{
    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;

    for (int i = threadID; i < size; i += numThreads)
    {
        a[i] = ComplexScale(ComplexMul(a[i], b[i]), scale);
    }
}

static __global__ void UpSample(Complex *signal, Complex *input, int nchan, int nsamp) {

    // For every chan (blockIdx.x) add nchan-1 zeros
    const int input_chan = blockIdx.x ;
    
    int inputID = 0;
    int outputID = 0;
    for (int sampnum = threadIdx.x ; sampnum < nsamp; sampnum=sampnum+blockDim.x) {

        inputID = input_chan*nchan + sampnum;
        outputID = blockIdx.x * (nsamp * nchan) + (sampnum * nchan);
        signal[outputID] = input[inputID];
        for (int zeros=1; zeros < nchan; zeros++) {
            outputID = input_chan*(nsamp * nchan) + (sampnum * nchan) + zeros;
            signal[outputID].x = 0.0;
            signal[outputID].y = 0.0;
        }
    }
}
static __global__ void ElementWiseMultiply(Complex *signal, Complex *filter, int nelements) {

    const int offset = blockIdx.x * nelements;

    for (int entry = threadIdx.x ; entry<nelements; entry=entry+blockDim.x){
        int inputID = entry + offset;
        int outputID = inputID;
        signal[outputID] = ComplexMul(signal[inputID],filter[entry]);
    }
}

static __global__ void ShiftandReduce(Complex *signal, int offset, int shift, int nelements) {

    const int element = blockDim.x ;
   
    int output_col = element+shift;
    if (output_col > nelements) {
       output_col = output_col - nelements;
    }

    int inputID = offset+element;
    int outputID = offset + output_col + nelements;

    signal[outputID] = ComplexAdd(signal[inputID],signal[outputID]);

} 

