#include "utils.h"

void getDeviceDimensions(int *nDevices)
{
    /* We need to know how many devices are available and its functionality. */

    printf("Querying system for device information --\n");
    cudaGetDeviceCount(nDevices); // get CUDA to count GPUs

    for (int i = 0; i < *nDevices; i++)
    {
        struct cudaDeviceProp prop; // create struct to store device info
        cudaGetDeviceProperties(&prop, i); // populate prop for this device
        printf("    Device number:                       %d\n", *nDevices-1);
        printf("    Device name:                         %s\n", prop.name);
        printf("    Total global memory available (MB):  %f\n", prop.totalGlobalMem/1e6);
        printf("    Max grid size (# blocks):           (%d, %d, %d)\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
        printf("    Max number of threads per block:    (%d, %d, %d)\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
    }
}

void requiredMemory(int size, int ntiles, int *niter, int *blockSize)
{
    /* Estimate how much device/host memory will be required per iteration.

       The device kernal consists of:
        
       az array = size * sizeof(double) - array of azimuth coords
       za array = size * sizeof(double) - array of zenith coords
       xpos = ntiles * sizeof(float) - x position of tiles from array center
       ypos = ntiles * sizeof(float) - y position of tiles from array center
       zpos = ntiles * sizeof(float) - z position of tiles from array center
       af_array = size * sizeof(cuDoubleComplex) - array factor

       plus we also pass:
       2 integers
       1 double

       and compute in the kernal:
       2 integers
       6 doubles
       1 cuDoubleComplex */

    int iters = 1, nDevices = 0;
    double tfrac = 0.9;
    long double azzaMem = 0, tileposMem = 0, afMem = 0, otherMem = 0, reqMem = 0, devMem = 0;
    size_t freeMem, totMem;
    cudaError_t res;
    struct cudaDeviceProp prop; // create struct to store device info

    // get info about avilable devices
    getDeviceDimensions(&nDevices);
    printf("Number of devices on system: %d\n", nDevices);
    printf("Using device: %d\n",                nDevices-1);
    cudaGetDeviceProperties(&prop, 0); 

    // check how much FREE memory is available
    res = cudaMemGetInfo(&freeMem, &totMem);
    if (res == cudaSuccess)
    {
        printf("Free device memory: %.2f MB\n", (double)freeMem/1.0e6);
    }
    else
    {
        printf("%s\n", cudaGetErrorString(res));
    }
    
    // get device max. threads per block
    *blockSize = prop.maxThreadsDim[0];

    // define the array sizes that will go onto the device
    azzaMem    = 2 * (size/1.0e6)   * sizeof(double);          // az and za arrays
    tileposMem = 3 * (ntiles/1.0e6) * sizeof(float);           // x,y,z positions of all tiles
    afMem      =     (size/1.0e6)   * sizeof(cuDoubleComplex); // "array factor" array
    // misc. memory requirments (likely inconsequential)

    otherMem   = (7 * sizeof(double) + 
                  4 * sizeof(int) + 
                  sizeof(cuDoubleComplex) + 
                  sizeof(wavenums)) / 1.0e6;
    
    reqMem = azzaMem + tileposMem + afMem + otherMem; // total required memory in MB
    devMem = (double)freeMem/1.0e6; // available memory in MB

    printf("Memory required for:\n");
    printf("    Az,ZA arrays: %Lf MB\n",   azzaMem);
    printf("    tile positions: %Lf MB\n", tileposMem);
    printf("    array fator: %Lf MB\n",    afMem);
    printf("    intermediate: %Lf MB\n",   otherMem);
    printf("Total: %Lf MB\n",              reqMem);

    if (reqMem < 0)
    {
        fprintf(stderr, "Negative required memory (%Lf)!! Aborting.\n", reqMem);
        exit(1);
    }
    else if ((tfrac*devMem) <= reqMem)
    {
        fprintf(stderr, "Arrays will not fit on device!\n");
        fprintf(stderr, "   total device memory required = %Lf MB\n",  reqMem);
        fprintf(stderr, "   useable device memory        = %Lf MB\n",  tfrac*devMem);
        fprintf(stderr, "       (useable fraction = %.2f of total)\n", tfrac);

        iters = (int)ceil(tfrac*reqMem / devMem)+1;
        printf("Will split task into: %d iterations (approx. %.2Lf MB per iteration)\n", iters, (reqMem/iters));
    }
    printf("\n");

    *niter = iters;
}
