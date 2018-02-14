#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <slalib.h>
#include <fitsio.h>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuComplex.h>

#define PI (acos(-1.0))         // Ensures PI is defined on all systems
#define RAD2DEG (180.0 / PI)
#define DEG2RAD (PI / 180.0)
#define SOL (299792458.0)       // Speed of light
#define KB (1.38064852e-23)     // Boltzmann's constant

#define MWA_LAT (-26.703319)    // Array latitude, degrees North
#define MWA_LON (116.67081)     // Array longitude, degrees East
#define MWA_HGT (377.827)       // Array elevation above sea level, in meters

// Macro to use gpuAssert function
#define gpuErrchk(ans) {gpuAssert((ans), __FILE__, __LINE__);}


/* struct to hold all the wavenumbers for each (Az,ZA) */
typedef struct wavenums_t
{
    double kx;
    double ky;
    double kz;
} wavenums;

/* struct to hold the target Azimuth and Zenith angle (in radians) */
typedef struct tazza_t
{
    double az;
    double za;
} tazza;



void getDeviceDimensions(int *nDevices);

void requiredMemory(int size, int ntiles, int *niter, int *blockSize);

void gpuAssert(cudaError_t code, const char *file, int line, bool abort);

void utc2mjd(char *utc_str, double *intmjd, double *fracmjd);

void mjd2lst(double mjd, double *lst);

#endif

