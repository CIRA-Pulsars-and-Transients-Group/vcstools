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


void gpuAssert(cudaError_t code, const char *file, int line, bool abort)
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


void utc2mjd(char *utc_str, double *intmjd, double *fracmjd)
{
    /* Convert a UTC string (YYYY-MM-DDThh:mm:ss.ss) into MJD in radians.
     * Accepts a stc string and pointers to the integer and fractional MJD values. */
    int year, month, day, hour, min, sec, jflag;

    sscanf(utc_str,"%d-%d-%dT%d:%d:%d", &year, &month, &day, &hour, &min, &sec);
    //fprintf(stderr,"Parsed date : yr %d, month %d, day %d, hour %d, min %d, sec %f\n", year, month, day, hour, min, sec);

    slaCaldj(year, month, day, intmjd, &jflag);
    if (jflag != 0) 
    {
        fprintf(stderr,"Failed to calculate MJD\n");
    }
    *fracmjd = (hour + (min/60.0) + (sec/3600.0))/24.0;
}


void mjd2lst(double mjd, double *lst)
{
    /* Greenwich Mean Sidereal Time to LMST
     * east longitude in hours at the epoch of the MJD */
    double lmst;
    lmst = slaRanorm(slaGmst(mjd) + MWA_LON*DEG2RAD);
    *lst = lmst;
}


void calcWaveNumber(double lambda, double az, double za, wavenums *p_wn)
{
    /* Calculate the 3D wavenumbers for a given wavelength (lambda) from the direction (az,za).
     * Accepts wavelength (m), azimuth (rad) and zenith angle (rad) and a pointer to a wavenums struct to populate.*/
    double a, ast, phi;

    a   = 2 * PI / lambda;
    ast = a * sin(za);
    phi = PI/2 - az;

    /* 
     * the standard equations are:
     *      a = 2 * pi / lambda
     *      kx = a * sin(theta) * cos(phi)
     *      ky = a * sin(theta) * sin(phi)
     *      kz = a * cos(theta)
     * this is assuming that the coordinates (theta,phi) are defined in 
     * the convention from Sutinjo et al. 2015, where
     *      theta = za
     *      phi = pi/2 - az
     * i.e. the azimuth is measured clockwise from East (standard for antenna theory, offset from astronomy)
     */

    p_wn->kx = ast * cos(phi); 
    p_wn->ky = ast * sin(phi); 
    p_wn->kz = a   * cos(za);   
}

void calcTargetAZZA(char *ra_hhmmss, char *dec_ddmmss, char *time_utc, tazza *p_tazza)
{
    int ra_ih, ra_im, ra_j;
    int dec_id, dec_im, dec_j;
    int sign;
    double ra_rad, ra_fs, ha_rad;
    double dec_rad, dec_fs;
    double az, el;
    double mjd, intmjd, fracmjd, lmst;
    double pr = 0.0, pd = 0.0, px = 0.0, rv = 0.0, eq = 2000.0;
    double ra_ap = 0.0, dec_ap = 0.0;
    char id_str[20];

    // Read RA string into hours, minutes and seconds
    sscanf(ra_hhmmss, "%d:%d:%lf", &ra_ih, &ra_im, &ra_fs);

    // Read Dec string into degrees, arcmin and arsec (extra steps for sign, '+' or '-')
    sscanf(dec_ddmmss, "%s:%d:%lf", id_str, &dec_im, &dec_fs);
    sign = (id_str[0] == '-' ? -1 : 1); // Check sign of dec

    sscanf(dec_ddmmss, "%d:%d:%lf", &dec_id, &dec_im, &dec_fs); // Assign values
    dec_id = dec_id * sign; // Ensure correct sign


    // Convert angles to radians
    slaCtf2r(ra_ih, ra_im, ra_fs, &ra_rad, &ra_j); // RA 
    slaDaf2r(dec_id, dec_im, dec_fs, &dec_rad, &dec_j); // Dec

    if (ra_j != 0) 
    {
        fprintf(stderr, "CRITICAL: Error parsing %s as hhmmss\nslalib error code: j=%d\n", ra_hhmmss, ra_j);
        fprintf(stderr, "          ih = %d, im = %d, fs = %f\n", ra_ih, ra_im, ra_fs);
        exit(EXIT_FAILURE);
    }
    if (dec_j != 0) 
    {
        fprintf(stderr, "CRITICAL: Error parsing %s as ddmmss\nslalib error code: j=%d\n", dec_ddmmss, dec_j);
        fprintf(stderr, "          ih = %d, im = %d, fs = %f\n", dec_id, dec_im, dec_fs);
        exit(EXIT_FAILURE);
    }

    // Convert UTC to MJD
    utc2mjd(time_utc, &intmjd, &fracmjd);
    mjd = intmjd + fracmjd;
    mjd2lst(mjd, &lmst);

    // Get apparent RA and Dec of target
    slaMap(ra_rad, dec_rad, pr, pd, px, rv, eq, mjd, &ra_ap, &dec_ap);
    printf("RA = %.4f  RA_app = %.4f  DEC = %.4f  DEC_app = %.4f\n", ra_rad, ra_ap, dec_rad, dec_ap);

    // Use RA and LST to get HA
    ha_rad = slaRanorm(lmst - ra_ap);

    // Convert (HA, Dec) to (az, el)
    slaDe2h(ha_rad, dec_rad, MWA_LAT*DEG2RAD, &az, &el);

    printf("Az = %.4f  ZA = %.4f\n", az, PI/2-el);
    p_tazza->az = az;
    p_tazza->za = (PI/2) - el;
}

