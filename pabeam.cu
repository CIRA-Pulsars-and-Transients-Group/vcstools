#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <slalib.h>
#include <fitsio.h>
#include <complex.h>
// CUDA specific includes
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuComplex.h>

#define PI (acos(-1.0))         //Ensures PI is defined on all systems
#define RAD2DEG (180.0 / PI)
#define DEG2RAD (PI / 180.0)
#define SOL (299792458.0)       //Speed of light
#define KB (1.38064852e-23)     //Boltzmann's constant

#define MWA_LAT (-26.703319)    //Array latitude, degrees North
#define MWA_LON (116.67081)     //Array longitude, degrees East
#define MWA_HGT (377.827)       //Array elevation above sea level, in meters

// in column major format, a matrix is indexed to an array by
#define IDX2C(i,j,ld) (((j) * (ld)) + (i))
// here ld is leading dimension of the matrix 
// (in this case, should be number of rows)

/* struct to hold all the wavenumbers for each (Az,ZA) */
typedef struct wavenums_t
{
    double kx;
    double ky;
    double kz;
} wavenums; // can just refer to this struct as type wavenums

/* struct to hold the target Azimuth and Zenith angle (in radians) */
typedef struct tazza_t
{
    double az;
    double za;
} tazza;


/* Define all the function prototypes */
void usage();
void utc2mjd(char *utc_str, double *intmjd, double *fracmjd);
void mjd2lst(double mjd, double *lst);
void calcWaveNumber(double lambda, double az, double za, wavenums *p_wn);
void calcTargetAZZA(char *ra_hhmmss, char *dec_ddmmss, char *time_utc, tazza *p_tazza);
int getNumTiles(const char *metafits);
void getTilePositions(const char *metafits, int ninput,\
                        float *n_pols, float *e_pols, float *h_pols,\
                        float *n_tile, float *e_tile, float *h_tile);
int getFlaggedTiles(const char *badfile, int *badtiles);
void removeFlaggedTiles(float *n_tile, float *e_tile, float *h_tile,\
                            int *badtiles, int nbad, int nelements);
void requiredMemory(int size, int ntiles, cudaDeviceProp *dprop, int *niter);
void getDeviceDimensions(int *nDevices);



void usage()
{
    printf("pabeam -- computes the array factor that represents the MWA tied-array beam for a given observation and pointing\n");
    printf("syntax:\n\
            pabeam -f <frequency (Hz)> -r <ra in hh:mm:ss> -d <dec in dd:mm:ss> -t <utc time string> -m <metafits file> -b <RTS flagged_tiles.txt file>\n\n");
    printf("Options:\n\
            -f observing frequncy, in Hz\n\
            -e radiation efficiency (default: 1.0)\n\
            -r target RA (J2000), in hh:mm:ss.ss format\n\
            -d target DEC (J2000), in dd:mm:ss.ss format\n\
            -t UTC time to evaluate, in format YYYY-MM-DDThh:mm:ss\n\
            -m metafits files for the observation\n\
            -b RTS flagged_tiles.txt file from calibration\n\
            -x Azimuth grid resolution element (default: 1.0)\n\
            -y Zenith angle grid resolution element (default: 1.0)\n");
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

    a = 2 * PI / lambda;
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
    p_wn->kz = ast;   
}

double complex calcArrayFactor(int ntiles, float *xp, float *yp, float *zp, wavenums *p_wn)
{
    /* CPU version to calc array factor for a given wavevector */
    double kx,ky,kz;
    double ph;
    double complex af = 0.0+0.0*I;

    kx = p_wn->kx;
    ky = p_wn->ky;
    kz = p_wn->kz;
    for (int j = 0; j < ntiles; j++)
    {
        ph = (kx * xp[j]) + (ky * yp[j]) + (kz * zp[j]);
        af = af + (cos(ph) + 1.0*I*sin(ph));
    }
    af = af/ntiles;

    return af;
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
    double pr=0, pd=0, px=0, rv=0, eq=2000, ra_ap=0, dec_ap=0; // all for conversion to apparent RA/DEC
    char id_str[20];

    // read ra string into hours, minutes and seconds
    sscanf(ra_hhmmss, "%d:%d:%lf", &ra_ih, &ra_im, &ra_fs);

    //read dec string into degrees, arcmin and arsec (extra steps for sign, '+' or '-')
    sscanf(dec_ddmmss, "%s:%d:%lf", id_str, &dec_im, &dec_fs);
    sign = (id_str[0] == '-' ? -1 : 1); // check sign of dec

    sscanf(dec_ddmmss, "%d:%d:%lf", &dec_id, &dec_im, &dec_fs); // assign values
    dec_id = dec_id * sign; // ensure correct sign


    // convert angles to radians
    slaCtf2r(ra_ih, ra_im, ra_fs, &ra_rad, &ra_j); //right ascension
    slaDaf2r(dec_id, dec_im, dec_fs, &dec_rad, &dec_j); //declination

    if (ra_j != 0) 
    {
        fprintf(stderr,"Error parsing %s as hhmmss\nslalib error code: j=%d\n", ra_hhmmss, ra_j);
        fprintf(stderr,"ih = %d, im = %d, fs = %f\n", ra_ih, ra_im, ra_fs);
        exit(-1);
    }
    if (dec_j != 0) 
    {
        fprintf(stderr,"Error parsing %s as ddmmss\nslalib error code: j=%d\n", dec_ddmmss, dec_j);
        fprintf(stderr,"ih = %d, im = %d, fs = %f\n", dec_id, dec_im, dec_fs);
        exit(-1);
    }

    // convert UTC to MJD
    utc2mjd(time_utc, &intmjd, &fracmjd);
    mjd = intmjd + fracmjd;
    mjd2lst(mjd, &lmst);

    // get apparent RA and Dec of target
    slaMap(ra_rad, dec_rad, pr, pd, px, rv, eq, mjd, &ra_ap, &dec_ap);
    printf("RA = %.4f  RA_app = %.4f  DEC = %.4f  DEC_app = %.4f\n", ra_rad, ra_ap, dec_rad, dec_ap);

    // use RA and LST to get HA
    ha_rad = slaRanorm(lmst - ra_ap);

    // convert (HA, Dec) to (az, el)
    slaDe2h(ha_rad, dec_rad, MWA_LAT*DEG2RAD, &az, &el);

    printf("Az = %.4f  ZA = %.4f\n", az, PI/2-el);
    p_tazza->az = az;
    p_tazza->za = (PI/2) - el;
}

int getNumTiles(const char *metafits)
{
    /* Figure out the number of tiles based on the information in the metafits.

       NOTE: we get warnings from this in compilation because the library functions
       expect char characters, but conversion from string literals to chars is bad.
       It works, but we get warnings... */

    fitsfile *fptr=NULL;
    int status=0;
    size_t ninput=0;

    fits_open_file(&fptr, metafits, READONLY, &status); // open metafits file
    fits_movnam_hdu(fptr, BINARY_TBL, "TILEDATA", 0, &status); // move to TILEDATA HDU
    if (status != 0)
    {
        fprintf(stderr,"Error: Failed to move to TILEDATA HDU\n");
        exit(-1);
    }

    fits_read_key(fptr, TINT, "NAXIS2", &ninput, NULL, &status); // read how many tiles are included
    if (status != 0)
    {
        fprintf(stderr,"Error: Failed to read size of binary table in TILEDATA\n");
        exit(-1);
    }
    fits_close_file(fptr, &status);
    return ninput;
}

void getTilePositions(const char *metafits, int ninput,\
        float *n_pols, float *e_pols, float *h_pols,\
        float *n_tile, float *e_tile, float *h_tile)
{
    /* Get the tile positions from the metafits file.
       Accepts the metafits file name, 
       number of items to read (i.e. 2x number of tiles, 1 per polarisation),
       the array to fill with the polarisation locations, and
       the array to fill with the tile locations (every second element of *_pols).
       
       NOTE: we get warnings from this in compilation because the library functions
       expect char characters, but conversion from string literals to chars is bad.
       It works, but we get warnings... */

    fitsfile *fptr=NULL;
    int status=0, anynull=0;
    int colnum=0;


    fits_open_file(&fptr, metafits, READONLY, &status); // open metafits file
    fits_movnam_hdu(fptr, BINARY_TBL, "TILEDATA", 0, &status); // move to TILEDATA HDU
    if (status != 0) 
    {
        fprintf(stderr,"Error: Failed to move to TILEDATA HDU\n");
        exit(-1);
    }

    fits_get_colnum(fptr, 1, "North", &colnum, &status); // get north coordinates of tiles
    fits_read_col_flt(fptr, colnum, 1, 1, ninput, 0.0, n_pols, &anynull, &status);
    if (status != 0)
    {
        fprintf(stderr,"Error: Failed to read  N coord in metafile\n");
        exit(-1);
    }

    fits_get_colnum(fptr, 1, "East", &colnum, &status); // get east coordinates of tiles
    fits_read_col_flt(fptr, colnum, 1, 1, ninput, 0.0, e_pols, &anynull, &status);
    if (status != 0)
    {
        fprintf(stderr,"Error: Failed to read E coord in metafile\n");
        exit(-1);
    }

    fits_get_colnum(fptr, 1, "Height", &colnum, &status); // get height a.s.l. of tiles
    fits_read_col_flt(fptr, colnum, 1, 1, ninput, 0.0, h_pols, &anynull, &status);
    if (status != 0)
    {
        fprintf(stderr,"Error: Failed to read H coord in metafile\n");
        exit(-1);
    }
    fits_close_file(fptr, &status);

    // populate the tile arrays with every second element of the pol arrays
    for (int i = 0; i < ninput; i+=2)
    {
        n_tile[i/2] = n_pols[i];
        e_tile[i/2] = e_pols[i];
        h_tile[i/2] = h_pols[i];
    }

    // convert heights into height above array center
    for (int i = 0; i < ninput/2; i++)
    {
        h_tile[i] = h_tile[i] - MWA_HGT;
    }
}

int getFlaggedTiles(const char *badfile, int *badtiles)
{
    /* Open the flagged tiles file, read into an array and count how many lines are read.
       Update the array pointer and return number of elements to read from that array
       (as it's initialised to be able to hold every tile) */
    
    FILE *fp;
    int tile=0, i=0;
    int nlines=0;

    fp = fopen(badfile,"r");
    if (fp == NULL)
    {
        fprintf(stderr,"Error opening flagged tiles file.\n");
        exit(-1);
    }

    while(fscanf(fp, "%d\n", &tile) > 0)
    {
        printf("    bad tile: %d\n",tile);
        badtiles[i++] = tile;
        nlines++;
    }

    fclose(fp);
    return nlines;
}

void removeFlaggedTiles(float *n_tile, float *e_tile, float *h_tile,\
        int *badtiles, int nbad, int nelements)
{
    /* Get rid of the bad/flagged tiles from the array. Basically just
       shift all the indexes around so that whenever a bad tile is there,
       it's data is over-written. We end up with an array of the same size,
       but the last nbad elements are all identical (and can be ignored). */

    int counter=0,bidx=0;

    for (int b=0; b < nbad; b++)
    {
        // for each bad tile index in badtiles
        bidx = badtiles[b];
        for (int i=(bidx-counter); i < nelements-1; i++)
        {
            // shift each element in tile positions to the left by one
            // excluding the last element
            n_tile[i] = n_tile[i+1];
            e_tile[i] = e_tile[i+1];
            h_tile[i] = h_tile[i+1];
        }
        // array shifted left one, but the bad indexes refer to original tile positions
        // so we need to move the bad index to the left by one, too
        counter++;
    }
}


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

    int iters=1;
    int nDevices=0;
    long double azzaMem=0, tileposMem=0, afMem=0, otherMem=0, reqMem=0, devMem=0;
    double tfrac = 0.9;
    size_t freeMem, totMem;
    cudaError_t res;
    struct cudaDeviceProp prop; // create struct to store device info

    // get info about avilable devices
    getDeviceDimensions(&nDevices);
    printf("Number of devices on system: %d\n", nDevices);
    printf("Using device: %d\n", nDevices-1);
    cudaGetDeviceProperties(&prop, 0); // populate prop for this device
    
    res = cudaMemGetInfo(&freeMem, &totMem);
    if (res == cudaSuccess)
    {
        printf("Free device memory: %.2f MB\n", (double)freeMem/1.0e6);
    }
    else
    {
        printf("%s\n", cudaGetErrorString(res));
    }

    *blockSize = prop.maxThreadsDim[0];

    azzaMem = 2 * (size/1.0e6) * sizeof(double);
    tileposMem = 3 * (ntiles/1.0e6) * sizeof(float);
    afMem = (size/1.0e6) * sizeof(cuDoubleComplex);
    otherMem = (7 * sizeof(double) + 4 * sizeof(int) + sizeof(cuDoubleComplex) + sizeof(wavenums))/1.0e6;
    
    reqMem = azzaMem + tileposMem + afMem + otherMem;
    devMem = (double)freeMem/1.0e6;

    printf("Memory required for:\n");
    printf("    Az,ZA arrays: %Lf MB\n",azzaMem);
    printf("    tile positions: %Lf MB\n",tileposMem);
    printf("    array fator: %Lf MB\n",afMem);
    printf("    intermediate: %Lf MB\n",otherMem);
    printf("Total: %Lf MB\n",reqMem);

    if (reqMem < 0)
    {
        fprintf(stderr, "Negative required memory - integer overflow? Probably requesting resolution elements that are too small!!\n");
        exit(1);
    }
    else if ((tfrac*devMem) <= reqMem)
    {
        fprintf(stderr, "Arrays will not fit on device!\n");
        fprintf(stderr, "   total device memory required = %Lf MB\n", reqMem);
        fprintf(stderr, "   useable device memory        = %Lf MB\n", tfrac*devMem);
        fprintf(stderr, "       (useable fraction = %.2f of total)\n", tfrac);

        iters = (int)ceil(tfrac*reqMem / devMem)+1;
        printf("Will split task into: %d iterations (approx. %.2Lf MB per iteration)\n", iters, (reqMem/iters));
    }
    printf("\n");

    *niter = iters;
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
#define gpuErrchk(ans) {gpuAssert((ans), __FILE__, __LINE__);}


//__global__ void getWavevectors(int nel, double a, double *za, double *az, wavenums *p_twn, wavenums *p_wn)
//{
//    /* Kernal to compute the wavevectors for a list of theta/phi.
//       nel   = total number of elements in each array
//       a     = amplitude factor (i.e. 2pi/lambda)
//       za    = array of zenith angles
//       az    = array of azimuths
//       p_twn = target wavenumber struct
//       p_wn  = array of structs to contain kx, ky and kz
//
//       The standard equations are:
//       a = 2 * pi / lambda
//       kx = a * sin(theta) * cos(phi)
//       ky = a * sin(theta) * sin(phi)
//       kz = a * cos(theta)
//
//       Assuming that (theta,phi) are in the convention from Sutinjo et al. 2015:
//       i.e. phi = pi/2 - az   AND   theta = za
//
//       The azimuth is measured clockwise from East (standard for antenna theory, offset from astronomy)
//
//       Will produce (k - k_target) which is the wavenumber required for the array factor calculation.
//    */
//    int index = blockIdx.x * blockDim.x + threadIdx.x;
//    int stride = blockDim.x * gridDim.x;
//    double ast, phi;
//
//    for (int i = index; i < nel; i+= stride)
//    {
//        ast = a * sin(za[i]);
//        phi = PI/2 - az[i];
//        p_wn[i].kx = ast * cos(phi) - p_twn->kx; 
//        p_wn[i].ky = ast * sin(phi) - p_twn->ky;
//        p_wn[i].kz = ast - p_twn->kz;
//    } 
//    __syncthreads();
//}
//
//
//__global__ void getArrayFactor(int nel, int ntiles, float *xp, float *yp, float *zp, wavenums *p_wn, cuDoubleComplex *af)
//{
//    /* Kernal to compute the array factor from the tiles positions and wave numbers.
//       nel    = total number of elements in final array
//       ntiles =  number of tiles used to form tied-array beam
//       xp     = array of tile x-positions (East)
//       yp     = array of tile y-positions (North)
//       zp     = array of tile z-positions (above array centre)
//       p_wn   = array of wavenumber structs for each pixel in af
//       af     = array containing the complex valued array factor
//    */
//    int index = blockIdx.x * blockDim.x + threadIdx.x;
//    int stride = blockDim.x * gridDim.x;
//    double ph;
//    double kx, ky, kz;
//    cuDoubleComplex n = make_cuDoubleComplex(ntiles,0);
//
//    for (int i = index; i < nel; i+=stride)
//    {
//        kx = p_wn[i].kx;
//        ky = p_wn[i].ky;
//        kz = p_wn[i].kz;
//        af[i] = make_cuDoubleComplex(0,0);
//        for (int j = 0; j < ntiles; j++)
//        {
//            ph = (kx * xp[j]) + (ky * yp[j]) + (kz * zp[j]);
//            af[i] = cuCadd(af[i], make_cuDoubleComplex(cos(ph), sin(ph)));         
//        }
//        af[i] = cuCdiv(af[i], n);
//    }
//    __syncthreads();
//}
//

__global__ void calcArrayFactor(int nel, int ntiles, double a,
                                double *za, double *az, 
                                float *xp, float *yp, float *zp, 
                                wavenums *p_twn, 
                                cuDoubleComplex *af)
{
    /* Kernal which takes (az,za) coordinates and tile positions to 
       compute the array factor. This minimises data transfer between 
       host and device.
   
       #####################
       # General variables #
       #####################
       nel = total number of elements in each array

       ##########################
       # Wavenumber computation #
       ##########################
       a     = amplitude factor (i.e. 2pi/lambda)
       za    = array of zenith angles
       az    = array of azimuths
       p_twn = target wavenumber struct

       NOTE: The standard equations are:
             a = 2 * pi / lambda
             kx = a * sin(theta) * cos(phi)
             ky = a * sin(theta) * sin(phi)
             kz = a * cos(theta)

             where:
             lambda is the observing wavelength, and
             assuming that (theta,phi) are in the convention from Sutinjo et al. 2015:
             phi = pi/2 - az   AND   theta = za

             The azimuth is measured clockwise from East (standard for antenna theory, offset from astronomy)
       
       ############################
       # Array factor computation #
       ############################
       nel    = total number of elements in final array
       ntiles =  number of tiles used to form tied-array beam
       xp     = array of tile x-positions (East)
       yp     = array of tile y-positions (North)
       zp     = array of tile z-positions (above array centre)
       p_wn   = array of wavenumber structs for each pixel in af
       af     = array containing the complex valued array factor

       NOTE: The analytical form for this is:
             
                f(theta,phi;tza,taz) = (1/ntiles) * sum{n=1,n=ntiles}( conj(psi_n(tza,taz)) * psi_n(theta,phi) )

             where:
             (taz,tza) are the target source azimuth and zenith angle
             (theta,phi) are the azimuth and zenith angle pixels over which we evalute, theta=[0,90], phi=[0,360)
             
             psi_n is the planar wave front as detected by the "nth" tile, defined as
                
                psi_n(theta,phi) = exp[ (2*pi*I/lambda) * (x_n*k_x + y_n*k_y + z_n*k_z) ]
             
             where x_n is the x-coordinate of tile n (similarly for y_n, z_n), with k_x, k_y, k_z and lambda as defined above.      
    */

    // set up CUDA thread indexing
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    
    // other intermediate variables
    double ast=0, phi=0;
    double ph=0;
    double kx=0, ky=0, kz=0;
    cuDoubleComplex n = make_cuDoubleComplex(ntiles,0);

    
    for (int i = index; i < nel; i += stride)
    {
        // pre-calculate coefficients/transforms
        ast = a * sin(za[i]);
        phi = PI/2 - az[i];

        // calculate (k - k_target)
        kx = ast * cos(phi) - p_twn->kx; 
        ky = ast * sin(phi) - p_twn->ky;
        kz = ast - p_twn->kz;

        // initialise this pixel's array factor value
        af[i] = make_cuDoubleComplex(0,0);
        
        // calculate array factor contribution from each tile and sum
        for (int j = 0; j < ntiles; j++)
        {
            ph = (kx * xp[j]) + (ky * yp[j]) + (kz * zp[j]);
            af[i] = cuCadd(af[i], make_cuDoubleComplex(cos(ph), sin(ph)));         
        }

        // normalise the array factor
        af[i] = cuCdiv(af[i], n);
    }
    __syncthreads();

}

void calcArrayFactorCPU(int nel, int ntiles, double a,
                                double *za, double *az, 
                                float *xp, float *yp, float *zp, 
                                wavenums *p_twn, 
                                double complex *af)
{
    /* CPU, serial version of the kernal */ 
    
    // other intermediate variables
    double ast=0, phi=0;
    double ph=0;
    double kx=0, ky=0, kz=0;
    double complex n = ntiles+I*0;

    //printf("    entering outer loop\n");
    for (int i = 0; i < nel; i ++)
    {
        // pre-calculate coefficients/transforms
        //printf("    precompute\n");
        ast = a * sin(za[i]);
        phi = PI/2 - az[i];

        // calculate (k - k_target)
        //printf("    calc k-ktarg\n");
        kx = ast * cos(phi) - p_twn->kx; 
        ky = ast * sin(phi) - p_twn->ky;
        kz = ast - p_twn->kz;

        // initialise this pixel's array factor value
        af[i] = 0.0 + I*0.0;
        
        // calculate array factor contribution from each tile and sum
        //printf("    entering tile sum loop\n");
        for (int j = 0; j < ntiles; j++)
        {
            ph = (kx * xp[j]) + (ky * yp[j]) + (kz * zp[j]);
            af[i] += cos(ph) + I*sin(ph);         
        }

        // normalise the array factor
        //printf("    normalising\n");
        af[i] /= n;
    }
}






int main(int argc, char *argv[])
{
    char *ra="00:00:00.00";
    char *dec="-26:00:00.00";
    char *time="2017-06-19T11:13:00";
    char *metafits=NULL;
    char *flagfile=NULL;
    int c=0;
    int ntiles=0;
    int blockSize,numBlocks;
    int niter=1, iter=0;
    int n_az, n_za;
    long int size;
    double freq=0, lambda=0;
    double az_step=1.0, za_step=1.0, eta=1.0;
    //double ph=0.0, omega_A=0.0, af_max=-1.0, eff_area=0.0, gain=0.0;
    tazza target;
    wavenums target_wn;

    /* Parse options */
    if (argc > 1)
    {
        while ((c = getopt(argc, argv, "f:e:r:d:t:m:b:x:y:")) != -1)
        {
            switch(c)
            {
                case 'f':
                    freq = atof(optarg);
                    lambda = SOL/freq;
                    break;
                case 'e': 
                    eta = atof(optarg);
                    break;
                case 'r':
                    ra = strdup(optarg);
                    break;
                case 'd': 
                    dec = strdup(optarg);
                    break;
                case 't':
                    time = strdup(optarg);
                    break;
                case 'm':
                    metafits = strdup(optarg);
                    break;
                case 'b':
                    flagfile = strdup(optarg);
                    break;
                case 'x':
                    az_step = atof(optarg);
                    break;
                case 'y':
                    za_step = atof(optarg);
                    break;
                default:
                    usage();
                    exit(1);
            }
        }
    }

    if (argc == 1)
    {
        usage();
        exit(1);
    }


    // calculate target az,za and wavevector 
    printf("Getting target (Az,ZA)\n");
    calcTargetAZZA(ra, dec, time, &target);
    printf("Computing wavenumbers towards target\n");
    calcWaveNumber(lambda, target.az, target.za, &target_wn);
    printf("    kx = %f    ky = %f    kz = %f\n", target_wn.kx, target_wn.ky, target_wn.kz);
    printf("\n");

    // get the number of tiles in array
    printf("Determining number of tiles from metafits\n");
    ntiles = getNumTiles(metafits); // returns 2x the number of tiles, 1 per pol.
    ntiles = ntiles / 2;
    printf("    number of tiles: %d\n",ntiles);

    // allocate dynamic memory for intermediate tile position arrays
    // (probably don't need to check as this should be ~MB scale)
    float *N_pols = (float *)malloc(2 * ntiles * sizeof(float));
    float *E_pols = (float *)malloc(2 * ntiles * sizeof(float));
    float *H_pols = (float *)malloc(2 * ntiles * sizeof(float));
    // allocate dynamic memory for tile positions
    float *N_tile = (float *)malloc(ntiles * sizeof(float));
    float *E_tile = (float *)malloc(ntiles * sizeof(float));
    float *H_tile = (float *)malloc(ntiles * sizeof(float));
    printf("Getting tile positions\n");
    getTilePositions(metafits, 2*ntiles,\
            N_pols, E_pols, H_pols,\
            N_tile, E_tile, H_tile);
    free(N_pols);
    free(E_pols);
    free(H_pols);

    // have to remove tiles from the flagged tiles list.
    // each element in the list is the index of a tile that needs to be removed.
    printf("Getting flagged tiles\n");
    int *flagged_tiles = (int *)malloc(ntiles * sizeof(int));
    int ntoread;
    ntoread = getFlaggedTiles(flagfile, flagged_tiles);
    int flagged[ntoread];

    for (int i=0; i < ntoread; i++)
    {
        flagged[i] = flagged_tiles[i];
    }
    printf("Removing %d flagged tiles\n",ntoread);
    printf("Tiles remaining: %d\n",ntiles-ntoread);
    removeFlaggedTiles(N_tile, E_tile, H_tile, flagged, ntoread, ntiles);
    free(flagged_tiles);
    printf("\n");

    // but, the last ntoread elements are pointless 
    // so now we can allocate static memory for the final list of positions
    ntiles = ntiles - ntoread;
    float xpos[ntiles], ypos[ntiles], zpos[ntiles];

    for (int i=0; i<ntiles; i++)
    {
        // x = East, y = North, z = Height
        xpos[i] = E_tile[i];
        ypos[i] = N_tile[i];
        zpos[i] = H_tile[i];
    }
    free(N_tile);
    free(E_tile);
    free(H_tile);


    // determine number of az/za elements from specified pixel size
    n_az = (int)(360.0/az_step);
    n_za = (int)(90.0/za_step)+1; // +1 because we want to evalute at 90deg too!
    size = n_az * n_za;
    printf("Number of az steps [0,360): %d\n", n_az); // step from [0,360) - 360 will double count the 0 values
    printf("Number of za steps [0,90] : %d\n", n_za); // step from [0,90]
    printf("Total number of elements to compute: %ld\n", size);
    niter = 1; // how many times do I need to split the problem up?
    printf("\n");
 
    // figure out how many iterations are required (being conservative)
    // and the device properties (as a consequence)
    requiredMemory(size, ntiles, &niter, &blockSize);

    /* We now have the relevant array configuration and target source information 
       needed to calculate the array factor. The best way is to split it up into 
       managable chunks (depending on the device capabilities). */

      // construct arrays for computation on host
    double *az_array, *za_array;
    
    // allocate memory on host and check
    // azimuth vector
    az_array = (double *)calloc(size, sizeof(double));
    if (!az_array)
    {
        fprintf(stderr,"Host memory allocation failed (allocate az_array)\n");
        return EXIT_FAILURE;
    }
    // zenith vector
    za_array = (double *)calloc(size, sizeof(double));
    if (!za_array)
    {
        fprintf(stderr,"Host memory allocation failed (allocate za_array)\n");
        return EXIT_FAILURE;
    }

    // populate the host vectors::
    // TODO: this is currently the most memory intensive part: ~20GB at 0.01x0.01 resolution
    //       maybe we want to move this initilisation part into the iteration loop
    //       which will then make the arrays MUCH smaller - need to figure out how to 
    //       populate correctly then...
    printf("Initialising az, za and result matrix\n");
    // want arrays to be something like:
    // az = [0 0 0 0 ... 1 1 1 1 ...]
    // za = [0 1 2 3 ... 0 1 2 3 ...]
    int cc=0;
    int i=0;
    do
    {
        for (int j=0; j<n_za; j++)
        {
            za_array[cc+j] = j * za_step * DEG2RAD;
            az_array[cc+j] = i * az_step * DEG2RAD;
        }
        cc += n_za;
        i++;
    } while(cc < size);
     
//    FILE *fptest;
//    fptest = fopen("azza.dat","w");
//    for (int i=0; i<size; i++)
//    {
//        fprintf(fptest, "%f %f\n", az_array[i]*RAD2DEG, za_array[i]*RAD2DEG);
//    }
//    fclose(fptest);

    // construct arrays for device computation
    double *d_az_array, *d_za_array;
    float *d_xpos, *d_ypos, *d_zpos;
    //wavenums *wn_array, *d_wn_array, *d_twn;
    wavenums *d_twn;
    int itersize;
    double af_max = -1, omega_A = 0.0;
    double *subAz, *subZA;
    //wavenums *tmpsubwn; 
    cuDoubleComplex *af_array, *d_af_array; // these will be used within the iteration loop
    
    int iter_n_az = (int)floor(size / niter);
    int iter_n_za = (int)floor(size / niter);
    int az_idx1, az_idx2, za_idx1, za_idx2; 

    printf("iter n_az, n_za = %d %d\n", iter_n_az, iter_n_za);


//    subAz = (double *)malloc(itersize * sizeof(double));
//    if (!subAz)
//    {
//        fprintf(stderr,"Host memory allocation failed (allocate subAz)\n");
//        return EXIT_FAILURE;
//    }
//    subZA = (double *)malloc(itersize * sizeof(double));
//    if (!subZA)
//    {
//        fprintf(stderr,"Host memory allocation failed (allocate subZA)\n");
//        return EXIT_FAILURE;
//    }
//    wn_array = (wavenums *)malloc(itersize * sizeof(wavenums));
//    if (!wn_array)
//    {
//        fprintf(stderr,"Host memory allocation failed (allocate wn_array)\n");
//        return EXIT_FAILURE;
//    }
//    tmpsubwn = (wavenums *)malloc(sizeof(wavenums));
//
//    af_array = (cuDoubleComplex *)malloc(itersize * sizeof(cuDoubleComplex));
//    if (!af_array)
//    {
//        fprintf(stderr,"Host memory allocation failed (allocate af_array)\n");
//        return EXIT_FAILURE;
//    }
    
    // before we get to the real computation, better open a file ready for writing
    int obsid;
    char output[100];
    sscanf(metafits, "%d%*s", &obsid);
    printf("Will output beam pattern to:\n");
    printf("    %d_%.2fMHz_%s.dat\n", obsid, freq/1.0e6, time);
    sprintf(output, "%d_%.2fMHz_%s.dat", obsid, freq/1.0e6, time);
    
    FILE *fp;
    fp = fopen(output,"w");  // open the file to write
    fprintf(fp, "Az\tZA\tP\n"); // and write the header info

    /* This is the primary loop which does the calculations */
    printf("%d az , %d za per iteration\n", iter_n_az, iter_n_za);
    for (iter = 0; iter < niter; iter++)
    {  
        printf("==== Iteration %d ====\n", iter);
        // figure out this iteration size, then allocate memory
        if (iter != niter-1)
        {
            itersize = iter_n_az; // = iter_n_za

            az_idx1 = iter * iter_n_az;
            az_idx2 = (iter+1) * iter_n_az;
            
            za_idx1 = iter * iter_n_za;
            za_idx2 = (iter+1) * iter_n_za;
        }
        else
        {
            /* If this is the last iteration, construct 
               iter_n_az/za such that it includes what ever
               is left over to compute.
               
               Should be ok in terms of memory because we made
               the number of iterations was computed on a 
               conservative estimate of the device memory. */
            
            iter_n_za = size - (iter * iter_n_za);
            iter_n_az = size - (iter * iter_n_az);
            itersize =  iter_n_az; // = iter_n_za

            az_idx1 = iter * iter_n_az;
            az_idx2 = az_idx1 + itersize - 1;

            za_idx1 = iter * iter_n_za;
            za_idx2 = za_idx1 + itersize - 1;
        }

        printf("iteration size = %d\n", itersize);
        printf("n az: %d  n za: %d\n",iter_n_az,iter_n_za); 
        
        subAz = (double *)malloc(iter_n_az * sizeof(double));
        if (!subAz)
        {
            fprintf(stderr,"Host memory allocation failed (allocate subAz)\n");
            return EXIT_FAILURE;
        }
        subZA = (double *)malloc(iter_n_za * sizeof(double));
        if (!subZA)
        {
            fprintf(stderr,"Host memory allocation failed (allocate subZA)\n");
            return EXIT_FAILURE;
        }
        af_array = (cuDoubleComplex *)malloc(iter_n_az * sizeof(cuDoubleComplex));
        if (!af_array)
        {
            fprintf(stderr,"Host memory allocation failed (allocate af_array)\n");
            return EXIT_FAILURE;
        }

        numBlocks = (itersize + blockSize - 1) / blockSize; // number of blocks required 

        printf("azimuth idx: %d - %d\n", az_idx1, az_idx2);
        printf("zenith  idx: %d - %d\n", za_idx1, za_idx2);
        printf("Number of GPU blocks used: %d\n", numBlocks);
        
        //printf("testing:: omega_A = %f\n",omega_A);
        
        // place subset of az/za array into subAz/subZA
        for (int i=0; i<itersize; i++)
        {
            subAz[i] = az_array[i+az_idx1];
            subZA[i] = za_array[i+za_idx1];
            //wn_array[i].kx = 0.0;
            //wn_array[i].ky = 0.0;
            //wn_array[i].kz = 0.0;
            af_array[i] = make_cuDoubleComplex(0,0);

        }

        //printf("First/Last Az element for this iteration: %f %f\n", subAz[0]*RAD2DEG, subAz[itersize-1]*RAD2DEG);
        //printf("First/Last ZA element for this iteration: %f %f\n", subZA[0]*RAD2DEG, subZA[itersize-1]*RAD2DEG);

        // allocate memory on device
        gpuErrchk( cudaMalloc((void**)&d_az_array, itersize * sizeof(*az_array)));
        gpuErrchk( cudaMalloc((void**)&d_za_array, itersize * sizeof(*za_array)));
        gpuErrchk( cudaMalloc((void**)&d_twn, sizeof(wavenums)));
        gpuErrchk( cudaMalloc((void**)&d_xpos, ntiles * sizeof(*xpos)));
        gpuErrchk( cudaMalloc((void**)&d_ypos, ntiles * sizeof(*ypos)));
        gpuErrchk( cudaMalloc((void**)&d_zpos, ntiles * sizeof(*zpos)));
        gpuErrchk( cudaMalloc((void**)&d_af_array, itersize * sizeof(*af_array)));


        // copy arrays onto device
        gpuErrchk( cudaMemcpy(d_az_array, subAz, itersize * sizeof(*subAz), cudaMemcpyHostToDevice));
        gpuErrchk( cudaMemcpy(d_za_array, subZA, itersize * sizeof(*subZA), cudaMemcpyHostToDevice));
        gpuErrchk( cudaMemcpy(d_twn, &target_wn, sizeof(wavenums), cudaMemcpyHostToDevice));

        // copy the array factor vector to device
        gpuErrchk( cudaMemcpy(d_af_array, af_array, itersize * sizeof(*af_array), cudaMemcpyHostToDevice));
        
        // copy over tile position arrays to device
        gpuErrchk( cudaMemcpy(d_xpos, xpos, ntiles * sizeof(*xpos), cudaMemcpyHostToDevice));
        gpuErrchk( cudaMemcpy(d_ypos, ypos, ntiles * sizeof(*ypos), cudaMemcpyHostToDevice));
        gpuErrchk( cudaMemcpy(d_zpos, zpos, ntiles * sizeof(*zpos), cudaMemcpyHostToDevice));

        printf("Launching kernal to compute array factor\n");
        calcArrayFactor<<<numBlocks, blockSize>>>(itersize, ntiles, 2*PI/lambda, d_za_array, d_az_array, d_xpos, d_ypos, d_zpos, d_twn, d_af_array);

        // copy relevant memory back to host
        gpuErrchk( cudaMemcpy(af_array, d_af_array, itersize * sizeof(*af_array), cudaMemcpyDeviceToHost));
        printf("==== Done ====\n");


        // test the CPU equivalent to make sure we get the same numbers
        //printf("    comparing CPU to GPU:\n");
        //printf("    real: %f imag: %f abs: %f\n", af_array[12].x, af_array[12].y, cuCabs(af_array[12]));
        //printf("initialise cpu af array\n");
        //double complex *aftmp;
        //aftmp = (double complex *)malloc(itersize * sizeof(*aftmp));
        //printf("inttialise elements\n");
        //for (int i=0; i<itersize; i++)
        //{
        //    aftmp[i] = 0+I*0;
        //}
        //printf("call function\n");
        //calcArrayFactorCPU(itersize, ntiles, 2*PI/lambda, subZA, subAz, xpos, ypos, zpos, &target_wn, aftmp);
        //printf("    real: %f imag: %f abs: %f\n", creal(aftmp[12]), cimag(aftmp[12]), cabs(aftmp[12]));

        // cool, we're done with the GPU computation
        printf("Freeing device memory\n");
        gpuErrchk( cudaFree(d_twn));
        gpuErrchk( cudaFree(d_xpos));
        gpuErrchk( cudaFree(d_ypos));
        gpuErrchk( cudaFree(d_zpos));
        gpuErrchk( cudaFree(d_af_array));
        gpuErrchk( cudaFree(d_az_array));
        gpuErrchk( cudaFree(d_za_array));

        // test the CPU equivalent to make sure we get the same numbers
        //printf("Calculation array factor power (|af|^2)\n");
        //printf("    comparing CPU to GPU:\n");
        //double tmp1 = pow(cuCabs(af_array[12]), 2);
        //double tmp2 = pow(cabs(aftmp[12]), 2);
        //printf("    gpu power: %f\n", tmp1);
        //printf("    cpu power: %f\n", tmp2);


        /* Write the output to a file */
        //fprintf(fp, "iteration %d\n",iter);
        double af_power = 0.0;
        //double cpu_power = 0.0;
        printf("Writing to file...\n");
        for (int i=0; i<itersize; i++)
        {
            af_power = pow(cuCabs(af_array[i]), 2); // need to use cuCabs given af_array is of cuComplexDouble type
            //cpu_power = pow(cabs(aftmp[i]), 2);
            fprintf(fp, "%f\t%f\t%f\n", subAz[i]*RAD2DEG, subZA[i]*RAD2DEG, af_power);
            if (af_power > af_max) {af_max = af_power;}
            
            // integrate over sky
            omega_A = omega_A + sin(subZA[i]) * af_power * (za_step*DEG2RAD) * (az_step*DEG2RAD);
        }
        printf("Done -- freeing intermediate host memory\n");
        //free(aftmp);
        free(subAz);
        free(subZA);
        free(af_array);
    }
    printf("\n");
    printf("Closing file\n");
    fclose(fp); // close the file
    
    printf("Freeing host memory\n");
    free(az_array);
    free(za_array);

    // compute the gain and effective area from simulation
    double eff_area = 0.0, gain = 0.0;
    printf("Finished -- now computing relevant parameters:\n");
    eff_area = eta * pow(lambda, 2) * (4 * PI / omega_A);
    gain = (1.0e-26) * eff_area / (2 * KB);

    printf("    Array factor max:                 %f\n", af_max);
    printf("    Beam solid angle (sr):            %f\n", omega_A);
    printf("    Radiation efficiency:             %f\n", eta);
    printf("    Effective collecting area (m^2):  %.4f\n", eff_area);
    printf("    Effective array gain (K/Jy):      %.4f\n", gain);
   

    return 0;
}
