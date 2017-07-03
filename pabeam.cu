#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <slalib.h>
#include <fitsio.h>
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
     * Accepts the metafits file name, 
     * number of items to read (i.e. 2x number of tiles, 1 per polarisation),
     * the array to fill with the polarisation locations, and
     * the array to fill with the tile locations (every second element of *_pols) */
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
     * Update the array pointer and return number of elements to read from that array
     * (as it's initialised to be able to hold every tile) */
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
    /* We need to know how many devices are available and how to split the 
     * problem up to maximise the occupancy of each device. */
    printf("Querying system for device information --\n");
    cudaGetDeviceCount(nDevices); // get CUDA to count GPUs

    for (int i = 0; i < *nDevices; i++)
    {
        struct cudaDeviceProp prop; // create struct to store device info
        cudaGetDeviceProperties(&prop, i); // populate prop for this device
        printf("    Device number:                       %d\n",*nDevices-1);
        printf("    Total global memory available (MB):  %f\n",prop.totalGlobalMem/1e6);
        printf("    Max grid size (# blocks):           (%d, %d, %d)\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
        printf("    Max number of threads per block:    (%d, %d, %d)\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
    }
}


//static const char *_cudaGetErrorEnum(cublasStatus_t error)
//{
//    /* Function to convert cublas error to string */
//    switch (error)
//    {
//        case CUBLAS_STATUS_SUCCESS:
//            return "CUBLAS_STATUS_SUCCESS";
//
//        case CUBLAS_STATUS_NOT_INITIALIZED:
//            return "CUBLAS_STATUS_NOT_INITIALIZED";
//
//        case CUBLAS_STATUS_ALLOC_FAILED:
//            return "CUBLAS_STATUS_ALLOC_FAILED";
//
//        case CUBLAS_STATUS_INVALID_VALUE:
//            return "CUBLAS_STATUS_INVALID_VALUE";
//
//        case CUBLAS_STATUS_ARCH_MISMATCH:
//            return "CUBLAS_STATUS_ARCH_MISMATCH";
//
//        case CUBLAS_STATUS_MAPPING_ERROR:
//            return "CUBLAS_STATUS_MAPPING_ERROR";
//
//       case CUBLAS_STATUS_EXECUTION_FAILED:
//            return "CUBLAS_STATUS_EXECUTION_FAILED";
//
//        case CUBLAS_STATUS_INTERNAL_ERROR:
//            return "CUBLAS_STATUS_INTERNAL_ERROR";
//    }
//
//    return "<unknown>";
//}

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
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


__global__ void getWavevectors(int nel, double a, double *za, double *az, wavenums *p_twn, wavenums *p_wn)
{
    /* Kernal to compute the wavevectors for a list of theta/phi.
       nel   = total number of elements in each array
       a     = amplitude factor (i.e. 2pi/lambda)
       za    = array of zenith angles
       az    = array of azimuths
       p_twn = target wavenumber struct
       p_wn  = array of structs to contain kx, ky and kz

       The standard equations are:
       a = 2 * pi / lambda
       kx = a * sin(theta) * cos(phi)
       ky = a * sin(theta) * sin(phi)
       kz = a * cos(theta)

       Assuming that (theta,phi) are in the convention from Sutinjo et al. 2015:
       i.e. phi = pi/2 - az   AND   theta = za

       The azimuth is measured clockwise from East (standard for antenna theory, offset from astronomy)

       Will produce (k - k_target) which is the wavenumber required for the array factor calculation.
    */
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    double ast, phi;

    for (int i = index; i < nel; i+= stride)
    {
        ast = a * sin(za[i]);
        phi = PI/2 - az[i];
        p_wn[i].kx = ast * cos(phi) - p_twn->kx; 
        p_wn[i].ky = ast * sin(phi) - p_twn->ky;
        p_wn[i].kz = ast - p_twn->kz;
    } 
    __syncthreads();
}


__global__ void getArrayFactor(int nel, int ntiles, float *xp, float *yp, float *zp, wavenums *p_wn, cuDoubleComplex *af)
{
    /* Kernal to compute the array factor from the tiles positions and wave numbers.
       nel    = total number of elements in final array
       ntiles =  number of tiles used to form tied-array beam
       xp     = array of tile x-positions (East)
       yp     = array of tile y-positions (North)
       zp     = array of tile z-positions (above array centre)
       p_wn   = array of wavenumber structs for each pixel in af
       af     = array containing the complex valued array factor
    */
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    double ph;
    double kx, ky, kz;
    cuDoubleComplex n = make_cuDoubleComplex(ntiles,0);

    for (int i = index; i < nel; i+=stride)
    {
        kx = p_wn[i].kx;
        ky = p_wn[i].ky;
        kz = p_wn[i].kz;
        af[i] = make_cuDoubleComplex(0,0);
        for (int j = 0; j < ntiles; j++)
        {
            ph = (kx * xp[j]) + (ky * yp[j]) + (kz * zp[j]);
            af[i] = cuCadd(af[i], make_cuDoubleComplex(cos(ph), sin(ph)));         
        }
        af[i] = cuCdiv(af[i], n);
        __syncthreads();
    }
    __syncthreads();
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
    int nDevices=0,blockSize,numBlocks;
    int n_az, n_za, size, niter=1, iter=0;
    double reqMem=0.0, devMem=0.0;
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
    printf("    kx = %f    ky = %f    kz = %f\n",target_wn.kx,target_wn.ky,target_wn.kz);

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
    // each row in the list is the index of the tile that needs to be removed.
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
    n_za = (int)(90.0/za_step);
    size = n_az * n_za;
    printf("Number of az steps: %d\n",n_az); // step from [0,360)
    printf("Number of za steps: %d\n",n_za); // step from [0,90]
    niter = 1; // how many times do I need to split the problem up?


    // get info about avilable devices
    getDeviceDimensions(&nDevices);
    printf("Number of devices on system: %d\n",nDevices);
    printf("Using device: %d\n", nDevices-1);
    struct cudaDeviceProp prop; // create struct to store device info
    cudaGetDeviceProperties(&prop, 0); // populate prop for this device

    // estimate how much device memory will be required
    //      az array = n_az * n_za * sizeof(double) - array of azimuth coords
    //      za array = n_za * n_az * sizeof(double) - array of zenith coords
    //      wv array = n_az * n_za * sizeof(wavenums) - 3D wavenumbers (structs)
    //      xpos = ntiles * sizeof(double) - x position of tiles from array center
    //      ypos = ntiles * sizeof(double) - y position of tiles from array center
    //      zpos = ntiles * sizeof(double) - z position of tiles from array center
    //      af_array = n_az * n_za * sizeof(double) - array factor 
    //      
    // so in total, will need:
    //       3 * n_az * n_za * sizeof(double) + n_az * n_za * sizeof(wavenums) + 3 * ntiles * sizeof(double)
    reqMem = (double)(3 * size * sizeof(double)) + (size * sizeof(wavenums)) + (3 * ntiles * sizeof(double));
    devMem = (double)prop.totalGlobalMem;

    if (devMem <= reqMem)
    {
        fprintf(stderr, "Matrices will not fit on device!\n");
        fprintf(stderr, "   total device memory required  = %f MB\n", reqMem/1.0e6);
        fprintf(stderr, "   device global memory          = %f MB\n", devMem/1.0e6);

        niter = (int)ceil(reqMem / devMem)+1;
        printf("Will split task into: %d (%f MB per task)\n", niter, (reqMem/(1.0e6*niter)));
    }
  
    blockSize = prop.maxThreadsDim[0]; // number of threads available per block, in x-direction

   // printf("%d %d\n",blockSize,numBlocks);



    /* We now have the relevant array configuration and target source information 
       needed to calculate the array factor. The best way is to use a matrix approach
       and split it up into managable chunks (depending on the device capabilities). */

    int iter_n_az = (int)floor(n_az / niter);
    int iter_n_za = (int)floor(n_za / niter);
    int az_idx1, az_idx2, za_idx1, za_idx2;

    // construct arrays for computation on host
    double *az_array, *za_array;
    
    // allocate memory on host and check
    // azimuth vector
    az_array = (double *)malloc(size * sizeof(double));
    if (!az_array)
    {
        fprintf(stderr,"Host memory allocation failed (allocate az_array)\n");
        return EXIT_FAILURE;
    }
    // zenith vector
    za_array = (double *)malloc(size * sizeof(double));
    if (!za_array)
    {
        fprintf(stderr,"Host memory allocation failed (allocate za_array)\n");
        return EXIT_FAILURE;
    }

    // populate the host vectors::
    //      this is currently the most memory intensive part: ~20GB at 0.01x0.01 resolution
    //      maybe we want to move this initilisation part into the iteration loop
    //      which will then make the arrays MUCH smaller - need to figure out how to 
    //      populate correctly then...
    printf("Initialising az, za and result matrix\n");
    // want arrays to be something like:
    // az = [0 0 0 0 ... 1 1 1 1 ...]
    // za = [0 1 2 3 ... 0 1 2 3 ...]
    int cc=0;
    do
    {
        for (int j=0; j<n_za; j++)
        {
            za_array[cc+j] = j * za_step * DEG2RAD;
        }
        cc += n_za;
    } while(cc < size);
    
    cc=0;
    int i=0;
    do
    {
        for (int j=0; j<n_az; j++)
        {
            az_array[cc+j] = i * az_step * DEG2RAD;
        }
        cc += n_az;
        i++;
    } while(cc < size);


    // construct arrays for device computation
    double *d_az_array, *d_za_array;
    float *d_xpos, *d_ypos, *d_zpos;
    wavenums *wn_array, *d_wn_array, *d_twn;
    int itersize = iter_n_az * iter_n_za;
    double af_max = -1, omega_A = 0.0;
    double *subAz, *subZA;
    wavenums *tmpsubwn; 
    cuDoubleComplex *af_array, *d_af_array; // these will be used within the iteration loop

    subAz = (double *)malloc(itersize * sizeof(double));
    if (!subAz)
    {
        fprintf(stderr,"Host memory allocation failed (allocate subAz)\n");
        return EXIT_FAILURE;
    }
    subZA = (double *)malloc(itersize * sizeof(double));
    if (!subZA)
    {
        fprintf(stderr,"Host memory allocation failed (allocate subZA)\n");
        return EXIT_FAILURE;
    }
    wn_array = (wavenums *)malloc(itersize * sizeof(wavenums));
    if (!wn_array)
    {
        fprintf(stderr,"Host memory allocation failed (allocate wn_array)\n");
        return EXIT_FAILURE;
    }
    tmpsubwn = (wavenums *)malloc(sizeof(wavenums));

    af_array = (cuDoubleComplex *)malloc(itersize * sizeof(cuDoubleComplex));
    if (!af_array)
    {
        fprintf(stderr,"Host memory allocation failed (allocate af_array)\n");
        return EXIT_FAILURE;
    }
    for (int i=0; i<itersize; i++)
    {
        af_array[i] = make_cuDoubleComplex(0,0);
    }

    // before we get to the real computation, better open a file ready for writing
    int obsid;
    char output[100];
    sscanf(metafits, "%d%*s", &obsid);
    printf("Will output beam pattern to:\n");
    printf("    %d_%.2fMHz_%s.dat\n", obsid, freq/1.0e6, time);
    sprintf(output, "%d_%.2fMHz_%s.dat", obsid, freq/1.0e6, time);
    FILE *fp;
    fp = fopen(output,"w");
  
    /* This is the primary loop which does the calculations */
    printf("%d az , %d za per iteration\n", iter_n_az, iter_n_za);
    int stride = 0;
    for (iter = 0; iter < niter; iter++)
    {  
        numBlocks = ((iter_n_az-1) * iter_n_za + blockSize - 1) / blockSize; // number of blocks required 
        printf("==== Iteration %d ====\n", iter);
        az_idx1 = iter * iter_n_az;
        az_idx2 = (iter+1) * iter_n_az - 1; // -1 to zero based for output: e.g. 4000 elements = 0--3999 = 0--(iter_n_az-1)
        za_idx1 = iter * iter_n_za;
        za_idx2 = (iter+1) * iter_n_za - 1; // -1 to zero based for output: e.g. 1000 elements = 0--999 = 0--(iter_n_za-1)

        printf("azimuth idx: %d - %d\n", az_idx1, az_idx2);
        printf("zentih  idx: %d - %d\n", za_idx1, za_idx2);
        printf("Number of GPU blocks used: %d\n", numBlocks);
    
        // place subset of az/za array into subAz/subZA
        for (int i=0; i<itersize; i++)
        {
            subAz[i] = az_array[i+stride];
            subZA[i] = za_array[i+stride];
            wn_array[i].kx = 0.0;
            wn_array[i].ky = 0.0;
            wn_array[i].kz = 0.0;
        }
        stride += itersize;

        /* Calculate the wavenumbers for each pixel */
        // allocate memory on device
        gpuErrchk( cudaMalloc((void**)&d_az_array, itersize * sizeof(*az_array)));
        gpuErrchk( cudaMalloc((void**)&d_za_array, itersize * sizeof(*za_array)));
        gpuErrchk( cudaMalloc((void**)&d_wn_array, itersize * sizeof(*wn_array)));
        gpuErrchk( cudaMalloc((void**)&d_twn, sizeof(wavenums)));

        // copy arrays onto device
        gpuErrchk( cudaMemcpy(d_az_array, subAz, itersize * sizeof(*subAz), cudaMemcpyHostToDevice));
        gpuErrchk( cudaMemcpy(d_za_array, subZA, itersize * sizeof(*subZA), cudaMemcpyHostToDevice));
        gpuErrchk( cudaMemcpy(d_wn_array, wn_array, itersize * sizeof(*wn_array), cudaMemcpyHostToDevice));
        gpuErrchk( cudaMemcpy(d_twn, &target_wn, sizeof(wavenums), cudaMemcpyHostToDevice));

        // launch kernal to compute wavenumbers
        printf("Launching kernal to compute wave vectors\n");
        getWavevectors<<<numBlocks, blockSize>>>(itersize, 2*PI/lambda, d_za_array, d_az_array, d_twn, d_wn_array);
        cudaDeviceSynchronize(); // don't let host code move on until device is finished

        // copy relevant memory back to host
        gpuErrchk( cudaMemcpy(wn_array, d_wn_array, itersize * sizeof(*wn_array), cudaMemcpyDeviceToHost));

        // test the CPU equivalent to make sure we get the same numbers
        //printf("    comparing GPU and CPU output\n");
        //calcWaveNumber(lambda, subAz[12], subZA[12], tmpsubwn);
        //printf("    %f %f %f\n", wn_array[12].kx, wn_array[12].ky, wn_array[12].kz);
        //printf("    %f %f %f\n", tmpsubwn->kx-target_wn.kx, tmpsubwn->ky-target_wn.ky, tmpsubwn->kz-target_wn.kz);

        // done for now - free the memory on the device
        gpuErrchk( cudaFree(d_az_array));
        gpuErrchk( cudaFree(d_za_array));

        
        /* Calculate the array factor for each pixel */
        // allocate memory on device
        gpuErrchk( cudaMalloc((void**)&d_xpos, ntiles * sizeof(*xpos)));
        gpuErrchk( cudaMalloc((void**)&d_ypos, ntiles * sizeof(*ypos)));
        gpuErrchk( cudaMalloc((void**)&d_zpos, ntiles * sizeof(*zpos)));
        gpuErrchk( cudaMalloc((void**)&d_af_array, itersize * sizeof(*af_array)));

        // copy the calculated wavenumbers back onto the device
        gpuErrchk( cudaMemcpy(d_wn_array, wn_array, itersize * sizeof(*wn_array), cudaMemcpyHostToDevice));

        // copy the array factor vector to device
        gpuErrchk( cudaMemcpy(d_af_array, af_array, itersize * sizeof(*af_array), cudaMemcpyHostToDevice));
        
        // copy over tile position arrays to device
        gpuErrchk( cudaMemcpy(d_xpos, xpos, ntiles * sizeof(*xpos), cudaMemcpyHostToDevice));
        gpuErrchk( cudaMemcpy(d_ypos, ypos, ntiles * sizeof(*ypos), cudaMemcpyHostToDevice));
        gpuErrchk( cudaMemcpy(d_zpos, zpos, ntiles * sizeof(*zpos), cudaMemcpyHostToDevice));

        printf("Launching kernal to compute array factor\n");
        getArrayFactor<<<numBlocks, blockSize>>>(itersize, ntiles, d_xpos, d_ypos, d_zpos, d_wn_array, d_af_array);
        cudaDeviceSynchronize();

        // copy relevant memory back to host
        gpuErrchk( cudaMemcpy(af_array, d_af_array, itersize * sizeof(*af_array), cudaMemcpyDeviceToHost));

        //printf("real: %f imag: %f abs: %f\n", af_array[itersize/2].x, af_array[itersize/2].y, cuCabs(af_array[itersize/2]));

        // cool, we're done with the GPU computation
        gpuErrchk( cudaFree(d_twn));
        gpuErrchk( cudaFree(d_wn_array));
        gpuErrchk( cudaFree(d_xpos));
        gpuErrchk( cudaFree(d_ypos));
        gpuErrchk( cudaFree(d_zpos));
        gpuErrchk( cudaFree(d_af_array));

        /* Write the output to a file */
        double af_power = 0.0;
        fprintf(fp, "Az\tZA\tP\n");
        for (int i=0; i<itersize; i++)
        {
            af_power = pow(cuCabs(af_array[i]), 2);
            fprintf(fp, "%f\t%f\t%f\n", subAz[i], subZA[i], af_power);
            if (af_power > af_max) {af_max = af_power;}
            omega_A += sin(subZA[i]) * af_power * (za_step*DEG2RAD) * (az_step*DEG2RAD);
        }
    }
    fclose(fp);
    /* TODO: is the cudaMalloc/cudaFree calls impacting performance drastically? they are relatively expensive... 
            Actually it seems that the kernals themselves are the longest, and then the synchronisations
     */
    free(az_array);
    free(za_array);
    free(af_array);
    free(wn_array);
    free(subAz);
    free(subZA);
    free(tmpsubwn);

    double eff_area = 0.0, gain = 0.0;
    printf("Finished -- now computing relevant parameters:\n");
    eff_area = eta * pow(lambda, 2) * (4 * PI / omega_A);
    gain = (1e-26) * eff_area / (2 * KB);

    printf("   Array factor max:                 %f\n",af_max);
    printf("   Beam solid angle (sr):            %f\n",omega_A);
    printf("   Effective collecting area (m^2):  %.4f\n",eff_area);
    printf("   Effective array gain (K/Jy):      %.4f\n",gain);
   

    return 0;
}
