#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <slalib.h>
#include <fitsio.h>

#define PI (acos(-1.0))         //Ensures PI is defined on all systems
#define RAD2DEG (180.0 / PI)
#define DEG2RAD (PI / 180.0)
#define SOL (299792458.0)       //Speed of light
#define KB (1.38064852e-23)     //Boltzmann's constant

#define MWA_LAT (-26.703319)    //Array latitude, degrees North
#define MWA_LON (116.67081)     //Array longitude, degrees East
#define MWA_HGT (377.827)       //Array elevation above sea level, in meters

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
void calcWaveNumber(double lambda, double phi, double theta, wavenums *p_wn);
void calcTargetAZZA(char *ra_hhmmss, char *dec_ddmmss, char *time_utc, tazza *p_tazza);
int getNumTiles(char *metafits);
void getTilePositions(char *metafits, int ninput,\
        float *n_pols, float *e_pols, float *h_pols,\
        float *n_tile, float *e_tile, float *h_tile);
int getFlaggedTiles(char *badfile, int *badtiles);
void removeFlaggedTiles(float *n_tile, float *e_tile, float *h_tile,\
        int *badtiles, int nbad, int nelements);



void usage()
{
    printf("pabeam -- computes the array factor that represents the MWA tied-array beam for a given observation and pointing\n");
    printf("syntax:\n\
    pabeam -f <frequency (Hz)> -r <ra in hh:mm:ss> -d <dec in dd:mm:ss> -t <utc time string> -m <metafits file> -b <RTS flagged_tiles.txt file>\n\n\
    Options:\n\
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

void calcWaveNumber(double lambda, double phi, double theta, wavenums *p_wn)
{
    /* Calculate the 3D wavenumbers for a given wavelength (lambda) from the direction (az,za).
     * Accepts wavelength (m), azimuth (rad) and zenith angle (rad) and a pointer to a wavenums struct to populate.*/
    double a, ast;

    a = 2 * PI / lambda;
    ast = a * sin(theta);

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

int getNumTiles(char *metafits)
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

void getTilePositions(char *metafits, int ninput,\
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

int getFlaggedTiles(char *badfile, int *badtiles)
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

void removeFlaggedTiles(float *n_tile, float *e_tile, float *h_tile, int *badtiles, int nbad, int nelements)
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



int main(int argc, char *argv[])
{
    char *ra="00:00:00.00";
    char *dec="-26:00:00.00";
    char *time="2017-06-19T11:13:00";
    char *metafits=NULL;
    char *flagfile=NULL;
    int c, ntiles=0;
    double freq=0, lambda=0;
    double az_step=1.0, za_step=1.0, eta=1.0;
    double ph, omega_A, af_max, eff_area, gain;
    double complex af;
    tazza target;
    wavenums wn, target_wn;


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

    

    // copy test values into appropriate variables
    //strcpy(ra,"05:34:31.97");
    //strcpy(dec,"+22:00:52.06");
    //strcpy(time,"2014-11-07T16:53:20");
    //strcpy(metafits,"1099414416_metafits_ppds.fits");
    //strcpy(flagfile,"flagged_tiles.txt");

    printf("Getting target (Az,ZA)\n");
    calcTargetAZZA(ra, dec, time, &target);
    printf("Computing wavenumbers towards target\n");
    calcWaveNumber(lambda, PI/2-target.az, target.za, &target_wn);

    // Have to get tile positions and remove the tiles that were flagged
    printf("Determining number of tiles from metafits\n");
    ntiles = getNumTiles(metafits); // returns 2x the number of tiles, 1 per pol.
    ntiles = ntiles / 2;
    printf("    number of tiles: %d\n",ntiles); 

    // allocate dynamic memory for intermediate tile position arrays
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
    printf("    removing %d flagged tiles\n",ntoread);
    removeFlaggedTiles(N_tile, E_tile, H_tile, flagged, ntoread, ntiles);

    // but, the last ntoread elements are pointless 
    // so now we can allocate static memory for the final list of positions
    float xpos[ntiles-ntoread], ypos[ntiles-ntoread], zpos[ntiles-ntoread];

    for (int i=0; i<(ntiles-ntoread); i++)
    {
        // x = East, y = North, z = Height
        xpos[i] = E_tile[i];
        ypos[i] = N_tile[i];
        zpos[i] = H_tile[i];
    }   
    free(N_tile);
    free(E_tile);
    free(H_tile);
  
     

    FILE *fp;
    ph = 0.0;
    omega_A = 0.0;
    eff_area = 0.0;
    af_max = -1.0;

    fp = fopen("test.dat","w");
    for (double az = 0.0; az < 360.0; az += az_step)
    {
        /* one loop is ok, but we'll want to vectorise the next parts... */ 
        for (double za = 0.0; za <= 90.0; za += za_step)
        {
            af = 0.0 + 0.0*I;
            calcWaveNumber(lambda, PI/2-(az*DEG2RAD), za*DEG2RAD, &wn);
            for (int i = 0; i < ntiles; i++)
            {
                // af = exp(i*k.r)exp(i*k'.r)* = exp(i*[k-k'].r)
                //    = cos([k-k'].r) + i*sin([k-k'].r)
                // where k' is the target pointing wavenumbers, thus
                // af is maximised when pointing directly at the target
                ph = (wn.kx - target_wn.kx) * xpos[i] +\
                     (wn.ky - target_wn.ky) * ypos[i] +\
                     (wn.kz - target_wn.kz) * zpos[i];
                af += cos(ph) + 1.0*I*sin(ph); // sum over all tiles
            }
            af = af / ntiles; // normalise it
            
            // keep array factor maximum up-to-date (booking keeping)
            if (pow(cabs(af),2) > af_max) {af_max = pow(cabs(af),2);}
            fprintf(fp, "%f\t%f\t%f\n", az*DEG2RAD, za*DEG2RAD, pow(cabs(af),2));

            // calculate this pixel's contribution to beam solid angle
            omega_A = omega_A + sin(za*DEG2RAD) * powl(cabs(af),2) * (az_step*DEG2RAD) * (za_step*DEG2RAD);

            printf("\rComputing array factor: %.1f%%",(az/360.0)*100); fflush(stdout);
        } 
    }
    printf("\n");
    printf("Finished -- now computing relevant parameters:\n");
    eff_area = eta * pow((SOL / freq),2) * (4 * PI / omega_A);
    gain = (1e-26) * eff_area / (2 * KB);

    printf("   Array factor max:                 %f\n",af_max);
    printf("   Beam solid angle (sr):            %f\n",omega_A);
    printf("   Effective collecting area (m^2):  %.4f\n",eff_area);
    printf("   Effective array gain (K/Jy):      %.4f\n",gain);
    

    return 0;
}
