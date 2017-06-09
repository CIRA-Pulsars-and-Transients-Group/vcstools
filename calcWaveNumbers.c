#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <slalib.h>

#define PI (acos(-1.0))
#define RAD2DEG (180.0 / PI)
#define DEG2RAD (PI / 180.0)

#define MWA_LAT -26.703319  //Array latitude, degrees North
#define MWA_LON 116.67081   //Array longitude, degrees East
#define MWA_HGT 377.827     //Array elevation above sea level, in meters

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
} tazza; // can just refer to this struct as type tazza


/* Define all the function prototypes */
void calcWaveNumber(double lambda, double az, double za, wavenums* p_wn);
void calcTargetAZZA(char* ra, char* dec);



void calcWaveNumber(double lambda, double az, double za, wavenums* p_wn)
{
    /* Calculate the 3D wavenumbers for a given wavelength (lambda) from the direction (az,za).
     * Accepts wavelength (m), azimuth (deg) and zenith angle (deg) and a pointer to a wavenums struct to populate.*/
    double a, a_sinTheta;

    a = 2 * PI / lambda;
    a_sinTheta = a * sin(za * DEG2RAD);

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
    p_wn->kx = a_sinTheta * cos(az); 
    p_wn->ky = a_sinTheta * sin(az); 
    p_wn->kz = a_sinTheta;   
}

void calcTargetAZZA(char* ra_hhmmss, char* dec_ddmmss)
{
    int ra_ih, ra_im, ra_j;
    int dec_id, dec_im, dec_j;
    int sign;
    double ra_rad, ra_fs;
    double dec_rad, dec_fs;
    char id_str[20];

    // read ra string into hours, minutes and seconds
    sscanf(ra_hhmmss,"%d:%d:%lf",&ra_ih,&ra_im,&ra_fs);

    //read dec string into degrees, arcmin and arsec (extra steps for sign, '+' or '-')
    sscanf(dec_ddmmss,"%s:%d:%lf",id_str,&dec_im,&dec_fs);
    if (id_str[0] == '-') 
    {
        sign = -1;
    }
    else 
    {
        sign = 1;
    }
    sscanf(dec_ddmmss,"%d:%d:%lf",&dec_id,&dec_im,&dec_fs); // assign values
    dec_id = dec_id * sign; // ensure correct sign

    
    // convert angles to radians
    slaCtf2r(ra_ih, ra_im, ra_fs, &ra_rad, &ra_j); //right ascension
    slaDaf2r(dec_id, dec_im, dec_fs, &dec_rad, &dec_j); //declination

    if (ra_j != 0) 
    {
        fprintf(stderr,"Error parsing %s as hhmmss\nslalib error code: j=%d\n",ra_hhmmss,ra_j);
        fprintf(stderr,"ih = %d, im = %d, fs = %lf\n", ra_ih, ra_im, ra_fs);
        exit(-1);
    }

    if (dec_j != 0) 
    {
        fprintf(stderr,"Error parsing %s as ddmmss\nslalib error code: j=%d\n",dec_ddmmss,dec_j);
        fprintf(stderr,"ih = %d, im = %d, fs = %lf\n", dec_id, dec_im, dec_fs);
        exit(-1);
    }

    printf("%f %f\n",ra_rad,dec_rad);
}



int main()
{
    char ra[64], dec[64];
    double lambda, az_step, za_step;
    wavenums wn; 

    lambda = 3.0e8/184.96e6;

    az_step = 1;
    za_step = 1;

    // copy RA and DEC coords into ra, dec variables
    strcpy(ra,"06:34:31.97");
    strcpy(dec,"+22:00:52.06");

    for (double az = 0.0; az < 360.0; az += az_step)
    {
        for (double za = 0.0; za < 90.0; za += za_step)
        {
            calcWaveNumber(lambda, az, za, &wn);
            printf("%f %f %f\n", wn.kx, wn.ky, wn.kz);
            calcTargetAZZA(ra, dec);
        }
    }

    return 0;
}
