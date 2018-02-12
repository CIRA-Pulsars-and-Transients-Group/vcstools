// Standard includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <slalib.h>
#include <fitsio.h>

// MWA tile beam
#include "FEE2016/beam2016implementation.h"
#include "FEE2016/mwa_beam_interface.h"
#include "FEE2016/system.h"
#include "H5Cpp.h"


#define PI (acos(-1.0))         // Ensures PI is defined on all systems
#define RAD2DEG (180.0 / PI)
#define DEG2RAD (PI / 180.0)
#define SOL (299792458.0)       // Speed of light
#define KB (1.38064852e-23)     // Boltzmann's constant

#define MWA_LAT (-26.703319)    // Array latitude, degrees North
#define MWA_LON (116.67081)     // Array longitude, degrees East
#define MWA_HGT (377.827)       // Array elevation above sea level, in meters


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

