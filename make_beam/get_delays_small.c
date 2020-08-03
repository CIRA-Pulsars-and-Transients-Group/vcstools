/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include "mycomplex.h"
#include "slalib.h"
#include "slamac.h"
#include "psrfits.h"
#include "fitsio.h"
#include <string.h>
#include "beam_common.h"
#include "mwa_hyperbeam.h"

/* make a connection to the MWA database and get the antenna positions.
 * Then: calculate the geometric delay to a source for each antenna
 *
 * Initially geocentric to match as best as possible the output of calc.
 *
 */


/* This now also reads in RTS calibration files and generates calibration matrices
 */

#define MWA_LAT -26.703319        // Array latitude. degrees North
#define MWA_LON 116.67081         // Array longitude. degrees East
#define MWA_HGT 377               // Array altitude. meters above sea level
#define MAXREQUEST 3000000
#define NCHAN   128
#define NANT    128
#define NPOL    2
#define VLIGHT 299792458.0        // speed of light. m/s
double arr_lat_rad=MWA_LAT*(M_PI/180.0),arr_lon_rad=MWA_LON*(M_PI/180.0),height=MWA_HGT;

/* these externals are needed for the mwac_utils library */
int nfrequency;
int npol;
int nstation;
//=====================//

double parse_dec( char* dec_ddmmss ) {
/* Parse a string containing a declination in dd:mm:ss format into
 * a double in units of degrees
 */

    int id=0, im=0, J=0, sign=0;
    double fs=0., dec_rad=0.;
    char id_str[16];

    sscanf(dec_ddmmss, "%s:%d:%lf", id_str, &im, &fs);

    if (id_str[0] == '-') {
        sign = -1;
    }
    else {
        sign = 1;
    }
    sscanf(dec_ddmmss, "%d:%d:%lf", &id, &im, &fs);
    id = id*sign;
    slaDaf2r(id, im, fs, &dec_rad, &J);

    if (J != 0) {
        fprintf(stderr,"Error parsing %s as dd:mm:ss - got %d:%d:%f -- error code %d\n",dec_ddmmss,id,im,fs,J);
        exit(EXIT_FAILURE);
    }

    return dec_rad*DR2D*sign;
}

double parse_ra( char* ra_hhmmss ) {
/* Parse a string containing a right ascension in hh:mm:ss format into
 * a double in units of hours
 */

    int ih=0, im=0, J=0;
    double fs=0., ra_rad=0.;

    sscanf(ra_hhmmss, "%d:%d:%lf", &ih, &im, &fs);

    slaCtf2r(ih, im, fs, &ra_rad, &J);

    if (J != 0) { // slalib returned an error
        fprintf(stderr,"Error parsing %s as hhmmss\nslalib error code: j=%d\n",ra_hhmmss,J);
        fprintf(stderr,"ih = %d, im = %d, fs = %lf\n", ih, im, fs);
        exit(EXIT_FAILURE);
    }

    return ra_rad*DR2H;
}

/*********************************
 convert coords in local topocentric East, North, Height units to
 'local' XYZ units. Local means Z point north, X points through the equator from the geocenter
 along the local meridian and Y is East.
 This is like the absolute system except that zero lon is now
 the local meridian rather than prime meridian.
 Latitude is geodetic, in radian.
 This is what you want for constructing the local antenna positions in a UVFITS antenna table.
 **********************************/

void ENH2XYZ_local(double E,double N, double H, double lat, double *X, double *Y, double *Z) {
    double sl,cl;

    sl = sin(lat);
    cl = cos(lat);
    *X = -N*sl + H*cl;
    *Y = E;
    *Z = N*cl + H*sl;
}

void calcUVW(double ha,double dec,double x,double y,double z,double *u,double *v,double *w) {
    double sh,ch,sd,cd;

    sh = sin(ha); sd = sin(dec);
    ch = cos(ha); cd = cos(dec);
    *u  = sh*x + ch*y;
    *v  = -sd*ch*x + sd*sh*y + cd*z;
    *w  = cd*ch*x  - cd*sh*y + sd*z;
}



int calcEjones(ComplexDouble response[MAX_POLS], // pointer to 4-element (2x2) voltage gain Jones matrix
               const long freq, // observing freq (Hz)
               const float lat, // observing latitude (radians)
               const float az0, // azimuth & zenith angle of tile pointing
               const float za0, // zenith angle to sample
               const float az, // azimuth & zenith angle to sample
               const float za);


void mjd2lst(double mjd, double *lst) {

    // Greenwich Mean Sidereal Time to LMST
    // east longitude in hours at the epoch of the MJD
    double arr_lon_rad = MWA_LON * M_PI/180.0;
    double lmst = slaRanorm(slaGmst(mjd) + arr_lon_rad);

    *lst = lmst;
}

void utc2mjd(char *utc_str, double *intmjd, double *fracmjd) {

    int J=0;
    struct tm *utc;
    utc = calloc(1,sizeof(struct tm));

    sscanf(utc_str, "%d-%d-%dT%d:%d:%d",
            &utc->tm_year, &utc->tm_mon, &utc->tm_mday,
            &utc->tm_hour, &utc->tm_min, &utc->tm_sec);

    slaCaldj(utc->tm_year, utc->tm_mon, utc->tm_mday, intmjd, &J);

    if (J !=0) {
        fprintf(stderr,"Failed to calculate MJD\n");
    }
    *fracmjd = (utc->tm_hour + (utc->tm_min/60.0) + (utc->tm_sec/3600.0))/24.0;
    free(utc);
}

void get_delays(
        // an array of pointings [pointing][ra/dec][characters]
        char                   pointing_array[][2][64],
        int                    npointing, // number of pointings
        long int               frequency,
        struct                 calibration *cal,
        float                  samples_per_sec,
        char                  *beam_model,
        char                  *time_utc,
        double                 sec_offset,
        struct delays          delay_vals[],
        struct metafits_info  *mi,
        ComplexDouble       ****complex_weights_array,  // output: cmplx[npointing][ant][ch][pol]
        ComplexDouble       ****invJi )                 // output: invJi[ant][ch][pol][pol]
{

    int row;     // For counting through nstation*npol rows in the metafits file
    int ant;     // Antenna number
    int pol;     // Polarisation number
    int ch;      // Channel number
    int p1, p2;  // Counters for polarisation

    int conjugate = -1;
    int invert = -1;

    /* easy -- now the positions from the database */

    double phase;

    /* Calibration related defines */
    /* set these here for the library */
    nfrequency = NCHAN;
    nstation   = NANT;
    npol       = NPOL;

    double amp = 0;

    ComplexDouble Jref[NPOL*NPOL];            // Calibration Direction
    ComplexDouble E[NPOL*NPOL];               // Model Jones in Desired Direction
    ComplexDouble G[NPOL*NPOL];               // Coarse channel DI Gain
    ComplexDouble Gf[NPOL*NPOL];              // Fine channel DI Gain
    ComplexDouble Ji[NPOL*NPOL];              // Gain in Desired Direction

    ComplexDouble  **M  = (ComplexDouble ** ) calloc(NANT, sizeof(ComplexDouble * )); // Gain in direction of Calibration
    ComplexDouble ***Jf = (ComplexDouble ***) calloc(NANT, sizeof(ComplexDouble **)); // Fitted bandpass solutions

    for (ant = 0; ant < NANT; ant++) {
        M[ant]  = (ComplexDouble * ) calloc(NPOL*NPOL,  sizeof(ComplexDouble));
        Jf[ant] = (ComplexDouble **) calloc(cal->nchan, sizeof(ComplexDouble *));
        for (ch = 0; ch < cal->nchan; ch++) { // Only need as many channels as used in calibration solution
            Jf[ant][ch] = (ComplexDouble *) calloc(NPOL*NPOL, sizeof(ComplexDouble));
        }
    }

    // Choose a reference tile
    int refinp = 84; // Tile012
    double N_ref = mi->N_array[refinp];
    double E_ref = mi->E_array[refinp];
    double H_ref = mi->H_array[refinp];

    double intmjd;
    double fracmjd;
    double lmst;
    double mjd;
    double pr=0, pd=0, px=0, rv=0, eq=2000, ra_ap=0, dec_ap=0;
    double mean_ra, mean_dec, ha;
    double app_ha_rad, app_dec_rad;
    double az,el;

    double unit_N;
    double unit_E;
    double unit_H;
    int    n;

    int *order = (int *)malloc( NANT*sizeof(int) );

    double dec_degs;
    double ra_hours;
    long int freq_ch;
    int cal_chan;

    double cable;

    double integer_phase;
    double X,Y,Z,u,v,w;
    double geometry, delay_time, delay_samples, cycles_per_sample;

    double Fnorm;
    // Read in the Jones matrices for this (coarse) channel, if requested
    ComplexDouble invJref[4];
    if (cal->cal_type == RTS || cal->cal_type == RTS_BANDPASS) {

        read_rts_file(M, Jref, &amp, cal->filename); // Read in the RTS DIJones file
        inv2x2(Jref, invJref);

        if  (cal->cal_type == RTS_BANDPASS) {

            read_bandpass_file(              // Read in the RTS Bandpass file
                    NULL,                    // Output: measured Jones matrices (Jm[ant][ch][pol,pol])
                    Jf,                      // Output: fitted Jones matrices   (Jf[ant][ch][pol,pol])
                    cal->chan_width,         // Input:  channel width of one column in file (in Hz)
                    cal->nchan,              // Input:  (max) number of channels in one file (=128/(chan_width/10000))
                    NANT,                    // Input:  (max) number of antennas in one file (=128)
                    cal->bandpass_filename); // Input:  name of bandpass file

        }

    }
    else if (cal->cal_type == OFFRINGA) {

        // Find the ordering of antennas in Offringa solutions from metafits file
        for (n = 0; n < NANT; n++) {
            order[mi->antenna_num[n*2]] = n;
        }
        read_offringa_gains_file(M, NANT, cal->offr_chan_num, cal->filename, order);
        free(order);

        // Just make Jref (and invJref) the identity matrix since they are already
        // incorporated into Offringa's calibration solutions.
        Jref[0] = CMaked( 1.0, 0.0 );
        Jref[1] = CMaked( 0.0, 0.0 );
        Jref[2] = CMaked( 0.0, 0.0 );
        Jref[3] = CMaked( 1.0, 0.0 );
        inv2x2(Jref, invJref);
    }

    /* get mjd */
    utc2mjd(time_utc, &intmjd, &fracmjd);

    /* get requested Az/El from command line */

    mjd = intmjd + fracmjd;
    mjd += (sec_offset+0.5)/86400.0;
    mjd2lst(mjd, &lmst);

    for ( int p = 0; p < npointing; p++ )
    {

        dec_degs = parse_dec( pointing_array[p][1] );
        ra_hours = parse_ra( pointing_array[p][0] );

        /* for the look direction <not the tile> */

        mean_ra = ra_hours * DH2R;
        mean_dec = dec_degs * DD2R;

        slaMap(mean_ra, mean_dec, pr, pd, px, rv, eq, mjd, &ra_ap, &dec_ap);

        // Lets go mean to apparent precess from J2000.0 to EPOCH of date.

        ha = slaRanorm(lmst-ra_ap)*DR2H;

        /* now HA/Dec to Az/El */

        app_ha_rad = ha * DH2R;
        app_dec_rad = dec_ap;

        slaDe2h(app_ha_rad, dec_ap, MWA_LAT*DD2R, &az, &el);

        /* now we need the direction cosines */

        unit_N = cos(el) * cos(az);
        unit_E = cos(el) * sin(az);
        unit_H = sin(el);
        // Everything from this point on is frequency-dependent
        for (ch = 0; ch < NCHAN; ch++) {

            // Calculating direction-dependent matrices
            freq_ch = frequency + ch*mi->chan_width;    // The frequency of this fine channel
            cal_chan = 0;
            if (cal->cal_type == RTS_BANDPASS) {
                cal_chan = ch*mi->chan_width / cal->chan_width;  // The corresponding "calibration channel number"
                if (cal_chan >= cal->nchan) {                        // Just check that the channel number is reasonable
                    fprintf(stderr, "Error: \"calibration channel\" %d cannot be ", cal_chan);
                    fprintf(stderr, ">= than total number of channels %d\n", cal->nchan);
                    exit(EXIT_FAILURE);
                }
            }
            if (strcmp(beam_model, "analytic") {
                // Analytic beam:
                calcEjones(E,                                 // pointer to 4-element (2x2) voltage gain Jones matrix
                        freq_ch,                              // observing freq of fine channel (Hz)
                        (MWA_LAT*DD2R),                       // observing latitude (radians)
                        mi->tile_pointing_az*DD2R,            // azimuth & zenith angle of tile pointing
                        (DPIBY2-(mi->tile_pointing_el*DD2R)),
                        az,                                   // azimuth & zenith angle of pencil beam
                        (DPIBY2-el));
            }
            else{
                // FEE2016 beam
                FEEBeam *beam = new_fee_beam( beam_model );
                int zenith_norm = 1; // TODO decide on 1 or 0 for zenith norm
                double *jones = calc_jones(beam, az, DPIBY2-el, freq_ch, mi->delays, mi->amps, zenith_norm);
                for (n = 0; n<NPOL*NPOL; n++){
                    E[n] = CMaked(jones[n*2], jones[n*2+1]);
                }
                free_fee_beam( beam );
                free(jones)
            }

            for (row=0; row < (int)(mi->ninput); row++) {

                // Get the antenna and polarisation number from the row
                ant = row / NPOL;
                pol = row % NPOL;

                mult2x2d(M[ant], invJref, G); // M x J^-1 = G (Forms the "coarse channel" DI gain)

                if (cal->cal_type == RTS_BANDPASS)
                    mult2x2d(G, Jf[ant][cal_chan], Gf); // G x Jf = Gf (Forms the "fine channel" DI gain)
                else
                    cp2x2(G, Gf); //Set the fine channel DI gain equal to the coarse channel DI gain

                mult2x2d(Gf, E, Ji); // the gain in the desired look direction

                // this automatically spots an RTS flagged tile
                if (fabs(CImagd(Ji[0])) < 0.0000001) {  // THIS TEST CAN PROBABLY BE IMPROVED
                    Ji[0] = CMaked( 0.0, 0.0 );
                    Ji[1] = CMaked( 0.0, 0.0 );
                    Ji[2] = CMaked( 0.0, 0.0 );
                    Ji[3] = CMaked( 0.0, 0.0 );
                }

                // Calculate the complex weights array
                if (complex_weights_array != NULL) {
                    if (mi->weights_array[row] != 0.0) {

                        cable = mi->cable_array[row] - mi->cable_array[refinp];
                        double El = mi->E_array[row];
                        double N = mi->N_array[row];
                        double H = mi->H_array[row];

                        ENH2XYZ_local(El,N,H, MWA_LAT*DD2R, &X, &Y, &Z);

                        calcUVW (app_ha_rad,app_dec_rad,X,Y,Z,&u,&v,&w);

                        // shift the origin of ENH to Antenna 0 and hoping the Far Field Assumption still applies ...

                        geometry = (El-E_ref)*unit_E + (N-N_ref)*unit_N + (H-H_ref)*unit_H ;
                        // double geometry = E*unit_E + N*unit_N + H*unit_H ;
                        // Above is just w as you should see from the check.

                        delay_time = (geometry + (invert*(cable)))/(VLIGHT);
                        delay_samples = delay_time * samples_per_sec;

                        // freq should be in cycles per sample and delay in samples
                        // which essentially means that the samples_per_sec cancels

                        // and we only need the fractional part of the turn
                        cycles_per_sample = (double)freq_ch/samples_per_sec;

                        phase = cycles_per_sample*delay_samples;
                        phase = modf(phase, &integer_phase);

                        phase = phase*2*M_PI*conjugate;

                        // Store result for later use
                        complex_weights_array[p][ant][ch][pol] =
                            CScld( CExpd( CMaked( 0.0, phase ) ), mi->weights_array[row] );

                    }
                    else {
                        complex_weights_array[p][ant][ch][pol] = CMaked( mi->weights_array[row], 0.0 ); // i.e. = 0.0
                    }
                }

                // Now, calculate the inverse Jones matrix
                if (invJi != NULL) {
                    if (pol == 0) {
                        conj2x2( Ji, Ji ); // The RTS conjugates the sky so beware
                        Fnorm = norm2x2( Ji, Ji );

                        if (Fnorm != 0.0)
                            inv2x2S( Ji, invJi[ant][ch] );
                        else {
                            for (p1 = 0; p1 < NPOL;  p1++)
                            for (p2 = 0; p2 < NPOL;  p2++)
                                invJi[ant][ch][p1][p2] = CMaked( 0.0, 0.0 );
                        }
                    }
                }

            } // end loop through antenna/pol (row)
        } // end loop through fine channels (ch)

        // Populate a structure with some of the calculated values
        if (delay_vals != NULL) {

            delay_vals[p].mean_ra  = mean_ra;
            delay_vals[p].mean_dec = mean_dec;
            delay_vals[p].az       = az;
            delay_vals[p].el       = el;
            delay_vals[p].lmst     = lmst;
            delay_vals[p].fracmjd  = fracmjd;
            delay_vals[p].intmjd   = intmjd;

        }
    }

    // Free up dynamically allocated memory

    for (ant = 0; ant < NANT; ant++) {
        for (ch = 0; ch < cal->nchan; ch++)
            free(Jf[ant][ch]);
        free(Jf[ant]);
        free(M[ant]);
    }
    free(Jf);
    free(M);
}

int calcEjones(ComplexDouble response[MAX_POLS], // pointer to 4-element (2x2) voltage gain Jones matrix
               const long freq, // observing freq (Hz)
               const float lat, // observing latitude (radians)
               const float az0, // azimuth & zenith angle of tile pointing
               const float za0,
               const float az, // azimuth & zenith angle of pencil beam
               const float za) {

    const double c = VLIGHT;                    // speed of light (m/s)
    float sza, cza, saz, caz;                   // sin & cos of azimuth & zenith angle (pencil beam)
    float sza0, cza0, saz0, caz0;               // sin & cos of azimuth & zenith angle (tile pointing)
    float ground_plane, ha, dec, beam_ha, beam_dec;

    float dipl_e,  dipl_n,  dipl_z;  // Location of dipoles relative to
    float proj_e,  proj_n,  proj_z;  // Components of pencil beam   in local (E,N,Z) coordinates
    float proj0_e, proj0_n, proj0_z; // Components of tile pointing in local (E,N,Z) coordinates

    float rot[2 * N_COPOL];
    ComplexDouble PhaseShift, multiplier;
    int i, j; // For iterating over columns (i; East-West) and rows (j; North-South) of dipoles within tiles
    int n_cols = 4, n_rows = 4; // Each tile is a 4x4 grid of dipoles
    int result = 0;

    float lambda = c / freq;                    //
    float radperm = 2.0 * DPI / lambda;         // Wavenumber (m^-1}
    float dpl_sep = 1.10;                       // Separation between dipoles (m)
    float dpl_hgt = 0.3;                        // Height of dipoles (m)
    float n_dipoles = (float)(n_cols * n_rows); // 4x4 = 16 dipoles per tile

    const int scaling = 1; /* 0 -> no scaling; 1 -> scale to unity toward zenith; 2 -> scale to unity in look-dir */

    response[0] = CMaked( 0.0, 0.0 );
    response[1] = CMaked( 0.0, 0.0 );
    response[2] = CMaked( 0.0, 0.0 );
    response[3] = CMaked( 0.0, 0.0 );

    sza = sin(za);
    cza = cos(za);
    saz = sin(az);
    caz = cos(az);
    sza0 = sin(za0);
    cza0 = cos(za0);
    saz0 = sin(az0);
    caz0 = cos(az0);

    proj_e = sza * saz;
    proj_n = sza * caz;
    proj_z = cza;

    proj0_e = sza0 * saz0;
    proj0_n = sza0 * caz0;
    proj0_z = cza0;

    multiplier = CMaked( 0.0, R2C_SIGN * radperm );

    /* loop over dipoles */
    for (i = 0; i < n_cols; i++) {
        for (j = 0; j < n_rows; j++) {
            dipl_e = (i + 1 - 2.5) * dpl_sep;
            dipl_n = (j + 1 - 2.5) * dpl_sep;
            dipl_z = 0.0;
            PhaseShift = CExpd( CScld( multiplier,
                                  (dipl_e * (proj_e - proj0_e)
                                 + dipl_n * (proj_n - proj0_n)
                                 + dipl_z * (proj_z - proj0_z)) ) );
            // sum for p receptors
            response[0] = CAddd( response[0], PhaseShift );
            response[1] = CAddd( response[1], PhaseShift );
            // sum for q receptors
            response[2] = CAddd( response[2], PhaseShift );
            response[3] = CAddd( response[3], PhaseShift );
        }
    }
    // assuming that the beam centre is normal to the ground plane, the separation should be used instead of the za.
    // ground_plane = 2.0 * sin( 2.0*pi * dpl_hgt/lambda * cos(za) );

    //fprintf(stdout,"calib: dpl_hgt %f radperm %f za %f\n",dpl_hgt,radperm,za);

    ground_plane = 2.0 * sin(dpl_hgt * radperm * cos(za)) / n_dipoles;

    slaH2e(az0, DPIBY2 - za0, lat, &beam_ha, &beam_dec);

    // ground_plane = 2.0*sin(2.0*DPI*dpl_hgt/lambda*cos(slaDsep(beam_ha,beam_dec,ha,dec))) / n_dipoles;

    if (scaling == 1) {
        ground_plane /= (2.0 * sin(dpl_hgt * radperm));
    }
    else if (scaling == 2) {
        ground_plane /= (2.0 * sin( 2.0*DPI* dpl_hgt/lambda * cos(za) ));
    }

    slaH2e(az, DPIBY2 - za, lat, &ha, &dec);

    rot[0] = cos(lat) * cos(dec) + sin(lat) * sin(dec) * cos(ha);
    rot[1] = -sin(lat) * sin(ha);
    rot[2] = sin(dec) * sin(ha);
    rot[3] = cos(ha);

    //fprintf(stdout,"calib:HA is %f hours \n",ha*DR2H);
    // rot is the Jones matrix, response just contains the phases, so this should be an element-wise multiplication.
    for (i = 0; i < 4; i++)
        response[i] = CScld( response[i], rot[i] * ground_plane );
    //fprintf(stdout,"calib:HA is %f groundplane factor is %f\n",ha*DR2H,ground_plane);
    return (result);

} /* calcEjones */

