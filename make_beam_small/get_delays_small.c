#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include "mwac_utils.h"
#include "slalib.h"
#include "slamac.h"
#include "psrfits.h"
#include "fitsio.h"
#include <string.h>
#include "beamer_version.h"
#include "make_beam_small.h"

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

void printf_psrfits( struct psrfits *pf ) {
    fprintf(stdout, "\nPSRFITS:\n");
    fprintf(stdout, "Basename of output file     [%s]\n", pf->basefilename);
    fprintf(stdout, "Filename of output file     [%s]\n", pf->filename);
    fprintf(stdout, "CFITSIO file pointer        [%p]\n", pf->fptr);

    fprintf(stdout, "\nPSRFITS HDRINFO:\n");
    fprintf(stdout, "Obs mode                    [%s]\n", pf->hdr.obs_mode);
    fprintf(stdout, "Telescope                   [%s]\n", pf->hdr.telescope);
    fprintf(stdout, "Observer                    [%s]\n", pf->hdr.observer);
    fprintf(stdout, "Source                      [%s]\n", pf->hdr.source);
    fprintf(stdout, "Front End                   [%s]\n", pf->hdr.frontend);
    fprintf(stdout, "Back End                    [%s]\n", pf->hdr.backend);
    fprintf(stdout, "Project ID                  [%s]\n", pf->hdr.project_id);
    fprintf(stdout, "Date Obs                    [%s]\n", pf->hdr.date_obs);
    fprintf(stdout, "RA (string)                 [%s]\n", pf->hdr.ra_str);
    fprintf(stdout, "Dec (string)                [%s]\n", pf->hdr.dec_str);
    fprintf(stdout, "Pol recorded (LIN or CIRC)  [%s]\n", pf->hdr.poln_type);
    fprintf(stdout, "Order of pols               [%s]\n", pf->hdr.poln_order);
    fprintf(stdout, "Track mode                  [%s]\n", pf->hdr.track_mode);
    fprintf(stdout, "Cal mode                    [%s]\n", pf->hdr.cal_mode);
    fprintf(stdout, "Feed mode                   [%s]\n", pf->hdr.feed_mode);
    fprintf(stdout, "Start MJD                   [%Lf]\n", pf->hdr.MJD_epoch);
    fprintf(stdout, "Sample Time (s)             [%lf]\n", pf->hdr.dt);
    fprintf(stdout, "Centre Frequency (MHz)      [%lf]\n", pf->hdr.fctr);
    fprintf(stdout, "Orig freq spacing (MHz)     [%lf]\n", pf->hdr.orig_df);
    fprintf(stdout, "Freq spacing (MHz)          [%lf]\n", pf->hdr.df);
    fprintf(stdout, "Bandwidth (MHz)             [%lf]\n", pf->hdr.BW);
    fprintf(stdout, "RA (2000) (deg)             [%lf]\n", pf->hdr.ra2000);
    fprintf(stdout, "Dec (2000) (deg)            [%lf]\n", pf->hdr.dec2000);
    fprintf(stdout, "Azimuth (deg)               [%lf]\n", pf->hdr.azimuth);
    fprintf(stdout, "Zenith Angle (deg)          [%lf]\n", pf->hdr.zenith_ang);
    fprintf(stdout, "Beam FWHM (deg)             [%lf]\n", pf->hdr.beam_FWHM);

    fprintf(stdout, "Length of scan in this file [%lf]\n", pf->hdr.scanlen);
    fprintf(stdout, "Seconds past 00h LST        [%lf]\n", pf->hdr.start_lst);
    fprintf(stdout, "Seconds past 00h UTC        [%lf]\n", pf->hdr.start_sec);

    fprintf(stdout, "Start MJD (whole day)       [%d]\n", pf->hdr.start_day);
    fprintf(stdout, "Scan Number                 [%d]\n", pf->hdr.scan_number);
    fprintf(stdout, "Number of bits per sample   [%d]\n", pf->hdr.nbits);

    fprintf(stdout, "Number of Channels          [%d]\n", pf->hdr.nchan);
    fprintf(stdout, "Number of polarisations     [%d]\n", pf->hdr.npol);
    fprintf(stdout, "Number of spectra per row   [%d]\n", pf->hdr.nsblk);

    fprintf(stdout, "Summed Polarisations? [1/0] [%d]\n", pf->hdr.summed_polns);
    fprintf(stdout, "Receiver Polarisation       [%d]\n", pf->hdr.rcvr_polns);
    fprintf(stdout, "Offset Subint               [%d]\n", pf->hdr.offset_subint);
    fprintf(stdout, "Dwnsmpl fact in time        [%d]\n", pf->hdr.ds_time_fact);
    fprintf(stdout, "Dwnsmpl fact in freq        [%d]\n", pf->hdr.ds_freq_fact);
    fprintf(stdout, "Only Stokes I?              [%d]\n", pf->hdr.onlyI);

    fprintf(stdout, "\nPSRFITS SUBINT:\n");
    fprintf(stdout, "Length of subint (sec)      [%lf]\n", pf->sub.tsubint);
    fprintf(stdout, "Offset (sec)                [%lf]\n", pf->sub.offs);
    fprintf(stdout, "LST (sec)                   [%lf]\n", pf->sub.lst);
    fprintf(stdout, "RA (J2000) (deg)            [%lf]\n", pf->sub.ra);
    fprintf(stdout, "Dec (J2000) (deg)           [%lf]\n", pf->sub.dec);
    fprintf(stdout, "Gal. long. (deg)            [%lf]\n", pf->sub.glon);
    fprintf(stdout, "Gal. lat. (deg)             [%lf]\n", pf->sub.glat);
    fprintf(stdout, "Feed angle (deg)            [%lf]\n", pf->sub.feed_ang);
    fprintf(stdout, "Pos angle of feed (deg)     [%lf]\n", pf->sub.pos_ang);
    fprintf(stdout, "Parallactic angle           [%lf]\n", pf->sub.par_ang);
    fprintf(stdout, "Telescope azimuth           [%lf]\n", pf->sub.tel_az);
    fprintf(stdout, "Telescope zenith angle      [%lf]\n", pf->sub.tel_zen);
    fprintf(stdout, "Bytes per row of raw data   [%d]\n", pf->sub.bytes_per_subint);
    fprintf(stdout, "FITS data typecode          [%d]\n", pf->sub.FITS_typecode);

}

double parse_dec( char* dec_ddmmss ) {
/* Parse a string containing a declination in dd:mm:ss format into
 * a double in units of degrees
 */

    int id=0, im=0, J=0, sign=0;
    double fs=0., dec_rad=0.;
    char id_str[4];

    sscanf(dec_ddmmss,"%s:%d:%lf",id_str,&im,&fs);

    if (id_str[0] == '-') {
        sign = -1;
    }
    else {
        sign = 1;
    }
    sscanf(dec_ddmmss,"%d:%d:%lf",&id,&im,&fs);
    id = id*sign;
    slaDaf2r(id,im,fs,&dec_rad,&J);

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

    sscanf(ra_hhmmss,"%d:%d:%lf",&ih,&im,&fs);

    slaCtf2r(ih,im,fs,&ra_rad,&J);

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



int calcEjones(complex double response[MAX_POLS], // pointer to 4-element (2x2) voltage gain Jones matrix
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
        int coarse_chan,
        char *dec_ddmmss,
        char *ra_hhmmss,
        long int frequency,
        char *metafits,
        int get_offringa,
        int get_rts,
        char *DI_Jones_file,
        float samples_per_sec,
        long int chan_width,
        char *time_utc,
        double sec_offset,
        struct delays *delay_vals,
        complex double **complex_weights_array,  // output
        double *weights_array,
        complex double **invJi                   // output
        )
{
    
    fprintf(stdout, "* RUNNING GET_DELAYS *\n");
    int row;
    int i, j;

    // some defaults for testing
    
    double dec_degs = parse_dec( dec_ddmmss );
    double ra_hours = parse_ra( ra_hhmmss );

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

    // Jref is the reference Jones direction - for which the calibration was generated
    // Mref is measured Jones in Jref direction
    // G = DI gain formed by inv(Jref).Mref
    // E = model is desired direction
    // Ji = Jones in desired direction (formed by G.E)

    complex double *Jref = (complex double *) calloc(NPOL*NPOL, sizeof(complex double));   // Calibration Direction
    complex double *E    = (complex double *) calloc(NPOL*NPOL, sizeof(complex double));   // Model Jones in Desired Direction
    complex double **G   = (complex double **) calloc(NANT, sizeof(complex double *)); // Actual DI Gain
    complex double **M   = (complex double **) calloc(NANT, sizeof(complex double *)); // Gain in direction of Calibration
    complex double **Ji  = (complex double **) calloc(NANT, sizeof(complex double *)); // Gain in Desired Direction

    for (i = 0; i < NANT; i++) { //
        G[i] = (complex double *) malloc(NPOL * NPOL * sizeof(complex double)); //
        M[i] = (complex double *) malloc(NPOL * NPOL * sizeof(complex double)); //
        Ji[i] =(complex double *) malloc(NPOL * NPOL * sizeof(complex double)); //
        if (G[i] == NULL || M[i] == NULL || Ji[i] == NULL) { //
            fprintf(stderr, "malloc failed: G[i], M[i], J[i]\n"); //
            exit(1); //
        } //
    } //

    // ===================================================================================== //
    // Get actual tile pointing Az El from metafits file

    fitsfile *fptr=NULL;
    int status = 0;
    double tile_pointing_ra = 0.0;
    double tile_pointing_dec = 0.0;
    double tile_pointing_az;
    double tile_pointing_el;

    fits_open_file(&fptr,metafits,READONLY,&status);
    if (fptr == NULL) {
        fprintf( stderr, "Failed to open metafits file \"%s\"\n", metafits );
        exit(EXIT_FAILURE);
    }
    fits_read_key(fptr,TDOUBLE,"RA",&tile_pointing_ra,NULL,&status);
    fits_read_key(fptr,TDOUBLE,"DEC",&tile_pointing_dec,NULL,&status);
    fits_read_key(fptr,TDOUBLE,"AZIMUTH",&tile_pointing_az,NULL,&status);
    fits_read_key(fptr,TDOUBLE,"ALTITUDE",&tile_pointing_el,NULL,&status);


    fits_close_file(fptr,&status);

     if (status != 0) {
        fprintf(stderr,"Fits status set: failed to read az/alt, ra/dec fitskeysfrom the header\n");
        exit(1);
    }
    // ====================================================================================== //
    // Get tile coordinate information from metafits file

    /* open the metafits file */
    status = 0;
    int anynull = 0;
    fits_open_file(&fptr,metafits,READONLY,&status);
    
    fits_movnam_hdu(fptr, BINARY_TBL, "TILEDATA", 0, &status);

    if (status != 0) {
        fprintf(stderr,"Error:Failed to move to TILEDATA HDU\n");
        exit(-1);
    }

    size_t ninput = 0;
    fits_read_key(fptr,TINT,"NAXIS2",&ninput,NULL,&status);

    //fprintf(stdout,"Status: will read %zu inputs\n",ninput);

    if (status != 0){
        fprintf(stderr,"Error:Failed to read size of binary table in TILEDATA\n");
        exit(-1);
    }

    /* allocate arrays for tile positions */
    float *N_array = (float *) malloc(ninput*sizeof(float));
    float *E_array = (float *) malloc(ninput*sizeof(float));
    float *H_array = (float *) malloc(ninput*sizeof(float));
    char **tilenames = (char **) malloc(ninput*sizeof(char *));
    for (i = 0; i < (int)ninput; i++) {
        tilenames[i] = (char *) malloc(32*sizeof(char));
    }
    float *cable_array = (float *) malloc(ninput*sizeof(float));
    char *testval = (char *) malloc(1024);
    short int *antenna_num = (short int *)malloc(ninput*sizeof(short int));
    int colnum;
    
    /* read the columns */
    for (i=0; i < (int)ninput; i++) {

        fits_get_colnum(fptr, 1, "Length", &colnum, &status);
        if(fits_read_col_str(fptr,colnum,i+1,1,1,"0.0",&testval,&anynull,&status)) {

            fprintf(stderr,"Error:Failed to cable column  in metafile\n");
            exit(-1);
        }

        sscanf(testval,"EL_%f",&cable_array[i]);
        //fprintf(stdout,"Input %d Cable %f\n",i,cable_array[i]);
    }

    fits_get_colnum(fptr, 1, "TileName", &colnum, &status);
    if (status != 0) {
        status = 0;
        fits_get_colnum(fptr, 1, "Tile", &colnum, &status);
    }
    if (status != 0) {
        fprintf(stderr, "Could not find either column \"TileName\" or \"Tile\" in metafits file\n");
        exit(-1);
    }
    fits_read_col(fptr,TSTRING,colnum,1,1,ninput,NULL,tilenames,&anynull,&status);
    if (status != 0){
        fprintf(stderr,"Error:Failed to read Tile(Name) in metafile\n");
        exit(-1);
    }

    fits_get_colnum(fptr, 1, "North", &colnum, &status);
    fits_read_col_flt(fptr,colnum,1,1,ninput,0.0,N_array,&anynull,&status);
    if (status != 0){
        fprintf(stderr,"Error:Failed to read  N coord in metafile\n");
        exit(-1);
    }

    fits_get_colnum(fptr, 1, "East", &colnum, &status);
    fits_read_col_flt(fptr,colnum,1,1,ninput,0.0,E_array,&anynull,&status);
    if (status != 0){
        fprintf(stderr,"Error:Failed to read E coord in metafile\n");
        exit(-1);
    }

    fits_get_colnum(fptr, 1, "Height", &colnum, &status);
    fits_read_col_flt(fptr,colnum,1,1,ninput,0.0,H_array,&anynull,&status);

    if (status != 0){
        fprintf(stderr,"Error:Failed to read H coord in metafile\n");
        exit(-1);
    }

    fits_get_colnum(fptr, 1, "Antenna", &colnum, &status);
    fits_read_col_sht(fptr,colnum,1,1,ninput,0.0,antenna_num,&anynull,&status);

    if (status != 0){
        fprintf(stderr,"Error:Failed to read field \"Antenna\" in metafile\n");
        exit(-1);
    }

    fits_close_file(fptr,&status);
    int refinp = 84; // Tile012
    double E_ref = E_array[refinp];
    double N_ref = N_array[refinp];
    double H_ref = H_array[refinp];

    // END (get tile coordinate information from metafits file)
    // ====================================================================================== //

    
    double intmjd;
    double fracmjd;
    double lmst;
    double mjd;
    double pr=0,pd=0,px=0,rv=0,eq=2000,ra_ap=0,dec_ap=0;
    double mean_ra,mean_dec,ha;
    double app_ha_rad,app_dec_rad;
    double az,el;
    
    double unit_N;
    double unit_E;
    double unit_H;
    int    n;

    // Read in the Jones matrices for this (coarse) channel, if requested
    complex double invJref[4];
    if (get_rts) {
        read_rts_file(M, Jref, NANT, &amp, DI_Jones_file);
        inv2x2(Jref, invJref);
    }
    else if (get_offringa) {
        // Find the ordering of antennas in Offringa solutions from metafits file
        int *order = (int *)malloc(NANT*sizeof(int));
        for (n = 0; n < NANT; n++) {
            order[antenna_num[n*2]] = n;
        }
        read_offringa_gains_file(M, NANT, coarse_chan, DI_Jones_file, order);
        //read_offringa_gains_file(M, NANT, coarse_chan, DI_Jones_file, NULL);
        free(order);
        // Just make Jref (and invJref) the identity matrix since they are already
        // incorporated into Offringa's calibration solutions.
        Jref[0] = 1 + I*0;
        Jref[1] = 0 + I*0;
        Jref[2] = 0 + I*0;
        Jref[3] = 1 + I*0;
        inv2x2(Jref, invJref);
    }

    /* get mjd */

    utc2mjd(time_utc, &intmjd, &fracmjd);

    /* get requested Az/El from command line */
  
    mjd = intmjd + fracmjd;
    mjd += (sec_offset+0.5)/86400.0;
    mjd2lst(mjd, &lmst);

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
    
    if (get_rts || get_offringa) {
        //fprintf(stdout, "Calculating direction-dependent matrices\n");
        double Fnorm;
        calcEjones(E, // pointer to 4-element (2x2) voltage gain Jones matrix
                   frequency, // observing freq (Hz)
                   (MWA_LAT*DD2R), // observing latitude (radians)
                   tile_pointing_az*DD2R, // azimuth & zenith angle of tile pointing
                   (DPIBY2-(tile_pointing_el*DD2R)), // zenith angle to sample
                   az, // azimuth & zenith angle to sample
                   (DPIBY2-el));
        for (i=0; i < 4;i++) {
            //fprintf(stdout,"calib:Jones Jref[%d] %f %f: Delay Jref[%d] %f %f\n",i,creal(Jref[i]),cimag(Jref[i]),i,creal(E[i]),cimag(E[i]));
            //fprintf(stdout,"calib:ratio RTS/Delay [%d]  %f %f \n",i,creal(Jref[i])/creal(E[i]),cimag(Jref[i])/cimag(E[i]));
        }
        for (i=0;i<NANT;i++){
            mult2x2d(M[i],invJref,G[i]); // forms the DI gain
            mult2x2d(G[i],E,Ji[i]); // the gain in the desired look direction
            
            // this automatically spots an RTS flagged tile
            Fnorm = 0;
            for (j=0; j < 4;j++) {
                Fnorm += (double) Ji[j][i] * conj(Ji[j][i]);
            }
            if (fabs(cimag(Ji[i][0])) < 0.0000001) {
                Ji[i][0] = 0.0 + I*0;
                Ji[i][1] = 0.0 + I*0;
                Ji[i][2] = 0.0 + I*0;
                Ji[i][3] = 0.0 + I*0;

                G[i][0] = 0.0 + I*0;
                G[i][1] = 0.0 + I*0;
                G[i][2] = 0.0 + I*0;
                G[i][3] = 0.0 + I*0;

            }

        }
    }
        
    
    /* for the tile <not the look direction> */

    int ch=0;

    for (row=0; row < (int)ninput; row++) {

        if (weights_array[row] != 0.0) {

            double cable = cable_array[row]-cable_array[refinp];
            double E = E_array[row];
            double N = N_array[row];
            double H = H_array[row];

            double integer_phase;
            double X,Y,Z,u,v,w;

            ENH2XYZ_local(E,N,H, MWA_LAT*DD2R, &X, &Y, &Z);

            calcUVW (app_ha_rad,app_dec_rad,X,Y,Z,&u,&v,&w);

            // shift the origin of ENH to Antenna 0 and hoping the Far Field Assumption still applies ...

            double geometry = (E-E_ref)*unit_E + (N-N_ref)*unit_N + (H-H_ref)*unit_H ;
            // double geometry = E*unit_E + N*unit_N + H*unit_H ;
            // Above is just w as you should see from the check.

            double delay_time = (geometry + (invert*(cable)))/(VLIGHT);
            double delay_samples = delay_time * samples_per_sec;

            for (ch=0; ch < NCHAN;ch++) {
                long int freq_ch = frequency + ch*chan_width;

                // freq should be in cycles per sample and delay in samples
                // which essentially means that the samples_per_sec cancels

                // and we only need the fractional part of the turn
                double cycles_per_sample = (double)freq_ch/samples_per_sec;

                phase = cycles_per_sample*delay_samples;
                phase = modf(phase, &integer_phase);

                phase = phase*2*M_PI*conjugate;

                // Store result for later use
                complex_weights_array[row][ch] = weights_array[row]*cexp(I*phase);

            }
        }
        else {
            for (ch=0; ch < NCHAN; ch++)
                complex_weights_array[row][ch] = weights_array[row]; // i.e. =0.0
        }

        // Now, calculate the inverse Jones matrix
        if (row % NPOL == 0) {

            int station = row / NPOL;
            double Fnorm;

            conj2x2( Ji[station], Ji[station] ); // The RTS conjugates the sky so beware
            Fnorm = norm2x2( Ji[station], Ji[station] );

            if (Fnorm != 0.0)
                inv2x2( Ji[station], invJi[station] );
            else
                for (i = 0; i < 4; i++)
                    invJi[station][i] = 0.0 + I*0.0;
        }
        
        
    }

    // Populate a structure with some of the calculated values
    if (delay_vals != NULL) {
        
        delay_vals->mean_ra  = mean_ra;
        delay_vals->mean_dec = mean_dec;
        delay_vals->az       = az;
        delay_vals->el       = el;
        delay_vals->lmst     = lmst;
        delay_vals->fracmjd  = fracmjd;
        delay_vals->intmjd   = intmjd;

    }
    

    // Free up memory

    for (i = 0; i < NANT; i++) { //
        free(G[i]);
        free(M[i]);
        free(Ji[i]);
    } //
    free(Jref);
    free(E);
    free(G);
    free(M);
    free(Ji);

    free(N_array);
    free(E_array);
    free(H_array);
    for (i = 0; i < (int)ninput; i++)
        free(tilenames[i]);
    free(tilenames);
    free(cable_array);
    free(testval);
    free(antenna_num);

    fprintf(stdout, "* EXITING GET_DELAYS *\n");
}

int calcEjones(complex double response[MAX_POLS], // pointer to 4-element (2x2) voltage gain Jones matrix
               const long freq, // observing freq (Hz)
               const float lat, // observing latitude (radians)
               const float az0, // azimuth & zenith angle of tile pointing
               const float za0, // zenith angle to sample
               const float az, // azimuth & zenith angle to sample
               const float za) {
    
    const double c = VEL_LIGHT;
    float sza, cza, saz, caz, sza0, cza0, saz0, caz0;
    float ground_plane, ha, dec, beam_ha, beam_dec;
    float dipl_e, dipl_n, dipl_z, proj_e, proj_n, proj_z, proj0_e, proj0_n,
    proj0_z;
    float rot[2 * N_COPOL];
    complex double PhaseShift, multiplier;
    int i, j, n_cols = 4, n_rows = 4, result = 0;
    
    float lambda = c / freq;
    float radperm = 2.0 * DPI / lambda;
    float dpl_sep = 1.10, dpl_hgt = 0.3, n_dipoles = (float) n_cols * n_rows;
    
    const int scaling = 1; /* 0 -> no scaling; 1 -> scale to unity toward zenith; 2 -> scale to unity in look-dir */
    
    response[0] = 0.0;
    response[1] = 0.0;
    response[2] = 0.0;
    response[3] = 0.0;
    
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
    
    n_cols = 4;
    n_rows = 4;
    
    multiplier = R2C_SIGN * I * radperm;
    
    /* loop over dipoles */
    for (i = 0; i < n_cols; i++) {
        for (j = 0; j < n_rows; j++) {
            dipl_e = (i + 1 - 2.5) * dpl_sep;
            dipl_n = (j + 1 - 2.5) * dpl_sep;
            dipl_z = 0.0;
            PhaseShift = cexp(
                              multiplier
                              * (dipl_e * (proj_e - proj0_e)
                                 + dipl_n * (proj_n - proj0_n)
                                 + dipl_z * (proj_z - proj0_z)));
            // sum for p receptors
            response[0] += PhaseShift;
            response[1] += PhaseShift;
            // sum for q receptors
            response[2] += PhaseShift;
            response[3] += PhaseShift;
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
    response[0] *= rot[0] * ground_plane;
    response[1] *= rot[1] * ground_plane;
    response[2] *= rot[2] * ground_plane;
    response[3] *= rot[3] * ground_plane;
    //fprintf(stdout,"calib:HA is %f groundplane factor is %f\n",ha*DR2H,ground_plane);
    return (result);
    
} /* calcEjones */

