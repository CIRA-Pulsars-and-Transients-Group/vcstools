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
#define NANT    128
#define VLIGHT 299792458.0        // speed of light. m/s
double arr_lat_rad=MWA_LAT*(M_PI/180.0),arr_lon_rad=MWA_LON*(M_PI/180.0),height=MWA_HGT;

int verbose;

/* these externals are needed for the mwac_utils library */
int nfrequency;
int npol;
int nstation;
//=====================//

void printf_psrfits( struct psrfits *pf ) {
    fprintf(stdout, "PSRFITS:\n");
    fprintf(stdout, "Basename of output file     [%s]\n", pf->basefilename);
    fprintf(stdout, "Filename of output file     [%s]\n", pf->filename);

    fprintf(stdout, "PSRFITS HDRINFO:\n");
    fprintf(stdout, "Obs mode                    [%s]\n", pf->hdr.obs_mode);
    fprintf(stdout, "Telescope                   [%s]\n", pf->hdr.telescope);
    fprintf(stdout, "Observer                    [%s]\n", pf->hdr.observer);
    fprintf(stdout, "Source                      [%s]\n", pf->hdr.source);
    fprintf(stdout, "Front End                   [%s]\n", pf->hdr.frontend);
    fprintf(stdout, "Back End                    [%s]\n", pf->hdr.backend);
    fprintf(stdout, "project_id                  [%s]\n", pf->hdr.project_id);
    fprintf(stdout, "Date Obs                    [%s]\n", pf->hdr.date_obs);
    fprintf(stdout, "Sample Time (s)             [%lf]\n", pf->hdr.dt);
    fprintf(stdout, "Centre Frequency (MHz)      [%lf]\n", pf->hdr.fctr);
    fprintf(stdout, "Bandwidth (MHz)             [%lf]\n", pf->hdr.BW);
    fprintf(stdout, "RA (2000) (deg)             [%lf]\n", pf->hdr.ra2000);
    fprintf(stdout, "Dec (2000) (deg)            [%lf]\n", pf->hdr.dec2000);
    fprintf(stdout, "Azimuth (deg)               [%lf]\n", pf->hdr.azimuth);
    fprintf(stdout, "Zenith Angle (deg)          [%lf]\n", pf->hdr.zenith_ang);
    fprintf(stdout, "Beam FWHM (deg)             [%lf]\n", pf->hdr.beam_FWHM);
    fprintf(stdout, "Seconds past 00h LST        [%lf]\n", pf->hdr.start_lst);
    fprintf(stdout, "Seconds past 00h UTC        [%lf]\n", pf->hdr.start_sec);
    fprintf(stdout, "Start MJD (whole day)       [%d]\n", pf->hdr.start_day);
    fprintf(stdout, "Scan Number                 [%d]\n", pf->hdr.scan_number);
    fprintf(stdout, "Receiver Polarisation       [%d]\n", pf->hdr.rcvr_polns);
    fprintf(stdout, "Summed Polarisations? [1/0] [%d]\n", pf->hdr.summed_polns);
    fprintf(stdout, "Number of polarisations     [%d]\n", pf->hdr.npol);
    fprintf(stdout, "Length of scan in this file [%lf]\n", pf->hdr.scanlen);
    fprintf(stdout, "Offset Subint               [%d]\n", pf->hdr.offset_subint);
    fprintf(stdout, "Number of Channels          [%d]\n", pf->hdr.nchan);
    fprintf(stdout, "Number of bits per sample   [%d]\n", pf->hdr.nbits);
    fprintf(stdout, "Number of spectra per row   [%d]\n", pf->hdr.nsblk);
    fprintf(stdout, "Start MJD                   [%Lf]\n", pf->hdr.MJD_epoch);
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

    // we need the eq of the equinoxes - as this might shift our LST by +/- a second

    // double gast = slaEqeqx(mjd) + gmst; // GAST in radians
    
    // lst (at epoch of date) is therefore
    
    /// double last = gast + (lon_hours*DH2R); // last in radians
    
    // last = last * DR2H;
    
    // if (last > 24.0) {
       // last = last - 24.0;
    //}
    //else if (last < 0.0) {
      //  last = 24.0 + last;
    //}


    *lst = lmst;
}

void utc2mjd(char *utc_str, double *intmjd, double *fracmjd) {
    
    extern int verbose;
    int J=0;
    struct tm *utc;
    utc = calloc(1,sizeof(struct tm));
    
    //fprintf(stderr,"Parsing UTC from %s\n",opts->utc_str);
    
    sscanf(utc_str,"%d-%d-%dT%d:%d:%d",&utc->tm_year,&utc->tm_mon,&utc->tm_mday,&utc->tm_hour,&utc->tm_min,&utc->tm_sec);
    if (verbose)
        fprintf(stderr,"yr %d, mon %d, day %d, hour %d, min %d, sec %d\n",utc->tm_year, utc->tm_mon, utc->tm_mday, utc->tm_hour, utc->tm_min,utc->tm_sec);
    
    slaCaldj(utc->tm_year,utc->tm_mon,utc->tm_mday,intmjd,&J);
    
    if (J !=0) {
        fprintf(stderr,"Failed to calculate MJD\n");
    }
    *fracmjd = (utc->tm_hour + (utc->tm_min/60.0) + (utc->tm_sec/3600.0))/24.0;
    free(utc);
}


/*
void usage() {
    
    fprintf(stdout,"get_delays - interrogates the database for telescope positions and cable lengths and returns the geometric delay and delay rate in seconds (and seconds per second");
    fprintf(stdout,"syntax:\n\
            get_delays -z <utc time string> -o obsid -r <ra in hh:mm:ss> -d <dec in dd:mm:ss>\n \
            \nOther options include:\n \
            \t -a <output file directory> directory to put the jones.txt & flags.txt & phases.txt\n \
            \t -C <coarse channel number> used in conjunction with -O. Default = '0' \n \
            \t -e number of low channels to skip\n \
            \t -f <middle of the first frequency channel in Hz> \n \
            \t -G Switch off Geometry [Expert] \n \
            \t -m <metafits file> for this obsID \n \
            \t -n number of channels \n \
            \t -o <obsID> \n \
            \t -O <DI Jones file from Offringa's tools> Jones matrix input \n \
            \t -p create a psrfits header for this obs\n \
            \t -R <DI Jones file from the RTS> Jones matrix input \n \
            \t -s samples per second \n \
            \t -t <input number in correlator product order>\n \
            \t -v <1 == verbose> \n\
            \t -V print the version number and exit \n\
            \t -w channel width (in Hz) \n \
            ");
    
    fprintf(stdout,"Some notes: What do we need to know \n \
            \n \
            1) The antenna positions\n \
            2) The cable lengths\n \
            3) The source position\n \
            4) The time\n \
            \n \
            We will obtain 1 & 2 via a database call. This needs a time in order to work.\n \
            \n \
            Therefore we first require a time. We can separate this into a sidereal time (or one where we can get to LST quickly) \
            and a GPS time for which we can get antenna information.\n \
            \n \
            So arguments are at the very least:\n \
            \n \
            reftime=\"gps time/obsid for reference and antenna ocations + flags\"\n \
            calctime=\"UTC time for calculation - good to at least the second\" in ISO 8601 e.g. 2014-03-18T19:51:54\n \
            source position =\"ra/dec of the look direction\" \n \
            \n I've decided that we should output phase of the complex weight. \
            This will need to be combined with an amplitude weighting - but this requires get_weights ....\n");
    
    fprintf(stdout,"Some notes on the Jones matrix information\n \
            \n \
            Essentially we load in the Direction Indepenendent Calibration files generated by the RTS and optionally the bandpass solutions. \
            These are used to form the intrumental Jones matrix for the look direction specified on the command line. \
            The outer product of the Jones matrices with themselves (J x J*) is then formed and inverted. \
            This are output in a file /tmp/jones.txt\n\n");
    ;
    
    
}
*/

void get_delays(
        int coarse_chan,
        char *dec_ddmmss,
        char *ra_hhmmss,
        int edge,
        long int frequency,
        char *metafits,
        int nchan,
        char *obsid,
        int get_offringa,
        int get_rts,
        char *DI_Jones_file,
        int get_psrfits,
        int write_psrfits,
        float samples_per_sec,
        int verbose,
        long int chan_width,
        char *time_utc,
        double sec_offset,
        char *phases_filename,            // soon-deprecated output
        char *flag_filename,              // soon-deprecated output
        char *jones_filename,             // soon-deprecated output
        char *psrfits_filename,           // soon-deprecated output
        double **weights_array,           // output
        double ***phases_array,           // output
        double ***complex_weights_array,  // output
        complex double ***invJi           // output
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
    int write_files = 1;
    
    FILE *phase_file = NULL;
    FILE *flag_file = NULL;
    FILE *jones_file = NULL;
    FILE *psrfits_file = NULL;
    int secs = 0;

    /*
    int c;
    if (argc > 1) {
        
        while ((c = getopt(argc, argv, "b:C:d:e:f:hm:n:o:O:pr:R:s:vVw:z:")) != -1) {
            switch(c) {
                case 'b':
                    nsecs = atoi(optarg);
                    break;
                case 'C':
                    coarse_chan = atoi(optarg);
                    break;
                case 'd':
                {
                    int id=0,im=0,J=0,sign=0;
                    double fs=0.,dec_rad=0.;
                    char id_str[4];
                    char *dec_ddmmss = strdup(optarg);
                    
                    sscanf(dec_ddmmss,"%s:%d:%lf",id_str,&im,&fs);
                    
                    if (id_str[0] == '-') {
                        sign = -1;
                    }
                    else {
                        sign = 1;
                    }
                    sscanf(dec_ddmmss,"%d:%d:%lf",&id,&im,&fs);
                    id=id*sign;
                    slaDaf2r(id,im,fs,&dec_rad,&J);
                    
                    if (J==0) {
                        dec_degs = dec_rad*DR2D*sign;
                    }
                    else {
                        fprintf(stderr,"Error parsing %s as dd:mm:ss - got %d:%d:%f -- error code %d\n",dec_ddmmss,id,im,fs,J);
                        usage();
                        exit(-1);
                    }
                    break;
                }
                case 'e':
                    edge = atoi(optarg);
                    break;
                case 'f':
                    frequency = atol(optarg);
                    break;
                case 'h':
                    usage();
                    exit(0);
                    break;
                case 'm':
                    metafits = strdup(optarg);
                    break;
                case 'n':
                    nchan = atoi(optarg);
                    break;
                case 'o':
                    obsid = strdup(optarg);
                    break;
                case 'O':
                    get_offringa = 1;
                    get_rts = 0;
                    DI_Jones_file = strdup(optarg);
                    break;
                case 'p':
                    get_psrfits=1;
                    write_psrfits=1;
                    break;
                case 'R':
                    get_rts = 1;
                    get_offringa = 0;
                    DI_Jones_file = strdup(optarg);
                    break;
                case 'r':
                {
                    int ih=0,im=0,J=0;
                    double fs=0.,ra_rad=0.;
                    
                    char *ra_hhmmss = strdup(optarg);
                    
                    sscanf(ra_hhmmss,"%d:%d:%lf",&ih,&im,&fs);
                    
                    slaCtf2r(ih,im,fs,&ra_rad,&J);
                    
                    if (J==0) {
                        ra_hours = ra_rad*DR2H;
                    }
                    else {
                        fprintf(stderr,"Error parsing %s as hhmmss\nslalib error code: j=%d\n",ra_hhmmss,J);
                        fprintf(stderr,"ih = %d, im = %d, fs = %lf\n", ih, im, fs);
                        usage();
                        exit(-1);
                    }
                    break;
                }
                case 's':
                    samples_per_sec = atof(optarg);
                    break;
                case 'v':
                    verbose = 1;
                    break;
                case 'V':
                    printf("%s\n", GET_DELAYS_VERSION);
                    exit(0);
                    break;
                case 'w':
                    chan_width = atol(optarg);
                    break;
                case 'z':
                    time_utc = strdup(optarg);
                    break;
                default:
                    usage();
                    exit(1);
            }
        }
    }
    
    if (argc == 1) {
        usage();
        exit(1);
    }
    */

    /* easy -- now the positions from the database */
    
    if (write_files) {
        
        phase_file = fopen(phases_filename, "w");
        flag_file  = fopen(flag_filename,   "w");
        jones_file = fopen(jones_filename,  "w");

        if (phase_file == NULL || flag_file == NULL || jones_file == NULL) {
            fprintf(stderr, "Failed to open %s or %s or %s\n", phases_filename, flag_filename, jones_filename );
            exit(EXIT_FAILURE);
        }
    }

    double *phase = calloc(nchan,sizeof(double));

    /* Calibration related defines */
    /* set these here for the library */
    nfrequency = nchan;
    nstation = NANT;
    npol = 2;
    
    
    double amp = 0;
    
    // Jref is the reference Jones direction - for which the calibration was generated
    // Mref is measured Jones in Jref direction
    // G = DI gain formed by inv(Jref).Mref
    // E = model is desired direction
    // Ji = Jones in desired direction (formed by G.E)
    
    complex double *Jref = (complex double *) calloc(npol*npol, sizeof(complex double)); // Calibration Direction
    complex double *E = (complex double *) calloc(npol*npol, sizeof(complex double)); // Model Jones	in Desired Direction	//
    complex double **G = (complex double **) calloc(nstation, sizeof(complex double *)); // Actual DI Gain			//
    complex double **M = (complex double **) calloc(nstation, sizeof(complex double *)); // Gain in direction of Calibration
    complex double **Ji = (complex double **) calloc(nstation, sizeof(complex double *)); // Gain in Desired Direction ..... da da da dum.....

    for (i = 0; i < nstation; i++) { //
        G[i] = (complex double *) malloc(npol * npol * sizeof(complex double)); //
        M[i] = (complex double *) malloc(npol * npol * sizeof(complex double)); //
        Ji[i] =(complex double *) malloc(npol * npol * sizeof(complex double)); //
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
    int *flag_array = (int *)malloc(ninput*sizeof(int));
    short int *antenna_num = (short int *)malloc(ninput*sizeof(short int));
    int colnum;
    
    /* read the columns */
    fits_get_colnum(fptr, 1, "Flag", &colnum, &status);
    fits_read_col_int(fptr,colnum,1,1,ninput,0.0,flag_array,&anynull,&status);
    if (status != 0){
        fprintf(stderr,"Error:Failed to read flags column in metafile\n");
        exit(-1);
    }
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
        read_rts_file(M, Jref, nstation, &amp, DI_Jones_file);
        inv2x2(Jref, invJref);
    }
    else if (get_offringa) {
        // Find the ordering of antennas in Offringa solutions from metafits file
        int *order = (int *)malloc(nstation*sizeof(int));
        for (n = 0; n < nstation; n++) {
            order[antenna_num[n*2]] = n;
        }
        read_offringa_gains_file(M, nstation, coarse_chan, DI_Jones_file, order);
        //read_offringa_gains_file(M, nstation, coarse_chan, DI_Jones_file, NULL);
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
  
    
    utc2mjd(time_utc,&intmjd,&fracmjd);
    

    //fprintf(stdout,"calib:start intmjd=%.1f fracmjd=%f\n",intmjd,fracmjd);
    
    /**/
    
    /* get requested Az/El from command line */
  
    mjd = intmjd + fracmjd;

    mjd = mjd + (sec_offset+0.5)/86400.0;

    mjd2lst(mjd,&lmst);

    //fprintf(stdout,"calib: current lmst (radian) = %f\n",lmst);
    // last is in epoch of MJD (
    
    /* for the look direction <not the tile> */
    
    mean_ra = ra_hours * DH2R;
    mean_dec = dec_degs * DD2R;

    //fprintf(stdout,"calib:Mean (HA RA /Dec) HA %lf RA %lf (radian)  Apparent Dec %f (radian) \n",ha*DH2R,mean_ra,mean_dec);
    slaMap(mean_ra,mean_dec,pr,pd,px,rv,eq,mjd,&ra_ap,&dec_ap);
    
    // Lets go mean to apparent precess from J2000.0 to EPOCH of date.
    //
    
    ha = slaRanorm(lmst-ra_ap)*DR2H;

    
    //fprintf(stdout,"calib:Apparent (look direction/precessed HA RA /Dec) HA %lf RA %lf hrs  Apparent Dec %f degrees \n",ha,ra_ap*DR2H,dec_ap*DR2D);
    //fprintf(stdout,"calib:Apparent (look direction/precessed HA RA /Dec) HA %lf RA %lf (radian)  Apparent Dec %f (radian) \n",ha*DH2R,ra_ap,dec_ap);
    //fprintf(stdout,"calib:LMST %lf (radian) \n",lmst);

    /* now HA/Dec to Az/El */
    
    app_ha_rad = ha * DH2R;
    app_dec_rad = dec_ap;
  
    
    slaDe2h(app_ha_rad,dec_ap,MWA_LAT*DD2R,&az,&el);
    

    //fprintf(stdout,"calib:Look direction Azimuth %lf (deg)  Elevation %lf (deg) \n",az*DR2D,el*DR2D);
    
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
        for (i=0;i<nstation;i++){
            mult2x2d(M[i],invJref,G[i]); // forms the DI gain
            mult2x2d(G[i],E,Ji[i]); // the gain in the desired look direction
            
            for (j=0; j < 4;j++) {
                //fprintf(stdout,"calib:Jones Mi[%d] %f %f: Delay Ji[%d] %f %f\n",i,creal(M[i][j]),cimag(M[i][j]),i,creal(Ji[i][j]),cimag(Ji[i][j]));
                //fprintf(stdout,"calib:ratio Mi/Ji [%d]  %f %f \n",i,creal(M[i][j])/creal(Ji[i][j]),cimag(M[i][j])/cimag(Ji[i][j]));
            }
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
    

    //fprintf(stdout,"calib: Tile position (degrees) (Az:%f,El:%f,RA:%f,Dec:%f) -- Not precessing -- assuming fixed Az-El \n", tile_pointing_az,tile_pointing_el,tile_pointing_ra,tile_pointing_dec );
     //fprintf(stdout,"calib: Tile position (radian) (Az:%f,El:%f,RA:%f,Dec:%f) -- Not precessing -- assuming fixed Az-El \n", tile_pointing_az*DD2R,tile_pointing_el*DD2R,tile_pointing_ra*DD2R,tile_pointing_dec*DD2R);
    //fprintf(stdout,"calib: Requested  Look direction (degrees) (Az:%f,El:%f,RA:%f,Dec:%f) -- After precession \n", az*DR2D,el*DR2D,ra_ap*DR2D,dec_ap*DR2D);
    //fprintf(stdout,"calib: Requested  Look direction (degrees) (Az:%f,El:%f,RA:%f,Dec:%f) -- After precession \n", az,el,ra_ap,dec_ap);

    
    
    int ch=0;

    for (row=0; row < (int)ninput; row++) {


        double cable = cable_array[row]-cable_array[refinp];
        double E = E_array[row];
        double N = N_array[row];
        double H = H_array[row];


        double integer_phase;

        int flag = flag_array[row];

        double X,Y,Z,u,v,w;
        
        ENH2XYZ_local(E,N,H, MWA_LAT*DD2R, &X, &Y, &Z);

        //fprintf(stdout,"calib: Antenna %d: HA %f Dec %f --  X: %f Y: %f Z: %f\n",row,app_ha_rad, app_dec_rad,X,Y,Z);
        calcUVW (app_ha_rad,app_dec_rad,X,Y,Z,&u,&v,&w);

        // shift the origin of ENH to Antenna 0 and hoping the Far Field Assumption still applies ...



        double geometry = (E-E_ref)*unit_E + (N-N_ref)*unit_N + (H-H_ref)*unit_H ;
        // double geometry = E*unit_E + N*unit_N + H*unit_H ;
        // Above is just w as you should see from the check.

        double delay_time = (geometry + (invert*(cable)))/(VLIGHT);
        double delay_samples = delay_time * samples_per_sec;

        //fprintf(stdout,"Antenna %d, E %f, N %f, H %f\n",row,E,N,H);
        //fprintf(stdout,"Distance from reference, E-E_ref %f, N-N_ref %f, H-N_ref %f\n",E-E_ref,N-N_ref,H-H_ref);

        //fprintf(stdout,"Look direction, E %f, N %f, H %f\n",unit_E,unit_N,unit_H);
        //fprintf(stdout,"calib:geom: %f w: %f cable(-cable_ref): %f time (s):%g (samples):%g \n",geometry, w, cable, delay_time,delay_samples);
        //fprintf(stdout,"calib:geom: u %f v %f w %f\n",u,v,w);// we have to get this amount of delay into the data




        for (ch=0;ch<nchan;ch++) {
            long int freq_ch = frequency + (edge+ch)*chan_width;

            // freq should be in cycles per sample and delay in samples
            // which essentially means that the samples_per_sec cancels

            // and we only need the fractional part of the turn
            double cycles_per_sample = (double)freq_ch/samples_per_sec;

            phase[ch] = cycles_per_sample*delay_samples;
            phase[ch] = modf(phase[ch], &integer_phase);

            phase[ch] = phase[ch]*2*M_PI*conjugate;

            if (ch == 0) {
                //fprintf(stdout,"Comp:ch %d Freq (Cycles/s) %ld\n",ch,freq_ch);
                //fprintf(stdout,"Comp:ch %d Freq (Cycles/sample) %lf\n",ch,(double)freq_ch/samples_per_sec);
                //fprintf(stdout,"Comp:Geo: %f Cable %f (total (s)) %g:Phase (raw) %f Phase (sample) %f\n",geometry,cable,(geometry+cable)/VLIGHT,phase[ch],phase[ch]/samples_per_sec);
            }

            if (phase_file != NULL) {
                fprintf(phase_file,"%lf\n",phase[ch]);
            }
            
        }
        if (secs == 0) {
            
            if (flag_file != NULL) {
                if (flag > 0) {
                    fprintf(flag_file,"0.0\n");
                }
                else {
                    fprintf(flag_file,"1.0\n");
                }
            }
        }
        if (row%npol == 0) {
            if (jones_file != NULL) {
                for (i=0;i<4;i++){
                    fprintf(jones_file,"%f %f ",creal(Ji[row/npol][i]), cimag(Ji[row/npol][i]));
                    //fprintf(jones_file,"%f %f ",creal(G[row/npol][i]), cimag(G[row/npol][i]));
                }
                fprintf(jones_file,"\n");
            }
        }
        
        
    }

    if (verbose)
        puts("==========================");
    
    if (phase_file != NULL)
        fclose(phase_file);
    if (flag_file !=NULL)
        fclose(flag_file);
    if (jones_file != NULL)
        fclose(jones_file);
      /* ========= Generate a FITS HEADER ==========*/
    if (get_psrfits) {
        
        
        struct psrfits pf;
        
        status=0;

        fits_open_file(&fptr,metafits,READONLY,&status);
        fits_read_key(fptr,TSTRING,"PROJECT",pf.hdr.project_id,NULL,&status);
        fits_close_file(fptr,&status);

        strcpy(pf.basefilename, "/tmp/obsfile");
        
        // Now set values for our hdrinfo structure
        strcpy(pf.hdr.obs_mode,"SEARCH");
        pf.hdr.scanlen = 1.0; // in sec
        strcpy(pf.hdr.observer, "MWA User");
        strcpy(pf.hdr.telescope, "MWA");
        strncpy(pf.hdr.source,obsid,23);
        
        
        strcpy(pf.hdr.frontend, "MWA-RECVR");
        char backend[24];
        snprintf(backend, 24*sizeof(char), "GD-%s-MB-%s-U-%s", GET_DELAYS_VERSION, MAKE_BEAM_VERSION, UTILS_VERSION);
        strcpy(pf.hdr.backend, backend);
        
        /* Now let us finally get the time right */
        
        strcpy(pf.hdr.date_obs, time_utc);
        strcpy(pf.hdr.poln_type, "LIN");
        strcpy(pf.hdr.track_mode, "TRACK");
        strcpy(pf.hdr.cal_mode, "OFF");
        strcpy(pf.hdr.feed_mode, "FA");
        
        pf.hdr.dt   = 1.0/samples_per_sec;                                   // sample rate (s)
        pf.hdr.fctr = (frequency + (edge+(nchan/2.0))*chan_width)/1.0e6;     // frequency (MHz)
        pf.hdr.BW   = (nchan*chan_width)/1.0e6;

        pf.hdr.ra2000    = mean_ra  * DR2D;
        pf.hdr.dec2000   = mean_dec * DR2D;

        dec2hms(pf.hdr.ra_str, pf.hdr.ra2000/15.0, 0);
        dec2hms(pf.hdr.dec_str, pf.hdr.dec2000, 1);

        pf.hdr.azimuth    = az*DR2D;
        pf.hdr.zenith_ang = 90.0 - (el*DR2D);

        pf.hdr.beam_FWHM = 0.25;
        pf.hdr.start_lst = lmst * 60.0 * 60.0;        // Local Apparent Sidereal Time in seconds
        pf.hdr.start_sec = roundf(fracmjd*86400.0);   // this will always be a whole second - so I'm rounding.
        pf.hdr.start_day = intmjd;

        pf.hdr.scan_number   = 1;
        pf.hdr.rcvr_polns    = 2;
        pf.hdr.summed_polns  = 0;
        pf.hdr.offset_subint = 0;

        pf.hdr.nchan      = nchan;
        pf.hdr.df         = chan_width/1.0e6;
        pf.hdr.orig_nchan = pf.hdr.nchan;
        pf.hdr.orig_df    = pf.hdr.df;
        pf.hdr.nbits      = 8;
        pf.hdr.npol       = 2;
        pf.hdr.nsblk      = samples_per_sec;  // block is always 1 second of data
        pf.hdr.MJD_epoch  = intmjd + fracmjd;
        
        pf.hdr.ds_freq_fact = 1;
        pf.hdr.ds_time_fact = 1;

        // some things that we are unlikely to change
        pf.hdr.fd_hand  = 0;
        pf.hdr.fd_sang  = 0.0;
        pf.hdr.fd_xyph  = 0.0;
        pf.hdr.be_phase = 0.0;
        pf.hdr.chan_dm  = 0.0;

        // Now set values for our subint structure
        pf.tot_rows     = 0;
        pf.sub.tsubint  = roundf(pf.hdr.nsblk * pf.hdr.dt);
        pf.sub.offs     = roundf(pf.tot_rows * pf.sub.tsubint) + 0.5*pf.sub.tsubint;
        pf.sub.lst      = pf.hdr.start_lst;
        pf.sub.ra       = pf.hdr.ra2000;
        pf.sub.dec      = pf.hdr.dec2000;
        slaEqgal(pf.hdr.ra2000*DD2R, pf.hdr.dec2000*DD2R,
                 &pf.sub.glon, &pf.sub.glat);
        pf.sub.glon    *= DR2D;
        pf.sub.glat    *= DR2D;
        pf.sub.feed_ang = 0.0;
        pf.sub.pos_ang  = 0.0;
        pf.sub.par_ang  = 0.0;
        pf.sub.tel_az   = pf.hdr.azimuth;
        pf.sub.tel_zen  = pf.hdr.zenith_ang;
        pf.sub.bytes_per_subint = (pf.hdr.nbits * pf.hdr.nchan *
                                   pf.hdr.npol * pf.hdr.nsblk) / 8;
        pf.sub.FITS_typecode = TBYTE;  // 11 = byte
        
        // Print out values to stdout
        printf_psrfits( &pf );
        
        /* now copy this to a file that can be read by the beamformer */
        if (write_psrfits) {
            psrfits_file = fopen( psrfits_filename, "w" );
            if (psrfits_file != NULL) {
                fwrite((void *) &pf, sizeof(pf),1,psrfits_file);
            }
            fclose(psrfits_file);
        }
    }
    

    // Free up memory
    free(phase);
    

    for (i = 0; i < nstation; i++) { //
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
    for (i = 0; i < ninput; i++)
        free(tilenames[i]);
    free(tilenames);
    free(cable_array);
    free(testval);
    free(flag_array);
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

