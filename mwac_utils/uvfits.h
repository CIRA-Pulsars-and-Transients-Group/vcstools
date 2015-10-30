/* structure definitions and utilities for reading & writing
 * UV FITS files in C. Randall Wayth. Sep, 2006
 */

#include <time.h>

#define SIZE_ANT_NAME 8
#define SIZE_SOURCE_NAME 16
#define SIZE_POL_TYPE 4

/* constants for WGS84 Geoid model for the Earth */
#define EARTH_RAD_WGS84 6378137.0 /* meters */
#define E_SQUARED 6.69437999014e-3

typedef struct  _ant_table {
  char name[SIZE_ANT_NAME+1];
  double xyz_pos[3];			/* position relative to array centre  */
  float  xyz_deriv[3];
  /* skip orbital params for now */
  int  station_num;
  int  mount_type;
  double axis_offset[3];
  float pol_angleA;
  float pol_calA;
  float pol_angleB;
  float pol_calB;
  char pol_typeA[SIZE_POL_TYPE];
  char pol_typeB[SIZE_POL_TYPE];
} ant_table;

typedef struct _source_table {
  char name[SIZE_SOURCE_NAME+1];
  int  id;
  int  qual;
  char calcode[4];
  int  freq_id;
  double ra;    /* decimal hours */
  double dec;   /* decimal degrees */
} source_table;

typedef struct _array_table {
  int   n_ant;
  double xyz_pos[3];    /* the X,Y,Z coord of the array in conventional radio astronomy units */
  char  name[16];
} array_data;


typedef struct _uvdata {
  int  n_pol;
  int  pol_type;        /* index to signal what kind of data: 1: Stokes, -1: Circ, -5: Linear */
  int  n_freq;
  int  n_vis;           /* number of sets of visibilities, one for each time instant */
  float cent_freq;
  float freq_delta;
  double *date;         /* Julian date. array the size of n_vis */
  int  *n_baselines;    /* array the size of n_vis. number of baselines for each scan. */
  source_table *source; /* a pointer to a source table */
  ant_table *antennas;  /* a pointer to an array of ant_tables the size of the number of antennas */
  array_data *array;    /* a pointer to an array struct */
  float **visdata;      /* array the size of n_vis whose elements point to arrays of visibiliites
                           the size of n_freq*n_pol*n_baselines complex floats */
  float **weightdata;   /* weight data for visibiliites. Same data ordering as above, just float,
                           not complex */
  float **baseline;     /* same ordering again. encoded baseline using Miriad encoding convention */
  double **u;           /* arry the size of n_vis whose elements point to arrays of uvw
                           data the size of n_baselines */
  double **v;
  double **w;
} uvdata;

/* public function prototypes */
int writeUVFITS(char *fname, uvdata *data);
int readUVFITS(char *fname, uvdata **data);
void printUVData(uvdata *data, FILE *fp);
void JD_to_Cal(double jd, int *year, int *month, int *day);
void JD_get_GSTIA0(double jd, double *GSTIA0);
void Cal_to_JD(int year, int month, int day, double *jd);
void uvfitsSetDebugLevel(int in_debug);
void freeUVFITSdata(uvdata *data);
void printAntennaData(uvdata *data,FILE *fp);
void EncodeBaseline(int b1, int b2, float *result);
void DecodeBaseline(float blcode, int *b1, int *b2);
void Geodetic2XYZ(double lat_rad, double lon_rad, double height_meters, double *X, double *Y, double *Z);
void ENH2XYZ_absolute(double E,double N, double H, double lat_rad, double lon_rad, double *X, double *Y, double *Z);
void ENH2XYZ_local(double E,double N, double H, double lat, double *X, double *Y, double *Z);

