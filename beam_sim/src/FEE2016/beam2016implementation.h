#ifndef _GET_FF_H__
#define _GET_FF_H__

/**
 * C++ implementation of Full Embeded Element beam model for MWA based on beam_full_EE.py
 * script and Sokolowski et al (2017) paper
 * Implemented by Marcin Sokolowski (May 2017) - marcin.sokolowski@curtin.edu.au
 */

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <boost/math/special_functions/legendre.hpp>

#include <H5Cpp.h>

// default 
#define DEFAULT_H5_FILE "mwa_full_embedded_element_pattern.h5"
#define DEFAULT_H5_FILE_PATH "mwapy/data/"
#define N_ANT_COUNT 16

/*
  Structure for Jones matrix :
  j00 j01
  j10 j11
*/ 
class JonesMatrix 
{
public :
   std::complex<double> j00; 
   std::complex<double> j01; 
   std::complex<double> j10; 
   std::complex<double> j11; 

   JonesMatrix( double j00_r=0.00, double j01_r=0.00, double j10_r=0.00, double j11_r=0.00 ) : j00(j00_r,0), j01(j01_r,0), j10(j10_r,0), j11(j11_r,0){
   }      
   
   // to be replaced with proper std initialisation :
   static void zeros( std::vector< std::vector<JonesMatrix> >& jones, int x_size, int y_size );
   
   void Print(const char* name,double az_deg=0,double za_deg=0);   
};


class Beam2016Implementation
{
public :
   Beam2016Implementation( const double* delays, const double* amps );      
   ~Beam2016Implementation();
   

   //-------------------------------------------------------------------- Calculation of Jones matrix ------------------------------------------------------------------------------------------------------------------
   // Calculate jones matrix in a specified direction for a given frequency, delays and amplitudes. Zenith normalisation can also be disabled - but by default is enabled :
   // This function will  :
   //    - calculate coefficients of spherical waves (SPW) if required (meaning if frequency, delays or amplitudes are different than last set of coefficients was calculated for)
   //    - calculate electric fields (equations 4 and 5 in Sokolowski et al paper)
   //    - normalise Jones matrix to one at zenith at the same frequency (if required by parameter):
   // 
   // Only thse the functions should be used for external calls. The rest is protected from external use.
   // 
   // INPUT : (az_deg,za_deg) - azimuth and zenith angles [in degrees]
   //         freq_hz_param   - frequency in Hz
   //         delays          - delays in beamformer steps 
   //         amps            - amplitudes 
   //         bZenithNorm     - normalise to zenith (>0) or not (<=0)
   // OUTPUT : Jones matrix (normalised or not - depending on the parameter bZenithNorm )
   JonesMatrix CalcJones( double az_deg, double za_deg, int freq_hz_param, bool bZenithNorm=true );
   JonesMatrix CalcJones( double az_deg, double za_deg, int freq_hz_param, const double* delays, const double* amps, bool bZenithNorm=true );

   // Calculation of Jones matrix for an image passed in the arrays azimuth and zenith angles maps. 
   // This function calls the single direction one (above) for all pixel in the input image.
   // INPUT : 
   //       2D array of azimuth angles in degrees
   //       2D array of zenith angles in degrees
   //       freq_hz_param   - frequency in Hz
   //       delays          - delays in beamformer steps 
   //       amps            - amplutudes          
   //       bZenithNorm     - normalise to zenith (>0) or not (<=0)
   // OUTPUT :
   //       2D array of JonesMatrix for each pixel 
   void CalcJonesArray( std::vector< std::vector<double> >& azim_arr, std::vector< std::vector<double> >& za_arr, std::vector< std::vector<JonesMatrix> >& jones,
                   int freq_hz_param, const double* delays=NULL, const double* amps=NULL, bool bZenithNorm=true );
protected :

   // Calculation of Jones matrix for a single pointing direction (internal function):
   // INPUT : (az_rad,za_rad) - azimuth and zenith in [radians]
   JonesMatrix CalcJonesDirect( double az_rad, double za_rad );

   //-------------------------------------------------------------------- Calculation of Jones matrix components for given polarisation (eq. 4 and 5 in the Sokolowski et al (2017) paper --------------------------------
   // Internal function to calculate Jones matrix componenets for a given polarisation (pol). The are calculated as electric field vectors E_theta_mn (eq.4) and E_phi_mn (eq.5 in the Sokolowski et al 2017 paper).
   // INPUT :
   //    phi,theta - FEKO convention coordinates (phi=90-azim), theta=za in radians already
   //    Q1_accum, Q2_accum, M_accum, N_accum, MabsM, Nmax - coefficients calculated earlier in CalcModes function for given frequency, delays and amplitudes 
   // OUTPUT :
   //    Jones matrix filled with components for given polarisation (pol input parameter)
   void CalcSigmas( double phi, double theta, 
                    std::vector< std::complex<double> >& Q1_accum, std::vector< std::complex<double> >& Q2_accum, 
                    std::vector<double>& M_accum, std::vector<double>& N_accum, std::vector<double>& MabsM, double Nmax, std::vector<double>& Cmn,
                    char pol,
                    JonesMatrix& jones_matrix );
   
   //---------------------------------------------------- Calculations and variables for spherical waves coefficients (see equations 3-6 in the Sokolowski et al paper) ----------------------------------------------------
   // Coefficients of spherical waves (SPW) - see equations 3-6 in the Sokolowski et al paper
   // 1 refers to s=1 (transverse electric - TE modes) - eq.5 
   // 2 refers to s=2 (transverse magnetic - TM modes) - eq.5 
   // These coefficients are calculated ones for a given frequency, delays, amplitudes so these parameters are stored to 
   // only calculate new coefficients when they change.
   // X polarisation :
   std::vector< std::complex<double> > Q1_accum_X;
   std::vector< std::complex<double> > Q2_accum_X;
   std::vector<double> M_accum_X;
   std::vector<double> N_accum_X;
   std::vector<double> MabsM_X; // precalculated m/abs(m) to make it once for all pointings
   double Nmax_X;          // maximum N coefficient for Y (=max(N_accum_X)) - to avoid relaculations 
   std::vector<double> Cmn_X; // coefficient under sumation in equation 3 for X pol.
    

   // Y polarisation :
   std::vector< std::complex<double> > Q1_accum_Y;
   std::vector< std::complex<double> > Q2_accum_Y;
   std::vector<double> M_accum_Y;
   std::vector<double> N_accum_Y;
   std::vector<double> MabsM_Y; // precalculated m/abs(m) to make it once for all pointings
   double Nmax_Y;          // maximum N coefficient for Y (=max(N_accum_Y)) - to avoid relaculations
   std::vector<double> Cmn_Y; // coefficient under sumation in equation 3 for Y pol
   
   std::vector< std::complex<double> > j_power_cache; // powers of j calculated once 
   
   // Information on last modes parameters - not to recalculate the same again and again !
   int     m_CalcModesLastFreqHz;
   double* m_CalcModesLastDelays;
   double* m_CalcModesLastAmps;
   
   // function comparing current parameters : frequency, delays and amplitudes with those previously used to calculate spherical waves coefficients (stored in the 3 variables above)
   int IsCalcModesRequired( int freq_hz, int n_ant, const double* delays, const double* amps );
   
   // Calculation of modes Q1 and Q2 and coefficients N and M and some derived variables (MabsM_X,MabsM_Y,Nmax_X and Nmax_Y) to make it once for a given pointing and 
   // then re-use for many different (azim,za) angles:
   
   // function calculating coefficients for X and Y and storing parameters frequency, delays and amplitudes 
   void CalcModes( int freq_hz, size_t n_ant, const double* delays, const double* amps );
   
   // function calculating all coefficients Q1, Q2, N, M and derived MabsM, Nmax for a given polarisation ("X" or "Y") - perhaps enum should be used here 
   double CalcModes( int freq_hz, size_t n_ant, const double* delays, const double* amp, char pol,
                     std::vector< std::complex<double> >& Q1_accum, std::vector< std::complex<double> >& Q2_accum,
                     std::vector<double>& M_accum, std::vector<double>& N_accum, std::vector<double>& MabsM, std::vector<double>& Cmn   );
   // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

   // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   // Calculation of normalisation matrix :
   void CalcZenithNormMatrix( int freq_hz );
   JonesMatrix m_NormJones;
   int         m_NormFreqHz;
   // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

   // Some static/constant variables like : default delays and amps:
   static const double m_DefaultDelays[];  // at zenith {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
   static const double m_DefaultAmps[];    // just {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
   static const double _delayStep;          // delay step in MWA beamformer (in picoseconds)
	double _delays[16], _amps[16];

   // parameters:
   static int m_VerbLevel;  // debug level
   
   // ----------------------------------------------------------------
   // Checking freuquences included in the H5 file (stored in vector m_freq_list) :
   int has_freq(int freq_hz);
   int find_closest_freq(int freq_hz);         

   //----------------------------------------------------------------------------------- HDF5 File interface and data structures for H5 data -----------------------------------------------------------------------------------
   // Interface to HDF5 file format and structures to store H5 data 
   // Functions for reading H5 file and its datasets :
   // Read data from H5 file, file name is specified in the object constructor 
   int Read();      
   
   // number of antennas simulated (as read from the H5 file) :   
   // the code can be made to automatically detect simulated number of antennas, just need to make some arrays dynamic size and 
   // change places with N_ANT_COUNT define used 
   int    m_AntennaCount;                     

   // Read dataset_name from H5 file 
   void ReadDataSet( const char* dataset_name, std::vector< std::vector<double> >& out_vector );
   
   // function for iteration call to H5Ovisit function :
   static herr_t list_obj_iterate(hid_t loc_id, const char *name, const H5O_info_t *info, void *operator_data);
      
      
   std::string _h5filename;   // H5 File name 
   H5::H5File* m_pH5File;   
   
   // Data structures for H5 file data :
   std::vector<std::string> m_obj_list;  // list of datasets in H5 file 
   std::vector<int>    m_freq_list; // list of simulated frequencies 
   std::vector< std::vector<double> > m_Modes;  // data in Modes DataSet 

   //------------------------------------------------------------------------------------------------------ maths functions and wrappers ---------------------------------------------------------------------------------------
   // local power calculations
   // WARNING : please do not remove - the program works 2x slower when std::power is used instead !!!
   std::complex<double> power_complex( std::complex<double> val, int n );

   // Calculations of Legendre polynomials :
   void lpmv( std::vector<double>& output, int n, double x );
   
   // OUTPUT : returns list of Legendre polynomial values calculated up to order nmax :
   int P1sin( int nmax, double theta, std::vector<double>& p1sin_out, std::vector<double>& p1_out );
      
   
   // factorial calculation :
   double factorial_wrapper_base( unsigned n );
   double factorial_wrapper( unsigned n ); // uses precalculated factorials stored in m_Factorial
   void cache_factorial( unsigned max_n );                     
   // precalculated factorials :
   std::vector<double> m_Factorial;   


   //----------------------------------------------------------------------------------- debuging / logging functions ----------------------------------------------------------------------------------------------------------------
   void print( std::vector< std::complex<double> >& Q1, std::vector< std::complex<double> >& Q1_accum, std::vector< std::complex<double> >& Q2, std::vector< std::complex<double> >& Q2_accum, int my_len_half, 
               std::vector<double>& N_accum, std::vector<double>& M_accum, double Nmax );
   void print( std::vector<double>& arr, const char* name, int force=0 );
   void print( std::vector<int>& arr, const char* name, int force=0 );
   void print( std::vector< std::vector<double> >& arr, const char* name, int force=0 );
};


#endif
