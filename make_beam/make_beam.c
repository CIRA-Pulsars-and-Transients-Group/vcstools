/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

// TODO: Remove superfluous #includes
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include "slalib.h"
#include "slamac.h"
#include "ascii_header.h"
#include "mwa_header.h"
#include <glob.h>
#include <fcntl.h>
#include <assert.h>
#include "beam_common.h"
#include "beam_psrfits.h"
#include "beam_vdif.h"
#include "make_beam.h"
#include "vdifio.h"
#include "filter.h"
#include "psrfits.h"
#include "mycomplex.h"
#include "form_beam.h"
#include <omp.h>
#include <cuda_runtime.h>

#ifdef HAVE_CUDA

#include <cuda_runtime.h>
#include "ipfb.h"
#define NOW  ((double)clock()/(double)CLOCKS_PER_SEC)

#else

#define NOW  (omp_get_wtime())

#endif

int main(int argc, char **argv)
{
    #ifndef HAVE_CUDA
    // Initialise FFTW with OpenMP
    fftw_init_threads();
    fftw_plan_with_nthreads( omp_get_max_threads() );
    #endif

    // A place to hold the beamformer settings
    struct make_beam_opts opts;

    /* Set default beamformer settings */

    // Variables for required options
    opts.obsid       = NULL; // The observation ID
    opts.begin       = 0;    // GPS time -- when to start beamforming
    opts.end         = 0;    // GPS time -- when to stop beamforming
    opts.time_utc    = NULL; // utc time string "yyyy-mm-ddThh:mm:ss"
    opts.pointings   = NULL; // list of pointings "dd:mm:ss_hh:mm:ss,dd:mm:ss_hh:mm:ss"
    opts.datadir     = NULL; // The path to where the recombined data live
    opts.metafits    = NULL; // filename of the metafits file
    opts.rec_channel = NULL; // 0 - 255 receiver 1.28MHz channel
    opts.frequency   = 0;    // = rec_channel expressed in Hz
    opts.beam_model  = NULL; // HDF5 file containing FEE 2016 beam model, or 'analytic'

    // Variables for MWA/VCS configuration
    opts.nstation      = 128;    // The number of antennas
    opts.nchan         = 128;    // The number of fine channels (per coarse channel)
    opts.chan_width    = 10000;  // The bandwidth of an individual fine chanel (Hz)
    opts.sample_rate   = 10000;  // The VCS sample rate (Hz)
    opts.custom_flags  = NULL;   // Use custom list for flagging antennas

    // Output options
    opts.out_incoh     = 0;  // Default = PSRFITS (incoherent) output turned OFF
    opts.out_coh       = 0;  // Default = PSRFITS (coherent)   output turned OFF
    opts.out_vdif      = 0;  // Default = VDIF                 output turned OFF
    opts.out_bf        = 1;  // Default = beamform all (non-flagged) antennas
    opts.out_ant       = 0;  // The antenna number (0-127) to write out if out_bf = 0
    opts.synth_filter  = NULL;
    opts.out_summed    = 0;  // Default = output only Stokes I output turned OFF

    // Variables for calibration settings
    opts.cal.filename          = NULL;
    opts.cal.bandpass_filename = NULL;
    opts.cal.chan_width        = 40000;
    opts.cal.nchan             = 0;
    opts.cal.cal_type          = NO_CALIBRATION;
    opts.cal.offr_chan_num     = 0;

    // Parse command line arguments
    make_beam_parse_cmdline( argc, argv, &opts );

    // Create "shorthand" variables for options that are used frequently
    int nstation             = opts.nstation;
    int nchan                = opts.nchan;
    const int npol           = 2;   // (X,Y)
    int outpol_coh           = 4;  // (I,Q,U,V)
    if ( opts.out_summed )
        outpol_coh           = 1;  // (I)
    const int outpol_incoh   = 1;  // ("I")

    float vgain = 1.0; // This is re-calculated every second for the VDIF output

    // Start counting time from here (i.e. after parsing the command line)
    double begintime = NOW;
    fprintf( stderr, "[%f]  Starting %s with GPU acceleration\n", NOW-begintime, argv[0] );

    // Calculate the number of files
    int nfiles = opts.end - opts.begin + 1;
    if (nfiles <= 0) {
        fprintf(stderr, "Cannot beamform on %d files (between %lu and %lu)\n", nfiles, opts.begin, opts.end);
        exit(EXIT_FAILURE);
    }

    // Parse input pointings
    int max_npointing = 120; // Could be more
    char RAs[max_npointing][64];
    char DECs[max_npointing][64];
    int npointing = sscanf( opts.pointings, 
            "%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,]," , 
                            RAs[0],  DECs[0],  RAs[1],  DECs[1],  RAs[2],  DECs[2],
                            RAs[3],  DECs[3],  RAs[4],  DECs[4],  RAs[5],  DECs[5],
                            RAs[6],  DECs[6],  RAs[7],  DECs[7],  RAs[8],  DECs[8],
                            RAs[9],  DECs[9],  RAs[10], DECs[10], RAs[11], DECs[11],
                            RAs[12], DECs[12], RAs[13], DECs[13], RAs[14], DECs[14],
                            RAs[15], DECs[15], RAs[16], DECs[16], RAs[17], DECs[17],
                            RAs[18], DECs[18], RAs[19], DECs[19], RAs[20], DECs[20],
                            RAs[21], DECs[21], RAs[22], DECs[22], RAs[23], DECs[23],
                            RAs[24], DECs[24], RAs[25], DECs[25], RAs[26], DECs[26],
                            RAs[27], DECs[27], RAs[28], DECs[28], RAs[29], DECs[29],
                            RAs[30], DECs[30], RAs[31], DECs[31], RAs[32], DECs[32],
                            RAs[33], DECs[33], RAs[34], DECs[34], RAs[35], DECs[35],
                            RAs[36], DECs[36], RAs[37], DECs[37], RAs[38], DECs[38],
                            RAs[39], DECs[39], RAs[40], DECs[40], RAs[41], DECs[41],
                            RAs[42], DECs[42], RAs[43], DECs[43], RAs[44], DECs[44],
                            RAs[45], DECs[45], RAs[46], DECs[46], RAs[47], DECs[47],
                            RAs[48], DECs[48], RAs[49], DECs[49], RAs[50], DECs[50],
                            RAs[51], DECs[51], RAs[52], DECs[52], RAs[53], DECs[53],
                            RAs[54], DECs[54], RAs[55], DECs[55], RAs[56], DECs[56],
                            RAs[57], DECs[57], RAs[58], DECs[58], RAs[59], DECs[59],
                            RAs[60], DECs[60], RAs[61], DECs[61], RAs[62], DECs[62],
                            RAs[63], DECs[63], RAs[64], DECs[64], RAs[65], DECs[65],
                            RAs[66], DECs[66], RAs[67], DECs[67], RAs[68], DECs[68],
                            RAs[69], DECs[69], RAs[70], DECs[70], RAs[71], DECs[71],
                            RAs[72], DECs[72], RAs[73], DECs[73], RAs[74], DECs[74],
                            RAs[75], DECs[75], RAs[76], DECs[76], RAs[77], DECs[77],
                            RAs[78], DECs[78], RAs[79], DECs[79], RAs[80], DECs[80],
                            RAs[81], DECs[81], RAs[82], DECs[82], RAs[83], DECs[83],
                            RAs[84], DECs[84], RAs[85], DECs[85], RAs[86], DECs[86],
                            RAs[87], DECs[87], RAs[88], DECs[88], RAs[89], DECs[89],
                            RAs[90], DECs[90], RAs[91], DECs[91], RAs[92], DECs[92],
                            RAs[93], DECs[93], RAs[94], DECs[94], RAs[95], DECs[95],
                            RAs[96], DECs[96], RAs[97], DECs[97], RAs[98], DECs[98],
                            RAs[99], DECs[99], RAs[100], DECs[100], RAs[101], DECs[101],
                            RAs[102], DECs[102], RAs[103], DECs[103], RAs[104], DECs[104],
                            RAs[105], DECs[105], RAs[106], DECs[106], RAs[107], DECs[107],
                            RAs[108], DECs[108], RAs[109], DECs[109], RAs[110], DECs[110],
                            RAs[111], DECs[111], RAs[112], DECs[112], RAs[113], DECs[113],
                            RAs[114], DECs[114], RAs[115], DECs[115], RAs[116], DECs[116],
                            RAs[117], DECs[117], RAs[118], DECs[118], RAs[119], DECs[119] );

    if (npointing%2 == 1)
    {
        fprintf(stderr, "Number of RAs do not equal the number of Decs given. Exiting\n");
        fprintf(stderr, "npointings : %d\n", npointing);
        fprintf(stderr, "RAs[0] : %s\n", RAs[0]);
        fprintf(stderr, "DECs[0] : %s\n", DECs[0]);
        exit(0);
    }
    else
        npointing /= 2; // converting from number of RAs and DECs to number of pointings

    char pointing_array[npointing][2][64];
    int p;
    for ( p = 0; p < npointing; p++) 
    {
       strcpy( pointing_array[p][0], RAs[p] );
       strcpy( pointing_array[p][1], DECs[p] );
       fprintf(stderr, "[%f]  Pointing Num: %i  RA: %s  Dec: %s\n", NOW-begintime,
                             p, pointing_array[p][0], pointing_array[p][1]);
    }

    // Allocate memory
    char **filenames = create_filenames( &opts );
    ComplexDouble  ****complex_weights_array = create_complex_weights( npointing, nstation, nchan, npol );
    ComplexDouble  ****invJi                 = create_invJi( nstation, nchan, npol );
    ComplexDouble  ****detected_beam         = create_detected_beam( npointing, 2*opts.sample_rate, nchan, npol );

    // Read in info from metafits file
    fprintf( stderr, "[%f]  Reading in metafits file information from %s\n", NOW-begintime, opts.metafits);
    struct metafits_info mi;
    get_metafits_info( opts.metafits, &mi, opts.chan_width );

    // If using bandpass calibration solutions, calculate number of expected bandpass channels
    if (opts.cal.cal_type == RTS_BANDPASS)
        opts.cal.nchan = (nchan * opts.chan_width) / opts.cal.chan_width;

    // If a custom flag file has been provided, use that instead of the metafits flags
    int i;
    if (opts.custom_flags != NULL)
    {
        // Reset the weights to 1
        for (i = 0; i < nstation*npol; i++)
            mi.weights_array[i] = 1.0;

        // Open custom flag file for reading
        FILE *flagfile = fopen( opts.custom_flags, "r" );
        if (flagfile == NULL)
        {
            fprintf( stderr, "error: couldn't open flag file \"%s\" for "
                             "reading\n", opts.custom_flags );
            exit(EXIT_FAILURE);
        }

        // Read in flags from file
        int nitems;
        int flag, ant;
        while (!feof(flagfile))
        {
            // Read in next item
            nitems = fscanf( flagfile, "%d", &ant );
            if (nitems != 1 && !feof(flagfile))
            {
                fprintf( stderr, "error: couldn't parse flag file \"%s\"\n",
                        opts.custom_flags );
                exit(EXIT_FAILURE);
            }

            // Flag both polarisations of the antenna in question
            flag = ant*2;
            mi.weights_array[flag]   = 0.0;
            mi.weights_array[flag+1] = 0.0;
        }

        // Close file
        fclose( flagfile );
    }

    // Issue warnings if any antennas are being used which are flagged in the metafits file
    for (i = 0; i < nstation*npol; i++)
    {
        if (mi.weights_array[i] != 0.0 &&
            mi.flag_array[i]    != 0.0)
        {
            fprintf( stderr, "warning: antenna %3d, pol %d is included even "
                             "though it is flagged in the metafits file\n",
                             i / npol,
                             i % npol );
        }
    }
    
    double wgt_sum = 0;
    for (i = 0; i < nstation*npol; i++)
        wgt_sum += mi.weights_array[i];
    double invw = 1.0/wgt_sum;

    // Run get_delays to populate the delay_vals struct
    fprintf( stderr, "[%f]  Setting up output header information\n", NOW-begintime);
    struct delays delay_vals[npointing];
    get_delays(
            pointing_array,     // an array of pointings [pointing][ra/dec][characters]
            npointing,          // number of pointings
            opts.frequency,     // middle of the first frequency channel in Hz
            &opts.cal,          // struct holding info about calibration
            opts.sample_rate,   // = 10000 samples per sec
            opts.beam_model,    // name of beam model file
            opts.time_utc,      // utc time string
            0.0,                // seconds offset from time_utc at which to calculate delays
            delay_vals,        // Populate psrfits header info
            &mi,                // Struct containing info from metafits file
            NULL,               // complex weights array (ignore this time)
            NULL                // invJi array           (ignore this time)
    );

    // Create structures for holding header information
    struct psrfits  *pf;
    struct psrfits  *pf_incoh;
    pf = (struct psrfits *)malloc(npointing * sizeof(struct psrfits));
    pf_incoh = (struct psrfits *)malloc(1 * sizeof(struct psrfits));
    vdif_header     vhdr;
    struct vdifinfo *vf;
    vf = (struct vdifinfo *)malloc(npointing * sizeof(struct vdifinfo));


    // Create structures for the PFB filter coefficients
    int ntaps, fil_size = 0;
    double *coeffs = NULL;

    // If no synthesis filter was explicitly chosen, choose the LSQ12 filter
    if (!opts.synth_filter)  opts.synth_filter = strdup("LSQ12");
    if (strcmp( opts.synth_filter, "LSQ12" ) == 0)
    {
        ntaps = 12;
        fil_size = ntaps * nchan; // = 12 * 128 = 1536
        coeffs = (double *)malloc( fil_size * sizeof(double) );
        double tmp_coeffs[] = LSQ12_FILTER_COEFFS; // I'll have to change the way these coefficients are stored
                                                   // in order to avoid this cumbersome loading procedure
        for (i = 0; i < fil_size; i++)
            coeffs[i] = tmp_coeffs[i];
    }
    else if (strcmp( opts.synth_filter, "MIRROR" ) == 0)
    {
        ntaps = 12;
        fil_size = ntaps * nchan; // = 12 * 128 = 1536
        coeffs = (double *)malloc( fil_size * sizeof(double) );
        double tmp_coeffs[] = MIRROR_FILTER_COEFFS;
        for (i = 0; i < fil_size; i++)
            coeffs[i] = tmp_coeffs[i];
    }
    else
    {
        fprintf( stderr, "error: unrecognised synthesis filter: %s\n",
                opts.synth_filter );
        exit(EXIT_FAILURE);
    }
    ComplexDouble *twiddles = roots_of_unity( nchan );

    // Adjust by the scaling that was introduced by the forward PFB,
    // along with any other scaling that I, Lord and Master of the inverse
    // PFB, feels is appropriate.
    double approx_filter_scale = 15.0/7.2; // 7.2 = 16384/117964.8
    for (i = 0; i < fil_size; i++)
        coeffs[i] *= approx_filter_scale;

    // Populate the relevant header structs
    populate_psrfits_header( pf,       opts.metafits, opts.obsid,
            opts.time_utc, opts.sample_rate, opts.frequency, nchan,
            opts.chan_width,outpol_coh, opts.rec_channel, delay_vals,
            mi, npointing, 1 );
    populate_psrfits_header( pf_incoh, opts.metafits, opts.obsid,
            opts.time_utc, opts.sample_rate, opts.frequency, nchan,
            opts.chan_width, outpol_incoh, opts.rec_channel, delay_vals,
            mi, 1, 0 );

    populate_vdif_header( vf, &vhdr, opts.metafits, opts.obsid,
            opts.time_utc, opts.sample_rate, opts.frequency, nchan,
            opts.chan_width, opts.rec_channel, delay_vals, npointing );

    // To run asynchronously we require two memory allocations for each data 
    // set so multiple parts of the memory can be worked on at once.
    // We control this by changing the pointer to alternate between
    // the two memory allocations
    
    // Create array for holding the raw data
    int bytes_per_file = opts.sample_rate * nstation * npol * nchan;

    //cudaMallocHost( (void**)&data, bytes_per_file * sizeof(uint8_t) );
    uint8_t *data = (uint8_t *)malloc( bytes_per_file * sizeof(uint8_t) );
    assert(data);

    // Create output buffer arrays
    float *data_buffer_coh    = NULL;
    float *data_buffer_incoh  = NULL;
    float *data_buffer_vdif   = NULL;

    data_buffer_coh   = create_pinned_data_buffer_psrfits( npointing * nchan *
                                                           outpol_coh * pf[0].hdr.nsblk );
    data_buffer_incoh = create_pinned_data_buffer_psrfits( nchan * outpol_incoh *
                                                           pf_incoh[0].hdr.nsblk );
    data_buffer_vdif  = create_pinned_data_buffer_vdif( vf->sizeof_buffer *
                                                        npointing );

    /* Allocate host and device memory for the use of the cu_form_beam function */
    // Declaring pointers to the structs so the memory can be alternated
    struct gpu_formbeam_arrays gf;
    struct gpu_ipfb_arrays gi;
    int nchunk = 10;
    malloc_formbeam( &gf, opts.sample_rate, nstation, nchan, npol, nchunk,
                     outpol_coh, outpol_incoh, npointing, NOW-begintime );

    if (opts.out_vdif)
    {
        malloc_ipfb( &gi, ntaps, opts.sample_rate, nchan, npol, fil_size, npointing );
        cu_load_filter( coeffs, twiddles, &gi, nchan );
    }

    // Set up parrel streams
    cudaStream_t streams[npointing];
    for ( p = 0; p < npointing; p++ )
        cudaStreamCreate(&(streams[p])) ;

    fprintf( stderr, "[%f]  **BEGINNING BEAMFORMING**\n", NOW-begintime);

    int offset;
    unsigned int s;
    int ch, pol;

    // Set up timing for each section
    long read_time[nfiles], delay_time[nfiles], calc_time[nfiles], write_time[nfiles][npointing];
    int file_no;

    for (file_no = 0; file_no < nfiles; file_no++)
    {
        // Read in data from next file
        clock_t start = clock();
        fprintf( stderr, "[%f] [%d/%d] Reading in data from %s \n", NOW-begintime,
                file_no+1, nfiles, filenames[file_no]);
        read_data( filenames[file_no], data, bytes_per_file  );
        read_time[file_no] = clock() - start;

        // Get the next second's worth of phases / jones matrices, if needed
        start = clock();
        fprintf( stderr, "[%f] [%d/%d] Calculating delays\n", NOW-begintime,
                                file_no+1, nfiles );
        get_delays(
                pointing_array,     // an array of pointings [pointing][ra/dec][characters]
                npointing,          // number of pointings
                opts.frequency,         // middle of the first frequency channel in Hz
                &opts.cal,              // struct holding info about calibration
                opts.sample_rate,       // = 10000 samples per sec
                opts.beam_model,        // name of beam model file
                opts.time_utc,          // utc time string
                (double)file_no,        // seconds offset from time_utc at which to calculate delays
                NULL,                   // Don't update delay_vals
                &mi,                    // Struct containing info from metafits file
                complex_weights_array,  // complex weights array (answer will be output here)
                invJi );                // invJi array           (answer will be output here)
        delay_time[file_no] = clock() - start;


        fprintf( stderr, "[%f] [%d/%d] Calculating beam\n", NOW-begintime,
                                file_no+1, nfiles);
        start = clock();

        if (!opts.out_bf) // don't beamform, but only procoess one ant/pol combination
        {
            // Populate the detected_beam, data_buffer_coh, and data_buffer_incoh arrays
            // detected_beam = [2*opts.sample_rate][nchan][npol] = [20000][128][2] (ComplexDouble)
            if (file_no % 2 == 0)
                offset = 0;
            else
                offset = opts.sample_rate;

            for (p   = 0; p   < npointing;        p++  )
            for (s   = 0; s   < opts.sample_rate; s++  )
            for (ch  = 0; ch  < nchan           ; ch++ )
            for (pol = 0; pol < npol            ; pol++)
            {
                detected_beam[p][s+offset][ch][pol] = UCMPLX4_TO_CMPLX_FLT(data[D_IDX(s,ch,opts.out_ant,pol,nchan)]);
                detected_beam[p][s+offset][ch][pol] = CMuld(detected_beam[p][s+offset][ch][pol], CMaked(128.0, 0.0));
            }
        }
        else // beamform (the default mode)
        {
            cu_form_beam( data, &opts, complex_weights_array, invJi, file_no,
                    npointing, nstation, nchan, npol, outpol_coh, invw, &gf,
                    detected_beam, data_buffer_coh, data_buffer_incoh,
                    streams, opts.out_incoh, nchunk );
        }

        // Invert the PFB, if requested
        if (opts.out_vdif)
        {
            fprintf( stderr, "[%f] [%d/%d] Inverting the PFB (full)\n",
                            NOW-begintime, file_no+1, nfiles);
            cu_invert_pfb_ord( detected_beam, file_no, npointing,
                               opts.sample_rate, nchan, npol, vf->sizeof_buffer,
                               &gi, data_buffer_vdif );
        }
        calc_time[file_no] = clock() - start;


        // Write out for each pointing
        for ( p = 0; p < npointing; p++)
        {
            start = clock();
            fprintf( stderr, "[%f] [%d/%d] [%d/%d] Writing data to file(s)\n",
                    NOW-begintime, file_no+1, nfiles, p+1, npointing );

            if (opts.out_coh)
                psrfits_write_second( &pf[p], data_buffer_coh, nchan,
                                      outpol_coh, p );
            if (opts.out_incoh && p == 0)
                psrfits_write_second( &pf_incoh[p], data_buffer_incoh,
                                      nchan, outpol_incoh, p );
            if (opts.out_vdif)
                vdif_write_second( &vf[p], &vhdr,
                                   data_buffer_vdif + p * vf->sizeof_buffer,
                                   &vgain );
            write_time[file_no][p] = clock() - start;
        }
    }
    
    // Calculate total processing times
    float read_sum = 0, delay_sum = 0, calc_sum = 0, write_sum = 0;
    for (file_no = 0; file_no < nfiles; file_no++)
    {
        read_sum  += (float) read_time[file_no];
        delay_sum += (float) delay_time[file_no];
        calc_sum  += (float) calc_time[file_no];
        for ( p = 0; p < npointing; p++)
            write_sum += (float) write_time[file_no][p];
    }
    float read_mean, delay_mean, calc_mean, write_mean;
    read_mean  = read_sum  / nfiles;
    delay_mean = delay_sum / nfiles / npointing;
    calc_mean  = calc_sum  / nfiles / npointing;
    write_mean = write_sum / nfiles / npointing;

    // Calculate the standard deviations
    float read_std = 0, delay_std = 0, calc_std = 0, write_std = 0;
    for (file_no = 0; file_no < nfiles; file_no++)
    {
        read_std  += pow((float)read_time[file_no]  - read_mean,  2);
        delay_std += pow((float)delay_time[file_no] - delay_mean / npointing, 2);
        calc_std  += pow((float)calc_time[file_no]  - calc_mean / npointing,  2);
        for ( p = 0; p < npointing; p++)
            write_std += pow((float)write_time[file_no][p] - write_mean / npointing, 2);
    }
    read_std  = sqrt( read_std  / nfiles );
    delay_std = sqrt( delay_std / nfiles / npointing );
    calc_std  = sqrt( calc_std  / nfiles / npointing );
    write_std = sqrt( write_std / nfiles / npointing );


    fprintf( stderr, "[%f]  **FINISHED BEAMFORMING**\n", NOW-begintime);
    fprintf( stderr, "[%f]  Total read  processing time: %9.3f s\n", 
                     NOW-begintime, read_sum / CLOCKS_PER_SEC);
    fprintf( stderr, "[%f]  Mean  read  processing time: %9.3f +\\- %8.3f s\n", 
                     NOW-begintime, read_mean / CLOCKS_PER_SEC, read_std / CLOCKS_PER_SEC);
    fprintf( stderr, "[%f]  Total delay processing time: %9.3f s\n", 
                     NOW-begintime, delay_sum / CLOCKS_PER_SEC);
    fprintf( stderr, "[%f]  Mean  delay processing time: %9.3f +\\- %8.3f s\n", 
                     NOW-begintime, delay_mean / CLOCKS_PER_SEC, delay_std / CLOCKS_PER_SEC);
    fprintf( stderr, "[%f]  Total calc  processing time: %9.3f s\n", 
                     NOW-begintime, calc_sum / CLOCKS_PER_SEC);
    fprintf( stderr, "[%f]  Mean  calc  processing time: %9.3f +\\- %8.3f s\n", 
                     NOW-begintime, calc_mean / CLOCKS_PER_SEC, calc_std / CLOCKS_PER_SEC);
    fprintf( stderr, "[%f]  Total write processing time: %9.3f s\n", 
                     NOW-begintime, write_sum  * npointing / CLOCKS_PER_SEC);
    fprintf( stderr, "[%f]  Mean  write processing time: %9.3f +\\- %8.3f s\n", 
                     NOW-begintime, write_mean / CLOCKS_PER_SEC, write_std / CLOCKS_PER_SEC);
    
    
    fprintf( stderr, "[%f]  Starting clean-up\n", NOW-begintime);

    // Free up memory
    destroy_filenames( filenames, &opts );
    destroy_complex_weights( complex_weights_array, npointing, nstation, nchan );
    destroy_invJi( invJi, nstation, nchan, npol );
    destroy_detected_beam( detected_beam, npointing, 2*opts.sample_rate, nchan );

    free( twiddles );
    free( coeffs );

    destroy_metafits_info( &mi );
    //free( data_buffer_coh    );
    //free( data_buffer_incoh  );
    //free( data_buffer_vdif   );
    cudaFreeHost( data_buffer_coh   );
    cudaFreeHost( data_buffer_incoh );
    cudaFreeHost( data_buffer_vdif  );
    cudaFreeHost( data );

    free( opts.obsid        );
    free( opts.time_utc     );
    free( opts.pointings    );
    free( opts.datadir      );
    free( opts.metafits     );
    free( opts.rec_channel  );
    free( opts.cal.filename );
    free( opts.custom_flags );
    free( opts.synth_filter );
    free( opts.beam_model   );

    if (opts.out_incoh)
    {
        free( pf_incoh[0].sub.data        );
        free( pf_incoh[0].sub.dat_freqs   );
        free( pf_incoh[0].sub.dat_weights );
        free( pf_incoh[0].sub.dat_offsets );
        free( pf_incoh[0].sub.dat_scales  );
    }
    for (p = 0; p < npointing; p++)
    {
        if (opts.out_coh)
        {
            free( pf[p].sub.data        );
            free( pf[p].sub.dat_freqs   );
            free( pf[p].sub.dat_weights );
            free( pf[p].sub.dat_offsets );
            free( pf[p].sub.dat_scales  );
        }
        if (opts.out_vdif)
        {
            free( vf[p].b_scales  );
            free( vf[p].b_offsets );
        }
    }
    free_formbeam( &gf );
    if (opts.out_vdif)
    {
        free_ipfb( &gi );
    }
    #ifndef HAVE_CUDA
    // Clean up FFTW OpenMP
    fftw_cleanup_threads();
    #endif

    return 0;
}


void usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "usage: make_beam [OPTIONS]\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "REQUIRED OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-o, --obsid=GPSTIME       ");
    fprintf(stderr, "Observation ID (GPS seconds).\n");
    fprintf(stderr, "\t-b, --begin=GPSTIME       ");
    fprintf(stderr, "Begin time of observation, in GPS seconds\n");
    fprintf(stderr, "\t-e, --end=GPSTIME         ");
    fprintf(stderr, "End time of observation, in GPS seconds\n");
    fprintf(stderr, "\t-z, --utc-time=UTCTIME    ");
    fprintf(stderr, "The UTC time that corresponds to the GPS time given by the -b\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "option. UTCTIME must have the format: yyyy-mm-ddThh:mm:ss\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-P, --pointings=hh:mm:ss.s_dd:mm:ss.s,hh:mm:ss.s_dd:mm:ss.s...\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "Right ascension and declinations of multiple pointings\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-d, --data-location=PATH  ");
    fprintf(stderr, "PATH is the directory containing the recombined data\n");
    fprintf(stderr, "\t-m, --metafits-file=FILE  ");
    fprintf(stderr, "FILE is the metafits file pertaining to the OBSID given by the\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr,  "-o option\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-H, --beam-model=FILE     ");
    fprintf(stderr, "The hdf5 FILE containing the full_embedded_element_pattern used\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "for calculating the FEE 2016 beam model. If FILE = 'analytic',\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "the original analytic beam is used instead.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-f, --coarse-chan=N       ");
    fprintf(stderr, "Absolute coarse channel number (0-255)\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "OUTPUT OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-i, --incoh                ");
    fprintf(stderr, "Turn on incoherent PSRFITS beam output.                             ");
    fprintf(stderr, "[default: OFF]\n");
    fprintf(stderr, "\t-p, --psrfits              ");
    fprintf(stderr, "Turn on coherent PSRFITS output (will be turned on if none of\n");
    fprintf(stderr, "\t                           ");
    fprintf(stderr, "-i, -p, -u, -v are chosen).                                         ");
    fprintf(stderr, "[default: OFF]\n");
    fprintf(stderr, "\t-v, --vdif                 ");
    fprintf(stderr, "Turn on VDIF output with upsampling                                 ");
    fprintf(stderr, "[default: OFF]\n");
    fprintf(stderr, "\t-s, --summed               ");
    fprintf(stderr, "Turn on summed polarisations of the coherent output (only Stokes I) ");
    fprintf(stderr, "[default: OFF]\n");
    fprintf(stderr, "\t-A, --antpol=ant           ");
    fprintf(stderr, "Do not beamform. Instead, only operate on the specified ant\n");
    fprintf(stderr, "\t                           ");
    fprintf(stderr, "stream (0-127)\n" );
    fprintf(stderr, "\t-S, --synth_filter=filter  ");
    fprintf(stderr, "Apply the named filter during high-time resolution synthesis.    ");
    fprintf(stderr, "[default: LSQ12]\n");
    fprintf(stderr, "\t                           ");
    fprintf(stderr, "filter can be MIRROR or LSQ12.\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "MWA/VCS CONFIGURATION OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-a, --antennas=N          ");
    fprintf(stderr, "The number of antennas in the array. For MWA Phase 2, N=128.     ");
    fprintf(stderr, "[default: 128]\n");
    fprintf(stderr, "\t-n, --num-fine-chans=N    ");
    fprintf(stderr, "The number of fine channels per coarse channel.                  ");
    fprintf(stderr, "[default: 128]\n");
    fprintf(stderr, "\t-w, --fine-chan-width=N   ");
    fprintf(stderr, "The bandwidth of an individual fine channel (Hz).                ");
    fprintf(stderr, "[default: 10000]\n");
    fprintf(stderr, "\t-r, --sample-rate=N       ");
    fprintf(stderr, "The VCS sample rate, in Hz. (The sample rate given in the meta-  ");
    fprintf(stderr, "[default: 10000]\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "metafits file matches the correlator settings at the time of\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "the observation, which is not necessarily the same as that of\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "the VCS. Hence the necessity of this option.)\n");
    fprintf(stderr, "\t-F, --custom-flags=file   ");
    fprintf(stderr, "Flag the antennas listed in file instead of those flagged in the ");
    fprintf(stderr, "[default: none]\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "metafits file given by the -m option.\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "CALIBRATION OPTIONS (RTS)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-J, --dijones-file=PATH   ");
    fprintf(stderr, "The direction-independent Jones matrix file that is output from\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "the RTS. Using this option instructs the beamformer to use the\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "RTS-generated calibration solution. Either -J or -O must be\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "supplied. If both are supplied the one that comes last will\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "override the former.\n");
    fprintf(stderr, "\t-B, --bandpass-file=PATH  ");
    fprintf(stderr, "The bandpass file that is output from the RTS. If this option\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "is given, the RTS calibration solution will be applied to each\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "fine channel. If -J is supplied but -B is not, then the coarse\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "channel solution will be applied to ALL fine channels\n");
    fprintf(stderr, "\t-W, --rts-chan-width      ");
    fprintf(stderr, "RTS calibration channel bandwidth (Hz)                           ");
    fprintf(stderr, "[default: 40000]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "CALIBRATION OPTIONS (OFFRINGA)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-O, --offringa-file=PATH  ");
    fprintf(stderr, "The calibration solution file that is output from the tools\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "made by Andre Offringa. Using this option instructs the beam-\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "former to use the Offringa-style calibration solution. Either\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "-J or -O must be supplied. If both are supplied the one that\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "comes last will override the former.\n");
    fprintf(stderr, "\t-C, --offringa-chan=N     ");
    fprintf(stderr, "The zero-offset position of the coarse channel solution in the   ");
    fprintf(stderr, "[default: 0]\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "calibration file given by the -O option.\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "OTHER OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-h, --help                ");
    fprintf(stderr, "Print this help and exit\n");
    fprintf(stderr, "\t-V, --version             ");
    fprintf(stderr, "Print version number and exit\n");
    fprintf(stderr, "\n");
}



void make_beam_parse_cmdline(
        int argc, char **argv, struct make_beam_opts *opts )
{
    if (argc > 1) {

        int c;
        while (1) {

            static struct option long_options[] = {
                {"obsid",           required_argument, 0, 'o'},
                {"begin",           required_argument, 0, 'b'},
                {"end",             required_argument, 0, 'e'},
                {"incoh",           no_argument,       0, 'i'},
                {"psrfits",         no_argument,       0, 'p'},
                {"vdif",            no_argument,       0, 'v'},
                {"summed",          no_argument,       0, 's'},
                {"synth_filter",    required_argument, 0, 'S'},
                {"beam-model",      required_argument, 0, 'H'},
                {"antpol",          required_argument, 0, 'A'},
                {"utc-time",        required_argument, 0, 'z'},
                {"pointings",       required_argument, 0, 'P'},
                {"data-location",   required_argument, 0, 'd'},
                {"metafits-file",   required_argument, 0, 'm'},
                {"coarse-chan",     required_argument, 0, 'f'},
                {"antennas",        required_argument, 0, 'a'},
                {"num-fine-chans",  required_argument, 0, 'n'},
                {"fine-chan-width", required_argument, 0, 'w'},
                {"sample-rate",     required_argument, 0, 'r'},
                {"custom-flags",    required_argument, 0, 'F'},
                {"dijones-file",    required_argument, 0, 'J'},
                {"bandpass-file",   required_argument, 0, 'B'},
                {"rts-chan-width",  required_argument, 0, 'W'},
                {"offringa-file",   required_argument, 0, 'O'},
                {"offringa-chan",   required_argument, 0, 'C'},
                {"help",            required_argument, 0, 'h'},
                {"version",         required_argument, 0, 'V'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv,
                             "a:A:b:B:C:d:e:f:F:hH:iJ:m:n:o:O:pP:r:sS:vVw:W:z:",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch(c) {

                case 'a':
                    opts->nstation = atoi(optarg);
                    break;
                case 'A':
                    opts->out_bf = 0; // Turn off normal beamforming
                    opts->out_ant = atoi(optarg); // 0-127
                    break;
                case 'b':
                    opts->begin = atol(optarg);
                    break;
                case 'B':
                    opts->cal.bandpass_filename = strdup(optarg);
                    opts->cal.cal_type = RTS_BANDPASS;
                    break;
                case 'C':
                    opts->cal.offr_chan_num = atoi(optarg);
                    break;
                case 'd':
                    opts->datadir = strdup(optarg);
                    break;
                case 'e':
                    opts->end = atol(optarg);
                    break;
                case 'f':
                    opts->rec_channel = strdup(optarg);
                    // The base frequency of the coarse channel in Hz
                    opts->frequency = atoi(optarg) * 1.28e6 - 640e3;
                    break;
                case 'F':
                    opts->custom_flags = strdup(optarg);
                    break;
                case 'h':
                    usage();
                    exit(0);
                    break;
                case 'H':
                    opts->beam_model = strdup(optarg);
                    break;
                case 'i':
                    opts->out_incoh = 1;
                    break;
                case 'J':
                    opts->cal.filename = strdup(optarg);
                    if (opts->cal.cal_type != RTS_BANDPASS)
                        opts->cal.cal_type = RTS;
                    break;
                case 'm':
                    opts->metafits = strdup(optarg);
                    break;
                case 'n':
                    opts->nchan = atoi(optarg);
                    break;
                case 'o':
                    opts->obsid = strdup(optarg);
                    break;
                case 'O':
                    opts->cal.filename = strdup(optarg);
                    opts->cal.cal_type = OFFRINGA;
                    break;
                case 'p':
                    opts->out_coh = 1;
                    break;
                case 'P':
                    opts->pointings = strdup(optarg);
                    break;
                case 'r':
                    opts->sample_rate = atoi(optarg);
                    break;
                case 'S':
                    opts->synth_filter = strdup(optarg);
                    break;
                case 's':
                    opts->out_summed = 1;
                    break;
                case 'v':
                    opts->out_vdif = 1;
                    break;
                case 'V':
                    fprintf( stderr, "MWA Beamformer %s\n", VERSION_BEAMFORMER);
                    exit(0);
                    break;
                case 'w':
                    opts->chan_width = atoi(optarg);
                    break;
                case 'W':
                    opts->cal.chan_width = atoi(optarg);
                    break;
                case 'z':
                    opts->time_utc = strdup(optarg);
                    break;
                default:
                    fprintf(stderr, "error: make_beam_parse_cmdline: "
                                    "unrecognised option '%s'\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
            }
        }
    }
    else {
        usage();
        exit(EXIT_FAILURE);
    }


    // Check that all the required options were supplied
    assert( opts->obsid        != NULL );
    assert( opts->begin        != 0    );
    assert( opts->end          != 0    );
    assert( opts->time_utc     != NULL );
    assert( opts->pointings    != NULL );
    assert( opts->datadir      != NULL );
    assert( opts->metafits     != NULL );
    assert( opts->rec_channel  != NULL );
    assert( opts->cal.cal_type != NO_CALIBRATION );
    assert( opts->beam_model   != NULL );

    // If neither -i, -p, nor -v were chosen, set -p by default
    if ( !opts->out_incoh && !opts->out_coh && !opts->out_vdif )
    {
        opts->out_coh = 1;
    }
}



char **create_filenames( struct make_beam_opts *opts )
{
    // Calculate the number of files
    int nfiles = opts->end - opts->begin + 1;
    if (nfiles <= 0) {
        fprintf( stderr, "Cannot beamform on %d files (between %lu and %lu)\n",
                 nfiles, opts->begin, opts->end);
        exit(EXIT_FAILURE);
    }
    // Allocate memory for the file name list
    char **filenames = NULL;
    filenames = (char **)malloc( nfiles*sizeof(char *) );

    // Allocate memory and write filenames
    int second;
    unsigned long int timestamp;
    for (second = 0; second < nfiles; second++) {
        timestamp = second + opts->begin;
        filenames[second] = (char *)malloc( MAX_COMMAND_LENGTH*sizeof(char) );
        sprintf( filenames[second], "%s/%s_%ld_ch%s.dat",
                 opts->datadir, opts->obsid, timestamp, opts->rec_channel );
    }

    return filenames;
}

void destroy_filenames( char **filenames, struct make_beam_opts *opts )
{
    int nfiles = opts->end - opts->begin + 1;
    int second;
    for (second = 0; second < nfiles; second++)
        free( filenames[second] );
    free( filenames );
}


ComplexDouble ****create_complex_weights( int npointing, int nstation, int nchan, int npol )
// Allocate memory for complex weights matrices
{
    int p, ant, ch; // Loop variables
    ComplexDouble ****array;
    
    array = (ComplexDouble ****)malloc( npointing * sizeof(ComplexDouble ***) );
    
    for (p = 0; p < npointing; p++)
    {
        array[p] = (ComplexDouble ***)malloc( nstation * sizeof(ComplexDouble **) );

        for (ant = 0; ant < nstation; ant++)
        {
            array[p][ant] = (ComplexDouble **)malloc( nchan * sizeof(ComplexDouble *) );

            for (ch = 0; ch < nchan; ch++)
                array[p][ant][ch] = (ComplexDouble *)malloc( npol * sizeof(ComplexDouble) );
        }
    }
    return array;
}


void destroy_complex_weights( ComplexDouble ****array, int npointing, int nstation, int nchan )
{
    int p, ant, ch;
    for (p = 0; p < npointing; p++)
    {
        for (ant = 0; ant < nstation; ant++)
        {
            for (ch = 0; ch < nchan; ch++)
                free( array[p][ant][ch] );

            free( array[p][ant] );
        }
        free( array[p] );
    }
    free( array );
}

ComplexDouble ****create_invJi( int nstation, int nchan, int npol )
// Allocate memory for (inverse) Jones matrices
{
    int ant, pol, ch; // Loop variables
    ComplexDouble ****invJi;
    invJi = (ComplexDouble ****)malloc( nstation * sizeof(ComplexDouble ***) );

    for (ant = 0; ant < nstation; ant++)
    {
        invJi[ant] =(ComplexDouble ***)malloc( nchan * sizeof(ComplexDouble **) );

        for (ch = 0; ch < nchan; ch++)
        {
            invJi[ant][ch] = (ComplexDouble **)malloc( npol * sizeof(ComplexDouble *) );

            for (pol = 0; pol < npol; pol++)
                invJi[ant][ch][pol] = (ComplexDouble *)malloc( npol * sizeof(ComplexDouble) );
        }
    }
    return invJi;
}


void destroy_invJi( ComplexDouble ****array, int nstation, int nchan, int npol )
{
    int ant, ch, pol;
    for (ant = 0; ant < nstation; ant++)
    {
        for (ch = 0; ch < nchan; ch++)
        {
            for (pol = 0; pol < npol; pol++)
                free( array[ant][ch][pol] );

            free( array[ant][ch] );
        }
        free( array[ant] );
    }
    free( array );
}


ComplexDouble ****create_detected_beam( int npointing, int nsamples, int nchan, int npol )
// Allocate memory for complex weights matrices
{
    int p, s, ch; // Loop variables
    ComplexDouble ****array;
    
    array = (ComplexDouble ****)malloc( npointing * sizeof(ComplexDouble ***) );
    for (p = 0; p < npointing; p++) 
    {
        array[p] = (ComplexDouble ***)malloc( nsamples * sizeof(ComplexDouble **) );

        for (s = 0; s < nsamples; s++)
        {
            array[p][s] = (ComplexDouble **)malloc( nchan * sizeof(ComplexDouble *) );

            for (ch = 0; ch < nchan; ch++)
                array[p][s][ch] = (ComplexDouble *)malloc( npol * sizeof(ComplexDouble) );
        }
    }
    return array;
}

void destroy_detected_beam( ComplexDouble ****array, int npointing, int nsamples, int nchan )
{
    int p, s, ch;
    for (p = 0; p < npointing; p++)    
    {
        for (s = 0; s < nsamples; s++)
        {
            for (ch = 0; ch < nchan; ch++)
                free( array[p][s][ch] );

            free( array[p][s] );
        }

        free( array[p] );
    }

    free( array );
}

float *create_data_buffer_psrfits( size_t size )
{
    float *ptr = (float *)malloc( size * sizeof(float) );
    return ptr;
}


float *create_data_buffer_vdif( size_t size )
{
    float *ptr  = (float *)malloc( size * sizeof(float) );
    return ptr;
}

