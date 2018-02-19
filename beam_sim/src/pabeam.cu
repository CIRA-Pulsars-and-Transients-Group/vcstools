#include "pabeam.h"


void usage()
{
    printf("pabeam_gpu --- computes the array factor that represents the naturally weighted synthesised MWA beam (tied-array/coherent beam) for a given configuration\n");
    printf("syntax:\n");
    printf("    pabeam -f <frequency in Hz> -r <ra in hh:mm:ss> -d <dec in dd:mm:ss> -t <UTC in ISOT format> -m <metafits file> -b <RTS flagged_tiles.txt file> [-e] [-x] [-y] [-g]\n\n");
    printf("Options:\n");
    printf("    -h this help\n");
    printf("    -f observing frequency, in Hz\n");
    printf("    -r target RA (J2000), in hh:mm:ss.ss format\n");
    printf("    -d target DEC (J2000), in dd:mm:ss.ss format\n");
    printf("    -t UTC time to evaluate, in format YYYY-MM-DDThh:mm:ss\n");
    printf("    -m metafits file for the observation\n");
    printf("    -b RTS flagged_tiles.txt file from calibration\n");
    printf("    -e radiation efficiency (if unsure, use 1.0)\n");
    printf("    -x Azimuth grid resolution element (>= 0.01)\n");
    printf("    -y Zenith angle grid resolution element (>=0.01)\n");
    printf("    -g Calculate and apply the FEE2016 tile beam with the given \"gridpoint\" number\n");
}


int getNumTiles(const char *metafits)
{
    /* Figure out the number of tiles based on the information in the metafits.

       NOTE: we get warnings from this in compilation because the library functions
       expect char characters, but conversion from string literals to chars is bad.
       It works, but we get warnings... 
     */

    fitsfile *fptr = NULL;
    int status = 0;
    size_t ninput = 0;
    char tiledata[] = "TILEDATA", naxis2[] = "NAXIS2";

    fits_open_file(&fptr, metafits, READONLY, &status); // Open metafits file
    if (status != 0)
    {
        fprintf(stderr, "CRITICAL: Error - Failed to open metafits file for reading\n");
        exit(EXIT_FAILURE);
    }

    fits_movnam_hdu(fptr, BINARY_TBL, tiledata, 0, &status); // Move to TILEDATA HDU
    if (status != 0)
    {
        fprintf(stderr, "CRITICAL: Error - Failed to move to TILEDATA HDU\n");
        exit(EXIT_FAILURE);
    }

    fits_read_key(fptr, TINT, naxis2, &ninput, NULL, &status); // Read how many tiles are included
    if (status != 0)
    {
        fprintf(stderr, "CRITICAL: Error - Failed to read size of binary table in TILEDATA\n");
        exit(EXIT_FAILURE);
    }

    fits_close_file(fptr, &status);
    if (status != 0)
    {
        fprintf(stderr, "CRITICAL: Error - Failed to close metafits file\n");
        exit(EXIT_FAILURE);
    }


    return ninput;
}


void getTilePositions(const char *metafits, int ninput,\
        float *n_pols, float *e_pols, float *h_pols,\
        float *n_tile, float *e_tile, float *h_tile)
{
    /* Get the tile positions from the metafits file.
       Accepts the metafits file name, 
       number of items to read (i.e. 2x number of tiles, 1 per polarisation),
       the array to fill with the polarisation locations, and
       the array to fill with the tile locations (every second element of *_pols).
       
       NOTE: we get warnings from this in compilation because the library functions
       expect char characters, but conversion from string literals to chars is bad.
       It works, but we get warnings... */

    fitsfile *fptr = NULL;
    int status = 0, anynull = 0;
    int colnum = 0;
    int i; // Loop counter
    char north[] = "North", east[] = "East", height[] = "Height", tiledata[] = "TILEDATA";


    fits_open_file(&fptr, metafits, READONLY, &status); // Open metafits file
    if (status != 0)
    {
        fprintf(stderr, "CRITICAL: Error - Failed to open metafits file for reading\n");
        exit(EXIT_FAILURE);
    }


    fits_movnam_hdu(fptr, BINARY_TBL, tiledata, 0, &status); // Move to TILEDATA HDU
    if (status != 0) 
    {
        fprintf(stderr, "CRITICAL: Error - Failed to move to TILEDATA HDU\n");
        exit(EXIT_FAILURE);
    }

    fits_get_colnum(fptr, 1, north, &colnum, &status); // Get north coordinates of tiles
    if (status != 0)
    {
        fprintf(stderr, "CRITICAL: Error - Failed to locate tile North coordinates column in metafits file\n");
        exit(EXIT_FAILURE);
    }

    fits_read_col_flt(fptr, colnum, 1, 1, ninput, 0.0, n_pols, &anynull, &status);
    if (status != 0)
    {
        fprintf(stderr, "CRITICAL: Error - Failed to read North coordinates from metafits file\n");
        exit(EXIT_FAILURE);
    }

    fits_get_colnum(fptr, 1, east, &colnum, &status); // Get east coordinates of tiles
    if (status != 0)
    {
        fprintf(stderr, "CRITICAL: Error - Failed to locate tile East coordinates column in metafits file\n");
        exit(EXIT_FAILURE);
    }

   fits_read_col_flt(fptr, colnum, 1, 1, ninput, 0.0, e_pols, &anynull, &status);
    if (status != 0)
    {
        fprintf(stderr, "CRITICAL: Error - Failed to read East coordinates in metafits file\n");
        exit(EXIT_FAILURE);
    }

    fits_get_colnum(fptr, 1, height, &colnum, &status); // Get height a.s.l. of tiles
    if (status != 0)
    {
        fprintf(stderr, "CRITICAL: Error - Failed to locate tile Height coordinates column in metafits file\n");
        exit(EXIT_FAILURE);
    }

    fits_read_col_flt(fptr, colnum, 1, 1, ninput, 0.0, h_pols, &anynull, &status);
    if (status != 0)
    {
        fprintf(stderr, "CRITICAL: Error - Failed to read H coord in metafile\n");
        exit(EXIT_FAILURE);
    }
    fits_close_file(fptr, &status);

    // Populate the tile arrays with every second element of the polarisation arrays
    for (i = 0; i < ninput; i+=2)
    {
        n_tile[i/2] = n_pols[i];
        e_tile[i/2] = e_pols[i];
        h_tile[i/2] = h_pols[i];
    }

    // Convert heights into height above array center
    for (i = 0; i < ninput/2; i++)
    {
        h_tile[i] -= MWA_HGT;
    }
}


int getFlaggedTiles(const char *badfile, int *badtiles)
{
    /* Open the flagged tiles file, read into an array and count how many lines are read.
       Update the array pointer and return number of elements to read from that array
       (as it's initialised to be able to hold every tile) 
     */
    
    FILE *fp;
    int tile = 0;
    int i = 0; // Loop counter

    fp = fopen(badfile, "r");

    if (fp == NULL)
    {
        fprintf(stderr, "CRITICAL: Error - Failed to open tile flagging file %s\n", badfile);
        exit(EXIT_FAILURE);
    }

    while(fscanf(fp, "%d\n", &tile) > 0)
    {
        printf("    bad tile: %d\n",tile);
        
        badtiles[i] = tile;
        
        i++;
    }

    fclose(fp);

    return i;
}


void removeFlaggedTiles(float *n_tile, float *e_tile, float *h_tile,
                        int *badtiles, int nbad, int nelements)
{
    /* Get rid of the bad/flagged tiles from the array. Basically just
       shift all the indexes around so that whenever a bad tile is there,
       it's data is over-written. We end up with an array of the same size,
       but the last nbad elements are all identical (and can be ignored). 
     */

    int counter = 0, bidx = 0;
    int b, i; // Loop counter

    for (b = 0; b < nbad; b++)
    {
        // For each bad tile index in badtiles
        bidx = badtiles[b];
        for (i = (bidx-counter); i < nelements-1; i++)
        {
            // Shift each element in tile positions to the left by one
            // excluding the last element
            n_tile[i] = n_tile[i+1];
            e_tile[i] = e_tile[i+1];
            h_tile[i] = h_tile[i+1];
        }
        // Array is now shifted left one place, but the bad indexes refer to original tile positions
        // so we need to move the bad index to the left by one, too
        counter++;
    }
}


int main(int argc, char *argv[])
{
    char *ra = NULL;
    char *dec = NULL;
    char *time = NULL;
    char *metafits = NULL;
    char *flagfile = NULL;
    int c = 0; // Arguments counter
    int i, j, k; // Loop counters
    double freq = 0.0, lambda = 0.0;
    double az_step = 1.0, za_step = 1.0, eta = 1.0;
    int nThreads = 1, nBlocks = 1;
    int use_tile_beam = 0, gridpoint = 0;

    // Parse options
    if (argc > 1)
    {
        while ((c = getopt(argc, argv, "h:f:e:r:d:t:m:b:x:y:g:")) != -1)
        {
            switch(c)
            {
                case 'h':
                    usage();
                    exit(EXIT_SUCCESS);
                case 'f':
                    freq   = atof(optarg);
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
                    if (az_step < 0.01)
                    {
                        fprintf(stderr, "WARNING: Using smaller that 0.01 degrees hasn't been tested...\n");
                        fprintf(stderr, "         Overriding to 0.01 degrees\n");
                        //usage();
                        //exit(EXIT_FAILURE);
                        az_step = 0.01;
                    }
                    break;
                case 'y':
                    za_step = atof(optarg);
                    if (az_step < 0.01)
                    {
                        fprintf(stderr, "WARNING: Using smaller that 0.01 degrees hasn't been tested...\n");
                        fprintf(stderr, "         Overriding to 0.01 degrees\n");
                        //usage();
                        //exit(EXIT_FAILURE);
                        za_step = 0.01;
                    }
                    break;
                case 'g':
                    use_tile_beam = 1;
                    gridpoint = atoi(optarg);
                    break;
                default:
                    usage();
                    exit(EXIT_SUCCESS);
            }
        }
    }

    if (argc == 1)
    {
        usage();
        exit(EXIT_SUCCESS);
    }
    
    // Let user know that using the FEE2016 tile beam model will slow down the simulation
    if (use_tile_beam == 1)
    {
        printf("Using FEE2016 tile beam model\n");
        printf("This will slow down the computation significantly...\n");
        printf("    grid point number provided: %d\n", gridpoint);
    }
    else
    {
        printf("Not using FEE2016 tile beam model\n");
        printf("Only computing array factor power...\n");
    }

    // Calculate target az,za and wavevector
    tazza    target;
    wavenums target_wn;

    printf("Getting target (Az,ZA)\n");
    calcTargetAZZA(ra, dec, time, &target);
    printf("Computing wavenumbers towards target\n");
    calcWaveNumber(lambda, target.az, target.za, &target_wn);
    printf("    kx = %f    ky = %f    kz = %f\n", target_wn.kx, target_wn.ky, target_wn.kz); 

    // Get the number of tiles in array
    int ntiles = 0;

    printf("Determining number of tiles from metafits file\n");
    ntiles = getNumTiles(metafits); // returns 2x the number of tiles, 1 per pol.
    ntiles = ntiles / 2;
    printf("    number of tiles: %d\n",ntiles);

    // Allocate dynamic memory for intermediate tile position arrays
    // (probably don't need to check as this should be ~MB scale)
    float *N_pols = (float *)malloc(2 * ntiles * sizeof(float));
    float *E_pols = (float *)malloc(2 * ntiles * sizeof(float));
    float *H_pols = (float *)malloc(2 * ntiles * sizeof(float));
    // Allocate dynamic memory for tile positions
    float *N_tile = (float *)malloc(ntiles * sizeof(float));
    float *E_tile = (float *)malloc(ntiles * sizeof(float));
    float *H_tile = (float *)malloc(ntiles * sizeof(float));
    printf("Getting tile positions\n");
    getTilePositions(metafits, 2*ntiles,
                     N_pols, E_pols, H_pols, 
                     N_tile, E_tile, H_tile);
    free(N_pols);
    free(E_pols);
    free(H_pols);

    // We have to remove tiles from the flagged tiles list.
    // Each element in the list is the index of a tile that needs to be removed.
    printf("Getting flagged tiles\n");
    int *flagged_tiles = (int *)malloc(ntiles * sizeof(int));
    int ntoread;

    ntoread = getFlaggedTiles(flagfile, flagged_tiles);
    
    int flagged[ntoread];

    for (i = 0; i < ntoread; i++)
    {
        flagged[i] = flagged_tiles[i];
    }
    
    printf("Removing %d flagged tiles\n", ntoread);
    printf("Tiles remaining: %d\n",       ntiles-ntoread);
    
    removeFlaggedTiles(N_tile, E_tile, H_tile, flagged, ntoread, ntiles);
    free(flagged_tiles);
    printf("\n");

    // But, the last ntoread elements are pointless,
    // so now we can allocate static memory for the final list of positions
    ntiles -= ntoread;
    float xpos[ntiles], ypos[ntiles], zpos[ntiles];

    for (i = 0; i < ntiles; i++)
    {
        // x = East, y = North, z = Height
        xpos[i] = E_tile[i];
        ypos[i] = N_tile[i];
        zpos[i] = H_tile[i];
    }
    free(N_tile);
    free(E_tile);
    free(H_tile);


    // Determine number of az/za elements from specified pixel size
    int niter = 1;
    int n_az, n_za;
    long int size;

    n_az = (int)(360.0/az_step);
    n_za = (int)(90.0/za_step)+1; // +1 because we want to evalute at 90deg too!
    size = n_az * n_za;
    printf("Number of az steps [0,360): %d\n",           n_az); // step from [0,360) - 360 will double count the 0 values
    printf("Number of za steps [0,90] : %d\n",           n_za); // step from [0,90]
    printf("Total number of elements to compute: %ld\n", size);
    niter = 1; // How many times do I need to split the problem up?
    printf("\n");
 
    // Figure out how many iterations are required (being conservative)
    // and the device properties (as a consequence)
    requiredMemory(size, ntiles, &niter, &nThreads);


    /* We now have the relevant array configuration and target source information 
       needed to calculate the array factor. The best way is to split it up into 
       managable chunks (depending on the device capabilities). 
     */

    // Construct arrays for computation on host
    double *az_array, *za_array;
    
    // Allocate memory and initialise to zeros on host and check
    az_array = (double *)calloc(size, sizeof(double)); // azimuth vector
    if (!az_array)
    {
        fprintf(stderr, "CRITCAL: Error - Host memory allocation failed (allocate az_array)\n");
        return EXIT_FAILURE;
    }

    za_array = (double *)calloc(size, sizeof(double)); // zenith vector
    if (!za_array)
    {
        fprintf(stderr, "CRITICAL: Error - Host memory allocation failed (allocate za_array)\n");
        return EXIT_FAILURE;
    }

    // Populate the host vectors:
    // TODO: this is currently the most memory intensive part on host.
    //       maybe we want to move this initilisation part into the iteration loop
    //       which will then make the arrays smaller --
    //           need to figure out how to populate correctly then...
    printf("Initialising (az, za) arrays\n");
    // We want arrays to be something like:
    // az = [0 0 0 0 ... 1 1 1 1 ...]
    // za = [0 1 2 3 ... 0 1 2 3 ...]
    
    k = 0; // Will be incremented by iteration size
    i = 0;
    do
    {
        for (j = 0; j < n_za; j++)
        {
            za_array[k+j] = j * za_step * DEG2RAD;
            az_array[k+j] = i * az_step * DEG2RAD;
        }
        k += n_za;
        i++;
    } while(k < size);
    printf("Done\n");


    // Construct arrays for device computation
    double *d_az_array, *d_za_array, *subAz, *subZA;
    double af_max = -1.0, omega_A = 0.0;
    cuDoubleComplex *af_array, *d_af_array;
    float *d_xpos, *d_ypos, *d_zpos;
    int itersize, az_idx1, az_idx2, za_idx1, za_idx2; 
    int iter_n_az = (int)floor(size / niter);
    int iter_n_za = (int)floor(size / niter);

    // Before we get to the real computation, better open a file ready for writing
    int obsid;
    char output[100];
    sscanf(metafits, "%d%*s", &obsid);
    printf("Will output beam pattern to:\n");
    printf("    %d_%.2fMHz_%s.dat\n", obsid, freq/1.0e6, time);
    sprintf(output, "%d_%.2fMHz_%s.dat", obsid, freq/1.0e6, time);
    
    FILE *fp;
    fp = fopen(output, "w");  // open the file to write

    // This is the primary loop which does the calculations
    printf("%d az , %d za per iteration\n", iter_n_az, iter_n_za);
    for (int iter = 0; iter < niter; iter++)
    {  
        printf("==== Iteration %d ====\n", iter);
        // Figure out this iteration size, then allocate memory
        if (iter != niter-1)
        {
            itersize = iter_n_az; // = iter_n_za
            az_idx1  = iter * iter_n_az;
            az_idx2  = (iter+1) * iter_n_az;
            za_idx1  = az_idx1;
            //za_idx1  = iter * iter_n_za;
            za_idx2  = az_idx2;
            //za_idx2  = (iter+1) * iter_n_za;
        }
        else
        {
            /* If this is the last iteration, construct 
               iter_n_az/za such that it includes what ever
               is left over to compute.
               
               Should be ok in terms of memory because we made
               the number of iterations was computed on a 
               conservative estimate of the device memory. 
             */
            
            iter_n_za = size - (iter * iter_n_za);
            iter_n_az = size - (iter * iter_n_az);
            itersize  = iter_n_az; // = iter_n_za

            az_idx1 = iter * iter_n_az;
            az_idx2 = az_idx1 + itersize - 1;

            za_idx1 = iter * iter_n_za;
            za_idx2 = za_idx1 + itersize - 1;
        }

        printf("# az: %d  # za: %d\n", iter_n_az, iter_n_za); 
        
        // Allocate interation array(s) memory
        subAz = (double *)malloc(iter_n_az * sizeof(double));
        if (!subAz)
        {
            fprintf(stderr, "CRITICAL: Error - Host memory allocation failed (allocate subAz)\n");
            return EXIT_FAILURE;
        }

        subZA = (double *)malloc(iter_n_za * sizeof(double));
        if (!subZA)
        {
            fprintf(stderr, "CRITICAL: Error - Host memory allocation failed (allocate subZA)\n");
            return EXIT_FAILURE;
        }

        af_array = (cuDoubleComplex *)malloc(iter_n_az * sizeof(cuDoubleComplex));
        if (!af_array)
        {
            fprintf(stderr, "CRITICAL: Error - Host memory allocation failed (allocate af_array)\n");
            return EXIT_FAILURE;
        }

        // Number of blocks required 
        nBlocks = (itersize + nThreads - 1) / nThreads; 
        

        printf("azimuth idx: %d - %d\n", az_idx1, az_idx2);
        printf("zenith  idx: %d - %d\n", za_idx1, za_idx2);
        printf("Number of GPU blocks and threads used:\n");
        printf("    # blocks = %d\n",             nBlocks);
        printf("    # threads = %d\n",           nThreads);

        
        // Place subset of az/za array into subAz/subZA
        for (i = 0; i < itersize; i++)
        {
            subAz[i]    = az_array[i+az_idx1];
            subZA[i]    = za_array[i+za_idx1];
            af_array[i] = make_cuDoubleComplex(0,0);
        }

        // Allocate memory on device
        gpuErrchk( cudaMalloc( (void **)&d_az_array, itersize * sizeof(*az_array) ));
        gpuErrchk( cudaMalloc( (void **)&d_za_array, itersize * sizeof(*za_array) ));
        gpuErrchk( cudaMalloc( (void **)&d_xpos,     ntiles   * sizeof(*xpos)     ));
        gpuErrchk( cudaMalloc( (void **)&d_ypos,     ntiles   * sizeof(*ypos)     ));
        gpuErrchk( cudaMalloc( (void **)&d_zpos,     ntiles   * sizeof(*zpos)     ));
        gpuErrchk( cudaMalloc( (void **)&d_af_array, itersize * sizeof(*af_array) ));


        // Copy arrays onto device
        gpuErrchk( cudaMemcpy( d_az_array, subAz, itersize * sizeof(*subAz), cudaMemcpyHostToDevice ));
        gpuErrchk( cudaMemcpy( d_za_array, subZA, itersize * sizeof(*subZA), cudaMemcpyHostToDevice ));

        // Copy the array factor vector to device
        gpuErrchk( cudaMemcpy( d_af_array, af_array, itersize * sizeof(*af_array), cudaMemcpyHostToDevice ));
        
        // Copy over tile position arrays to device
        gpuErrchk( cudaMemcpy( d_xpos, xpos, ntiles * sizeof(*xpos), cudaMemcpyHostToDevice ));
        gpuErrchk( cudaMemcpy( d_ypos, ypos, ntiles * sizeof(*ypos), cudaMemcpyHostToDevice ));
        gpuErrchk( cudaMemcpy( d_zpos, zpos, ntiles * sizeof(*zpos), cudaMemcpyHostToDevice ));

        printf("Launching kernal to compute array factor\n");
        calcArrayFactor<<<nBlocks, nThreads>>>(itersize, ntiles, 2*PI/lambda, 
                                               d_za_array, d_az_array, 
                                               d_xpos, d_ypos, d_zpos, 
                                               target_wn.kx, target_wn.ky, target_wn.kz, 
                                               d_af_array);
        cudaDeviceSynchronize();

        // Copy relevant memory back to host
        gpuErrchk( cudaMemcpy(af_array, d_af_array, itersize * sizeof(*af_array), cudaMemcpyDeviceToHost));
        printf("==== Done ====\n");

        // Cool, we're done with the GPU computation
        printf("Freeing device memory\n");
        gpuErrchk( cudaFree(d_xpos));
        gpuErrchk( cudaFree(d_ypos));
        gpuErrchk( cudaFree(d_zpos));
        gpuErrchk( cudaFree(d_af_array));
        gpuErrchk( cudaFree(d_az_array));
        gpuErrchk( cudaFree(d_za_array));

        // Write the output to a file
        double af_power = 0.0, tile_power = 1.0, total_power = 1.0; 
      
        printf("Done.\n");
        printf("Writing to file...\n");
        for (i = 0; i < itersize; i++)
        {
            af_power = pow(cuCabs(af_array[i]), 2); // need to use cuCabs given af_array is of cuComplexDouble type
            
            // Compute the tile beam power
            if (use_tile_beam == 1)
            {
                tile_power = CalcMWABeam(subAz[i]-PI/2, subZA[i], freq, 'X', gridpoint, 1);
            }
            else
            {
                tile_power = 1.0;
            }
            
            total_power = af_power * tile_power;

            //if (i % 10000 == 0) {printf("\rWriting element %d/%d", i, itersize); fflush(stdout);}

            fprintf(fp, "%f\t%f\t%f\n", subAz[i]*RAD2DEG, subZA[i]*RAD2DEG, total_power);
            
            // Check what the maximum array factor power is (should be ~1)
            if (af_power > af_max)
            {
                af_max = af_power;
            }
            
            // Integrate over sky
            omega_A += sin(subZA[i]) * af_power * (za_step*DEG2RAD) * (az_step*DEG2RAD);
        } // End power evaluation and file writing
        
        printf("\nDone -- freeing intermediate host memory\n");
        free(subAz);
        free(subZA);
        free(af_array);
    }

    printf("\n");
    printf("Closing file\n");
    fclose(fp); // Close the file
    
    printf("Freeing host memory\n");
    free(az_array);
    free(za_array);

    // Compute the gain and effective area from simulation
    double eff_area = 0.0, gain = 0.0;
    printf("Finished -- now computing relevant parameters:\n");
    eff_area = eta * pow(lambda, 2) * (4 * PI / omega_A);
    gain = (1.0e-26) * eff_area / (2 * KB);

    printf("    Array factor max:                 %f\n",   af_max);
    printf("    Beam solid angle (sr):            %f\n",   omega_A);
    printf("    Radiation efficiency:             %f\n",   eta);
    printf("    Effective collecting area (m^2):  %.4f\n", eff_area);
    printf("    Effective array gain (K/Jy):      %.4f\n", gain);

    return 0;
}
