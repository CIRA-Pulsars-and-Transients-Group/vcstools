/* --------------------------- header secton ----------------------------*/
#include <iostream>     /* yes we are moving to C++ */
#include <pthread.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <sys/socket.h> /* for socket(), bind(), and connect() */
#include <sys/un.h>
#include <sys/shm.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <arpa/inet.h>  /* for sockaddr_in and inet_ntoa() */
#include <errno.h>
#include "acquire_data.h"
#include "ringbuffer.h"
#include "buffer_sizes.h"
#include "packet.h"
#include <complex.h>
#include <unistd.h>
#include "fitsio.h"
#include "fourbit.h"
#include <fcntl.h>

/* -------------------------------- Correlator ------------------------------ */
#include "xgpu.h"
#include "cuda.h"
#include "cuda_runtime.h"
/* -------------------------------- Correlator Utils ------------------------------ */
#include "corr_utils.h"
/* -------------------------------- Beamformer ------------------------------ */
#include "run_beamer.h"
/* -------------------------------- End Beamformer ------------------------------ */

#define checkCudaError() do {                           \
cudaError_t error = cudaGetLastError();             \
if (error != cudaSuccess) {                         \
fprintf(stderr, "(CUDA) %s", cudaGetErrorString(error));  \
fprintf(stderr, " (" __FILE__ ":%d)\n", __LINE__);                \
return XGPU_CUDA_ERROR;                                           \
}                                                   \
} while (0)

/*
 
 IMPORTANT
 
 Data ordering for input vectors is (running from slowest to fastest)
 [time][channel][station][polarization][complexity]
 
 Output matrix has ordering
 [channel][station][station][polarization][polarization][complexity]
 
 We there is a wrinkle in that the station order out of the PFB and thus
 the correlator is *not* the same as the input order.
 
 Please see antenna_mapping.h
 
 
 */

ringbuf_t ring;




/*--------------------------------------------------------------------------*/

void usage()
{
    std::cout << "offline_correlator: a light-weight correlator for the MWA. "
              << "Takes a NCHAN of data from stdin and correlates as per the "
              << "parameters of the linked xGPU library" << std::endl;
    std::cout << "offline_correlator from VCS Tools v" << VERSION_BEAMFORMER
              << std::endl << std::endl;
    std::cout << "offline_correlator: -o <obsid> -s <time_t> -f nchan"
              << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -r <dump_rate> how many correlator dumps per second [1]"
              << std::endl;
    std::cout << " -n <number of channels to average> how many adjacent "
              << "channels to average " << std::endl;
    std::cout << " -i <number of correlator dumps to average> how many "
              << "correlator dumps to average " << std::endl;
    std::cout << "It will take data from stdin. In this case you need to "
              << "give it the start second of the dataset and the associated "
              << "obsid." << std::endl;
}

void *manager(void *context) {
    
    // this just manages the data ingest
    
    
    volatile manager_t *config = (manager_t *) context;
    char *raw_buffer = (char *) malloc(config->ring->bufsize);
    FILE *input = stdin;

    if (raw_buffer == NULL) {
        std::cerr << "error: raw data buffer on start" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    std::cout << "Building lookup... ";
    build_eight_bit_lookup();
    std::cout << "Ready" << std::endl;;
   
    if (config->infile)
    {
        input= fdopen( config->infile, "r" );
    }
 
    while(1)
    {
        // we are getting data from stdin
        // lets just fill up a buffer
        
        int ninputs = config->nstation * config->npol;
        int nchan = config->nfrequency;
        int ntime = config->ntime;
        int ndim = config->ndim;
        int edge = config->edge;
        int nbit = config->nbit;
        int dumps_per_second = config->dumps_per_sec;
        
        
        char *buf = NULL;
        
        while (buf == NULL) {
            buf = get_buffer_to_fill_sync(config->ring); //gets the current buffer
            
        }
        get_buffer_status(config->ring);
        // now need to fill it with 1 seconds worth of input data
        // which is now a variable size as it depends on the number of edge channels that were removed.
        // becuase it is stored on disk as fourbit numbers
        
        size_t nread = 0;
        size_t to_read = (ntime * (nchan-2*edge) * ninputs * ndim * nbit)/8;
        std::cout << "Attempting to read in "<< to_read << " bytes" << std::endl;
        char *raw_buffer_ptr = raw_buffer;
        
        if (nbit == 4) {
            nread = fread(raw_buffer_ptr,1,to_read,input);
            // check fred return status ... just in case
            if (nread != to_read) {
                std::cerr << "error: incomplete read on STDIN (" << nread << " of " << to_read
                          << "). Likely EOD" << std::endl;
                config->ring->EOD = 1;
                break;
            }
            else  {
                
                /* four to eight bit expansion */
                //
                size_t samps = 0;
                size_t timestep = 0;
                size_t chanstep = 0;
                
                
                int16_t *current_raw_ptr = (int16_t *) raw_buffer;
                int8_t *current_out_ptr = (int8_t *) buf;
                size_t samps_per_chan =ninputs*ndim; //NINPUTS*NDIM
                
                
                
                while (timestep < ntime) { // NCHAN*NINPUTS*NDIM*NTIME
                    
                    chanstep = 0;
                    
                    while (chanstep < nchan) {
                        
                        samps = 0;
                        while (samps < samps_per_chan) {
                            if ((chanstep < edge ) || chanstep >= (nchan-edge)) {
                                
                                current_out_ptr[0] = 0;
                                current_out_ptr[1] = 0;
                                current_out_ptr[2] = 0;
                                current_out_ptr[3] = 0;
                            }
                            else {
                                expand_4bit((uint16_t *) current_raw_ptr, (int8_t *) current_out_ptr);
                                current_raw_ptr++; // move 16 bits or 4 samples
                            }
                            current_out_ptr = current_out_ptr + 4; // mover 4 samples
                            samps=samps+4;
                        }
                        chanstep=chanstep+1;
                        
                    }
                    timestep = timestep+1;
                }
            }
            /* done 4 to 8bit expansion -- and droppped in the edges*/
            
        } else {
            nread = fread(buf,1,to_read,stdin);
            // check fred return status ... just in case
            if (nread != to_read)
            {
                std::cerr << "error: incomplete read on STDIN (" << nread << " of " << to_read
                          << "). Likely EOD" << std::endl;
            }
        }
        mark_buffer_filled(config->ring); // this marks the buffer full
        
    }
    if (config->infile) {
	fclose(input);
    }
    free(raw_buffer);

    return NULL;
}

int main(int argc, char **argv) {
    
    XGPUInfo xgpu_info;
    int xgpu_error = 0;
    
    char *buf = 0x0;

    pthread_t buffer_handler;
    extern int buffer_handler_arg;
    
    char *obsid=NULL;
    char *in_file=NULL;    
    
    /* picked up from the inbound header */
    
    struct tm start_utctime;
    struct tm current_utctime;
    /* picked up from the commandline */
    
    time_t starttime = -1;
    
    /* constructed from the tm struct */
    
    time_t start_time_t = 0;
    time_t current_time_t = 0;
    time_t incremented_time_t = 0;
    
    timeval clock1,clock2,clock3;
    double elapsed = 0.0;
    
    
    char file_time[128];
    char dump_filename[128];
    unsigned int npol = 2, nstation = 128, nfrequency = 128, ntime = 10000;
    
    int dumps_per_second = 1; //correlator output dumps per second;
    int chan_to_aver = 1; // number of channels to combine on output
    int dumps_to_aver = 1; // number of correlator dumps to combine on output
    
    int offline = 0;
    
    int edge = 0;
    int nbit = 4;
    int coarse_chan = -1; // only set in the header if this is >= 0
    
    int arg = 0;
    
    while ((arg = getopt(argc, argv, "b:c:d:e:f:hi:n:o:r:s:V")) != -1) {
        
        switch (arg) {
            case 'b':
                nbit = atoi(optarg);
                break;
            case 'c':
                coarse_chan = atoi(optarg);
                break;
            case 'd':
                in_file = strdup(optarg);
                break;
            case 'e':
                edge = atoi(optarg);
                break;
            case 'f':
                // number of channels to correlate per coarse
                nfrequency = atoi(optarg);
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
            case 'i':
                // correlator dumps to sum
                dumps_to_aver=atoi(optarg);
                break;
            case 'n':
                // number of channels to sum
                chan_to_aver=atoi(optarg);
                break;
            case 'o':
                offline = 1;
                obsid = strdup(optarg);
                break;
            case 'r':
                // correlator dump rate
                dumps_per_second = atoi(optarg);
                break;
            case 's':
                starttime = (time_t) atol(optarg);
                break;
            case 'V':
                std::cout << "offline_correlator from VCS Tools v"
                          << VERSION_BEAMFORMER << std::endl;
                exit(EXIT_SUCCESS);
                break;
        }
    }
    
    if (argc == 1)
    {
        usage();
        exit(EXIT_FAILURE);
    }
    
    manager_t the_manager; // dropped the volatile
    
    the_manager.shutdown=0;
    the_manager.offline=offline;
    the_manager.integrate=dumps_to_aver;
    the_manager.chan_to_aver=chan_to_aver;
    the_manager.edge = edge;
    the_manager.nbit = nbit;
    the_manager.coarse_chan = coarse_chan;
    the_manager.nstation = 128;
    the_manager.nfrequency = nfrequency;
    the_manager.ndim = 2;
    the_manager.npol = 2;
    the_manager.dumps_per_sec = dumps_per_second;
    the_manager.infile = 0;

    if (in_file != NULL) {
        // we have an input file
        if ((the_manager.infile = open(in_file,O_RDONLY)) == -1)
        {
            std::cerr << "error: input (" << in_file << ") file selected but cannot be opened" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if (starttime < 0)
    {
        usage();
        std::cerr << "error: offline mode selected but no starttime on command line" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    if (obsid < 0)
    {
        usage();
        std::cerr << "error: offline mode selected but no obsid on command line" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    
    
    /*
     First define the input ring buffers, as a thow back to earlier code this
     ringbuffer technology is the same technologu employed by some parts of the media conversion code
     it is not the same buffer technology employed for the output
     a buffer should be the size of an input buffer
     */
    
    /* lets allocate some buffers
     * need pinned memory
     */
    // Get sizing info from library
    
    xgpuInfo(&xgpu_info);
    if (npol != xgpu_info.npol)
    {
        std::cerr << "error: fatal missmatch between XGPU library and requested npol XGPU: "
                  << xgpu_info.npol << ", REQUESTED: " << npol << std::endl;
        exit(EXIT_FAILURE);
    }
    if (nstation != xgpu_info.nstation)
    {
        std::cerr << "error: fatal missmatch between XGPU library and requested nstation XGPU: "
                  << xgpu_info.nstation << ", REQUESTED: " << nstation << std::endl;
        exit(EXIT_FAILURE);
    }
    if (nfrequency != xgpu_info.nfrequency)
    {
        std::cerr << "error: fatal missmatch between XGPU library and requested channels XGPU: "
                  << xgpu_info.nfrequency << ", REQUESTED: " << nfrequency << std::endl;
        exit(EXIT_FAILURE);
    }
    ntime = xgpu_info.ntime;
    
    the_manager.ntime = ntime;
    size_t full_matLength = nfrequency * nstation * nstation * npol * npol;
    size_t full_size = dumps_per_second * full_matLength * sizeof(Complex);
    size_t baseLength = nfrequency;
    
    size_t ring_bufsz = xgpu_info.vecLength * sizeof(ComplexInput);
    
    char **cuda_buffers = (char **) calloc ((RING_NBUFS+1),sizeof(char*));
    
    size_t numbytes =(((ring_bufsz)+4095)/4096)*4096; // page size and page aligned
    
    for (int i = 0; i <= RING_NBUFS; i++)
    {
        std::cout << "allocating buffer " << i << " of " << numbytes << std::endl;
        
        cuda_buffers[i] = (char *) valloc(numbytes);
        
        if (cuda_buffers[i] == NULL) {
            std::cerr << "error: failed to allocate buffer " << i << std::endl;
            exit(EXIT_FAILURE);
        }
        
        if ((i>0) && (i<RING_NBUFS)) { // xgpu_Init will register the first buffer but not the rest
            cudaHostRegister(cuda_buffers[i],numbytes,0);
            checkCudaError();
        }
        
    }
    
    std::cout << "assigning buffers" << std::endl;
    if (assign_ring_buffers(RING_NBUFS,ring_bufsz,cuda_buffers,&ring) < 0)
    {
        std::cerr << "Failed to ASSIGN ringbuffer" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    Complex *full_matrix_h = NULL;
    Complex *baseline_h    = NULL;
    
    /*
     * the beamformer results. Format ... 8 bit int
     
     rember this is a total intensity sum.
     
     */
#ifdef RUN_BEAMER
    uint8_t *beam_h = NULL;
    uint8_t *beam_d = NULL;
    int8_t *data_d = NULL;
    size_t timesteps_per_call = TIMESTEPS_PER_CALL;
    // how much data is processed per beamformer call
    size_t step_data_size = xgpu_info.nfrequency * xgpu_info.nstation * xgpu_info.npol * 2 * timesteps_per_call; // complex
    
    size_t step_results_size = xgpu_info.nfrequency*timesteps_per_call*sizeof(uint8_t);
    size_t beam_size = xgpu_info.nfrequency*xgpu_info.ntime*dumps_per_second*sizeof(uint8_t);
    /* FIXME: please check return codes here */
    
    /* ntime is the number of timesamples per GPU call for the correlator */
    /* dumps_per_second is the number of correlator dumps there are every second */
    
    /* full seconds worth of output beam */
    
    beam_h = (uint8_t *) malloc(beam_size);
    
    /* the input data for the beamformer will simple be offset into the buffer */
    /* But we will need to assign the device memory - this only need to be a number of timesteps equal to a pipelength
     */
    
    cudaMalloc(&beam_d,step_results_size*sizeof(uint8_t)); //
    cudaMalloc(&data_d,step_data_size*sizeof(int8_t)); //
#endif
    
    full_matrix_h = (Complex *) malloc(full_size);
    
    baseline_h = (Complex *) malloc(baseLength * sizeof(Complex));
    
    the_manager.ring = &ring;
    
    the_manager.nfrequency = nfrequency;
    the_manager.nstation = nstation;
    the_manager.npol = npol;
    
    if (pthread_create(&buffer_handler, NULL, manager,
                       (void *) &the_manager)) {
        std::cerr << "error: could not launch manager thread" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    std::cout << "Launched manager thread" << std::endl;
    
    uint64_t blockSize = 0;
    int hdu_num = 0;
    
    // The PRIMARY header + image + padding
    uint64_t n_visibilities = ((uint64_t) xgpu_info.nbaseline *4) ;
    
    while (hdu_num < dumps_per_second) {
        blockSize = blockSize + 2880; // header
        blockSize = blockSize +  (n_visibilities * (uint64_t) xgpu_info.nfrequency * sizeof(Complex));; // sizeof a data cube
        
        int remainder = (blockSize%2880); // pad out to the end of the HDU
        blockSize = blockSize + (2880 - remainder);
        hdu_num++; // hdu increment
    }
    
    assert(!(blockSize%2880));
    
    std::cout << "Correlating " << nstation << " stations, with " << npol
              << " signals, with " << nfrequency << " channels" << std::endl;
    
    // allocate the GPU X-engine memory
    XGPUContext context;
    context.array_h = (ComplexInput *) cuda_buffers[0]; // already asssigned this above
    context.matrix_h = NULL; // we are letting the xgpu library configure tis memory - should make sure it is big enough
    context.array_len = xgpu_info.vecLength;
    context.matrix_len = xgpu_info.matLength;
    
    xgpu_error = xgpuInit(&context,0);
    
    if(xgpu_error)
    {
        std::cerr << "error: xgpuInit returned error code " << xgpu_error << std::endl;
        xgpuFree(&context);
        return xgpu_error;
    }
    
    Complex *cuda_matrix_h = context.matrix_h;
    
    int dumps_integrated = 0;
    
    
    tzset();
    extern long timezone;
    
    while (1) {
        
        
        strptime((const char *) the_manager.start_obs_UTC,"%Y-%m-%d-%H:%M:%S",&start_utctime);
        
        start_time_t = starttime;
        
        if (start_time_t != current_time_t)
        {
            std::cout << "start_time is " << the_manager.start_obs_UTC << ": "
                      << "decodes to: " << start_time_t << ", "
                      << "current is: " << current_time_t << std::endl;
            std::cout << "start_time is not current (" << start_time_t << ": "
                      << current_time_t << "): restart detected" << std::endl;
            std::cout << "Integrate " << the_manager.integrate << ": "
                      << "Chan to aver " << the_manager.chan_to_aver << std::endl;
            /* there has been a restart therefor the start time in the header is different to the expexted
             * start time*/
            current_time_t = start_time_t;
            incremented_time_t = current_time_t;
        }
        
        
        int x_done = 0; // how much of the second have we done
        
        // we only pass on a full second to the FITs builder. There are dumps_per_second correlator integrations
        // however there is nothing stopping the internal integration time of the correlator being much less than that
        // We are now capturing this with the integrate flag.
        
        // This integration now corresponds to howmany internal cuda invocations make up an integration, dumps per second is how many correlator integrations there are per second
        
        while (x_done < dumps_per_second) { // how many correlator dumps per second
            
            
            
            buf = NULL;
            static int count = 0; // just a check to see if the buffers are taking too long to drain,
            while (buf == NULL) {
                
                buf = wait_for_buffer(&ring); // the only way this returns is if there is a full buffer to read/or EOD/or overrun
                
                if (ring.EOD) // this can be set and still there can be data in the ring
                {
                    std::cout << "NOTICE: EOD on input buffer" << std::endl;
                    if (buffer_EOD(&ring)== 0) {
                        std::cout << "NOTICE:: EOD on input buffer - but ring not yet empty : no reset yet" << std::endl;
                        count++;
                        sleep(1);
                        if (count > 5) {
                            std::cerr << "warning: waited > 5s for buffer to drain -- forcing reset" << std::endl;
                            reset_ring_buffers(&ring);
                            buf = NULL;
                            goto SHUTDOWN;
                        }
                    }
                    else if (buffer_EOD(&ring) == 1)
                    {
                        std::cout << "NOTICE:: EOD on input buffer drained - reset" << std::endl;
                        reset_ring_buffers(&ring);
                        buf = NULL;
                        goto SHUTDOWN;
                    }
                }
                else if (ring.overrun)
                {
                    std::cerr << "error: OVERRUN hard reset" << std::endl;
                    get_buffer_status(&ring);
                    reset_ring_buffers(&ring);
                    goto SHUTDOWN;
                }
                
                // get_buffer_status(&ring);
                
                
            } // get buffer
            
            /* --------------------------- Run The Correlator --------------------------- */
            
            /* We now should have PIPE_LENGTH chunks of NTIME_PIPE samples in the pinned
             memory - we should be free to run the correlator now
             */
            
            
            gettimeofday(&clock1,NULL);
            /*
             * this is the point where we increase our dump time. I am going to offset into the buffer the correct amount for the
             * sub second and process as normal - this should restrict the complication to here only
             *
             */
            
            
            xgpu_error = 0 ;
            
            if ((the_manager.integrate == 0) || dumps_integrated == 0) { // first internal integration (on the GPU)
                
               xgpu_error = xgpuClearDeviceIntegrationBuffer(&context);
                
                if(xgpu_error) {
                    std::cerr << "error: xgpuCudaXengine returned error code " << xgpu_error << std::endl;
                    xgpuFree(&context);
                    return xgpu_error;
                }
                dumps_integrated = 0;
                
            }
            
            
            context.array_h = (ComplexInput *) buf;
            // context.input_offset =  0; // there is no offset in this mode as the buffer only contains enough data for NTIME
            
            xgpu_error = xgpuCudaXengine(&context,SYNCOP_DUMP);
            
            checkCudaError();
            
            if(xgpu_error) {
                std::cerr << "error: xgpuCudaXengine returned error code " << xgpu_error << std::endl;
                xgpuFree(&context);
                return xgpu_error;
            }
            
            cudaThreadSynchronize();
            
            if (the_manager.integrate != 0) {
                dumps_integrated++; // internal cuda counter
            }
            else {
                dumps_integrated = 0;
            }
 
            std::cout << "GPU X-Engine done (" << dumps_integrated << ":"
                      << the_manager.integrate << ")" << std::endl;
            gettimeofday(&clock2,NULL);
            elapsed = (clock2.tv_sec - clock1.tv_sec) * 1000.0;      // sec to ms
            elapsed += (clock2.tv_usec - clock1.tv_usec) / 1000.0;   // us to ms
            
            std::cout << "Correlator/Beamformer took " << elapsed << " milliseconds" << std::endl;
            
            //
            
            // Lets copy the cube out into another area (it is only small)
            //
            //
            // Launch the beamformer
            
            
            
           
            
#ifdef RUN_BEAMER
            
            int time_step = 0;
            int steps = 10000;
            int8_t *data = NULL;
            uint8_t *results = NULL;
            data = (int8_t *) buf;
            results = beam_h;
            
            while(time_step < steps) {
                run_beamer(data_d,data,beam_d,results,step_data_size,step_results_size);
                data =  data + step_data_size;
                results =  results + step_results_size;
                time_step=time_step+TIMESTEPS_PER_CALL;
            }
            
#endif
            
            if (dumps_integrated == the_manager.integrate) {
                
              
                
                Complex *ptr = full_matrix_h + (x_done * context.matrix_len);
                memcpy(ptr,cuda_matrix_h,(context.matrix_len * sizeof(Complex)));
                
                xgpu_error = xgpuClearDeviceIntegrationBuffer(&context);
                
                if(xgpu_error) {
                    std::cerr << "error: xgpuCudaXengine returned error code " << xgpu_error << std::endl;
                    xgpuFree(&context);
                    return xgpu_error;
                }
                dumps_integrated = 0;
                x_done++; // external sub-int count - we must copy out the accumulation and reset the integration buffer
                
            }
            mark_buffer_empty(&ring); // GPU input buffer mark clear
            
        } // done 1 seconds worth
        /* -------------------------- Write out the Product -------------------- */
        
        
        blockSize=0;
        hdu_num = 0;
        while (hdu_num < dumps_per_second) {
            
        
        
            blockSize += 2880; // header
            blockSize += n_visibilities * (uint64_t) xgpu_info.nfrequency *
                         sizeof(Complex) / the_manager.chan_to_aver; // sizeof a data cube
            
            int remainder = (blockSize%2880); // pad out to the end of the HDU
            blockSize += 2880 - remainder;
            hdu_num++; // hdu increment
        }
        
        assert(!(blockSize % 2880));
        
        if (the_manager.integrate != 0) {
            std::cout << "Integrated " << dumps_integrated << " invocations" << std::endl;
            std::cout << "blockSize " << (int)blockSize << std::endl;
        }
        
        char *outbuffer = (char *) malloc(blockSize);
        
        buildFITSBuffer( xgpu_info,
                         full_matrix_h,
                         blockSize,
                         (void *)outbuffer,
                         incremented_time_t,
                         dumps_per_second,
                         &the_manager );
        
        
        std::cout << "FITS file built" << std::endl;
        
        gmtime_r(&incremented_time_t,&current_utctime);
        
        std::cout << "Buffer time set " << asctime(&current_utctime) << std::endl;
        
        strftime(file_time,15,"%Y%m%d%H%M%S",&current_utctime);
        
        sprintf( dump_filename, "/%s_%s_gpubox%02d_00.fits",
                 obsid, file_time, coarse_chan );
        
        FILE *outf = fopen( dump_filename, "w" );
        
        if (outf != NULL) {
            
            fwrite(outbuffer, blockSize, 1, outf);
            fclose(outf);
            std::cout << "Last cube dumped" << std::endl;
            
        }
        else {
            
            std::cerr << "error: failed to open dump file " << strerror(errno)
                      << std::endl;
            
        }
        
        free(outbuffer);
        
        
        
        incremented_time_t++; // increment time by a second;
        
        gettimeofday(&clock3,NULL);
        elapsed = (clock3.tv_sec - clock2.tv_sec) * 1000.0;      // sec to ms
        elapsed += (clock3.tv_usec - clock2.tv_usec) / 1000.0;   // us to ms
        std::cout << "Data output (FITS building etc took a further "
                  << elapsed << " milliseconds" << std::endl;
        
        elapsed = (clock3.tv_sec - clock1.tv_sec) * 1000.0;      // sec to ms
        elapsed += (clock3.tv_usec - clock1.tv_usec) / 1000.0;   // us to ms
        std::cout << "Total processing took " << elapsed
                  << " milliseconds" << std::endl;
        
        
        
    }
    /* ------------------------ Do the Book-keeping ---------------------- */
SHUTDOWN:
    sleep(2);
    for (int i = 0; i < RING_NBUFS; i++) {
        cudaFreeHost(cuda_buffers[i]);
        cudaHostUnregister(cuda_buffers[i]); // device zero only
    }
    
    xgpuFree(&context);
#ifdef RUN_BEAMER
    cudaFree(beam_d); //
    cudaFree(data_d); //
#endif
    free(full_matrix_h);
    free(baseline_h);
    return xgpu_error;
    
}
