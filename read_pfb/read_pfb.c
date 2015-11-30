#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <assert.h>

#include "mwac_utils.h"
#include "antenna_mapping.h"


void usage() {
    
    fprintf(stderr,"read_pfb -n <nchan [128]> -s <nsteps [all]> -b {turns OFF binary} \n");
    fprintf(stderr,"input from stdin the contents of a gpu dump file\n");
    fprintf(stderr,"output to stdout the unpacked integers\n");
    fprintf(stderr,"note: transposes so whole antennas are output\n");
    fprintf(stderr,"note: reorders to PFB input order\n");
    fprintf(stderr,"\t -a nstation [128]\n");
    fprintf(stderr,"\t -n nchan [128]\n");
    fprintf(stderr,"\t -p npol [2]\n");
    fprintf(stderr,"\t -b [0 | 1 | 4 | 16] - nbytes of binary output 1byte int, 4byte float, 16 byte complex\n");
    fprintf(stderr,"\t -i <input file> -- defaults to stdin");
    fprintf(stderr,"\t -o <output file> -- defaults to stdout");
    fprintf(stderr,"\t -s nsteps [-1]: number of seconds -1 is ALL\n");
    fprintf(stderr,"\t -c count [10000]: number of samples per step\n");
    fprintf(stderr,"\t -8 : 8 bit input samples\n");
    fprintf(stderr,"\t -4 : 4 bit input samples\n");
    
}
int nstation;
int npol;
int nfrequency;

int main(int argc, char **argv) {
    
    int c;
    
    int nchan=128;
    int ncount=10000;
    int nstep=-1;
    int pol=0;
    
    int8_t *input_buffer = NULL;
    int8_t *output_buffer = NULL;
    int8_t *binary_buffer = NULL;

    int ascii = 0;
    int binary = 1;
    size_t input_nbit = 4;
    
    char *ifname,*ofname;
    int in_fd = 0; // stdin
    int out_fd = 1; // stdout
    FILE *fin = NULL;
    FILE *fout = NULL;
    mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
    nstation = 128;
    npol = 2;
    
    if (argc > 1) {
        
        while ((c = getopt(argc, argv, "a:b:c:hi:n:o:p:s:84")) != -1) {
            switch(c) {
                case 'a':
                    nstation = atoi(optarg);
                    break;
                case 'b':
                    binary = atoi(optarg);
                    if (binary == 0) {
                        ascii = 1;
                    }
                    else {
                        ascii = 0;
                    }
                    break;
                case 'c':
                    ncount = atoi(optarg);
                    break;
                case 'h':
                    usage();
                    exit(-1);
                    break;
                case 'i':
                    ifname = strdup(optarg);
                    in_fd = open(ifname,O_RDONLY);
                    if (in_fd < 0) {
                        perror("Failed to open input:");
                        exit(EXIT_FAILURE);
                    }
                    break;
                case 'n':
                    nchan = atoi(optarg);
                    break;
                case 'o':
                    
                    ofname = strdup(optarg);
                    out_fd = open(ofname,O_WRONLY|O_CREAT|O_TRUNC,mode);
                    if (out_fd < 0) {
                        perror("Failed to open output");
                        exit(EXIT_FAILURE);
                    }
                    break;
                case 'p':
                    npol = atoi(optarg);
                    break;
                case 's':
                    nstep = atoi(optarg);
                    break;
                case '8':
                    input_nbit = 8;
                    break;
                case '4':
                    input_nbit = 4;
                    break;
                default:
                    usage();
                    exit(EXIT_FAILURE);
            }
        }
    }
    // Initialise
    //
    fin = fdopen(in_fd,"r");
    
    if (fin == NULL) {
        perror ("Failed to fdopen the input");
        exit(EXIT_FAILURE);
    }
    
    fout = fdopen(out_fd,"w");
    
    if (fout == NULL) {
        perror ("Failed to fdopen the output");
        exit(EXIT_FAILURE);
    }
    
    nfrequency=nchan;
    npol=2;
    
    input_buffer = malloc((nstation*npol*nchan*2*input_nbit)/8);
    assert(input_buffer);
    output_buffer = malloc(nstation*npol*nchan*2*sizeof(int8_t));
    assert(output_buffer);
    if (binary != 0) {
        binary_buffer = malloc(nstation*npol*nchan*2*binary);
        assert(binary_buffer);
    }
    int index = 0;
    int chan = 0;
    int out_index = 0;
    int step = 0;
    

    size_t rtn = 0;
    size_t items_to_read = 0;
    size_t count = 0;
    
    uint8_t original;
    int8_t value;
    
    fill_mapping_matrix();

    size_t gulp = (nstation*npol*nchan*2*input_nbit)/8;
    int i = 0;
    for (i = 0; i < NINPUT; i++) {
        fprintf(stdout,"(READ_PFB::)Index %d is antenna %d\n",i,pfb_output_to_input[i]);
        fflush(stdout);
    }


    while (!feof(fin) && (nstep < 0 || step < nstep)) {
        
        for (count = 0 ; count < ncount ; count++) { 
            items_to_read = gulp;
            rtn = 0;
            rtn = fread(input_buffer,1,items_to_read,fin);
            
            if (feof(fin) || rtn != items_to_read) {
                if (feof(fin) && in_fd == 1) {
                    fprintf(stderr,"read_pfb: EOF on input on step %d read no.:%lu\n",step,count);
                    fclose(fin);
                    exit(1);
                }
                else if (feof(fin) && in_fd !=1) {
                    fprintf(stderr,"read_pfb finished (%lu returned on step %d read no. %lu):\n", rtn,step,count);
                    fclose(fin);
                    exit(0);
                    
                }
                
                fclose(fin);
                exit(1);
            } 
           
            // time then frequency runs slowest in this data block
            
            int8_t *inp_ptr = NULL;
            int map_index = 0;
            // #pragma omp parallel num_threads(4)
            // {
            // #pragma omp for private (chan,index,out_index,pol,inp_ptr,original,value)
            for (chan = 0 ; chan < nchan ; chan++) {
                
                inp_ptr = (int8_t *) input_buffer + ((chan*nstation*npol*2*input_nbit)/8);
                for (index = 0; index < nstation; index++) {
                    
                    for (pol =0 ;pol < npol; pol++) {
                        // this returns the desired output index (0 to 64)
                        map_index = pfb_output_to_input[index*npol + pol];
                        // output matrix index - this index reorders the
                        // stations AND transposes
                        // the order will now be [time][input][chan][complexity]
                        //
                        out_index = map_index*nchan*2 + chan*2;
                                              
                        // .. assigns the value to be the current de-ref input_ptr
                        if (input_nbit == 8) {
                            output_buffer[out_index] = (int8_t) *inp_ptr;
                            // increment both pointers to get to the imaginary part
                            inp_ptr++;
                            out_index++;
                            // assign imaginary part
                            output_buffer[out_index] = (int8_t) *inp_ptr;
                            inp_ptr++;
                            // increment the input pointer
                        }
                        else if (input_nbit == 4) {
                            original = (uint8_t) *inp_ptr;
                            original = original & 0xf; // first sample
                            if (original >= 0x8) {
                                value = original - 0x10;
                            }
                            else {
                                value = original;
                            }
                            
                            output_buffer[out_index] = (int8_t) (value & 0xff);
                            
                            //fprintf(stdout,"unpack:chan %d in %d map %d out %d original %u:%d \n",chan,index*npol+pol,map_index,out_index,original,output_buffer[out_index]);
                            //fflush(stdout);
                            
                            out_index++;
                            
                            original = (uint8_t) *inp_ptr ;
                            original = original >> 4;
                            original = original & 0xf; // second sample
                            if (original >= 0x8) {
                                value = original - 0x10;
                            }
                            else {
                                value = original;
                            }
                            output_buffer[out_index] = (int8_t) (value & 0xff);

                            //fprintf(stdout,"unpack:chan %d in %d map %d out %d original %u:%d \n",chan,index*npol+pol,map_index,out_index,original,output_buffer[out_index]);

                            //fflush(stdout);

                            
                            inp_ptr++;
                        }		    
                        
                    } // for all pol
                } // for all stations
            } // for all channels
            //}
            
            if (binary!=0) { // write out the whole block for this time sample
                int8_t *out_ptr;
                size_t out_size;
                out_ptr = output_buffer;
                out_size = 2*nstation*npol*nchan;
                
                if (binary == 1) {				
                    fwrite(out_ptr,sizeof(int8_t),out_size,fout);
                }
                if (binary == 4) {
                    int bindex=0;
                    float fval=0.0;
                    for (bindex=0;bindex<out_size;bindex++) {
                        fval = (float) *out_ptr;
                        memcpy(&binary_buffer[binary*bindex],&fval,binary);
                        out_ptr++;
                    }
                    fwrite(binary_buffer,binary,out_size,fout);
                }
                else if (binary == 16) {
                    int bindex=0;
                    complex float fval=0.0;
                    for (bindex=0;bindex<out_size/2;bindex++) {
                        fval = (double) *out_ptr + 0.0*I;
                        out_ptr++;
                        fval += (double) I * (double) *out_ptr;
                        out_ptr++;
                        memcpy(&binary_buffer[binary*bindex],&fval,binary);
                    }
                    fwrite(binary_buffer,binary,out_size/2,fout);
                }
            }
            else if (ascii == 1) {
                int8_t *sample= (int8_t *) output_buffer;
                int samps = 2*nstation*npol*nchan;
                int sampnum = 0;
                for(sampnum=0;sampnum < samps;sampnum=sampnum+2) {
                    fprintf(fout,"%d %d\n",*sample,*(sample+1));
                    sample=sample+2;
                }
            }
            
        } // next time step
        step++; // increment the steps
    }// still data?
    
    fprintf(stderr,"nstep: %d\n",step);
    free(input_buffer);
    free(output_buffer);
    free(binary_buffer);
}
