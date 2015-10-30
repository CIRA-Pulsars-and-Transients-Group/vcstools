//
//  compress.cpp
//  uvcompress
//
//  Created by Slava Kitaeff on 2/04/13.
//  Copyright (c) 2013 ICRAR/UWA. All rights reserved.
//

#include "compress.h"

char errmsg[80];

#define round(r) (r>0.0)?floor(r+0.5):ceil(r-0.5)

//
// fits_write_compressed()
//
// The function scales, decimates the data, and stores it in
// another FITS iamge file in compressed form.
//
// Arguments:
// *out - pointer to already opened destination FITS file
// *buff_in - pinter to the data array
// bscale - scaling factor used to force the precition
// comp - type of binary compression to use. Can be either RICE_1 or GZIP_1
// nelements - number of elements in array
// naxis & naxes[] - as per cfitsio library description
//
int fits_write_compressed(fitsfile *out, float *buff_in, LONGLONG nelements, int naxis, long *naxes, double bscale, int comp)
{
    int status = 0;    // returned status of FITS functions
    int *buff_out = (int*) malloc(sizeof(int) * nelements);
    if( buff_in == NULL ){
        cerr << "Err: Could not allocate memory.";
        return 1;
    }

    // round, decimate and convert
    for(int i=0; i<nelements; ++i)
        buff_out[i] = (int)round(bscale * buff_in[i]);
    
    //lets create HDU for encoded data
    fits_set_compression_type(out, comp, &status);
    PRINTERRMSG(status);
    
    fits_create_img(out, LONG_IMG, naxis, naxes, &status);
    PRINTERRMSG(status);
    
    fits_write_img(out, TINT, 1, nelements, buff_out, &status);
    PRINTERRMSG(status);

    // add the keys
    double bzero = 0;
    bscale = 1.0/bscale;
    fits_write_key(out, TDOUBLE, "BSCALE", &bscale, NULL, &status);
    fits_write_key(out, TDOUBLE, "BZERO", &bzero, NULL, &status);
    PRINTERRMSG(status);
    
    free(buff_out);

    return 0;
}

//
// Compress()
//
// The function reads the source FITS image data, and stores it in
// another FITS iamge file in compressed form. The function works with
// 32-bit floating-point data only.
//
// Arguments:
// *in & *out - pointers to already opened source and destination FITS files correspondently
// bscale - scaling factor used to force the precition
// comp - type of binary compression to use. Can be either RICE_1 or GZIP_1
//
int Compress(fitsfile *in, fitsfile *out, double bscale, int comp)
{
    long naxes[]={1,1,1,1,1,1,1,1,1};
    int naxis;
    int status = 0;    // returned status of FITS functions
    int bitpix;
    
    // loop through till the end of file
    while(status != END_OF_FILE){
        
        // get image dimensions and total number of pixels in image
        fits_get_img_param(in, 9, &bitpix, &naxis, naxes, &status);
        PRINTERRMSG(status);

        // try reading the whole image
        LONGLONG nelements = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4] * naxes[5] * naxes[6] * naxes[7] * naxes[8];
        
        float *buff_in = (float*) malloc(sizeof(float) * nelements);
        if( buff_in == NULL ){
            cerr << "Err: Could not allocate memory.";
            return 1;
        }
        
        float nulval = 0.0;
        fits_read_img(in, TFLOAT, 1, nelements, &nulval, buff_in, NULL, &status);
        PRINTERRMSG(status);
        
        fits_write_compressed(out, buff_in, nelements, naxis, naxes, bscale, comp);

        free(buff_in);
        
        // try next HDU
        fits_movrel_hdu(in, 1, NULL, &status);
    }
    
    return 0;
}


//
// Decompress()
//
// The function reads compressed source FITS image data, and stores it in
// another FITS iamge file in uncompressed form. The function works with
// DEC compressed data only.
//
// Arguments:
// *in & *out - pointers to already opened source and destination FITS files correspondently
//

int Decompress(fitsfile *in, fitsfile *out)
{
    int status = 0;    // returned status of FITS functions
    long naxes[]={1,1,1,1,1,1,1,1,1};
    int naxis;
    int bitpix;
    int hdutype;
    float nulval = 0.0;
    
    //loop through
    while(status != END_OF_FILE){

        //check if it's compressed image HDU
        fits_read_key(in, TLOGICAL, "ZIMAGE", &hdutype, NULL, &status);
        if( !hdutype ) {
            cerr << "Warrning: The data doesn't seem to be compressed. Nothing to be done.";
            return 1;
        }
        // get image dimensions and total number of pixels in image
        fits_get_img_param(in, 9, &bitpix, &naxis, naxes, &status);
        PRINTERRMSG(status);

        LONGLONG nelements = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4] * naxes[5] * naxes[6] * naxes[7] * naxes[8];
        
        float *buff = (float*) malloc(sizeof(float) * nelements);
        if( buff == NULL ){
            cerr << "Err: Could not allocate memory!";
            return 1;
        }
        
        fits_read_img(in, TFLOAT, 1, nelements, &nulval, buff, NULL, &status);
        PRINTERRMSG(status);
        
        fits_create_img(out, FLOAT_IMG, naxis, naxes, &status);
        PRINTERRMSG(status);
        
        fits_write_img(out, TFLOAT, 1, nelements, buff, &status);
        PRINTERRMSG(status);
        
        free(buff);
        
        //move to the next HDU
        fits_movrel_hdu(in, 1, NULL, &status);
    }
    
    return 0;
}