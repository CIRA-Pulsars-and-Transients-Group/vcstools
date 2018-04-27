#ifndef PABEAM_H
#define PABEAM_H

// Standard includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "utils.h"
#include "pabeam_kernal.h"


// MWA tile beam
#include "FEE2016/beam2016implementation.h"
#include "FEE2016/mwa_beam_interface.h"
#include "FEE2016/system.h"
#include "H5Cpp.h"


void usage();

int getNumTiles(const char *metafits);

void getTilePositions(const char *metafits, int ninput, 
                        float *n_pols, float *e_pols, float *h_pols,
                        float *n_tile, float *e_tile, float *h_tile);

int getFlaggedTiles(const char *badfile, int *badtiles);

void removeFlaggedTiles(float *n_tile, float *e_tile, float *h_tile, 
                        int *badtiles, int nbad, int nelements);


#endif
