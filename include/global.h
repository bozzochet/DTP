
// containing global variables

#ifndef GLOBAL_INCLUDE
#define GLOBAL_INCLUDE


#include "types.h"

//detector geometry
const int Nlayers = 10;
const int Nsquares = 8; //squares per side on a layer
const int Nrows = 2; //ladders rows on a single layer
const int Nlad = Nsquares*Nrows*Nlayers; //total number of ladders
const int Nstrips = 640; //strips per ladder
const double pitch = 0.015; //cm = 150 um
const double squareSide = Nstrips*pitch;

//energy released to generate a couple electron-hole in Si
const energy_t E_COUPLING = 8.145794455; //eV/couple
const charge_t FOND_CHARGE = 1.602e-19; //C


#endif
