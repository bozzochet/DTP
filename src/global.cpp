
#include "global.h"

const int Nlayers = 10;
const int Nsquares = 8; //squares per side on a layer
const int Nlad = Nsquares*2*Nlayers; //total number of ladders
const int Nstrips = 640; //strips per ladder
const double pitch = 0.015625;
const double squareSide = Nstrips*pitch;
