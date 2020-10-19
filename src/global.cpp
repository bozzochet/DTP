
#include "global.h"

const int Nlayers = 10;
const int Nsquares = 8; //squares per side
const int Nlad = Nsquares*2*Nlayers; //number of ladders
const double squareSide = 10;
const double pitch = 0.015625;
const int Nstrips = (int(squareSide/pitch)); //strips per ladder
const int jump = 2;
