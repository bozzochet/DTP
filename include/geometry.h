
// detector geometry parameters

#ifndef GEOMETRY_INCLUDE
#define GEOMETRY_INCLUDE


const int Nlayers = 10;
const int Nsquares = 8; //squares per side on a layer
const int Nrows = 2; //ladders rows on a single layer
const int Nladders = Nsquares*Nrows*Nlayers; //total number of ladders
const int Nstrips = 640; //strips per ladder
const double pitch = 0.015; //cm = 150 um
const double squareSide = Nstrips*pitch;


#endif
