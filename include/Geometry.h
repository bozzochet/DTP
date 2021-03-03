
/*****************************************************************
   useful in Digitization.cpp to make use of detector geometric
   parameters
 *****************************************************************/

#ifndef GEOMETRY_INCLUDE
#define GEOMETRY_INCLUDE


#include "physics.h"


struct Geometry
{
  int Nlayers = -1 ;
  int Nstrips = -1 ;
  int Nrows = -1 ;
  int Nsquares = -1 ;
  length_t pitch = -1 ;
  length_t thickness = -1 ;

  length_t squareSide = pitch * Nstrips;
  int Nladders = Nsquares * Nrows * Nlayers;
};


#endif //include guard
