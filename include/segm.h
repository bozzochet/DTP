
// general functions for time and position segmentations

#ifndef SEGM_INCLUDE
#define SEGM_INCLUDE


#include "geometry.h"
#include "vector2.h"


namespace one_every_N
{
  //active strips are one every jump
  void shareEnergy(vector2<double> &, const int &jump);
}


#endif
