
//simple struct to save simulation results

#ifndef MEASURE_INCLUDE
#define MEASURE_INCLUDE


#include "physics.h"


struct measure
{
  mytime_t time[2]   = {-9999, -9999};
  int      xy        = -9999;           //0 if pos is x, 1 if pos is y
  length_t position  = -9999;
  energy_t energy[2] = {-9999, -9999};
};


#endif //include guard
