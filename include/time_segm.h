
/*********************************************************************

  time segmentation functions and variables;

  there are implemented three types (A,B,C) of segm:

    A - measures are given considering groups of N strips as one

    B - measures are given considering i-th strip, the (i+N)-th strip
      etc... as one.

    C - measures are given considering Si pads placed on the backplate
      of the layer and oriented perpendicular to the strip direction

 ********************************************************************/

#ifndef TIME_SEGM_INCLUDE
#define TIME_SEGM_INCLUDE


#include "segm.h"
#include "vector2.h"


enum time_segm : const int {A, B, C};

namespace A_time_segm
{
  const int jump = 10;

  inline void shareEnergy(vector2<double> &energy)
  { one_every_N::shareEnergy(energy, jump); }
}

namespace B_time_segm
{
  const int jump = 10;
}

namespace C_time_segm
{
  // // // // // // //
}


#endif //include guard
