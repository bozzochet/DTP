
// position segmentation functions and variables

#ifndef POS_SEGM_INCLUDE
#define POS_SEGM_INCLUDE


#include "segm.h"
#include "vector2.h"


namespace pos_segm
{
  const int jump = 3;

  inline void shareEnergy(vector2<double> &energy)
  { one_every_N::shareEnergy(energy, jump); }
}


#endif
