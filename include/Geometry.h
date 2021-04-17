
/*****************************************************************
   useful in Digitization.cpp to make use of detector geometric
   parameters
 *****************************************************************/

#ifndef GEOMETRY_INCLUDE
#define GEOMETRY_INCLUDE

#include "TObject.h"
#include "physics.h"


class Geometry: public TObject {

public:
  length_t CaloSide = -1;
  length_t CaloStkGap = -1;
  int Nsquares = -1 ;
  int Nrows = -1 ;
  int Nlayers = -1 ;
  length_t LayerGap = -1;
  length_t PlaneGap = -1;  
  int Nstrips = -1 ;
  length_t pitch = -1 ;
  length_t thickness = -1 ;

  length_t squareSide = -1;
  int Nladders = -1;
  
  void ComputeDerived() {
    squareSide = pitch * Nstrips;
    Nladders = Nsquares * Nrows * Nlayers;
  }

  ClassDef(Geometry, 2)
};


#endif //include guard
