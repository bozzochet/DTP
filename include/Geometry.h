
/*****************************************************************
   useful in Digitization.cpp to make use of detector geometric
   parameters
 *****************************************************************/

#ifndef GEOMETRY_INCLUDE
#define GEOMETRY_INCLUDE

#include "TObject.h"
#include "physics.h"
#include <iostream>

class Geometry: public TObject {

public:

  // VARIABILI DTP
  length_t CaloSide = -1;
  int Nrows = -1 ;
  int Nsquares = -1 ;
  length_t PlaneGap = -1;
//////////////////////////////////////////////////////////////////////
  int Nstrips = -1 ;  //Queste due si possono mettere anche su comuni ma per ora le lasciamo qui
  length_t pitch = -1 ;
//////////////////////////////////////////////////////////////////////
  length_t squareSide = -1;
  int Nladders = -1;
  
  
  // VARIABILI SLA
  int NCaloPlanes = -1 ;
  length_t PlanesDistance = -1 ;
  length_t CubeX = -1 ;
  length_t CubeY = -1 ;
  length_t CubeZ = -1 ;
  int RowsOfCubes = -1 ;
  int ColumnsOfCubes = -1 ;
  int NScintillators = -1 ;
  length_t LayerGapScintillator = -1 ;

  // VARIABILI COMUNI
  int Nlayers = -1 ;
  length_t LayerGap = -1;
  length_t thickness = -1 ;
  length_t CaloStkGap = -1;  
  
    
  
  
  

  
  
  void ComputeDerived() {
    squareSide = pitch * Nstrips;
    Nladders = Nsquares * Nrows * Nlayers;
  }

  void Dump() const {
    std::cout << std::endl;
    std::cout <<"=================================" << std::endl;
    std::cout <<"Geometric parameters:       " << std::endl;
    std::cout <<"  Calo side:                " << CaloSide << std::endl;
    std::cout <<"  Calo-Stk gap:             " << LayerGap << std::endl;
    std::cout <<"  wafers per side:          " << Nsquares  << std::endl;
    std::cout <<"  ladders per 'column':     " << Nrows  << std::endl;
    std::cout <<"  layers:                   " << Nlayers << std::endl;
    std::cout <<"  gap between layers:       " << LayerGap << std::endl;
    std::cout <<"  gap between planes:       " << PlaneGap << std::endl;
    std::cout <<"  strips per ladder:        " << Nstrips << std::endl;
    std::cout <<"  implant pitch:            " << pitch  << std::endl;
    std::cout <<"  layers thickness:         " << thickness << std::endl;
    std::cout <<"  wafer side:               " << squareSide << std::endl;
    std::cout <<"  total # of ladders:       " << Nladders << std::endl;
    std::cout <<"  Calo planes:              " << NCaloPlanes << std::endl;
    std::cout <<"  planes distance:          " << PlanesDistance << std::endl;
    std::cout <<"  cube width:               " << CubeX << std::endl;
    std::cout <<"  cube height:              " << CubeY << std::endl;
    std::cout <<"  cube depth:               " << CubeZ << std::endl;
    std::cout <<"  rows of cubes:            " << RowsOfCubes << std::endl;
    std::cout <<"  columns of cubes:         " << ColumnsOfCubes << std::endl;
    std::cout <<"  number of scintillators:  " << NScintillators << std::endl;
    std::cout <<"  gap between scintillator: " << LayerGapScintillator << std::endl;

    std::cout <<"=================================" << std::endl;
  }
  
  ClassDef(Geometry, 2)
};


#endif //include guard
