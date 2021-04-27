/*
 * DetectorConstruction.h
 *
 *  Created on: 22 Feb 2017
 *      Authors: Junjing Wang, Ming Xu, Zheng Quan
 *      Adapted by: Nicola Mori
 */

#ifndef DETECTORCONSTRUCTION_HH_
#define DETECTORCONSTRUCTION_HH_

// GGS headers
#include "geometry/GGSVGeometryConstruction.h"

#include "G4SystemOfUnits.hh"

// Geant4 headers
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class DetectorMessenger;

class G4GenericMessenger;

class DetectorConstruction: public GGSVGeometryConstruction {
public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

  //bool ExportParameters();
  //const std::string GetVersion();

public:
  virtual G4VPhysicalVolume* Construct();
  //virtual void ConstructSDandField();
  bool ExportParameters();

public:
  void updateGeometry();

  /*! @brief Return the pointer to physical world */
  G4VPhysicalVolume* GetVolume() {
    return fPhysicalWorld;
  }
private:
  void DefineMaterials();

  G4bool fCheckOverlaps;

  DetectorMessenger* detMessenger;

  G4VPhysicalVolume* fPhysicalWorld;

  G4GenericMessenger *messenger_;


  //GEOMETRIC DEFAULT VALUES

  G4int Nlayers = 40;                          //silicon layers
  G4double thickness = 0.015 * cm;             //layers thickness
  G4double interdistance = 0.1512820513 * cm;  //distance between layers in cm
  G4double DistanceLayersCalo = 10 * cm;       //distance between layers and calorimeter
  G4int NCaloPlanes = 1;    //calorimeter planes
  G4double PlanesDistance = 3 * cm; //distance between calo planes
  G4double CubeX = 3 * cm;  //width of cube
  G4double CubeY = 3 * cm;  //height of cube
  G4double CubeZ = 3 * cm;  //depth of cube
  G4int rows_of_cubes = 3; //number of cubes in line
  G4int columns_of_cubes = 3; //number of cubes in column
  G4int NScintillators = 4; //number of bars of scintillators
  G4double LayerGapScintillator = 0.1 * cm; //depth gap betwenn scintillator layers in cm

};

#endif /* DETECTORCONSTRUCTION_HH_ */
