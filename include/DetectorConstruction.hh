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

  G4int Nlayers = 10;
  G4int Nsquares = 8;             //squares per side on layer
  G4int Nrows = 2;                //rows of ladders on a layer
  G4int Nstrips = 640;            //strips per ladder
  G4double pitch = 0.015 * cm;    //strips pitch in cm
  G4double thickness = 0.3 * mm;  //thickness of layers in mm

};

#endif /* DETECTORCONSTRUCTION_HH_ */
