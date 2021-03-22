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

  G4double CaloSide = 60;         // Calo side in cm
  G4double CaloStkGap = 2;        // STK-Calo ditance in cm
  G4int Nsquares = 8;             // squares per side on layer
  G4int Nrows = 2;                // number of ladders per 'column' (i.e. each ladders has Nsquares/Nrows wafers)
  G4int Nlayers = 10;             // number of layers (every two layers, i.e. a plane, one layer would be X and one would be Y)
  G4double LayerGap = 0.2;        // gap between two ladders of the same plane in cm
  G4double PlaneGap = 2;          // gap between planes in cm
  G4int Nstrips = 640;            //strips per ladder
  G4double pitch = 0.015 * cm;    //strips pitch in cm
  G4double thickness = 0.03 * cm;  //thickness of layers in cm

};

#endif /* DETECTORCONSTRUCTION_HH_ */
