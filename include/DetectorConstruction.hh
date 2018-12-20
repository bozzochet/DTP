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

// Geant4 headers
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class DetectorMessenger;

class DetectorConstruction: public GGSVGeometryConstruction {
public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

public:
  virtual G4VPhysicalVolume* Construct();
  // virtual void ConstructSDandField();

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
};

#endif /* DETECTORCONSTRUCTION_HH_ */
