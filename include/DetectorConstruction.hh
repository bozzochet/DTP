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

// VARIABILI DTP
  G4double CaloSide = 60 ;
  G4int Nrows = 1 ;
  G4int Nsquares = 1 ;
  G4double PlaneGap = 2 ; // Non so se questa ha bisogno di un altro valore

//////////////////////////////////////////////////////////////////////
  G4int Nstrips = 576 ;  //Queste due si possono mettere anche su comuni ma per ora le lasciamo qui
  G4double pitch = 0.015625 * cm ; //156.25 um
//////////////////////////////////////////////////////////////////////
  
  
  // VARIABILI SLA
  G4int NCaloPlanes = 1 ;
  G4double PlanesDistance = 3.0 * cm ;
  G4double CubeX = 3.0 * cm ;
  G4double CubeY = 3.0 * cm ;
  G4double CubeZ = 3.0 * cm ;
  G4int RowsOfCubes = 3 ;
  G4int ColumnsOfCubes = 3 ;
  G4int NScintillators = 4 ;
  G4double LayerGapScintillator = 0.0 * cm ;

  // VARIABILI COMUNI
  G4int Nlayers = 40 ;
  G4double thickness = 0.015 * cm ;
  G4double LayerGap = 0.1512820513 * cm ;
  G4double CaloStkGap = 10.0 * cm ;

};

#endif /* DETECTORCONSTRUCTION_HH_ */
