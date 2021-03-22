#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

// GGS headers
#include "geometry/pluginmanagers/GGSGeoPluginMacros.h"

// Geant4 headers
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Sphere.hh"
#include "G4Transform3D.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VisAttributes.hh"
#include "G4UImanager.hh"
#include "G4GenericMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Make DetectorConstruction a GGS plugin
GeometryPlugin(DetectorConstruction)

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
    : GGSVGeometryConstruction(), fCheckOverlaps(false), fPhysicalWorld(NULL) {
  DefineMaterials();
  detMessenger = new DetectorMessenger(this);

  messenger_ = new G4GenericMessenger(this, "/Detector/");
  
  // messenger_->DeclareProperty("CaloSide", CaloSide, "Calo side in cm").SetUnit("cm");
  // messenger_->DeclareProperty("CaloStkGap", CaloStkGap, "STK-Calo distance in cm").SetUnit("cm");
  // messenger_->DeclareProperty("Nsquares", Nsquares, "wafers per side on layer");
  // messenger_->DeclareProperty("Nrows", Nrows, "number of ladders per column");
  // messenger_->DeclareProperty("Nlayers", Nlayers, "layer number");
  // messenger_->DeclareProperty("LayerGap", LayerGap, "gap between 2 ladders on the same plane in cm").SetUnit("cm");
  // messenger_->DeclareProperty("PlaneGap", PlaneGap, "gap between planes in cm").SetUnit("cm");
  // messenger_->DeclareProperty("Nstrips", Nstrips, "strips per ladder");
  // messenger_->DeclareProperty("pitch", pitch, "strips pitch in cm").SetUnit("cm");
  // messenger_->DeclareProperty("thickness", thickness, "thickness of layers in cm").SetUnit("cm");
  
  messenger_->DeclarePropertyWithUnit("CaloSide", "cm", CaloSide, "Calo side in cm");
  messenger_->DeclarePropertyWithUnit("CaloStkGap", "cm", CaloStkGap, "STK-Calo distance in cm");
  messenger_->DeclareProperty("Nsquares", Nsquares, "wafers per side on layer");
  messenger_->DeclareProperty("Nrows", Nrows, "number of ladders per column");
  messenger_->DeclareProperty("Nlayers", Nlayers, "layer number");
  messenger_->DeclarePropertyWithUnit("LayerGap", "cm", LayerGap, "gap between 2 ladders on the same plane in cm");
  messenger_->DeclarePropertyWithUnit("PlaneGap", "cm", PlaneGap, "gap between planes in cm");
  messenger_->DeclareProperty("Nstrips", Nstrips, "strips per ladder");
  messenger_->DeclarePropertyWithUnit("pitch", "cm", pitch, "strips pitch in cm");
  messenger_->DeclarePropertyWithUnit("thickness", "cm", thickness, "thickness of layers in cm");

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() { delete detMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials() {
  // G4NistManager *man = G4NistManager::Instance();

  // G4bool isotopes = false;

  // G4Element *Si = man->FindOrBuildElement("Si", isotopes);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct() {
  G4cout << "begin of Detector Construction" << G4endl;

  //Run geometry configuration script before building
  if (_geoDataCard != "") {
    G4UImanager::GetUIpointer()->ApplyCommand(G4String("/control/execute " + _geoDataCard));
  }

  // Delete the messenger so that the commands for configuring the geometry won't be available anymore
  delete messenger_;
  messenger_ = NULL;

  
  //  G4double l = 0.0001 * mm;//0.1 um? MD: maybe just to avoid overlaps
  G4double l = 0.0;

  G4NistManager *nist = G4NistManager::Instance();
  G4Material *silicon = nist->FindOrBuildMaterial("G4_Si");
  G4Material *default_mat = nist->FindOrBuildMaterial("G4_Galactic");
#ifndef _NOCALO_
  G4Material *BGO = nist->FindOrBuildMaterial("G4_BGO");
#endif

  G4double world_diameter = 400. * cm;
  G4Sphere *solidWorld =
      new G4Sphere("World", 0., 0.5 * world_diameter, 0., 2 * pi, 0., pi); // For compatibility with VGM [V.F.]

  G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld,  // its solid
                                                    default_mat, // its material
                                                    "World");    // its name

  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

  fPhysicalWorld = new G4PVPlacement(0,                      // no rotation
                                     G4ThreeVector(0, 0, 0), // at (0,0,0)
                                     logicWorld,             // its logical volume
                                     "World",                // its name
                                     0,                      // its mother  volume
                                     false,                  // no boolean operation
                                     0,                      // copy number
                                     fCheckOverlaps);        // checking overlaps

  G4int N = Nsquares; //e.g 8
  G4double dim = pitch * Nstrips; //e.g. 9.6 cm
  G4double pad_x = N * dim; //e.g. 8*9.6 cm
  G4double pad_y = N * dim; //e.g. 8*9.6 cm
  G4double pad_z = thickness; //e.g 300 um

  //  printf("%f %f %f\n", pad_x, pad_y, pad_z);

  G4double StkDepth = Nlayers * (pad_z + LayerGap) + (Nlayers/2.0-1) * PlaneGap; // profondità tracker
  G4Box *padMother = new G4Box("pad", 0.5 * (pad_x + l), 0.5 * (pad_y + l), 0.5 * (StkDepth + l));
  G4LogicalVolume *padLogic = new G4LogicalVolume(padMother, default_mat, "pad");//sort of box containing the whole Stk
  new G4PVPlacement(0,                            // no rotation
   		    G4ThreeVector(0, 0, -StkDepth/2.0), // at (0,0,depth/2)
   		    // = its "center" is at depth/2
   		    // --> the top face is touching (0,0,0)
   		    padLogic,                     // its logical volume
   		    "pad",                        // its name
   		    logicWorld,                   // its mother  volume
   		    false,                        // no boolean operation
   		    0,                            // copy number
   		    fCheckOverlaps);              // checking overlaps

  G4Box *siLayer = new G4Box("siLayer", 0.5 * pad_x, 0.5 * pad_y, 0.5 * pad_z);
  G4LogicalVolume *siLayerLogic = new G4LogicalVolume(siLayer, silicon, "siLayer");
  
  G4Box *siLadder = new G4Box("siLadder", 0.5 * dim, 0.5 * pad_y, 0.5 * pad_z);
  G4LogicalVolume *siLadderLogic = new G4LogicalVolume(siLadder, silicon, "siLadder");
  
  for (int i = 0; i < N; i++) { //8 wafers per ladder but also 8 ladders side-by-side
    //    printf("%f, %f, %f\n", i * dim - 0.5 * (pad_x - dim), 0.0, 0.0),
    new G4PVPlacement(0,
  		      G4ThreeVector(i * dim - 0.5 * (pad_x - dim), 0.0, 0.0),
  		      siLadderLogic, //current
  		      "siLadder", //name
  		      siLayerLogic, //mother
  		      false, //many
  		      i, //copy number
  		      fCheckOverlaps);
  } //MD: non capisco perchè questi layer non vengono disegnati nella geometria Leonard
  // probabilmente perchè solamente gli N layer sono messi nella box che era messa nel World
  //  new G4PVReplica("siLadderp", siLadderLogic, siLayerLogic, kXAxis, N, dim);//should be alternative to the loop above
  
  G4Box *siSensor = new G4Box("siSensor", 0.5 * dim, 0.5 * dim, 0.5 * pad_z);
  G4LogicalVolume *siSensorLogic = new G4LogicalVolume(siSensor, silicon, "siSensor");  
  new G4PVReplica("siSensorp", siSensorLogic, siLadderLogic, kYAxis, N, dim);//the ladder is filled of wafers

  G4RotationMatrix* rotationMatrixX = new G4RotationMatrix();
  G4RotationMatrix* rotationMatrixY = new G4RotationMatrix();
  rotationMatrixY->rotateZ(90.*deg);
  
  G4double z = -StkDepth/2.0;
  for (int i = 0; i < Nlayers; i++) {
    G4RotationMatrix* rotationMatrix = rotationMatrixX;
    if (i%2) rotationMatrix = rotationMatrixY;
    new G4PVPlacement(rotationMatrix,
		      G4ThreeVector(0, 0, z),
		      siLayerLogic,
		      "siLayer",
		      padLogic,
		      false,
		      i,
		      fCheckOverlaps);
    if (!(i%2)) z += LayerGap;
    else z += PlaneGap;
    z += thickness;
  }
  
#ifndef _NOCALO_
  //  printf("%f\n", CaloSide);
  G4Box *calorimeter = new G4Box("calorimetro", CaloSide/2.0, CaloSide/2.0, CaloSide/2.0);
  G4LogicalVolume *calorimeterLogic = new G4LogicalVolume(calorimeter, BGO, "calorimeter");
  new G4PVPlacement(0, G4ThreeVector(0, 0, -1.0 * (CaloSide/2.0 + StkDepth + CaloStkGap)), calorimeterLogic, "calorimeter",
   		    logicWorld, false, 0, fCheckOverlaps);
#endif
  
  // always return the physical World
  return fPhysicalWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void DetectorConstruction::ConstructSDandField() {
// G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

// declare crystal as a MultiFunctionalDetector scorer
//
/*
 G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
 G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
 cryst->RegisterPrimitive(primitiv1);
 SetSensitiveDetector("Crystal",cryst);
 */
// G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
// G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
// cryst->RegisterPrimitive(primitiv1);
// SetSensitiveDetector("Crystal", cryst);
// }

void DetectorConstruction::updateGeometry() { G4RunManager::GetRunManager()->DefineWorldVolume(Construct()); }


bool DetectorConstruction::ExportParameters() {
  bool result = true;
  result = result && ExportIntParameter("CaloSide", CaloSide / cm);
  result = result && ExportIntParameter("CaloStkGap", CaloStkGap / cm);
  result = result && ExportIntParameter("Nsquares", Nsquares);
  result = result && ExportIntParameter("Nrows", Nrows);
  result = result && ExportIntParameter("Nlayers", Nlayers);
  result = result && ExportIntParameter("LayerGap", LayerGap / cm);
  result = result && ExportIntParameter("PlaneGap", PlaneGap / cm);
  result = result && ExportIntParameter("Nstrips", Nstrips);
  result = result && ExportRealParameter("pitch", pitch / cm);
  result = result && ExportRealParameter("thickness", thickness / cm);
  return result;
}
