#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

// GGS headers
#include "geometry/pluginmanagers/GGSGeoPluginMacros.h"

// Geant4 headers
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4UnionSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Make DetectorConstruction a GGS plugin
GeometryPlugin (DetectorConstruction)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() :
    GGSVGeometryConstruction(), fCheckOverlaps(false), fPhysicalWorld(NULL) {
  DefineMaterials();
  detMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {
  delete detMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials() {
  G4NistManager* man = G4NistManager::Instance();

  G4bool isotopes = false;

  G4Element* Si = man->FindOrBuildElement("Si", isotopes);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct() {
  G4cout << "begin of Detector Construction" << G4endl;
  G4double l = 0.0001 * mm;

  G4NistManager* nist = G4NistManager::Instance();
  G4Material* silicon = nist->FindOrBuildMaterial("G4_Si");
  G4Material* default_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* BGO = nist->FindOrBuildMaterial("G4_BGO");

  G4double world_diameter = 400. * cm;
  G4Sphere* solidWorld = new G4Sphere("World", 0., 0.5 * world_diameter, 0., 2*pi, 0., pi); //For compatibility with VGM [V.F.]

  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,          //its solid
      default_mat,         //its material
      "World");            //its name

  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

  fPhysicalWorld = new G4PVPlacement(0,                     //no rotation
      G4ThreeVector(),       //at (0,0,0)
      logicWorld,            //its logical volume
      "World",               //its name
      0,                     //its mother  volume
      false,                 //no boolean operation
      0,                     //copy number
      fCheckOverlaps);       // checking overlaps

  G4double pad_x = 20  * cm;
  G4double pad_y = 20  * cm;
  G4double pad_z = 0.3 * mm;

  G4double lp = 5*(2*pad_z+2*mm) + 4*(2*cm); //lunghezza pacchetto di silicio

  G4Box* padMother = new G4Box("pad", 0.5 * (pad_x + l), 0.5 * (pad_y + l), 0.5 * (lp + l));
  G4LogicalVolume* padLogic = new G4LogicalVolume(padMother, default_mat, "pad");
  new G4PVPlacement(0,                     //no rotation
    G4ThreeVector(),       //at (0,0,0)
    padLogic,              //its logical volume
    "pad",                 //its name
    logicWorld,            //its mother  volume
    false,                 //no boolean operation
    0,                     //copy number
    fCheckOverlaps);       // checking overlaps

  G4Box* siLayer = new G4Box("siLayer", 0.5 * pad_x, 0.5 * pad_y, 0.5 * pad_z);
  G4LogicalVolume* siLayerLogic = new G4LogicalVolume(siLayer, silicon, "siLayer");

G4int z=0;
for(int i=0; i<5; i++){
  z= i*(2*cm + 2*mm + 0.3*mm);
  new G4PVPlacement(0,G4ThreeVector(0, 0, z),siLayerLogic, "siLayer", padLogic, false, 2*i, fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0, 0, z+2*mm),siLayerLogic, "siLayer", padLogic, false, 2*i+1, fCheckOverlaps);
}

//  G4Box* calorimeter = new G4Box("calorimetro", 60*cm, 60*cm, 60*cm);
//  G4LogicalVolume* calorimeterLogic = new G4LogicalVolume(calorimeter, BGO, "calorimeter");
//  new G4PVPlacement(0,G4ThreeVector(0, 0, 60*cm + lp + 2*cm), calorimeterLogic, "calorimeter", logicWorld, false, 0, fCheckOverlaps);

  //always return the physical World
  return fPhysicalWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void DetectorConstruction::ConstructSDandField() {
  //G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // declare crystal as a MultiFunctionalDetector scorer
  //
  /*
   G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
   G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
   cryst->RegisterPrimitive(primitiv1);
   SetSensitiveDetector("Crystal",cryst);
   */
  //G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
  //G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
  //cryst->RegisterPrimitive(primitiv1);
  //SetSensitiveDetector("Crystal", cryst);
// }

void DetectorConstruction::updateGeometry() {
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}
