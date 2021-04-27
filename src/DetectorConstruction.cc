

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

// AGGIUNTI DOPO

#include "G4Element.hh"
#include "G4Material.hh"
#include "G4NistElementBuilder.hh"
#include "G4NistMaterialBuilder.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Make DetectorConstruction a GGS plugin
GeometryPlugin(DetectorConstruction)

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
    : GGSVGeometryConstruction(), fCheckOverlaps(false), fPhysicalWorld(NULL) {
  DefineMaterials();
  detMessenger = new DetectorMessenger(this);

  messenger_ = new G4GenericMessenger(this, "/Detector/");
  messenger_->DeclareProperty("Nlayers", Nlayers, "layers number");
  messenger_->DeclarePropertyWithUnit("thickness", "cm", thickness, "thickness of layers in cm");
  messenger_->DeclarePropertyWithUnit("interdistance", "cm", interdistance, "interdistance between layers in cm");
  messenger_->DeclarePropertyWithUnit("DistanceLayersCalo", "cm", DistanceLayersCalo, "distance bewteen layers and calorimeter in cm");
  messenger_->DeclareProperty("NCaloPlanes", NCaloPlanes, "number of planes of calorimeter");
  messenger_->DeclarePropertyWithUnit("PlanesDistance", "cm", PlanesDistance, "distance between planes of calorimeter in cm");
  messenger_->DeclarePropertyWithUnit("CubeX", "cm", CubeX, "width of LYSO cube in cm");
  messenger_->DeclarePropertyWithUnit("CubeY", "cm", CubeY, "height of LYSO cube in cm");
  messenger_->DeclarePropertyWithUnit("CubeZ", "cm", CubeZ, "depth of LYSO cube in cm");
  messenger_->DeclareProperty("rows_of_cubes", rows_of_cubes, "number of cubes in line");
  messenger_->DeclareProperty("columns_of_cubes", columns_of_cubes, "number of cubes in column");
  messenger_->DeclareProperty("NScintillators", NScintillators, "number of bars of scintillators");
  messenger_->DeclarePropertyWithUnit("LayerGapScintillator", "cm", LayerGapScintillator, "depth gap betwenn scintillator layers in cm");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() { delete detMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials() {
  G4NistManager *man = G4NistManager::Instance();

  G4bool isotopes = false;

  G4Element *Si = man->FindOrBuildElement("Si", isotopes);
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


  G4NistManager *nist = G4NistManager::Instance();
  G4Material *silicon = nist->FindOrBuildMaterial("G4_Si");
  G4Material *default_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material *BGO = nist->FindOrBuildMaterial("G4_BGO");
  G4Material *scintillator_material = nist -> FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  
  //LYSO ------------------------------------------------------------------------------------------------------------------------
  // https://geant4.web.cern.ch/sites/geant4.web.cern.ch/files/geant4/collaboration/working_groups/geometry/training/D1-Materials.pdf

  G4double a= 0. * g/mole;
  G4double z=0.;
  a=174.967 * g/mole;
  G4Element *elLu= new G4Element("Lutetium","Lu", z = 71., a);
  a=88.90585 * g/mole;
  G4Element *elY= new G4Element("Yttrium","Y", z = 39., a);
  a=28.0855 * g/mole;
  G4Element *elSi= new G4Element("Silicon","Si", z = 14., a);
  a=15.9994 * g/mole;
  G4Element *elO= new G4Element("Oxygen","O", z = 8., a);
  G4double density = 7.10 * g/cm3; //http://www.iem.cfmac.csic.es/departamentos/nuclear/fnexp/master-fusion/lab_master_MCsim.pdf
  G4Material *LYSO = new G4Material("LYSO", density,4);
  LYSO -> AddElement(elLu, 71.45 * perCent);
  LYSO -> AddElement(elY, 4.03 * perCent);
  LYSO -> AddElement(elSi, 6.37 * perCent);
  LYSO -> AddElement(elO, 18.15 * perCent);


  G4double world_diameter = 1000. * cm;
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


  //PARAMETRI LGAD
  G4double LGADSide = 9.0 * cm;
  G4double N = Nlayers;
  G4double LGADthickness = thickness;
  G4double inter = interdistance;
  


  //DISTANZA LGAD-CALORIMETRO
  G4double DistanceLC = DistanceLayersCalo;


  //PARAMETRI CALORIMETRO
  G4int planes = NCaloPlanes;
  G4double distanceCaloPlanes = PlanesDistance;
  G4double CaloSideX = CubeX;
  G4double CaloSideY = CubeY;
  G4double CaloSideZ = CubeZ;
  G4int LYSOrows = rows_of_cubes;
  G4int LYSOcolumns = columns_of_cubes;
  
  
  //PARAMETRI SCINTILLATORE
  G4int number_of_scintillators = NScintillators;
  G4double ScintillatorLayerGap = 0.1 * cm;
  G4double ScintillatorX = 9.6 * cm / number_of_scintillators; //width of scintillator
  G4double ScintillatorY = 0.5 * cm; //thickness of scintillator  
  G4double ScintillatorZ = 0.5 * cm; //length of scintillator
  G4double ScintillatorGap = 0.1 * cm; //gap between sheets of scintillators
  G4double ScintillatorGapZ = LayerGapScintillator;


  //DISTANZA PERCORSA NEL CORSO DELLA SIMULAZIONE
  G4double Distance = 0.0 * cm;


  G4Box *layerscintillator = new G4Box("layerscintillator", LGADSide/2.0, LGADSide/2.0, ScintillatorY/2.0);
  G4LogicalVolume *layerscintillatorLogic = new G4LogicalVolume(layerscintillator, scintillator_material, "layerscintillator");

  //LAYER SCINTILLATORE (DAVANTI AI LAYER)
  new G4PVPlacement(0,                            // no rotation
                    G4ThreeVector(0,0, ScintillatorY/2.0),
                    layerscintillatorLogic,       // its logical volume
                    "layerscintillator",          // its name
                    logicWorld,                   // its mother  volume
                    false,                        // no boolean operation
                    0,                            // copy number
                    fCheckOverlaps);              // checking overlaps
  
  Distance += ScintillatorY + ScintillatorGapZ;


  G4Box *LGADlayer = new G4Box("LGADLayer", LGADSide/2.0, LGADSide/2.0, LGADthickness/2.0);
  G4LogicalVolume *LGADLayerLogic = new G4LogicalVolume(LGADlayer, silicon, "LGADLayer"); 

  
  G4double save_first_layer = 0. * cm; 
  G4double save_last_layer = 0. * cm;


  for (int i=0; i<N; i++) {
    new G4PVPlacement(0,                          // no rotation
                    G4ThreeVector(0, 0, Distance + LGADthickness/2.0),
                    LGADLayerLogic,               // its logical volume
                    "LGADLayer",                  // its name
                    logicWorld,                   // its mother  volume
                    false,                        // no boolean operation
                    0,                            // copy number
                    fCheckOverlaps);              // checking overlaps

    if (i != N-1) { // Questa istruzione fa sì che una volta piazzato l'ultimo layer poi Distance non venga più modificata
      Distance += inter + LGADthickness/2.0;
    }
    if (i == 0) {
      save_first_layer = Distance + LGADthickness/2.0;
    }
    if (i == N-1) {
      save_last_layer = Distance + LGADthickness/2.0;
    }
  }

  // E' necessario inserire metà dello spessore del layer in maniera tale da trovarci sul bordo e distanziare in maniera corretta
  // l'ultimo layer e il calorimetro
  Distance += LGADthickness/2.0 + DistanceLC;

  //LAYER SCINTILLATORE (DAVANTI AL CALORIMETRO)
  new G4PVPlacement(0,                          // no rotation
                    G4ThreeVector(0,0, Distance - (ScintillatorGapZ + ScintillatorY/2.0)),
                    layerscintillatorLogic,                     // its logical volume
                    "layerscintillator",                        // its name
                    logicWorld,                   // its mother  volume
                    false,                        // no boolean operation
                    0,                            // copy number
                    fCheckOverlaps);              // checking overlaps


  //LYSO CUBE
  G4Box *LYSOcube = new G4Box("LYSOcube", CaloSideX/2.0, CaloSideY/2.0, CaloSideZ/2.0); //LYSO CUBE
  G4LogicalVolume *LYSOcubeLogic = new G4LogicalVolume(LYSOcube, LYSO, "LYSOcube");

  for (int k = 0; k < planes; k++) {
    if(LYSOrows % 2 != 0 ) { // Dispari
      for (int i = round(-LYSOcolumns/2); i <= round(LYSOcolumns/2); i++) {
        for (int j = round(-LYSOrows/2); j <= round(LYSOrows/2); j++) {
                    new G4PVPlacement(0,                            // no rotation
                        G4ThreeVector(i * CaloSideX, j * CaloSideY, Distance + CaloSideZ/2.0),
                        LYSOcubeLogic,                     // its logical volume
                        "LYSOcube",                        // its name
                        logicWorld,                        // its mother  volume
                        false,                             // no boolean operation
                        0,                                 // copy number
                        fCheckOverlaps);
          }
        }



    } else if (LYSOrows % 2 == 0) { // Pari
      for (int i = round(-LYSOcolumns/2); i <= round(LYSOcolumns/2); i++) { 
        for (int j = round(-LYSOrows/2); j < round(LYSOrows/2); j++) {
                    new G4PVPlacement(0,                            // no rotation
                        G4ThreeVector( i * CaloSideX, j * CaloSideY + CaloSideY/2.0, Distance + CaloSideZ/2.0), // at (0,0,0)
                        LYSOcubeLogic,                     // its logical volume
                        "LYSOcube",                        // its name
                        logicWorld,                        // its mother  volume
                        false,                             // no boolean operation
                        0,                                 // copy number
                        fCheckOverlaps);
        }
      }
     }
  }


  G4double longScintillatorZ =Distance + CaloSideZ - save_first_layer  + 1 * cm;
  
  G4Box *longscintillator = new G4Box("longscintillator", ScintillatorX/2.0, ScintillatorY/2.0, longScintillatorZ/2.0);
  G4LogicalVolume *longscintillatorLogic = new G4LogicalVolume(longscintillator, scintillator_material, "longscintillator");
  

  G4RotationMatrix *rotationMatrix = new G4RotationMatrix();
  rotationMatrix -> rotateZ(90. * deg);

  


  //BARRE SCINTILLATORE
  for (int i = 0; i<number_of_scintillators; i++) {     

    //sopra
    new G4PVPlacement(0,                          // no rotation
                    G4ThreeVector(LGADSide/2.0  - ScintillatorX/2.0 - i * (ScintillatorGap + ScintillatorX),
                    LGADSide/2.0 + ScintillatorY/2.0 + ScintillatorGap,
                    (Distance + CaloSideZ + (save_first_layer - LGADthickness/2.0))/2.0),
                    longscintillatorLogic,                     // its logical volume
                    "longscintillator",                        // its name
                    logicWorld,                   // its mother  volume
                    false,                        // no boolean operation
                    0,                            // copy number
                    fCheckOverlaps);              // checking overlaps
    
    //sotto
    new G4PVPlacement(0,                          // no rotation
                    G4ThreeVector(-LGADSide/2.0  + ScintillatorX/2.0 + i * (ScintillatorGap + ScintillatorX),
                    -LGADSide/2.0 - ScintillatorY/2.0 - ScintillatorGap,
                    (Distance + CaloSideZ + (save_first_layer - LGADthickness/2.0))/2.0),
                    longscintillatorLogic,                     // its logical volume
                    "longscintillator",                        // its name
                    logicWorld,                   // its mother  volume
                    false,                        // no boolean operation
                    0,                            // copy number
                    fCheckOverlaps);              // checking overlaps  
 
    
    //sinistra 
    new G4PVPlacement(rotationMatrix,                         
                    G4ThreeVector(LGADSide/2.0 + ScintillatorY/2.0 + ScintillatorGap,
                    -LGADSide/2.0  + ScintillatorX/2.0 + i * (ScintillatorGap + ScintillatorX),
                    (Distance + CaloSideZ + (save_first_layer - LGADthickness/2.0))/2.0),
                    longscintillatorLogic,                     // its logical volume
                    "longscintillator",                        // its name
                    logicWorld,                   // its mother  volume
                    false,                        // no boolean operation
                    0,                            // copy number
                    fCheckOverlaps);              // checking overlaps
    


    //destra 
    new G4PVPlacement(rotationMatrix,                         
                    G4ThreeVector(- LGADSide/2.0 - ScintillatorY/2.0 - ScintillatorGap,
                    LGADSide/2.0  - ScintillatorX/2.0 - i * (ScintillatorGap + ScintillatorX), 
                    (Distance + CaloSideZ + (save_first_layer - LGADthickness/2.0))/2.0),
                    longscintillatorLogic,                     // its logical volume
                    "longscintillator",                        // its name
                    logicWorld,                   // its mother  volume
                    false,                        // no boolean operation
                    0,                            // copy number
                    fCheckOverlaps);              // checking overlaps

 
    

  }

  //LAYER SCINTILLATORE (DIETRO AL CALORIMETRO)
  new G4PVPlacement(0,                          // no rotation
                    G4ThreeVector(0,0, Distance + ScintillatorGapZ + CaloSideZ + ScintillatorY/2.0),
                    layerscintillatorLogic,                     // its logical volume
                    "layerscintillator",                        // its name
                    logicWorld,                   // its mother  volume
                    false,                        // no boolean operation
                    0,                            // copy number
                    fCheckOverlaps);              // checking overlaps



  // always return the physical World
  return fPhysicalWorld;
}

void DetectorConstruction::updateGeometry() { G4RunManager::GetRunManager()->DefineWorldVolume(Construct()); }

// QUESTO E' IMPORTANTE SU Geometry.h
// RICORDATI DI AGGIORNARE CON LE VARIABILI DEL GEO.MAC
bool DetectorConstruction::ExportParameters() {
  bool result = true;
  result = result && ExportIntParameter("Nlayers", Nlayers);
  result = result && ExportRealParameter("thickness", thickness / cm);
  result = result && ExportRealParameter("interdistance", interdistance / cm);
  result = result && ExportRealParameter("DistanceLayersCalo", DistanceLayersCalo / cm);
  result = result && ExportIntParameter("NCaloPlanes", NCaloPlanes);
  result = result && ExportRealParameter("CubeX", CubeX / cm);
  result = result && ExportRealParameter("CubeY", CubeY / cm);
  result = result && ExportRealParameter("CubeZ", CubeZ / cm);
  result = result && ExportIntParameter("rows_of_cubes", rows_of_cubes);
  result = result && ExportIntParameter("columns_of_cubes", columns_of_cubes);
  result = result && ExportIntParameter("NScintillators", NScintillators);
  result = result && ExportRealParameter("LayerGapScintillator", LayerGapScintillator / cm);
  
  
  return result;
}

  /*
  G4double x = 0.0 * cm;
  G4double y = 0.0 * cm;


  
  for (int k = 0; k < planes; k++) {

      
    for (int i = round(-LYSOcolumns/2); i <= round(LYSOcolumns/2); i++) { // (int i=-1; i<LYSOcolumns-1; i++)
        for (int j = round(-LYSOrows/2); j <= round (LYSOrows/2); j++) { //(int j=-1; j<LYSOrows-1; j++)

          if (LYSOcolumns % 2 != 0)  {
            x = i * CaloSideX; //Dispari
          } else if (LYSOcolumns % 2 == 0) {
            x = i * CaloSideX + CaloSideX; //Pari
          }
          

          if (LYSOrows % 2 != 0)  {
            y = j * CaloSideY; //Dispari
          } else if (LYSOrows % 2 == 0) {
            y = j * CaloSideY + CaloSideY; //Pari
          }

          new G4PVPlacement(0,                            // no rotation
                    G4ThreeVector(x, y, Distance + CaloSideZ/2.0), // at (0,0,0)
                    LYSOcubeLogic,                     // its logical volume
                    "LYSOcube",                        // its name
                    logicWorld,                        // its mother  volume
                    false,                             // no boolean operation
                    0,                                 // copy number
                    fCheckOverlaps);
        }
      }

    if (k != planes - 1) { // Questa istruzione fa sì che Distance non venga modificata una volta piazzato l'ultimo piano
      Distance += CubeZ + distanceCaloPlanes;
    }
  }
  */