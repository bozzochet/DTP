

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

  messenger_->DeclarePropertyWithUnit("CaloSide", "cm", CaloSide, "DTP");
  messenger_->DeclareProperty("Nsquares", Nsquares, "DTP");
  messenger_->DeclareProperty("Nrows", Nrows, "DTP");
  messenger_->DeclarePropertyWithUnit("PlaneGap", "cm", PlaneGap, "DTP");
  messenger_->DeclareProperty("Nstrips", Nstrips, "DTP");
  messenger_->DeclareProperty("Nsquares", Nsquares, "DTP");
  messenger_->DeclarePropertyWithUnit("pitch", "cm", pitch, "DTP");


  messenger_->DeclareProperty("Nlayers", Nlayers, "layers number");
  messenger_->DeclarePropertyWithUnit("thickness", "cm", thickness, "thickness of layers in cm");
  messenger_->DeclarePropertyWithUnit("LayerGap", "cm", LayerGap, "LayerGap between layers in cm");
  messenger_->DeclarePropertyWithUnit("CaloStkGap", "cm", CaloStkGap, "distance bewteen layers and calorimeter in cm");
  messenger_->DeclareProperty("NCaloPlanes", NCaloPlanes, "number of planes of calorimeter");
  messenger_->DeclarePropertyWithUnit("PlanesDistance", "cm", PlanesDistance, "distance between planes of calorimeter in cm");
  messenger_->DeclarePropertyWithUnit("CubeX", "cm", CubeX, "width of LYSO cube in cm");
  messenger_->DeclarePropertyWithUnit("CubeY", "cm", CubeY, "height of LYSO cube in cm");
  messenger_->DeclarePropertyWithUnit("CubeZ", "cm", CubeZ, "depth of LYSO cube in cm");
  messenger_->DeclareProperty("RowsOfCubes", RowsOfCubes, "number of cubes in line");
  messenger_->DeclareProperty("ColumnsOfCubes", ColumnsOfCubes, "number of cubes in column");
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
  G4Material *scintillator_mat = nist -> FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  
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
  G4double LayerSide = 9.0 * cm;
  G4double N = Nlayers;
  G4double LGADthickness = thickness;
  G4double inter = LayerGap;
  


  //DISTANZA LGAD-CALORIMETRO
  G4double DistanceLC = CaloStkGap;


  //PARAMETRI CALORIMETRO
  G4int planes = NCaloPlanes;
  G4double distanceCaloPlanes = PlanesDistance;
  G4double CaloCubeSideX = CubeX;
  G4double CaloCubeSideY = CubeY;
  G4double CaloCubeSideZ = CubeZ;
  G4int CaloCubeRows = RowsOfCubes;
  G4int CaloCubeColumns = ColumnsOfCubes;
  
  
  //PARAMETRI SCINTILLATORE
  G4int NumberOfScintillatorsBars = NScintillators;
  G4double ScintillatorLayerGap = 0.1 * cm;
  G4double ScintillatorWidth = 9.6 * cm / NumberOfScintillatorsBars; //width of scintillator
  G4double ScintillatorThickness = 0.5 * cm; //thickness of scintillator  
  G4double ScintillatorGap = 0.1 * cm; //gap between sheets of scintillators
  G4double ScintillatorGapZ = LayerGapScintillator;


  //DISTANZA PERCORSA NEL CORSO DELLA SIMULAZIONE
  G4double Distance = 0.0 * cm;


  G4Box *layerscintillator = new G4Box("layerscintillator", LayerSide/2.0, LayerSide/2.0, ScintillatorThickness/2.0);
  G4LogicalVolume *layerscintillatorLogic = new G4LogicalVolume(layerscintillator, scintillator_mat, "layerscintillator");

  //LAYER SCINTILLATORE DAVANTI AI LAYER
  new G4PVPlacement(0,                            // no rotation
                    G4ThreeVector(0,0, ScintillatorThickness/2.0),
                    layerscintillatorLogic,       // its logical volume
                    "layerscintillator",          // its name
                    logicWorld,                   // its mother  volume
                    false,                        // no boolean operation
                    0,                            // copy number
                    fCheckOverlaps);              // checking overlaps
  
  Distance += ScintillatorThickness + ScintillatorGapZ;


  G4Box *siSensor = new G4Box("siSensor", LayerSide/2.0, LayerSide/2.0, LGADthickness/2.0);
  G4LogicalVolume *siSensorLogic = new G4LogicalVolume(siSensor, silicon, "siSensor"); 

  
  for (int i=0; i<N; i++) {
    new G4PVPlacement(0,                          // no rotation
                    G4ThreeVector(0, 0, Distance + LGADthickness/2.0),
                    siSensorLogic,               // its logical volume
                    "siSensor",                  // its name
                    logicWorld,                   // its mother  volume
                    false,                        // no boolean operation
                    0,                            // copy number
                    fCheckOverlaps);              // checking overlaps

    if (i != N-1) { // Questa istruzione fa sì che una volta piazzato l'ultimo layer poi Distance non venga più modificata
      Distance += inter + LGADthickness/2.0;
    }
  }

  Distance += LGADthickness/2.0 + DistanceLC;

  //LAYER SCINTILLATORE DAVANTI AL CALORIMETRO
  new G4PVPlacement(0,                          // no rotation
                    G4ThreeVector(0,0, Distance - (ScintillatorGapZ + ScintillatorThickness/2.0)),
                    layerscintillatorLogic,                     // its logical volume
                    "layerscintillator",                        // its name
                    logicWorld,                   // its mother  volume
                    false,                        // no boolean operation
                    0,                            // copy number
                    fCheckOverlaps);              // checking overlaps


  //LYSO CUBE
  G4Box *calorimeter = new G4Box("calorimeter", CaloCubeSideX/2.0, CaloCubeSideY/2.0, CaloCubeSideZ/2.0); //LYSO CUBE
  G4LogicalVolume *calorimeterLogic = new G4LogicalVolume(calorimeter, LYSO, "calorimeter");


  // In questo caso il posizionamento dei cubetti di LYSO è gestito da più for in maniera tale da differenziare 
  // righe pari o dispari (se pari vanno spostate di CaloCubeSideY/2.0)
  for (int k = 0; k < planes; k++) {
    if(CaloCubeRows % 2 != 0 ) { // Dispari
      for (int i = round(-CaloCubeColumns/2); i <= round(CaloCubeColumns/2); i++) {
        for (int j = round(-CaloCubeRows/2); j <= round(CaloCubeRows/2); j++) {
                    new G4PVPlacement(0,                            // no rotation
                        G4ThreeVector(i * CaloCubeSideX, j * CaloCubeSideY, Distance + CaloCubeSideZ/2.0),
                        calorimeterLogic,                     // its logical volume
                        "calorimeter",                        // its name
                        logicWorld,                        // its mother  volume
                        false,                             // no boolean operation
                        0,                                 // copy number
                        fCheckOverlaps);
          }
        }



    } else if (CaloCubeRows % 2 == 0) { // Pari
      for (int i = round(-CaloCubeColumns/2); i <= round(CaloCubeColumns/2); i++) { 
        for (int j = round(-CaloCubeRows/2); j < round(CaloCubeRows/2); j++) {
                    new G4PVPlacement(0,                            // no rotation
                        G4ThreeVector( i * CaloCubeSideX, j * CaloCubeSideY + CaloCubeSideY/2.0, Distance + CaloCubeSideZ/2.0),
                        calorimeterLogic,                     // its logical volume
                        "calorimeter",                        // its name
                        logicWorld,                        // its mother  volume
                        false,                             // no boolean operation
                        0,                                 // copy number
                        fCheckOverlaps);
        }
      }
     }
  }


  G4double ScintillatorBarsLength = Distance + CaloCubeSideZ + ScintillatorThickness + ScintillatorGapZ;
  
  G4Box *scintillatorBar = new G4Box("scintillatorBar", ScintillatorWidth/2.0, ScintillatorThickness/2.0, ScintillatorBarsLength/2.0);
  G4LogicalVolume *scintillatorBarLogic = new G4LogicalVolume(scintillatorBar, scintillator_mat, "scintillatorBar");
  

  G4RotationMatrix *rotationMatrix = new G4RotationMatrix();
  rotationMatrix -> rotateZ(90. * deg);

  


  //BARRE SCINTILLATORE
  //La lunghezza dei fogli sembra ok
  for (int i = 0; i < NumberOfScintillatorsBars; i++) {     

    //sopra
    new G4PVPlacement(0,                          // no rotation
                    G4ThreeVector(LayerSide/2.0  - ScintillatorWidth/2.0 - i * (ScintillatorGap + ScintillatorWidth), 
                    LayerSide/2.0 + ScintillatorThickness/2.0 + ScintillatorGap,
                    (Distance + CaloCubeSideZ + ScintillatorGapZ + ScintillatorThickness)/2.0),
                    scintillatorBarLogic,                     // its logical volume
                    "scintillatorBar",                        // its name
                    logicWorld,                   // its mother  volume
                    false,                        // no boolean operation
                    0,                            // copy number
                    fCheckOverlaps);              // checking overlaps
    
    //sotto
    new G4PVPlacement(0,                          // no rotation
                    G4ThreeVector(-LayerSide/2.0  + ScintillatorWidth/2.0 + i * (ScintillatorGap + ScintillatorWidth),
                    -LayerSide/2.0 - ScintillatorThickness/2.0 - ScintillatorGap,
                    (Distance + CaloCubeSideZ + ScintillatorGapZ + ScintillatorThickness)/2.0),
                    scintillatorBarLogic,                     // its logical volume
                    "scintillatorBar",                        // its name
                    logicWorld,                   // its mother  volume
                    false,                        // no boolean operation
                    0,                            // copy number
                    fCheckOverlaps);              // checking overlaps  
 
    
    //sinistra 
    new G4PVPlacement(rotationMatrix,                         
                    G4ThreeVector(LayerSide/2.0 + ScintillatorThickness/2.0 + ScintillatorGap,
                    -LayerSide/2.0  + ScintillatorWidth/2.0 + i * (ScintillatorGap + ScintillatorWidth),
                    (Distance + CaloCubeSideZ + ScintillatorGapZ + ScintillatorThickness)/2.0),
                    scintillatorBarLogic,                     // its logical volume
                    "scintillatorBar",                        // its name
                    logicWorld,                   // its mother  volume
                    false,                        // no boolean operation
                    0,                            // copy number
                    fCheckOverlaps);              // checking overlaps
    


    //destra 
    new G4PVPlacement(rotationMatrix,                         
                    G4ThreeVector(- LayerSide/2.0 - ScintillatorThickness/2.0 - ScintillatorGap,
                    LayerSide/2.0  - ScintillatorWidth/2.0 - i * (ScintillatorGap + ScintillatorWidth), 
                    (Distance + CaloCubeSideZ + ScintillatorGapZ + ScintillatorThickness)/2.0),
                    scintillatorBarLogic,                     // its logical volume
                    "scintillatorBar",                        // its name
                    logicWorld,                   // its mother  volume
                    false,                        // no boolean operation
                    0,                            // copy number
                    fCheckOverlaps);              // checking overlaps

 
    

  }

  //LAYER SCINTILLATORE DIETRO AL CALORIMETRO
  new G4PVPlacement(0,                          // no rotation
                    G4ThreeVector(0,0, Distance + ScintillatorGapZ + CaloCubeSideZ + ScintillatorThickness/2.0),
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


bool DetectorConstruction::ExportParameters() {
  bool result = true;
  result = result && ExportRealParameter("CaloSide", CaloSide / cm);
  result = result && ExportIntParameter("Nsquares", Nsquares);
  result = result && ExportIntParameter("Nrows", Nrows);
  result = result && ExportRealParameter("PlaneGap", PlaneGap / cm);
  result = result && ExportIntParameter("Nstrips", Nstrips);
  result = result && ExportRealParameter("pitch", pitch / cm);

  result = result && ExportIntParameter("Nlayers", Nlayers);
  result = result && ExportRealParameter("thickness", thickness / cm);
  result = result && ExportRealParameter("LayerGap", LayerGap / cm);
  result = result && ExportRealParameter("CaloStkGap", CaloStkGap / cm);
  result = result && ExportRealParameter("PlanesDistance", PlanesDistance / cm);
  result = result && ExportIntParameter("NCaloPlanes", NCaloPlanes);
  result = result && ExportRealParameter("CubeX", CubeX / cm);
  result = result && ExportRealParameter("CubeY", CubeY / cm);
  result = result && ExportRealParameter("CubeZ", CubeZ / cm);
  result = result && ExportIntParameter("RowsOfCubes", RowsOfCubes);
  result = result && ExportIntParameter("ColumnsOfCubes", ColumnsOfCubes);
  result = result && ExportIntParameter("NScintillators", NScintillators);
  result = result && ExportRealParameter("LayerGapScintillator", LayerGapScintillator / cm);
  
  
  return result;
}

  /*
  G4double x = 0.0 * cm;
  G4double y = 0.0 * cm;


  
  for (int k = 0; k < planes; k++) {

      
    for (int i = round(-CaloCubeColumns/2); i <= round(CaloCubeColumns/2); i++) { // (int i=-1; i<CaloCubeColumns-1; i++)
        for (int j = round(-CaloCubeRows/2); j <= round (CaloCubeRows/2); j++) { //(int j=-1; j<CaloCubeRows-1; j++)

          if (CaloCubeColumns % 2 != 0)  {
            x = i * CaloCubeSideX; //Dispari
          } else if (CaloCubeColumns % 2 == 0) {
            x = i * CaloCubeSideX + CaloCubeSideX; //Pari
          }
          

          if (CaloCubeRows % 2 != 0)  {
            y = j * CaloCubeSideY; //Dispari
          } else if (CaloCubeRows % 2 == 0) {
            y = j * CaloCubeSideY + CaloCubeSideY; //Pari
          }

          new G4PVPlacement(0,                            // no rotation
                    G4ThreeVector(x, y, Distance + CaloCubeSideZ/2.0), // at (0,0,0)
                    calorimeterLogic,                     // its logical volume
                    "calorimeter",                        // its name
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