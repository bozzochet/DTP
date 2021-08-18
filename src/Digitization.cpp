#include "TrCluster.hh"
#include "physics.h"
#include "info.h"
#include "measure.h"
#include "Geometry.h"
#include "TimeSegm.h"
#include "TimeSim.h"
#include "PosSim.h"

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGraph.h"

#include "montecarlo/readers/GGSTHadrIntReader.h"
#include "montecarlo/readers/GGSTHitsReader.h"
#include "montecarlo/readers/GGSTMCTruthReader.h"
#include "montecarlo/readers/GGSTRootReader.h"
#include "utils/GGSSmartLog.h"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

float eps = std::numeric_limits<float>::epsilon(); // Return the machine epsilon -> difference between 1.0 ancd the next 
// value representable by the floating point type T

#define COSTANT_FRACTION_FOR_TIMING 0.10 // Questo credo sia l'avanzamento temporale della simulazione ->
// -> avanza di 0.10s a iterazione?
//#define COSTANT_FRACTION_FOR_TIMING 0.15



// Includere? ------------------------------------------------------------------------------------------------
int findStrip
(double position, int Nsquares, int Nstrips, double pitch)
{
  int stripHit = (position+(Nsquares*Nstrips*pitch*0.5))/pitch;
  return stripHit % (int(Nstrips));
}

// Includere? ------------------------------------------------------------------------------------------------
int findLadder
(int ID, int layer, bool condition, int Nsquares)
{
  if(!condition) {
    int row = (ID - (pow(Nsquares,2)*layer))/Nsquares;
    bool leftOrRight = ID % (int(Nsquares)) > (Nsquares*0.5);
    return row + (Nsquares*leftOrRight) + (Nsquares*2*layer);
  }
  else
    return (ID % ((int(Nsquares))*2)) + (layer*Nsquares*2);
}


int fillGeoTree(TTree*, TDirectory*, Geometry*);

int fillEvTree(GGSTRootReader&, TTree*, TDirectory*, Geometry*);

//int fillCaloTree(GGSTRootReader&, TTree*, TDirectory*);

//int fillMeasTree(TTree*, TTree*, TDirectory*, Geometry*);



int main(int argc, char **argv)
{
  
  static const string routineName("Digitization::main");
  GGSSmartLog::verboseLevel = GGSSmartLog::INFO; // Print only INFO messages or more important

  if (argc<3) {
    COUT(ERROR) << "You need to pass at least an input file name and an output one" << ENDL;
    COUT(ERROR) << argv[0] << " [output file name] <input file name 1> <input file name 2> ..." << ENDL;
    return -1;
  }

  TString outFileName = argv[1];

  // Create the reader container and open the data file
  GGSTRootReader reader;

  int shift=2;
  for (int ii=shift; ii<argc; ii++) {
    TString inputFileName = argv[ii];
    //    printf("%d) %s\n", ii, argv[ii]);

    static GGSTFilesHandler* handler = NULL;
    handler = reader.Open(inputFileName, handler);
    
    if (!handler) {
      COUT(ERROR) << "Cannot open input file " << inputFileName << ENDL;
      return 1;
    }
  }

  //create output file
  TFile *outFile = TFile::Open(outFileName, "RECREATE");
  if (!outFile || outFile->IsZombie()) {
    COUT(ERROR) << "Cannot create digitization output file " << outFileName << ENDL;
    return 1;
  }

  // geometry


  const GGSTGeoParams *GEO = reader.GetGeoParams();
  Geometry *geo = new Geometry();

  geo->CaloSide = GEO->GetRealGeoParam("CaloSide");
  geo->LayerGap = GEO->GetRealGeoParam("LayerGap");
  geo->Nsquares = GEO->GetIntGeoParam("Nsquares");
  geo->Nrows = GEO->GetIntGeoParam("Nrows");
  geo->Nlayers = GEO->GetIntGeoParam("Nlayers");
  geo->CaloStkGap = GEO->GetRealGeoParam("CaloStkGap");
  geo->PlaneGap = GEO->GetRealGeoParam("PlaneGap");
  geo->Nstrips = GEO->GetIntGeoParam("Nstrips");
  geo->pitch = GEO->GetRealGeoParam("pitch");
  geo->thickness = GEO->GetRealGeoParam("thickness");
  geo->NCaloPlanes = GEO->GetIntGeoParam("NCaloPlanes");
  geo->PlanesDistance = GEO->GetRealGeoParam("PlanesDistance");
  geo->CubeX = GEO->GetRealGeoParam("CubeX");
  geo->CubeY = GEO->GetRealGeoParam("CubeY");
  geo->CubeZ = GEO->GetRealGeoParam("CubeZ");
  geo->RowsOfCubes = GEO->GetIntGeoParam("RowsOfCubes");
  geo->ColumnsOfCubes = GEO->GetIntGeoParam("ColumnsOfCubes");
  geo->NScintillators = GEO->GetIntGeoParam("NScintillators");
  geo->LayerGapScintillator = GEO->GetRealGeoParam("LayerGapScintillator");

  geo->ComputeDerived();
  geo->Dump();

  //trees

  outFile->cd();
  TTree *geo_tree = new TTree("geometry", "siSensorGeoParams");
  TTree *events_tree = new TTree("events", "siSensorHits");
  events_tree->SetAutoSave(-3000000); 
  events_tree->SetAutoFlush(-3000000); 


  // PER ORA ESCLUSI PER CAPIRE COME E SE FUNZIONA

  /*
  TTree *meas_tree = new TTree("measures", "siSensorMeasures");
  meas_tree->SetAutoSave(-3000000);
  meas_tree->SetAutoFlush(-3000000);
  TTree *calo_tree = new TTree("calorimeter", "calorimeterEdep");
  calo_tree->SetAutoSave(-3000000);
  calo_tree->SetAutoFlush(-3000000);
  */

 /*
  if
    (
     fillGeoTree(geo_tree, outFile, geo) != 0
     || fillEvTree(reader, events_tree, outFile, geo) != 0
     || fillCaloTree(reader, calo_tree, outFile) != 0
     || fillMeasTree(events_tree, meas_tree, outFile, geo) != 0
     )
    return 1;
  */
  // Save histograms
  outFile->cd();

  outFile->Close();
  delete outFile;


  COUT(INFO) <<ENDL;
  COUT(INFO) <<"Data saved in " <<outFileName <<ENDL;

  return 0;
}



int fillGeoTree(TTree *geo_tree, TDirectory* outFile, Geometry *geo)
{
  static const string routineName("Digitization::fillGeoTree");


  COUT(INFO) <<ENDL;
  COUT(INFO) <<"Saving geometric parameters..." <<ENDL;

  // i valori sono ottenuti da geo che è stata definita precedentemente a riga 111
  geo_tree->Branch("Geometry", &geo);
  
  outFile->cd();
  geo_tree->Fill();

  outFile->cd();
  geo_tree->Write();
  
  return 0;
}


// In teoria questo dovrebbe solamente valutare gli eventi che immagino siano le hit
int fillEvTree
(GGSTRootReader &reader, TTree *events_tree, TDirectory* outFile, Geometry *geo)
{
  static const string routineName("Digitization::fillEvTree"); //
  
  
  COUT(INFO) <<ENDL;
  COUT(INFO) <<"Saving MC truth..." <<ENDL;

  TClonesArray a("TrCluster", 200); // An array of clone objects
  // In particolare abbiamo che TrCluster viene clonato 200 volte?
  events_tree->Branch("Events", &a); // Viene creata una branch su events_tree e vengono copiati gli array di a


  // Create and retrieve the hits sub-reader
  GGSTHitsReader *hReader = reader.GetReader<GGSTHitsReader>();

  // Set which hit detectors are to be read
  // The name is the same of the sensitive logical volume name in the simulation
  hReader->SetDetector("siSensor", kTRUE); // Nome su SLA uguale a quello di DTP, qui non c'è bisogno di cambiare nulla

  // Retrieve the MC truth sub-reader
  GGSTMCTruthReader *mcReader = NULL;
  mcReader = reader.GetReader<GGSTMCTruthReader>();

  // Create and retrieve the hadronic interaction sub-reader
  GGSTHadrIntReader *hadrReader = reader.GetReader<GGSTHadrIntReader>();


  COUT(INFO) << "Begin loop over " << reader.GetEntries() << " events" << ENDL;

  std::clock_t start = std::clock();

  for (int iEv = 0; iEv < reader.GetEntries(); iEv++) {

    //print and update progress bar
    info::progress(start, iEv, reader.GetEntries());

    reader.GetEntry(iEv); // Reads all the data objects whose sub-readers have already been created
    GGSTPartHit *phit;
    GGSTIntHit *inthit;
    int nHits = hReader->GetNHits("siSensor"); // Number of hit siLayers for current event
    a.Clear();

    GGSTHadrIntInfo *intInfo = NULL;
    intInfo = hadrReader->GetInelastic();

    // Hits loop
    int ncluster = 0;
    
    for (int iHit = 0; iHit < nHits; iHit++) {
      inthit = (GGSTIntHit *)hReader->GetHit("siSensor", iHit);
      int nPHit = inthit->GetNPartHits(); // Questa riga salva le hit dei fotoni?
      int llayer = inthit->GetVolumeID() / (geo->Nsquares * geo->Nsquares); // Questa riga ci dice in quale layer ha convertito?
      //cout<<"layer : "<<llayer<<endl;
      //cout<<"square : "<<inthit->GetVolumeID()<<endl;
      
      for (int i = 0; i < nPHit; i++) {
        TrCluster *c = (TrCluster *)a.ConstructedAt(ncluster++);
        phit = (GGSTPartHit *)inthit->GetPartHit(i);
        //cout<<"Entry #"<<iHit+i+1<<endl;
	
	
        //fill with MC truth

        c->xy = llayer % 2 == 0; 
        c->pos[0] = 1e-2 * 0.5 * (phit->entrancePoint[0] + phit->exitPoint[0]); // x coordinate
        c->pos[1] = 1e-2 * 0.5 * (phit->entrancePoint[1] + phit->exitPoint[1]); // y coordinate
        c->pos[2] = 1e-2 * 0.5 * (phit->entrancePoint[2] + phit->exitPoint[2]); // z coordinate
        c->time = 1e-9 * phit->time; // Questo credo sia escludibile
        c->eDep = 1e+9 * phit->eDep; // Energia depositata
        //c->spRes = 0.00007;
        c->layer = llayer; // layer dove ha hittato?
        c->parID = phit->parentID;
        c->parPdg = phit->particlePdg;
        c->ID = inthit->GetVolumeID();
	
        //Find the nearest strip to the left
        c->strip = findStrip(c->pos[c->xy], geo->Nsquares, geo->Nstrips, 1e-2*geo->pitch);
	
        // Find the ladder it belongs
        c->ladder = findLadder(c->ID, c->layer, c->xy, geo->Nsquares);
	
        //Deposit energy and create cluster
        double thisPos = ((c->ladder%geo->Nsquares)*geo->squareSide) + (c->strip*(1e-2*geo->pitch)) - (geo->Nsquares*geo->squareSide*0.5);
        double fraction = (c->pos[c->xy]-thisPos)/(1e-2*geo->pitch);
	
        c->clust[0] = c->eDep * (1-fraction);
        c->clust[1] = c->eDep * (fraction);
	
	
        //errors
	
        int error = 0;
	
        if (
	        c->strip < 0 || c->strip >= geo->Nstrips || c->ladder < 0 || c->ladder >= geo->Nladders * geo->Nlayers
	    ) {
          COUT(ERROR) <<ENDL <<"ladder or strip out of range: lad:" << c->ladder << " str:" << c->strip << ENDL;  
          error = 1;
        }
	
        if ( c->pos[0] < -geo->Nsquares*0.5*geo->squareSide - eps || c->pos[0] > geo->Nsquares*0.5*geo->squareSide + eps
	    || c->pos[1] < -geo->Nsquares*0.5*geo->squareSide - eps || c->pos[1] > geo->Nsquares*0.5*geo->squareSide + eps )
	    {
	  
          COUT(ERROR) << ENDL << "pos x or y out of range:" << ENDL;
	        COUT(ERROR) << ENDL << "x:" << c->pos[0] << " (" << -geo->Nsquares*0.5*geo->squareSide << ", " << geo->Nsquares*0.5*geo->squareSide << ")" << ENDL;
	        COUT(ERROR) << ENDL << "y:" << c->pos[1] << " (" << -geo->Nsquares*0.5*geo->squareSide << ", " << geo->Nsquares*0.5*geo->squareSide << ")" << ENDL;
          error = 1;
        }

        if(error==1) return 1;

      } //for i
    } //for iHit

    outFile->cd();
    events_tree->Fill();
  } //for iEv

  outFile->cd();
  events_tree->Write();
  
  COUT(INFO) <<" ";
  info::elapsed_time(start);

  return 0;
}
