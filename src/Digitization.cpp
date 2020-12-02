

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



int findStrip
  (double position, int Nsquares, int Nstrips, double pitch)
{
  int stripHit = (position+(Nsquares*Nstrips*pitch*0.5))/pitch;
  return stripHit % (int(Nstrips));
}


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


int fillGeoTree(TTree*, Geometry*);

int fillEvTree(GGSTRootReader&, TTree*, Geometry*);

int fillCaloTree(GGSTRootReader&, TTree*);

int fillMeasTree(TTree*, TTree*, Geometry*);



int main(int argc, char **argv)
{

  static const string routineName("Digitization::main");
  GGSSmartLog::verboseLevel = GGSSmartLog::INFO; // Print only INFO messages or more important

  TString inputFileName = argv[1];
  TString outFileName = argv[2];

  // Create the reader container and open the data file
  GGSTRootReader reader;
  if (!(reader.Open(inputFileName))) {
    COUT(ERROR) << "Cannot open input file " << inputFileName << ENDL;
    return 1;
  }

  //create output file
  TFile *outFile = TFile::Open(outFileName, "RECREATE");
  if (!outFile || outFile->IsZombie()) {
    COUT(ERROR) << "Cannot create digitization output file " << outFileName << ENDL;
    return 1;
  }


  //geometry

  const GGSTGeoParams *GEO = reader.GetGeoParams();
  Geometry *geo = new Geometry;

  geo->Nlayers = GEO->GetIntGeoParam("Nlayers");
  geo->Nstrips = GEO->GetIntGeoParam("Nstrips");
  geo->Nrows = GEO->GetIntGeoParam("Nrows");
  geo->Nsquares = GEO->GetIntGeoParam("Nsquares");
  geo->pitch = 1e-2 * GEO->GetRealGeoParam("pitch");
  geo->thickness = 1e-3 * GEO->GetRealGeoParam("thickness");

  geo->squareSide = geo->pitch * ((double) geo->Nstrips);
  geo->Nladders = geo->Nsquares * geo->Nrows * geo->Nlayers;

  COUT(INFO) <<ENDL;
  COUT(INFO) <<"=================================" <<ENDL;
  COUT(INFO) <<"Geometric parameters:     " <<ENDL;
  COUT(INFO) <<"  layers:                 " <<geo->Nlayers <<ENDL;
  COUT(INFO) <<"  strips per ladder:      " <<geo->Nstrips <<ENDL;
  COUT(INFO) <<"  ladders rows per layer: " <<geo->Nrows  <<ENDL;
  COUT(INFO) <<"  squares per side:       " <<geo->Nsquares  <<ENDL;
  COUT(INFO) <<"  implant pitch:          " <<geo->pitch  <<ENDL;
  COUT(INFO) <<"  layers thickness:       " <<geo->thickness <<ENDL;
  COUT(INFO) <<"  squares side:           " <<geo->squareSide <<ENDL;
  COUT(INFO) <<"  ladders:                " <<geo->Nladders <<ENDL;
  COUT(INFO) <<"=================================" <<ENDL;


  //trees

  TTree *geo_tree = new TTree("geometry", "siSensorGeoParams");
  TTree *events_tree = new TTree("events", "siSensorHits");
  TTree *meas_tree = new TTree("measures", "siSensorMeasures");
  TTree *calo_tree = new TTree("calorimeter", "calorimeterEdep");

  if
  (
    fillGeoTree(geo_tree, geo) != 0
    || fillEvTree(reader, events_tree, geo) != 0
    || fillCaloTree(reader, calo_tree) != 0
    || fillMeasTree(events_tree, meas_tree, geo) != 0
  )
    return 1;

  // Save histograms
  outFile->cd();

  geo_tree->Write();
  events_tree->Write();
  calo_tree->Write();
  meas_tree->Write();

  outFile->Close();
  delete outFile;


  COUT(INFO) <<ENDL;
  COUT(INFO) <<"Data saved in " <<outFileName <<ENDL;

  return 0;
}


int fillGeoTree(TTree *geo_tree, Geometry *geo)
{
  static const string routineName("Digitization::fillGeoTree");


  COUT(INFO) <<ENDL;
  COUT(INFO) <<"Saving geometric parameters..." <<ENDL;

  geo_tree->Branch
  (
    "Geometry", &(geo->Nlayers),
    "Nlayers:Nstrips:Nrows:Nsquares:pitch/D:thickness/D:squareSide/D:Nladders"
  );

  geo_tree->Fill();

  return 0;
}


int fillEvTree
  (GGSTRootReader &reader, TTree *events_tree, Geometry *geo)
{
  static const string routineName("Digitization::fillEvTree");


  COUT(INFO) <<ENDL;
  COUT(INFO) <<"Saving MC truth..." <<ENDL;

  TClonesArray a("TrCluster", 200);
  events_tree->Branch("Events", &a);


  // Create and retrieve the hits sub-reader
  GGSTHitsReader *hReader = reader.GetReader<GGSTHitsReader>();

  // Set which hit detectors are to be read
  // The name is the same of the sensitive logical volume name in the simulation
  hReader->SetDetector("siSensor", kTRUE);

  // Retrieve the MC truth sub-reader
  GGSTMCTruthReader *mcReader = reader.GetReader<GGSTMCTruthReader>();

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

    GGSTHadrIntInfo *intInfo = hadrReader->GetInelastic();

    /*
    if (intInfo)
      cout << "\n  Inelastic interaction happened at z = " << intInfo->GetInteractionPoint()[2] << endl;
    if (hadrReader->GetNQuasiElastic() > 0) {
      for (int iqe = 0; iqe < hadrReader->GetNQuasiElastic(); iqe++) {
        GGSTHadrIntInfo *qintInfo = hadrReader->GetQuasiElastic(iqe);
        if (qintInfo)
          cout << "\n  QuasiElastic interaction happened at z = " << qintInfo->GetInteractionPoint()[2] << endl;
      }
    }
    */

    // Hits loop
    int ncluster = 0;

    for (int iHit = 0; iHit < nHits; iHit++) {
      inthit = hReader->GetHit("siSensor", iHit);
      int nPHit = inthit->GetNPartHits();
      int llayer = inthit->GetVolumeID() / (geo->Nsquares * geo->Nsquares);
      //cout<<"layer : "<<llayer<<endl;
      //cout<<"square : "<<inthit->GetVolumeID()<<endl;

      for (int i = 0; i < nPHit; i++) {
        TrCluster *c = (TrCluster *)a.ConstructedAt(ncluster++);
        phit = inthit->GetPartHit(i);
        //cout<<"Entry #"<<iHit+i+1<<endl;


        //fill with MC truth

        c->xy = llayer % 2 == 0;
        c->pos[0] = 1e-2 * 0.5 * (phit->entrancePoint[0] + phit->exitPoint[0]); // x coordinate
        c->pos[1] = 1e-2 * 0.5 * (phit->entrancePoint[1] + phit->exitPoint[1]); // y coordinate
        c->pos[2] = 1e-2 * 0.5 * (phit->entrancePoint[2] + phit->exitPoint[2]); // z coordinate
        c->time = 1e-9 * phit->time;
        c->eDep = 1e+9 * phit->eDep;
        //c->spRes = 0.00007;
        c->layer = llayer;
        c->parID = phit->parentID;
        c->parPdg = phit->particlePdg;
        c->ID = inthit->GetVolumeID();

        //Find the nearest strip to the left
        c->strip = findStrip(c->pos[c->xy], geo->Nsquares, geo->Nstrips, geo->pitch);

        // Find the ladder it belongs
        c->ladder = findLadder(c->ID, c->layer, c->xy, geo->Nsquares);

        //Deposit energy and create cluster
        double thisPos = ((c->ladder%geo->Nsquares)*geo->squareSide) + (c->strip*geo->pitch) - (geo->Nsquares*geo->squareSide*0.5);
        double fraction = (c->pos[c->xy]-thisPos)/geo->pitch;

        c->clust[0] = c->eDep * (1-fraction);
        c->clust[1] = c->eDep * (fraction);


        //errors

        int error = 0;

        if
        (
          c->strip < 0 || c->strip >= geo->Nstrips
          || c->ladder < 0 || c->ladder >= geo->Nladders * geo->Nlayers
        )
        {
          COUT(ERROR) <<ENDL <<"ladder or strip out of range: lad:"
            <<c->ladder <<" str:" <<c->strip <<ENDL;

          error = 1;
        }


        if
        (
          c->pos[0] < -geo->Nsquares*0.5*geo->squareSide
          || c->pos[0] > geo->Nsquares*0.5*geo->squareSide
          || c->pos[1] < -geo->Nsquares*0.5*geo->squareSide
          || c->pos[1] > geo->Nsquares*0.5*geo->squareSide
        )
        {
          COUT(ERROR) <<ENDL <<"pos x or y out of range: x:" <<c->pos[0]
            <<" y:" <<c->pos[1] <<ENDL;

          error = 1;
        }

        if(error==1) return 1;

      } //for i
    } //for iHit

    events_tree->Fill();
  } //for iEv

  COUT(INFO) <<" ";
  info::elapsed_time(start);

  return 0;
}


int fillCaloTree(GGSTRootReader &reader, TTree *calo_tree)
{
  static const string routineName("Digitization::fillCaloTree");


  COUT(INFO) <<ENDL;
  COUT(INFO) <<"Saving calorimeter energies..." <<ENDL;

  energy_t calo_energy = 0;
  calo_tree->Branch("Events", &calo_energy); //TTree can understand that calo_energy is a double?


  // Create and retrieve the hits sub-reader
  GGSTHitsReader *hReader = reader.GetReader<GGSTHitsReader>();

  // Set which hit detectors are to be read
  // The name is the same of the sensitive logical volume name in the simulation
  try
  {
    hReader->SetDetector("calorimeter", kTRUE);
  }
  catch(const std::runtime_error &e) //calorimeter not sensitive
  {
    COUT(INFO) <<"!!! calorimeter hits not available" <<ENDL;
    COUT(INFO) <<"If needed, setup sensitive calo in" <<ENDL;
    COUT(INFO) <<" macros/run.mac and run again the simulation" <<ENDL;
    return 0;
  }


  COUT(INFO) << "Begin loop over " << reader.GetEntries() << " events" << ENDL;

  std::clock_t start = std::clock();

  for (int iEv = 0; iEv < reader.GetEntries(); iEv++)
  {

    //print and update progress bar
    info::progress(start, iEv, reader.GetEntries());

    //reset energy for new event
    calo_energy = 0;

    reader.GetEntry(iEv); // Reads all the data objects whose sub-readers have already been created

    GGSTPartHit *phit;
    GGSTIntHit *inthit;

    int nHits = hReader->GetNHits("calorimeter"); // Number of hit siLayers for current event

    for (int iHit = 0; iHit < nHits; iHit++)
    {
      inthit = hReader->GetHit("calorimeter", iHit);
      int nPHit = inthit->GetNPartHits();

      for (int i = 0; i < nPHit; i++)
      {
        phit = inthit->GetPartHit(i);

        //add to total energy deposited in calo by event iEv
        calo_energy += 1e+9 * phit->eDep;

      } //for i
    } //for iHit

    calo_tree->Fill();

  } //for iEv

  COUT(INFO) <<" ";
  info::elapsed_time(start);

  return 0;
}


int fillMeasTree(TTree *events_tree, TTree *meas_tree, Geometry *geo)
{
  static const string routineName("Digitization::fillMeasTree");

  COUT(INFO) <<ENDL;
  COUT(INFO) <<"Saving measures..." <<ENDL;

  //set branch in which write measures
  measure meas;
  meas_tree->Branch
  (
    "Measures", &(meas.energy),
    "energy[2]/D:time[2]/D:position/D:xy/I"
  );

  //get MC truth
  TClonesArray *a = new TClonesArray("TrCluster", 200);
	events_tree->SetBranchAddress("Events", &a);


  TRandom3 *tr = new TRandom3(9298);

  TimeSegm *time_segm = new TimeSegm(geo, A, 1);
  TimeSim *time_sim = new TimeSim(time_segm, tr);

  PosSim *pos_sim = new PosSim(geo,1);


  COUT(INFO) << "Begin loop over " << events_tree->GetEntries() << " events" << ENDL;

  std::clock_t start = std::clock();

  for (int i = 0; i < events_tree->GetEntries(); i++)
  {

    //clear from previous events
    pos_sim->Reset();

    //print and update progress bar
    info::progress(start, i, events_tree->GetEntries());

		events_tree->GetEntry(i);

		for (int j = 0; j < a->GetEntries(); j++)
    {
      TrCluster *cl = (TrCluster*) a->At(j);

      // scan clusts
      for(int k = 0; k<2; ++k)
      {
        int strip = cl->strip;
        int ladder = cl->ladder;


        //clear from hits already processed
        time_segm->Reset();


        //go to next strip or stay?

        //hit on the last strip of the last ladder of the layer row
        if
        (
          k==1
          && cl->strip == geo->Nstrips-1
          && (cl->ladder+1) % geo->Nsquares == 0
        )
        {
          //other fraction of energy is lost: no strip on right
          meas.energy[k] = 0;
          meas.time[k] = -9999;
          continue;
        }

    		else if(k==1 && cl->strip==geo->Nstrips-1)
        {
          strip = 0;
          ++ladder;
        }

    		else if(k==1)
          ++strip;


        //set hit pos one time for clust (clust refers to same pos)
        if(k==0) pos_sim->SetHitPos(cl->layer, cl->pos[cl->xy]);


        //energy

        energy_t noise = time_sim->GetChargeNoise();

        meas.energy[k] = cl->clust[k]
          + noise / FOND_CHARGE * ENERGY_COUPLE ;


        if(meas.energy[k] > 0)
          pos_sim->DepositEnergy(ladder, strip, meas.energy[k]);

        else
        //deposit 0: a negative energy could affect pos calculation
          pos_sim->DepositEnergy(ladder, strip, 0);


        //time measures

        int gr = time_segm->SetHit(cl, noise);

        TGraph *charge = new TGraph();

        time_sim->GetChargeSignal(gr, charge);
        meas.time[k] = time_sim->GetMeas(charge, 0.15);

        delete charge;

		  } //for k


      //signal lost completely => unable to measure pos
      if(meas.energy[0] <= 0 && meas.energy[1] <= 0)
        meas.position = -9999;

      else
      {
        pos_sim->ShareEnergy(); //share energies between active strips
        meas.position = pos_sim->GetMeas();
      }

      meas.xy = cl->xy;

      /*********************
      * DEBUG:
      * data calculated here does not correspond with the ones
      * read in DataAnalysis. This happens not regularly it seems and
      * only for sme parameters; now it seems to affect position
      * measures pos_meas read in analysis is inf while in
      * Digitization was good (also well reconstructed);
      * probably there is some unexpected behaviuour in the use
      * of TTree object meas_tree
      **********************/

/*
      if(i==2 && j==4)
        COUT(INFO) <<ENDL <<"\n\tE: " <<cl->clust[0] <<" " <<cl->clust[1] <<"\n\tE + noise: " <<meas.energy[0] <<" " <<meas.energy[1] <<"\n\tt: " <<cl->time <<"\n\tt meas: " <<meas.time[0] <<" " <<meas.time[1] <<"\n\tpos: " <<cl->pos[cl->xy] <<"\n\tpos meas: " <<meas.position <<ENDL;
*/

      meas_tree->Fill();

    } //for j
  } //for i


  COUT(INFO) <<" ";
  info::elapsed_time(start);


  delete pos_sim;
  delete time_sim;
  delete time_segm;
  delete tr;

  return 0;
}
