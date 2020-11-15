

#include "TrCluster.hh"
#include "physics.h"
#include "progress.h"
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


int fillTree(GGSTRootReader&, TTree*, Geometry*);

int digitization(TTree*, Geometry*);



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

  COUT(INFO) <<"Geometric parameters:" <<ENDL;
  COUT(INFO) <<"layers:                 " <<geo->Nlayers <<ENDL;
  COUT(INFO) <<"strips per ladder:      " <<geo->Nstrips <<ENDL;
  COUT(INFO) <<"ladders rows per layer: " <<geo->Nrows  <<ENDL;
  COUT(INFO) <<"squares per side:       " <<geo->Nsquares  <<ENDL;
  COUT(INFO) <<"implant pitch:          " <<geo->pitch  <<ENDL;
  COUT(INFO) <<"layers thickness:       " <<geo->thickness <<ENDL;
  COUT(INFO) <<"squares side:           " <<geo->squareSide <<ENDL;
  COUT(INFO) <<"ladders:                " <<geo->Nladders <<ENDL;

  COUT(INFO) <<ENDL;


  //save data

  TTree *data = new TTree("Data", "siSensorHits");

  if(fillTree(reader, data, geo) != 0 || digitization(data, geo) != 0)
    return 1;

  // Save histograms
  outFile->cd();
  data->Write();

  outFile->Close();
  delete outFile;


  COUT(INFO) <<ENDL;
  COUT(INFO) <<"Data saved in " <<outFileName <<ENDL;

  return 0;
}



int fillTree(GGSTRootReader &reader, TTree *tree, Geometry *geo)
{
  static const string routineName("Digitization::fillTree");

  COUT(INFO) <<ENDL;
  COUT(INFO) <<"Saving MC truth..." <<ENDL;

  TClonesArray a("TrCluster", 200);
  tree->Branch("Events", &a);

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

  for (int iEv = 0; iEv < reader.GetEntries(); iEv++) {

    reader.GetEntry(iEv); // Reads all the data objects whose sub-readers have already been created
    GGSTPartHit *phit;
    GGSTIntHit *inthit;
    int nHits = hReader->GetNHits("siSensor"); // Number of hit siLayers for current event
    a.Clear();

    GGSTHadrIntInfo *intInfo = hadrReader->GetInelastic();
    if (intInfo)
      cout << "\n  Inelastic interaction happened at z = " << intInfo->GetInteractionPoint()[2] << endl;
    if (hadrReader->GetNQuasiElastic() > 0) {
      for (int iqe = 0; iqe < hadrReader->GetNQuasiElastic(); iqe++) {
        GGSTHadrIntInfo *qintInfo = hadrReader->GetQuasiElastic(iqe);
        if (qintInfo)
          cout << "\n  QuasiElastic interaction happened at z = " << qintInfo->GetInteractionPoint()[2] << endl;
      }
    }
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

    tree->Fill();
  } //for iEv

  return 0;
}



int digitization(TTree *tree, Geometry *geo)
{
  static const string routineName("Digitization::digitization");

  COUT(INFO) <<ENDL;
  COUT(INFO) <<"Digitization..." <<ENDL;

  //set branch in which write measures
  measure meas;
  TBranch *branch = tree->Branch
  (
    "Measures", &(meas.time),
    "time[2]/D:xy/I:position/D:energy[2]/D"
  );

  //get MC truth
  TClonesArray *a = new TClonesArray("TrCluster", 200);
	tree->SetBranchAddress("Events", &a);


  TRandom3 *tr = new TRandom3(9298);

  TimeSegm *time_segm = new TimeSegm(geo, A, 1);
  TimeSim *time_sim = new TimeSim(time_segm, tr);

  PosSim *pos_sim = new PosSim(geo,1);


  COUT(INFO) << "Begin loop over " << tree->GetEntries() << " events" << ENDL;

  std::clock_t start = std::clock();

  for (int i = 0; i < tree->GetEntries(); i++)
  {

    //print and update progress bar
    progress(std::clock() - start, i, tree->GetEntries());

		tree->GetEntry(i);

		for (int j = 0; j < a->GetEntries(); j++)
    {
      TrCluster *cl = (TrCluster*) a->At(j);

      pos_sim->Reset(); //clear previuos hit

      // scan clusts
      for(int k = 0; k<2; ++k)
      {
        int strip = cl->strip;
        int ladder = cl->ladder;

        if(k==1) //go to next strip or stay?
        {
          //hit on the last strip of the last ladder of the layer row
          if(cl->strip == geo->Nstrips-1 && (cl->ladder+1) % geo->Nsquares == 0)
          {
            // do nothing
            // => deposit energy on the same strip of clust[0]
          }

    			else if(cl->strip==geo->Nstrips-1)
          {
            strip = 0;
            ++ladder;
          }

    			else
            ++strip;

        } //if k==1


        //time detection simulation

        TGraph *charge = new TGraph();

        time_sim->GetChargeSignal
          (cl->time, cl->clust[k], charge, false);

        //eDep + noise
        meas.energy[k] = cl->clust[k]
          + time_sim->AddChargeNoise(charge) / FOND_CHARGE
            * ENERGY_COUPLE ;

        if(meas.energy[k] < 0) //signal lost due to noise
          meas.energy[k] = 0;

        meas.time[k] = time_sim->GetMeas(charge, 0.15);

        /* if GetMeas returns -9999 => unable to measure time because
         * of noise */

        /* if t_meas = -9999 it is saved anyway: in an analysis of time
         * measures they will be excluded */

        delete charge;


        //position simulation

        //set hit pos one time for clust (clust refers to same pos)
        if(k==0) pos_sim->SetHitPos(cl->layer, cl->pos[cl->xy]);

        pos_sim->DepositEnergy(ladder, strip, meas.energy[k]);

		  } //for k

      pos_sim->ShareEnergy();
      /* must share energy between active strips stored in pos_sim
       * before calling PosSim::GetMeas */

      length_t simPos = pos_sim->GetMeas();

      meas.xy = cl->xy;
      meas.position = simPos;

      branch->Fill();

    } //for j
  } //for i


  //digitization execution time

  int sec = (start - std::clock()) / CLOCKS_PER_SEC;

  if(sec != 0)
    COUT(INFO) <<"Digitization took ";

  if(sec >= 3600)
  {
    COUT(INFO) <<sec/3600 <<"h ";
    sec %= 3600;
  }

  if(sec >= 60)
  {
    COUT(INFO) <<sec/60 <<"min ";
    sec %= 60;
  }

  if(sec > 0)
    COUT(INFO) <<sec <<"s " <<ENDL;


  delete pos_sim;
  delete time_sim;
  delete time_segm;
  delete tr;

  return 0;
}
