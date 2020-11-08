/*
 * SimpleAnalysis.C
 *
 *  Created on: 23 feb 2017
 *      Author: Nicola Mori
 */

/* Analysis ROOT script for the output file created by SimpleRun.mac. */

#include "TrCluster.hh"
#include "physics.h"
#include "progress.h"
#include "Geometry.h"

#include "TFile.h"
#include "TString.h"
#include "TTree.h"

#include "montecarlo/readers/GGSTHadrIntReader.h"
#include "montecarlo/readers/GGSTHitsReader.h"
#include "montecarlo/readers/GGSTMCTruthReader.h"
#include "montecarlo/readers/GGSTRootReader.h"
#include "utils/GGSSmartLog.h"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

//global geometric parameters

Geometry GEO;

const int Nlayers = GEO.GetNlayers();
const int Nstrips = GEO.GetNstrips();
const int Nrows = GEO.GetNrows();
const int Nsquares = GEO.GetNsquares();
const int pitch = GEO.GetPitch();
const int squareSide = GEO.GetSquareSide();
const int Nladders = GEO.GetNladders();



int findStrip(double position) {
  int stripHit = (position+(Nsquares*squareSide*0.5))/pitch;
  return stripHit % (int(Nstrips));
}

int findLadder(int ID, int layer, bool condition) {
  if(!condition) {
    int row = (ID - (pow(Nsquares,2)*layer))/Nsquares;
    bool leftOrRight = ID % (int(Nsquares)) > (Nsquares*0.5);
    return row + (Nsquares*leftOrRight) + (Nsquares*2*layer);
    }
  else
    return (ID % ((int(Nsquares))*2)) + (layer*Nsquares*2);
}

int main(int argc, char **argv) {
  TString inputFileName = argv[1];
  TString outputFileName = argv[2];
  static const string routineName("simpleanalysis");

  GGSSmartLog::verboseLevel = GGSSmartLog::INFO; // Print only INFO messages or more important

  COUT(INFO) << "Begin analysis" << ENDL;

  // Create the reader container and open the data file
  GGSTRootReader reader;
  if (!(reader.Open(inputFileName))) {
    COUT(ERROR) << "Cannot open input file " << inputFileName << ENDL;
    return 1;
  }

  // Create the output file
  TFile *outFile = TFile::Open(outputFileName, "RECREATE");
  if (!outFile || outFile->IsZombie()) {
    COUT(ERROR) << "Cannot create output file " << outputFileName << ENDL;
    return 1;
  }

  // Create and retrieve the hits sub-reader
  GGSTHitsReader *hReader = reader.GetReader<GGSTHitsReader>();

  // Set which hit detectors are to be read
  // The name is the same of the sensitive logical volume name in the simulation
  hReader->SetDetector("siSensor", kTRUE);

  // Retrieve the MC truth sub-reader
  GGSTMCTruthReader *mcReader = reader.GetReader<GGSTMCTruthReader>();

  // Create and retrieve the hadronic interaction sub-reader
  GGSTHadrIntReader *hadrReader = reader.GetReader<GGSTHadrIntReader>();

  TTree *tree = new TTree("Tree", "siSensorHits");
  // TH1F *h = new TH1F("h", "test", 1000, 0, 5);

  COUT(INFO) << "Begin loop over " << reader.GetEntries() << " events" << ENDL;

	/*const GGSTGeoParams *geoParams = reader.GetGeoParams();

	int Nsquares = geoParams->GetRealGeoParam("Nsquares"); //squares per side
	int Nladders = Nsquares*5; //number of ladders
	int Nstrips = geoParams->GetRealGeoParam("Nstrips"); //strips per ladder
	double squareSide = geoParams->GetRealGeoParam("squareSide");*/

  TClonesArray a("TrCluster", 200);
  tree->Branch("Events", &a);

  for (int iEv = 0; iEv < reader.GetEntries(); iEv++) {

    //print and update progress bar;
    progress(iEv, reader.GetEntries());

    /* IMPORTANT: cout and printf MUST print strings beginning and
     * ending with a new line (i.e. "\n" or endl) to not overwrite
     * the bar */

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
		int llayer = inthit->GetVolumeID() / (Nsquares * Nsquares);
		//cout<<"layer : "<<llayer<<endl;
		//cout<<"square : "<<inthit->GetVolumeID()<<endl;



		for (int i = 0; i < nPHit; i++) {
			TrCluster *c = (TrCluster *)a.ConstructedAt(ncluster++);
			phit = inthit->GetPartHit(i);
      //cout<<"Entry #"<<iHit+i+1<<endl;

			c->segm = llayer % 2 == 0;
			c->pos[0] = 0.5 * (phit->entrancePoint[0] + phit->exitPoint[0]); // x coordinate
			c->pos[1] = 0.5 * (phit->entrancePoint[1] + phit->exitPoint[1]); // y coordinate
 			c->pos[2] = 0.5 * (phit->entrancePoint[2] + phit->exitPoint[2]); // z coordinate
			c->time = phit->time;
			c->eDep = phit->eDep;
			c->spRes = 0.00007;
			c->layer = llayer;
			c->parID = phit->parentID;
			c->parPdg = phit->particlePdg;
			c->ID = inthit->GetVolumeID();

			//Find the nearest strip to the left
			c->strip = findStrip(c->pos[c->segm]);


			// Find the ladder it belongs
      c->ladder = findLadder(c->ID, c->layer, c->segm);

      //Deposit energy and create cluster
      double thisPos = ((c->ladder%Nsquares)*squareSide) + (c->strip*pitch) - (Nsquares*squareSide*0.5);
      double fraction = (c->pos[c->segm]-thisPos)/pitch;

			c->clust[0] = c->eDep * (1-fraction);
			c->clust[1] = c->eDep * (fraction);

      }
    }

    tree->Fill();

  }

  cout <<endl <<endl;

  // Save histograms
  outFile->cd();
  tree->Write();
  outFile->Close();
  delete outFile;

}
