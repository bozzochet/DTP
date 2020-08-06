/*
 * SimpleAnalysis.C
 *
 *  Created on: 23 feb 2017
 *      Author: Nicola Mori
 */

/* Analysis ROOT script for the output file created by SimpleRun.mac. */

// GGS headers
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TrCluster.hh"
#include "montecarlo/readers/GGSTHadrIntReader.h"
#include "montecarlo/readers/GGSTHitsReader.h"
#include "montecarlo/readers/GGSTMCTruthReader.h"
#include "montecarlo/readers/GGSTRootReader.h"
#include "utils/GGSSmartLog.h"
#include <iostream>
#include <vector>
using namespace std;

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
	int Nlad = Nsquares*5; //number of ladders
	int Nstrips = geoParams->GetRealGeoParam("Nstrips"); //strips per ladder
	double squareSide = geoParams->GetRealGeoParam("squareSide");*/
	const int Nsquares = 8; //squares per side
	const int Nlad = Nsquares*2*10; //number of ladders
	const int Nstrips = 640; //strips per ladder
	const double squareSide = 10;
	const double lenght = squareSide/(double(Nstrips));

  TClonesArray a("TrCluster", 200);
  tree->Branch("Events", &a);

  for (int iEv = 0; iEv < reader.GetEntries(); iEv++) {

    reader.GetEntry(iEv); // Reads all the data objects whose sub-readers have already been created
    GGSTPartHit *phit;
    GGSTIntHit *inthit;
    int nHits = hReader->GetNHits("siSensor"); // Number of hit siLayers for current event
    a.Clear();

    GGSTHadrIntInfo *intInfo = hadrReader->GetInelastic();
    if (intInfo)
      cout << "  Inelastic interaction happened at z = " << intInfo->GetInteractionPoint()[2] << endl;
    if (hadrReader->GetNQuasiElastic() > 0) {
      for (int iqe = 0; iqe < hadrReader->GetNQuasiElastic(); iqe++) {
        GGSTHadrIntInfo *qintInfo = hadrReader->GetQuasiElastic(iqe);
        if (qintInfo)
          cout << "  QuasiElastic interaction happened at z = " << qintInfo->GetInteractionPoint()[2] << endl;
      }
    }
    // Hits loop
    int ncluster = 0;

	 for (int iHit = 0; iHit < nHits; iHit++) {
		inthit = hReader->GetHit("siSensor", iHit);
		int nPHit = inthit->GetNPartHits();
		int llayer = inthit->GetVolumeID() / (Nsquares * Nsquares);
		cout<<"layer : "<<llayer<<endl;
		cout<<"square : "<<inthit->GetVolumeID()<<endl;



		for (int i = 0; i < nPHit; i++) {
			TrCluster *c = (TrCluster *)a.ConstructedAt(ncluster++);
			phit = inthit->GetPartHit(i);
      cout<<"Entry #"<<(iHit+1)*(i+1)<<endl;
			c->segm = llayer % 2 == 0;
			c->pos[0] = 0.5 * (phit->entrancePoint[0] + phit->exitPoint[0]); // x coordinate
			c->pos[1] = 0.5 * (phit->entrancePoint[1] + phit->exitPoint[1]); // y coordinate
 			c->pos[2] = 0.5 * (phit->entrancePoint[2] + phit->exitPoint[2]); // z coordinate
			cout << c->pos[c->segm] <<endl;
			c->time = phit->time;
			c->eDep = phit->eDep;
			c->spRes = 0.00007;
			c->layer = llayer;
			c->parID = phit->parentID;
			c->parPdg = phit->particlePdg;
			c->ID = inthit->GetVolumeID();

			//Find the nearest strip to the left

		  int stripHit = (c->pos[c->segm]+((Nsquares*squareSide)/2))/lenght;
			int strip = stripHit % Nstrips;

			// Find the ladder it belongs

			int ladder, uladder;
			if(!c->segm) {
				uladder = (c->ID - Nsquares*Nsquares*c->layer)/Nsquares;
				if(c->ID % Nsquares > Nsquares/2)
					ladder= uladder + Nsquares + Nsquares*2*c->layer;
				else
					ladder= uladder + Nsquares*2*c->layer;
				}
			else
				ladder = (c->ID % (Nsquares*2)) + (c->layer*(Nsquares*2));
			if(strip == 639)
				ladder++;
			c->ladder = ladder;

			double fraction = remainder(c->pos[c->segm],lenght)/lenght;

			c->strip = strip;

			//Deposit energy and create cluster

			if(fraction < 0) {
				c->clust[0] = c->eDep * (1-abs(fraction));
				c->clust[1] = c->eDep * abs(fraction);
				}
			else {
				c->clust[0] = c->eDep * fraction;
				c->clust[1] = c->eDep * (1 - fraction);
      }
		}
  }

    tree->Fill();
  }

  //COUT(INFO) << "Event loop finished" << ENDL;

  // Save histograms
  outFile->cd();
  tree->Write();
  outFile->Close();
  delete outFile;

  //COUT(INFO) << "Analysis finished" << ENDL;
}
