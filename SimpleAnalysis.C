/*
 * SimpleAnalysis.C
 *
 *  Created on: 23 feb 2017
 *      Author: Nicola Mori
 */

/* Analysis ROOT script for the output file created by SimpleRun.mac. */

// GGS headers
#include "utils/GGSSmartLog.h"
#include <vector>
#include<iostream>
using namespace std;

struct TrCluster{
	double segm;
	double pos[2];
	double time;
	double eDep; 
	double spRes;
	double tRes;
};

void SimpleAnalysis(TString inputFileName, TString outputFileName) {
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
  hReader->SetDetector("siSensor", kTRUE); // The name is the same of the sensitive logical volume name in the simulation

  // Retrieve the MC truth sub-reader
  GGSTMCTruthReader *mcReader = reader.GetReader<GGSTMCTruthReader>();

  TTree *tree = new TTree("Tree","siSensorHits");
  //TH1F *h = new TH1F("h", "test", 1000, 0, 5);
  
COUT(INFO) << "Begin loop over " << reader.GetEntries() << " events" << ENDL;

int numberoflayers=10;
int N=5;
vector<vector<TrCluster> > v(numberoflayers, vector<TrCluster>());
TBranch *Events = tree->Branch("Events", "vector<vector<TrCluster> >", &v);

for (int iEv = 0; iEv < reader.GetEntries(); iEv++) {
	reader.GetEntry(iEv); // Reads all the data objects whose sub-readers have already been created
	GGSTPartHit* phit;
	GGSTIntHit* inthit;
	int nHits = hReader->GetNHits("siSensor"); //Number of hit siLayers for current event

	v.clear();

	// Hits loop
	for (int iHit = 0; iHit < nHits; iHit++) {
	TrCluster c;
	inthit  = hReader->GetHit("siSensor", iHit);
	phit  = inthit->GetPartHit(iHit);
	int layer = inthit->GetVolumeID()/(N*N);
	
	if(layer%2==0) c.segm= 0.5*(phit->entrancePoint[0]+phit->exitPoint[0]);
	else c.segm= 0.5*(phit->entrancePoint[1]+phit->exitPoint[1]);
	c.pos[0]= c.segm; //posizione lungo x o y
	c.pos[1]= 0.5*(phit->entrancePoint[2]+phit->exitPoint[2]); //posizione lungo z
	c.time= phit->time;
	c.eDep= phit->eDep;
	c.spRes= 0.00007;
	
	v[layer].push_back(c);
	}
	
	tree->Fill();
}



COUT(INFO) << "Event loop finished" << ENDL;

tree->Fill();

  // Save histograms
  outFile->cd();
  tree->Write();
  outFile->Close();
  delete outFile;

COUT(INFO) << "Analysis finished" << ENDL;

}
