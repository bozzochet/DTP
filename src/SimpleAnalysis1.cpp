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

  int N = 5;

  TClonesArray a("TrCluster", 200);
  tree->Branch("Events", &a);

  for (int iEv = 0; iEv < reader.GetEntries(); iEv++) {
    reader.GetEntry(iEv); // Reads all the data objects whose sub-readers have already been created
    GGSTPartHit *phit;
    GGSTIntHit *inthit;
    int nHits = hReader->GetNHits("siSensor"); // Number of hit siLayers for current event
    a.Clear();

    // std::cout << " ---------- " << std::endl;
    // std::cout << "  Event  " << iEv << std::endl;
    // std::cout << " ---------- " << std::endl;

    //    std::cout << nHits << " hits" << std::endl;
    GGSTHadrIntInfo *intInfo = hadrReader->GetInelastic();
    if (intInfo)
      std::cout << "  Inelastic interaction happened at z = " << intInfo->GetInteractionPoint()[2] << std::endl;
    if (hadrReader->GetNQuasiElastic() > 0) {
      for (int iqe = 0; iqe < hadrReader->GetNQuasiElastic(); iqe++) {
        GGSTHadrIntInfo *qintInfo = hadrReader->GetQuasiElastic(iqe);
        if (qintInfo)
          std::cout << "  QuasiElastic interaction happened at z = " << qintInfo->GetInteractionPoint()[2] << std::endl;
      }
    }
    // Hits loop
    int ncluster = 0;
    for (int iHit = 0; iHit < nHits; iHit++) {
      inthit = hReader->GetHit("siSensor", iHit);
      int nPHit = inthit->GetNPartHits();
      int llayer = inthit->GetVolumeID() / (N * N);

      for (int i = 0; i < nPHit; i++) {
        TrCluster *c = (TrCluster *)a.ConstructedAt(ncluster++);
        phit = inthit->GetPartHit(i);

        c->segm = llayer % 2 == 0;
        c->pos[0] = 0.5 * (phit->entrancePoint[c->segm] + phit->exitPoint[c->segm]); // posizione lungo x o y
        c->pos[1] = 0.5 * (phit->entrancePoint[2] + phit->exitPoint[2]);             // posizione lungo z
        c->time = phit->time;
        c->eDep = phit->eDep;
        c->spRes = 0.00007;
        c->layer = llayer;
        c->parID = phit->parentID;
        c->parPdg = phit->particlePdg;
      }
    }

    tree->Fill();
  }

  COUT(INFO) << "Event loop finished" << ENDL;

  // Save histograms
  outFile->cd();
  tree->Write();
  outFile->Close();
  delete outFile;

  COUT(INFO) << "Analysis finished" << ENDL;
}
