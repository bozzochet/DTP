/*
 * SimpleAnalysis.C
 *
 *  Created on: 23 feb 2017
 *      Author: Nicola Mori
 */

/* Analysis ROOT script for the output file created by SimpleRun.mac. */

/* ********** INSTRUCTIONS ***********
 * Before loading this script in the Root shell you must load the dictionaries
 * for the GGS analysis classes, e.g.:
 *
 *   gROOT->ProcessLine(".include /path/to/GGSinstallation/lib);
 *   gSystem->Load("libGGSReader.so");
 *   gSystem->Load("libGGSDataObjects.so");
 *
 */

// GGS headers
#include "utils/GGSSmartLog.h"

void SimpleAnalysis(TString inputFileName, TString outputFileName) {
  static const std::string routineName("simpleanalysis");

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
  hReader->SetDetector("siLayer", kTRUE); // The name is the same of the sensitive logical volume name in the simulation

  // Retrieve the MC truth sub-reader
  GGSTMCTruthReader *mcReader = reader.GetReader<GGSTMCTruthReader>();

  // Create and retrieve the hadronic interaction sub-reader
  GGSTHadrIntReader *hadrReader = reader.GetReader<GGSTHadrIntReader>();

  // Prepare the histograms
  TH1F *eDepHisto = new TH1F("eDepHisto", "Total energy deposit;E (MeV); Counts", 1000, 0, 15); // Total energy release
  TH1F *zIntHisto = new TH1F("zIntHisto", "Z coordinate of interaction point [cm]", 1000, -0.05, 0.05); // Interaction point

  // Event loop
  int nInt = 0;
  COUT(INFO) << "Begin loop over " << reader.GetEntries() << " events" << ENDL;
  for (int iEv = 0; iEv < reader.GetEntries(); iEv++) {
  // for (int iEv = 0; iEv < 1; iEv++) {
    reader.GetEntry(iEv); // Reads all the data objects whose sub-readers have already been created


    // Retrieve inelastic interaction information
    GGSTHadrIntInfo *intInfo = hadrReader->GetInelastic(); // Get info about the inelastic interaction of the primary particle
    // Check if the inelastic interaction actually happened for this event
    if (intInfo) {
      zIntHisto->Fill(intInfo->GetInteractionPoint()[2]);
      nInt++;
    }

    // Compute total energy release
    GGSTIntHit* thisHit;
    GGSTPartHit* thisPHit;
    int nHits = hReader->GetNHits("siLayer"); //Number of hit siLayers for current event
    float totEDep = 0.;
    // Hits loop
    for (int iHit = 0; iHit < nHits; iHit++) {
      thisHit  = hReader->GetHit("siLayer", iHit);
      totEDep += hReader->GetHit("siLayer", iHit)->eDep;
    }
    if(!intInfo) eDepHisto->Fill(1e3*totEDep);


  }



  COUT(INFO) << "Event loop finished" << ENDL;

  // Save histograms
  outFile->cd();
  eDepHisto->Write();
  outFile->Close();
  delete outFile;

  COUT(INFO) << "Analysis finished" << ENDL;

}
