#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TrCluster.hh"
#include "iostream"
#include <iostream>
#include <vector>
using namespace std;

int main(int argc, char **argv) {

  const int nLayer = 10;

  auto inFile = TFile::Open(argv[1]);

  TTree *events;
  inFile->GetObject("Tree", events);
  events->Print();
  TClonesArray *a = new TClonesArray("TrCluster", 200);
  events->SetBranchAddress("Events", &a);

  TFile *outFile = new TFile("histos.root", "recreate");

  TH1F *h = new TH1F("h", "h", 100, -0.5, 0.5);

  vector<vector<TrCluster>> v;

  for (int i = 0; i < events->GetEntries(); i++) {
    v.resize(nLayer);

    std::cout << " ------------ " << std::endl;
    std::cout << "Event " << i << std::endl;

    events->GetEntry(i);
    for (int j = 0; j < a->GetEntries(); j++) {
      TrCluster *cl = (TrCluster *)a->At(j);
      v[cl->layer].push_back(*cl);
    }

    // Analisi
    float tMean = 0;
    int _n = 0;
    for (auto il : v) {
      for (auto hit : il) {
        if (hit.parID != 0)
          continue;

        tMean += hit.time;
        _n++;
        // std::cout << hit.time << std::endl;
      }
    }
    tMean /= _n;

    for (auto il : v) {
      for (auto hit : il) {
        if (hit.parID != 0)
          continue;

        std::cout << hit.time - tMean << std::endl;
        h->Fill(hit.time - tMean);
      }
    }

    std::cout << " ------------ " << std::endl;

    v.clear();
  }

  outFile->WriteTObject(h);
  outFile->Close();
}
