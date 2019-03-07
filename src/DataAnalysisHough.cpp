#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TRandom3.h"
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

  TFile *outFile = new TFile("histos_el-ghosts.root", "recreate");
  TH1F *h1 = new TH1F("", "", 1000, 0, 10);
  TRandom3 *tr = new TRandom3();
  vector<vector<TrCluster>> v;

  for (int i = 0; i < events->GetEntries() - 1; i++) {
    v.resize(nLayer);

    events->GetEntry(i);
    for (int j = 0; j < a->GetEntries(); j++) {
      TrCluster *cl = (TrCluster *)a->At(j);
      v[cl->layer].push_back(*cl);
    }

    float tStart = v[0][0].time;
    float tMean = 0;
    int _n = 0;
    for (auto il : v) {
      for (auto hit : il) {
        if (tStart > hit.time)
          tStart = hit.time;

        if (hit.parID != 0)
          continue;

        tMean += hit.time;
        _n++;
      }
    }
    tMean /= _n;

    for (int il = 0; il < v.size(); il++) { // evento n
      for (int hit = 0; hit < v[il].size(); hit++) {
        if (v[9].size() > 5 && v[il][hit].eDep > 0.00001) {

          // cout<<	v[il][hit].parPdg <<endl;

          h1->Fill(log10(v[il][hit].time - tStart + 1));
        }
      }
    }

    v.clear();
  }

  outFile->WriteTObject(h1);
  outFile->Close();
}
