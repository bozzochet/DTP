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

  TFile *outFile = new TFile("histos-multi.root", "recreate");
  TCanvas *c1 = new TCanvas();
  TH1F *htotal = new TH1F("totale", "totale", 1000, 0, 4);
  TH1F *hprotons = new TH1F("protoni", "protoni", 1000, 0, 4);
  TH1F *hantip = new TH1F("antiprotoni", "antiprotoni", 1000, 0, 4);
  TH1F *hneutrons = new TH1F("neutroni", "neutroni", 1000, 0, 4);
  TH1F *hgamma = new TH1F("fotoni", "gamma", 1000, 0, 4);
  TH1F *hisotopes = new TH1F("isotopo", "isotopo", 1000, 0, 4);
  TH1F *helectron = new TH1F("elettroni", "elettroni", 1000, 0, 4);
  TH1F *hpositron = new TH1F("positroni", "positroni", 1000, 0, 4);
  TH1F *helectronmu = new TH1F("muoni", "muoni", 1000, 0, 4);
  TH1F *hpi = new TH1F("pi", "pi", 1000, 0, 4);
  TH1F *hk = new TH1F("k", "kaoni", 1000, 0, 4);
  TH1F *hprimari = new TH1F("primari", "primari", 1000, 0, 4);
  TH1F *hPrimEdep = new TH1F("PrimaryEdep", "edep", 500, 0, 0.001);
  TH1F *h6 = new TH1F("edep", "edep", 100, 0, 0.0001);
  TH1F *h8 = new TH1F("edep1", "edep", 100, 0, 0.001);
  TH1F *h9 = new TH1F("edep2", "edep", 100, 0, 0.001);
  TRandom3 *tr = new TRandom3();
  vector<vector<TrCluster>> v;

  for (int i = 0; i < events->GetEntries(); i++) {
    v.resize(nLayer);

    // cout << " ------------ " << endl;
    // cout << "Event " << i << endl;

    events->GetEntry(i);

    for (int j = 0; j < a->GetEntries(); j++) {
      TrCluster *cl = (TrCluster *)a->At(j);
      v[cl->layer].push_back(*cl);
    }

    // Analisi

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

    for (int il = 0; il < v.size(); il++) {
      for (int hit = 0; hit < v[il].size(); hit++) {
        // cout << v[il][hit].parID << " " << v[il][hit].parPdg << " " << v[il][hit].eDep << endl;

        // Get interacting events
        // if( v[9].size() < 6 ) continue;

        // 9 keV cut
        if (v[il][hit].eDep > 9e-6) {
          htotal->Fill(log10(v[il][hit].time - tStart + 1));
          if (v[il][hit].parID == 0) { // primario
            hprimari->Fill(log10(v[il][hit].time - tStart + 1));
            hPrimEdep->Fill(v[il][hit].eDep);
          } else { // secondari
            if (v[il][hit].parPdg == 2212) {
              hprotons->Fill(log10(v[il][hit].time - tStart + 1));
            } else if (v[il][hit].parPdg == -2212) {
              hantip->Fill(log10(v[il][hit].time - tStart + 1));
            } else if (v[il][hit].parPdg == 2112) {
              hneutrons->Fill(log10(v[il][hit].time - tStart + 1));
              h6->Fill(v[il][hit].eDep);
            } else if (v[il][hit].parPdg == 11) {
              helectron->Fill(log10(v[il][hit].time - tStart + 1));
              if (v[il][hit].time - tStart < 10000)
                h8->Fill(v[il][hit].eDep);
              if (v[il][hit].time - tStart >= 10000)
                h9->Fill(v[il][hit].eDep);
            } else if (v[il][hit].parPdg == -11) {
              hpositron->Fill(log10(v[il][hit].time - tStart + 1));
              if (v[il][hit].time - tStart < 10000)
                h8->Fill(v[il][hit].eDep);
              if (v[il][hit].time - tStart >= 10000)
                h9->Fill(v[il][hit].eDep);
            } else if (v[il][hit].parPdg == 211 || v[il][hit].parPdg == -211) {
              hpi->Fill(log10(v[il][hit].time - tStart + 1));
            } else if (v[il][hit].parPdg == 130 || v[il][hit].parPdg == 310 || v[il][hit].parPdg == 311 ||
                       v[il][hit].parPdg == 321 || v[il][hit].parPdg == -321) {
              hk->Fill(log10(v[il][hit].time - tStart + 1));
            } else if (v[il][hit].parPdg == 13 || v[il][hit].parPdg == -13) {
              helectronmu->Fill(log10(v[il][hit].time - tStart + 1));
            } else if (v[il][hit].parPdg == 22) {
              hgamma->Fill(log10(v[il][hit].time - tStart + 1));
            } else if (v[il][hit].parPdg > 1000000000) { // isotope
              hisotopes->Fill(log10(v[il][hit].time - tStart + 1));
              int code = v[il][hit].parPdg - 1000000000;
              int Z = (int)(code / 10000);
              int A = code - Z * 10000;
              A = (int)(A / 10);
              //						printf("%d -> Z=%d, A=%d\n", v[il][hit].parPdg, Z, A);
            } else {
              printf("%d\n", v[il][hit].parPdg);
            }
          }
        }
      }
    }

    v.clear();
  }

  htotal->SetLineColor(kBlack);
  htotal->Draw();
  hprotons->SetLineColor(kRed);
  hprotons->Draw("SAME");
  hantip->SetLineColor(kRed - 2);
  hantip->Draw("SAME");
  hneutrons->SetLineColor(kYellow);
  hneutrons->Draw("SAME");
  hgamma->SetLineColor(kGreen + 1);
  hgamma->Draw("SAME");
  hisotopes->SetLineColor(kRed + 2);
  hisotopes->Draw("SAME");
  helectron->SetLineColor(kBlue + 1);
  helectron->Draw("SAME");
  hpositron->SetLineColor(kCyan);
  hpositron->Draw("SAME");
  helectronmu->SetLineColor(kOrange + 1);
  helectronmu->Draw("SAME");
  hpi->SetLineColor(kGreen);
  hpi->Draw("SAME");
  hk->SetLineColor(kGray);
  hk->Draw("SAME");
  hprimari->SetLineColor(kViolet);
  hprimari->Draw("SAME");

  auto legend = new TLegend();
  legend->SetHeader("Legenda", "C");
  legend->AddEntry(htotal, "Total", "l");
  legend->AddEntry(hprotons, "Protons", "l");
  legend->AddEntry(hneutrons, "Neutrons", "l");
  legend->AddEntry(hgamma, "Gamma", "l");
  legend->AddEntry(hisotopes, "Isotopes", "l");
  legend->AddEntry(helectron, "Electrons", "l");
  legend->AddEntry(hpositron, "Positrons", "l");
  legend->AddEntry(helectronmu, "Muons", "l");
  legend->AddEntry(hpi, "Pions", "l");
  legend->AddEntry(hk, "Kaons", "l");
  legend->AddEntry(hprimari, "Primaries", "l");
  legend->Draw();

  outFile->WriteTObject(c1);

  outFile->WriteTObject(h6);

  TCanvas *c2 = new TCanvas();
  h8->SetLineColor(1);
  h8->Draw();
  h9->SetLineColor(2);
  h9->Draw("SAME");

  outFile->WriteTObject(c2);

  outFile->WriteTObject(htotal);
  outFile->WriteTObject(hprotons);
  outFile->WriteTObject(hneutrons);
  outFile->WriteTObject(hgamma);
  outFile->WriteTObject(hisotopes);
  outFile->WriteTObject(helectron);
  outFile->WriteTObject(hpositron);
  outFile->WriteTObject(helectronmu);
  outFile->WriteTObject(hpi);
  outFile->WriteTObject(hk);
  outFile->WriteTObject(hprimari);
  outFile->WriteTObject(hPrimEdep);

  outFile->Close();
}
