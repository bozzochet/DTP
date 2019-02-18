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
#include "TRandom3.h"
#include "TLegend.h"
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
//	TCanvas *c1 = new TCanvas();
	TH1F *h1 = new TH1F("", "", 1000, 0, 10);
/*	TH1F *h2 = new TH1F("protoni", "b", 1000, 0, 10);
	TH1F *h3 = new TH1F("neutroni", "c", 1000, 0, 10);
	TH1F *h4 = new TH1F("elettroni/positroni", "d", 1000, 0, 10);
	TH1F *h5 = new TH1F("pi", "e", 1000, 0, 10);
	TH1F *h7 = new TH1F("protoni primari", "f", 1000, 0, 10);
	TH1F *h6 = new TH1F("edep", "edep", 100, 0, 0.0001);
*/	TRandom3 *tr = new TRandom3();
	vector<vector<TrCluster>> v;

	for (int i = 0; i < events->GetEntries()-1; i++) {
		v.resize(nLayer);

		events->GetEntry(i);
		for (int j = 0; j < a->GetEntries(); j++) {
			TrCluster *cl = (TrCluster *)a->At(j);
			v[cl->layer].push_back(*cl);
		}

		events->GetEntry(i+1);
		for (int j = 0; j < a->GetEntries(); j++) {
			TrCluster *cl = (TrCluster *)a->At(j);
			v[cl->layer].push_back(*cl);
		}


		float tStart = v[0][0].time;
		float tMean = 0;
		int _n = 0;
		for (auto il : v) {
			for (auto hit : il) {
				if(tStart>hit.time) tStart = hit.time;
        
				if (hit.parID != 0) continue;

				tMean += hit.time;
				_n++;

			}
		}
		tMean /= _n;
//		cout<<"tempo minimo "<<tStart<<endl;

		for (int il = 0; il<v.size(); il++) {                                 //evento n
			for (int hit = 0; hit<v[il].size(); hit++) {
				if(v[9].size()>5 && v[il][hit].eDep>0.00001){
					
					//cout<<	v[il][hit].parPdg <<endl;					
					 
					h1->Fill(log(v[il][hit].time - tStart + 1)); 
/*					if (v[il][hit].parID == 0) h7->Fill(log(v[il][hit].time - tStart + 1));
					if(v[il][hit].parPdg == 2212 && v[il][hit].parID > 0) h2->Fill(log(v[il][hit].time - tStart + 1));
					if(v[il][hit].parPdg == 2112) { h3->Fill(log(v[il][hit].time - tStart + 1)); h6->Fill(v[il][hit].eDep);}
					if(v[il][hit].parPdg == 11 || v[il][hit].parPdg == -11 ) h4->Fill(log(v[il][hit].time - tStart + 1));
					if(v[il][hit].parPdg == -211) h5->Fill(log(v[il][hit].time - tStart + 1));
*/				}
			}
		}

		v.clear();
	}

/*
	h1->SetLineColor(1);
	h1->Draw();
	h2->SetLineColor(2);
	h2->Draw("SAME");
	h3->SetLineColor(5);
	h3->Draw("SAME");
	h4->SetLineColor(4);
	h4->Draw("SAME");
	h5->SetLineColor(3);
	h5->Draw("SAME");
	h7->SetLineColor(6);
	h7->Draw("SAME");

   auto legend = new TLegend();
   legend->SetHeader("Legenda","C"); 
   legend->AddEntry(h1,"totale","l");
   legend->AddEntry(h2,"protoni","l");
   legend->AddEntry(h3,"neutroni","l");
   legend->AddEntry(h4,"elettroni/positroni","l");
   legend->AddEntry(h5,"pioni","l");
   legend->AddEntry(h7,"protoni primari","l");
   legend->Draw();
*/
	outFile->WriteTObject(h1);
//	outFile->WriteTObject(h6);
	outFile->Close();
}
