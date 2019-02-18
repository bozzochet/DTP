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
	TH1F *h = new TH1F("disttemp", "disttemp", 100, -0.5, 0.5);
	TH1F *h1 = new TH1F("bt", "bt", 1000, 0, 2);
	TH1F *h2 = new TH1F("btls", "btls", 1000, 0, 2000);
	TH1F *h3 = new TH1F("btls_log", "btls_log", 1000, 0, 10);
	TH1F *h4 = new TH1F("bt_sim", "bt_sim", 1000, 0, 2);
	TH1F *h5 = new TH1F("bt_sim", "bt_sim", 1000, 0, 10);
	TRandom3 *tr = new TRandom3();
	vector<vector<TrCluster>> v;

	for (int i = 0; i < events->GetEntries(); i++) {
		v.resize(nLayer);

//		cout << " ------------ " << endl;
//		cout << "Event " << i << endl;

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
				if(tStart>hit.time) tStart = hit.time;
        
				if (hit.parID != 0) continue;

				tMean += hit.time;
				_n++;

			}
		}
		tMean /= _n;
//		cout<<"tempo minimo "<<tStart<<endl;

		for (int il = 0; il<v.size(); il++) {
			for (int hit = 0; hit<v[il].size(); hit++) {
				if(v[9].size()>5 && v[il][hit].eDep>0.00001){
					
//					cout<<	v[il][hit].parPdg <<endl;					
					 
					h1->Fill(v[il][hit].time - tStart); 
					h2->Fill(v[il][hit].time - tStart);
					h3->Fill(log(v[il][hit].time - tStart + 1));
					h4->Fill(tr->Gaus(v[il][hit].time - tStart,0.1));
					h5->Fill(tr->Gaus(v[il][hit].time - tStart,0.1));
				}
				if (v[il][hit].parID != 0) continue;

//				cout << v[il][hit].time - tMean << endl;
				h->Fill(v[il][hit].time - tMean); //ns
			}
		}

//		std::cout << " ------------ " << std::endl;

		v.clear();
	}

	outFile->WriteTObject(h);
	outFile->WriteTObject(h1);
	outFile->WriteTObject(h2);
	outFile->WriteTObject(h3);
	outFile->WriteTObject(h4);
	outFile->WriteTObject(h5);
	outFile->Close();
}
