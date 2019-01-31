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
	TH1F *h1 = new TH1F("h1", "h1", 1000, 0, 500);
	vector<vector<TrCluster>> v;

	for (int i = 0; i < events->GetEntries(); i++) {
		v.resize(nLayer);

		cout << " ------------ " << endl;
		cout << "Event " << i << endl;

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
		cout<<"tempo minimo "<<tStart<<endl;

		for (int il = 0; il<v.size(); il++) {
			for (int hit = 0; hit<v[il].size(); hit++) {
				if(v[9].size()>5) h1->Fill(v[il][hit].time - tStart);

				if (v[il][hit].parID != 0) continue;

				cout << v[il][hit].time - tMean << endl;
				h->Fill(v[il][hit].time - tMean); //ns
			}
		}

		std::cout << " ------------ " << std::endl;

		v.clear();
	}

	outFile->WriteTObject(h);
	outFile->WriteTObject(h1);
	outFile->Close();
}
