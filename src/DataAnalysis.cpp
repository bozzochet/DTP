#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include <TGraph.h>
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TrCluster.hh"
#include <TApplication.h>
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

int main(int argc, char **argv) {

	const int nLayer = 25;
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
	TH1F *h3 = new TH1F("btls_log", "btls_log", 1000, 0, 4);
	TH1F *h4 = new TH1F("bt_sim", "bt_sim", 1000, 0, 2);
	TH1F *h5 = new TH1F("btls_sim", "btls_sim", 1500, 0, 15);
	TH1F *h5cut = new TH1F("btls_sim_cut", "btls_sim_cut", 1500, 0, 15);
	TH1F *h5nopri = new TH1F("btls_sim_nopri", "btls_sim_nopri", 1500, 0, 15); 
	TH1F *h5nomip = new TH1F("btls_sim_nomip", "btls_sim_nomip", 1500, 0, 15);
	TH1F *hPrimEdep = new TH1F("PrimaryEdep", "edep", 500, 0, 0.001);
	TH1F *hPrimEdep0 = new TH1F("PrimaryEdep0", "edep0", 500, 0, 0.001);
	TH1F *hEdep = new TH1F("Edep", "edep", 500, 0, 0.001);
	TH1F *hEdep0 = new TH1F("Edep0", "edep0", 500, 0, 0.001);
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
	TH2D *realp= new TH2D("positions","positions",640*16,0,640,640*16,-40,40);
	TRandom3 *tr = new TRandom3();
	vector<vector<TrCluster>> v;
	

	int Nsquares = 8; //squares per side
	int Nlad = Nsquares*2*10; //number of ladders
	int Nstrips = 640; //strips per ladder
	double squareSide = 10;
	double lenght = squareSide/(double(Nstrips));
	double eDepSegm[Nlad][Nstrips];
	double PrimeDepSegm[Nlad][Nstrips];
	double posf[Nlad][Nstrips];
	

	//Array initialisation

	for (int i = 0; i < Nlad; i++) {
		for (int j = 0; j < Nstrips; j++) {
			eDepSegm[i][j] = 0;
			PrimeDepSegm[i][j] = 0;
			}
		}

	// Fill the eDep histograms (no segmentation)
	
	for (int i = 0; i < events->GetEntries(); i++) {
		
		v.resize(nLayer);
		
	  	//cout << " ------------ " << endl;
	  	//cout << "Event " << i << endl;
		
		events->GetEntry(i);
		if (a->GetEntries()>10) {
			printf("Event %d: %d hits\n", i, a->GetEntries());
	   		for (int j = 0; j < a->GetEntries(); j++) {
				TrCluster *cl = (TrCluster *)a->At(j);
				printf("%d) %d %f\n", j, cl->parID, cl->eDep);
				}
			}
		for (int j = 0; j < a->GetEntries(); j++) {

			cout<<endl<<"Entry #"<<j<<endl;
			TrCluster *cl = (TrCluster *)a->At(j);
			v[cl->layer].push_back(*cl);
			if (cl->parID == 0) hPrimEdep->Fill(cl->eDep); //primary
			if(cl->eDep > 9e-6) hEdep->Fill(cl->eDep); //total

			// Implement segmentation
			
			cout<<"ID = "<<cl->ID<<endl; //find the nearest strip to the left
			int stripHit = (cl->pos[cl->segm]+((Nsquares*squareSide)/2))/lenght;
			int strip = stripHit % Nstrips;
			cout<<"strip = "<<strip<<endl;

			int ladder, uladder;	//find the ladder
			if(!cl->segm) {
				uladder = (cl->ID - Nsquares*Nsquares*cl->layer)/Nsquares;
				if(cl->ID % Nsquares > Nsquares/2)
					ladder= uladder + Nsquares + Nsquares*2*cl->layer;
				else
					ladder= uladder + Nsquares*2*cl->layer;
				}
			else
				ladder = (cl->ID % (Nsquares*2)) + (cl->layer*(Nsquares*2));

			cout << "ladder = " << ladder << endl;

			// Fill the eDepSegm array

			double fraction = abs (remainder(cl->pos[cl->segm],lenght));
			fraction = fraction / lenght;			

			eDepSegm[ladder][strip] = eDepSegm[ladder][strip] + (cl->eDep * fraction);

			if (strip == 639)
				if (ladder % Nsquares != Nsquares - 1) 
					eDepSegm[ladder+1][0] = eDepSegm[ladder+1][0]+(cl->eDep*(1-fraction));
			else 
				eDepSegm[ladder][strip+1] = eDepSegm[ladder][strip+1]+(cl->eDep*(1-fraction));
			
			PrimeDepSegm[ladder][strip] = PrimeDepSegm[ladder][strip]+(cl->eDep*fraction);	

			if (cl->parID == 0) {
				if (ladder % Nsquares != Nsquares - 1) 
					PrimeDepSegm[ladder+1][0] = PrimeDepSegm[ladder+1][0]+(cl->eDep*(1-fraction));
				PrimeDepSegm[ladder][strip+1] = PrimeDepSegm[ladder][strip+1]+(cl->eDep*(1-fraction));
				}


			//Simulate position

			double fakep[2];
			fakep[cl->segm] = strip*lenght + 10*(ladder%Nsquares)-(Nsquares*5);
			int temp = ladder/Nsquares;
			if(temp % 2 == 0)
				fakep[!cl->segm] = (Nsquares*squareSide)/4;
			else
				fakep[!cl->segm] = -(Nsquares*squareSide)/4;
				
			cout<<"Simulated position = "<<fakep[cl->segm]<<endl;
			cout<<"Real position = "<<cl->pos[cl->segm]<<endl;
			posf[ladder][strip] = fakep[cl->segm];
			
			}
		}
	

	// Fill the eDep histograms (w/ segmentation)

	double insert[Nstrips*Nlad];
	double number[Nlad*Nstrips];
	long long int index=0;
	for (int i = 0; i < Nlad; i++) {
		int index = 0;
		for (int j = 0; j < Nstrips; j++) {
			if(eDepSegm[i][j] > 9e-6) {
				hEdep0->Fill(eDepSegm[i][j]); //total
				hPrimEdep0->Fill(eDepSegm[i][j]); //primary
			}
			
			realp->Fill(index,posf[i][j]);
			insert[j+(i*Nlad)] = posf[i][j];
			number[j+(i*Nlad)] = index;
			index++;			
		}
	}

	// Find tStart and tMean
	  
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
	
	//Analysis
	
	for (int il = 0; il<v.size(); il++) {
		for (int hit = 0; hit<v[il].size(); hit++) {
			//if(1) (
			//v[9].size()>5&&v[il][hit].eDep>0.00001 
			//10 keV 
			//cout<<v[il][hit].parPdg <<endl;			
			//eDep
			if (v[il][hit].parPdg == 2212)
				hprotons->Fill(log10(v[il][hit].time - tStart + 1));
			if (v[il][hit].parPdg == -2212)
				hantip->Fill(log10(v[il][hit].time - tStart + 1));
			if (v[il][hit].parPdg == 2112)
				hneutrons->Fill(log10(v[il][hit].time - tStart + 1));
			if (v[il][hit].parPdg == 11) {
				helectron->Fill(log10(v[il][hit].time - tStart + 1));
			if (v[il][hit].parPdg == 211 || v[il][hit].parPdg == -211)
				hpi->Fill(log10(v[il][hit].time - tStart + 1));
			if (v[il][hit].parPdg == 130 || v[il][hit].parPdg == 310 || v[il][hit].parPdg == 311 || v[il][hit].parPdg == 321 || v[il][hit].parPdg == -321)
				hk->Fill(log10(v[il][hit].time - tStart + 1));
			if (v[il][hit].parPdg == 13 || v[il][hit].parPdg == -13)
				helectronmu->Fill(log10(v[il][hit].time - tStart + 1));
			if (v[il][hit].parPdg == 22)
				hgamma->Fill(log10(v[il][hit].time - tStart + 1));
			if (v[il][hit].parPdg > 1000000000)
				hisotopes->Fill(log10(v[il][hit].time - tStart + 1));
			
			//time
			double smearedtime = tr->Gaus(v[il][hit].time - tStart,0.1);
			h1->Fill(v[il][hit].time - tStart);
			h2->Fill(v[il][hit].time - tStart);
			h3->Fill(log10(v[il][hit].time - tStart + 1));
			h4->Fill(smearedtime);
			h5->Fill(smearedtime);
			
			if (smearedtime<0.55) h5cut->Fill(smearedtime);
			if (v[il][hit].parID != 0) h5nopri->Fill(smearedtime);
			if (v[9].size()>5) h5nomip->Fill(smearedtime);
			if (v[il][hit].parID == 0) h->Fill(v[il][hit].time - tMean); //ns
			}
		}
	v.clear();

	outFile->WriteTObject(h);
	outFile->WriteTObject(h1);
	outFile->WriteTObject(h2);
	outFile->WriteTObject(h3);
	outFile->WriteTObject(h4);
	outFile->WriteTObject(h5);
	outFile->WriteTObject(h5cut);
	outFile->WriteTObject(h5nopri);
	outFile->WriteTObject(h5nomip);
	outFile->WriteTObject(hPrimEdep);
	outFile->WriteTObject(hEdep);
	outFile->WriteTObject(hPrimEdep0);
	outFile->WriteTObject(hEdep0);
	outFile->WriteTObject(hprotons);
	outFile->WriteTObject(hneutrons);
	outFile->WriteTObject(hgamma);
	outFile->WriteTObject(hisotopes);
	outFile->WriteTObject(helectron);
	outFile->WriteTObject(hpositron);
	outFile->WriteTObject(helectronmu);
	outFile->WriteTObject(hpi);
	outFile->WriteTObject(hk);
	outFile->WriteTObject(realp);
	outFile->Close();
	}
}


