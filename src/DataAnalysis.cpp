#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include <TGraph.h>
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TrCluster.hh"
#include <iostream>
#include <vector>
#include "montecarlo/readers/GGSTHadrIntReader.h"
#include "montecarlo/readers/GGSTHitsReader.h"
#include "montecarlo/readers/GGSTMCTruthReader.h"
#include "montecarlo/readers/GGSTRootReader.h"
#include "utils/GGSSmartLog.h"
using namespace std;

int main(int argc, char **argv) {

	const int nLayer = 26;
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
	TH1F *segmp = new TH1F("segmpositions", "segmpositions", 1000, -0.05, 0.05);
	//TH2D *realp = new TH2D("positions","positions",640*16,0,640,640*16,-40,40);
	TRandom3 *tr = new TRandom3();
	vector<vector<TrCluster>> v;

	/*GGSTRootReader reader;
	const GGSTGeoParams *geoParams = reader.GetGeoParams();

	const int Nsquares = geoParams->GetRealGeoParam("Nsquares"); //squares per side
	const int Nlad = Nsquares*5; //number of ladders
	const int Nstrips = geoParams->GetRealGeoParam("Nstrips"); //strips per ladder
	const double squareSide = geoParams->GetRealGeoParam("squareSide");
	const double lenght = squareSide/(double(Nstrips));
	double eDepSegm[Nlad][Nstrips];
	double PrimeDepSegm[Nlad][Nstrips];
	double posf[Nlad][Nstrips];*/

	const int Nsquares = 8; //squares per side
	const int Nlad = Nsquares*2*10; //number of ladders
	const int Nstrips = 640; //strips per ladder
	const double squareSide = 10;
	const double lenght = squareSide/(double(Nstrips));
	double eDepSegm[Nlad][Nstrips];
	double PrimeDepSegm[Nlad][Nstrips];
	double posf[Nlad][Nstrips];


	// eDepSegm array initialisation

	for (int i = 0; i < Nlad; i++) {
		for (int j = 0; j < Nstrips; j++) {
			double fluct = tr->Gaus(0,9e-6);
			eDepSegm[i][j] = fluct;
			PrimeDepSegm[i][j] = fluct;
			}
		}

	// Fill the eDep histograms and identify clusters

	for (int i = 0; i < events->GetEntries(); i++) {

		v.resize(nLayer);
		events->GetEntry(i);
		if (a->GetEntries()>10) {
			printf("Event %d: %d hits\n", i, a->GetEntries());
	   		for (int j = 0; j < a->GetEntries(); j++) {
				TrCluster *cl = (TrCluster *)a->At(j);
				printf("%d) %d %f\n", j, cl->parID, cl->eDep);
				}
			}
		for (int j = 0; j < a->GetEntries(); j++) {

			cout<<endl<<"Entry #"<<i+j<<endl;
			TrCluster *cl = (TrCluster *)a->At(j);
			v[cl->layer].push_back(*cl);
			if (cl->parID == 0) hPrimEdep->Fill(cl->eDep); //primary
			if(cl->eDep > 9e-6) hEdep->Fill(cl->eDep); //total

			//Fill the strips with the current cluster

			eDepSegm[cl->ladder][cl->strip] = cl->clust[0];
			if(cl->strip==639) {
				if(cl->ladder%Nsquares !=7)
					eDepSegm[cl->ladder+1][0] = cl->clust[1];
				}
			else
				eDepSegm[cl->ladder][cl->strip+1] = cl->clust[1];

			for (int ix = cl->layer * Nsquares * 2; ix < (Nsquares*2*(cl->layer+1) - 1); ix++) {
				for (int jx = 0; jx < Nstrips; jx++) {

				//Find the boundaries of the clusters

				if(eDepSegm[ix][jx] >= 27e-6) {
					cout<<"Analysing cluster\n";

					int j1 = jx-1;
					while(j1>=0 && eDepSegm[ix][j1] > 9e-6)
						j1--;

					int j2 = jx+1;
					while(j1<Nstrips && eDepSegm[ix][j2] > 9e-6)
						j2++;


					//Find the boundaries of the peaks

					for(int k = j1; k <= j2; k++) {
						if(eDepSegm[ix][k] >= 27e-6) {

						cout<<"Analysing peak with "<<eDepSegm[ix][k]<<" energy\n";
						int k1 = k-1;
						int k2 = k+1;
						while(k1!=j1 && eDepSegm[ix][k1] >= 27e-6)
							k1--;
						while(k2!=j2 && eDepSegm[ix][k2] >= 27e-6)
							k2++;

						//Find the peaks

						int kmax = k1;
						for(int khold = k1; khold < k2; khold++)
							if(eDepSegm[ix][khold]>=eDepSegm[ix][kmax])
								kmax = khold;

						double kPos = (kmax)*lenght + 10*(cl->ladder%Nsquares)-((Nsquares*squareSide)/2);
						cout<<"Simulated position: "<<kPos<<endl;
						cout<<"Real position: "<<cl->pos[cl->segm]<<endl;
						double kLeftPos = kPos - lenght;
						double kRightPos = kPos + lenght;
						if(eDepSegm[ix][kmax-1]>eDepSegm[ix][kmax+1]) {
							double realpos = ((kPos*eDepSegm[ix][kmax]) + (kLeftPos*eDepSegm[ix][kmax-1]))/(eDepSegm[ix][kmax] + eDepSegm[ix][kmax-1]);
							segmp->Fill(realpos-cl->pos[cl->segm]); }
						if(eDepSegm[ix][kmax+1]>eDepSegm[ix][kmax-1]) {
							double realpos = ((kPos*eDepSegm[ix][kmax]) + (kRightPos*eDepSegm[ix][kmax+1]))/(eDepSegm[ix][kmax] + eDepSegm[ix][kmax+1]);
							segmp->Fill(realpos-cl->pos[cl->segm]); }
						if(eDepSegm[ix][kmax-1]==eDepSegm[ix][kmax+1]) {
							segmp->Fill(kPos-cl->pos[cl->segm]); }
						k = k2;
								}
				jx = j2;
							}
						}

				}
		}
	}

			//Empty the strips

			for (int i = 0; i < Nlad; i++) {
				for (int j = 0; j < Nstrips; j++) {
					double fluct = tr->Gaus(0,9e-6);
					eDepSegm[i][j] = fluct;
					PrimeDepSegm[i][j] = fluct;
					}
				}
			}

	// Fill the eDep histograms (w/ segmentation but not time resolution)

	/*double insert[Nstrips*Nlad];
	double number[Nlad*Nstrips];
	long long int index=0;
	for (int i = 0; i < Nlad; i++) {
		int index = 0;
		for (int j = 0; j < Nstrips; j++) {
			if(eDepSegm[i][j] > 9e-6)
				hEdep0->Fill(eDepSegm[i][j]); //total
			if(PrimeDepSegm[i][j] > 9e-6)
				hPrimEdep0->Fill(PrimeDepSegm[i][j]); //primary

			// Fill the realp histogram

			//realp->Fill(index,posf[i][j]);
		}
	}*/

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

	//Particle identification

	for (int il = 0; il<v.size(); il++) {
		for (int hit = 0; hit<v[il].size(); hit++) {

			if (v[il][hit].parPdg == 2212)
				hprotons->Fill(log10(v[il][hit].time - tStart + 1));
			if (v[il][hit].parPdg == -2212)
				hantip->Fill(log10(v[il][hit].time - tStart + 1));
			if (v[il][hit].parPdg == 2112)
				hneutrons->Fill(log10(v[il][hit].time - tStart + 1));
			if (v[il][hit].parPdg == 11)
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
	//outFile->WriteTObject(realp);
	outFile->WriteTObject(segmp);
	outFile->Close();
	}
