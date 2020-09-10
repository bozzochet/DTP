#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include <TGraph.h>
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TrCluster.hh"
#include <iostream>
#include <vector>
#include "utils/GGSSmartLog.h"
using namespace std;

void stripReset(double array[][640], int dim1, int dim2) {
	TRandom3 *tr = new TRandom3();
	for (int ii = 0; ii < dim1; ii++) {
		for (int jj = 0; jj < dim2; jj++) {
			double fluct = tr->Gaus(0,9e-6);
			array[ii][jj] = fluct;
			}
		}
}

void hitReset(vector<double> array[]) {
	const int Nlayers = 10;
	for(int i = 0; i < Nlayers; i++) {
		array[i].clear();
	}
}

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

	TH1F *segmp = new TH1F("segmpositions", "segmpositions", 10000, -0.05, 0.05);

	TRandom3 *tr = new TRandom3();
	vector<vector<TrCluster>> v;

	/*

	Passing the parameters from DetectorConstruction.cc, not working yet

	GGSTRootReader reader;
	const GGSTGeoParams *geoParams = reader.GetGeoParams();

	const int Nsquares = geoParams->GetRealGeoParam("Nsquares"); //squares per side
	const int Nlad = Nsquares*5; //number of ladders
	const int Nstrips = geoParams->GetRealGeoParam("Nstrips"); //strips per ladder
	const double squareSide = geoParams->GetRealGeoParam("squareSide");
	const double pitch = squareSide/(double(Nstrips));
	double eDepSegm[Nlad][Nstrips];
	double PrimeDepSegm[Nlad][Nstrips];

	*/
	const int Nlayers = 10;
	const int Nsquares = 8; //squares per side
	const int Nlad = Nsquares*2*Nlayers; //number of ladders
	const int Nstrips = 640; //strips per ladder
	const double squareSide = 10;
	const double pitch = squareSide/(double(Nstrips));
	const int jump = 1;
	double eDepSegm[Nlad][Nstrips];
	vector<double> hitPos[Nlayers];


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

			// eDepSegm array initialisation
			stripReset(eDepSegm,Nlad,Nstrips);
			hitReset(hitPos);

		for (int j = 0; j < a->GetEntries(); j++) {

			cout<<endl<<"Entry #"<<j<<endl;
			TrCluster *cl = (TrCluster *)a->At(j);
			v[cl->layer].push_back(*cl);
			if (cl->parID == 0) hPrimEdep->Fill(cl->eDep); //primary
			if(cl->eDep > 9e-6) hEdep->Fill(cl->eDep); //total

			//Fill the strips with the current cluster

			eDepSegm[cl->ladder][cl->strip] += cl->clust[0];
			if(cl->strip==Nstrips-1)
				eDepSegm[cl->ladder+1][0] += cl->clust[1];
			else
				eDepSegm[cl->ladder][cl->strip+1] += cl->clust[1];

			//Fill the hit array

			hitPos[cl->layer].push_back(cl->pos[cl->segm]);

			// Sharing of the energy from non-active strips

			/*for (int ix = cl->layer * Nsquares * 2; ix < (Nsquares*2*(cl->layer+1) - 1); ix++) {
				for (int jx = 0; jx < Nstrips; jx+=jump) {

					if(jx > 0 && jx < Nstrips-1) {
						eDepSegm[ix][jx-1] += eDepSegm[ix][jx]*0.5;
						eDepSegm[ix][jx+1] += eDepSegm[ix][jx]*0.5;
					}
					if(jx==0)
						eDepSegm[ix][jx+1] += eDepSegm[ix][jx];
					if(jx==Nstrips-1)
						eDepSegm[ix][jx-1] += eDepSegm[ix][jx];

					eDepSegm[ix][jx] = 0;

					}
				}*/

			}

		for (int ix = 0; ix < Nlad; ix++) {
			for (int jx = 0; jx < Nstrips; jx++) {


			//Find the boundaries of the clusters

			if(eDepSegm[ix][jx] >= 27e-6) {
				cout<<"Analysing cluster\n";

				vector< pair<double,double>> strip;

				int i1 = ix;
				int j1 = jx-1;

				while(eDepSegm[i1][j1] > 9e-6){
					j1-=jump;
					if (j1<0) {
						i1--;
						j1 = Nstrips-1-j1;
						}
					}

				int i2 = ix;
				int j2 = jx+1;

				while(eDepSegm[i2][j2] > 9e-6) {
					j2+=jump;
					if (j1>=Nstrips) {
						i2++;
						j2 = j2-Nstrips+1;
						}
					}

				int clustSize = ((i2-i1)*Nstrips)+(j2-j1)+1;

				for(int ij = 0; ij<clustSize; ij+=jump) {

						double thisPos = ((i1%Nsquares)*squareSide) + (j1*pitch) - (Nsquares*squareSide*0.5);
						strip.push_back(make_pair(thisPos,eDepSegm[i1][j1]));

						if(j1>=Nstrips-1) {
							i1++;
							j1 = 0;
						}

						j1+=jump;
					}

				//Find the boundaries of the peaks

				for(int k = 0; k < strip.size(); k++) {
					if(strip[k].second >= 27e-6) {

					int k1 = k-1;
					int k2 = k+1;

					while(k1>0 && strip[k1].second >= 27e-6)
						k1--;
					while(k2<(strip.size()-1) && strip[k2].second >= 27e-6)
						k2++;

					//Find the peak

					int kMax = k1;
					for(int kHold = k1; kHold <= k2; kHold++)
						if(strip[kHold].second>strip[kMax].second)
							kMax = kHold;


					cout<<"Analysing peak of "<< strip[kMax].second <<"keV\n";

					//Find neighbour strip

					double kNext;

					if(strip[kMax+1].second > strip[kMax-1].second)
						kNext = kMax+1;
					else
						kNext = kMax-1;
					if(kMax == 0)
						kNext = 1;
					if (kMax == strip.size())
							kNext = strip.size()-1;



					double simPos = ((strip[kMax].first*strip[kMax].second) + (strip[kNext].first*strip[kNext].second)) / (strip[kMax].second + strip[kNext].second);
					int cLayer = ix/(Nsquares*2);

					for(int m = 0 ; m < hitPos[cLayer].size(); m++)
						segmp->Fill(simPos-hitPos[cLayer][m]);

					cout<<"Simulated position: "<<simPos<<endl;

					k = k2;
						}
					}
		ix = i2;
		jx = j2;
				}
			}
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
	outFile->WriteTObject(segmp);
	outFile->Close();
	}
