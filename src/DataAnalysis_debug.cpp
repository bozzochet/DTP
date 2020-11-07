

#include "DEBUG.h"

#include "globals_and_types.h"
#include "vector2.h"
#include "progress.h"
#include "Stopwatch.h"
#include "TimeSimulation.h"
#include "TrCluster.hh"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"
#include "TTreeReader.h"

#include "utils/GGSSmartLog.h"

#include <iostream>
#include <vector>

using namespace std;

void stripReset(vector2<double> &array) {
	for (int ii = 0; ii < array.size(); ii++)
		for (int jj = 0; jj < array[ii].size(); jj++)
			array[ii][jj] = 0;
}

void addNoise(vector2<double> &array, TRandom3* tr) {

	for (int ii = 0; ii < array.size(); ii++) {
		for (int jj = 0; jj < array[ii].size(); jj++) {
			double fluct = tr->Gaus(0, 9e-6);

			array[ii][jj] += fluct;
			}
		}
}

void hitReset(vector2<double> &array, int dim) {
		for(int i=0; i<array.size();i++) {
			array[i].clear();
			array[i].shrink_to_fit();
		}
}

void shareEnergy(vector2<double> &array, int jump) {
	vector<double> fill;
	for (int ix = 0; ix < array.size(); ix++) {
		for (int jx = 0; jx < array[0].size(); jx+=jump) {
			fill.clear();
			fill.shrink_to_fit();
			int i0 = ix;
			int j0 = jx;

			//Saving the energy from non-active strips

			for(int jp = 1; jp<jump; jp++) {

				if(j0!=array[0].size()-1) {
					j0++;
				}

				else {
					i0++;
					j0 = 0;
				}

				if(i0==array.size())
					break;

				if(array[i0][j0]!=0) {
					fill.push_back(array[i0][j0]);
					array[i0][j0] = 0;
					}
			}

			//Moving to right-side strip

			if(j0!=array.size()-1) {
				j0++;
			}
			else if( (i0+1) % Nsquares != 0) {
						i0++;
						j0 = 0;
					}
			else {
			/* no active strips between (ix,jx) and layer
			* row end: distribute the entire energy to
			* (ix,jx) strip
			*/
			i0 = ix;
			j0 = jx;
			}

        	//Distributing the energy to left and right-side strips

			for(int s = 0; s < fill.size(); s++) {
				array[ix][jx] += fill[s]/(2*(s+1));
				array[i0][j0] += fill[s]/(2*(fill.size()-s));
			}

		}
	}
}

int main(int argc, char **argv) {

  debug::start_debug();

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

  std::vector<TGraph*> vec_current;

	TRandom3 *tr = new TRandom3();
	tr->SetSeed(time(NULL));
	vector2<TrCluster> v;

	/*

	Passing the parameters from DetectorConstruction.cc, not working yet

	GGSTRootReader reader;
	const GGSTGeoParams *geoParams = reader.GetGeoParams();

	const int Nsquares = geoParams->GetRealGeoParam("Nsquares"); //squares per side
	const int Nladders = Nsquares*5; //number of ladders
	const int Nstrips = geoParams->GetRealGeoParam("Nstrips"); //strips per ladder
	const double squareSide = geoParams->GetRealGeoParam("squareSide");
	const double pitch = squareSide/(double(Nstrips));
	double eDepSegm[Nladders][Nstrips];
	double PrimeDepSegm[Nladders][Nstrips];

	*/

  const int jump = 3;
	vector2<double> eDepSegm(Nladders, vector<double>(Nstrips));
	vector2<double> hitPos(Nlayers);

  TimeSimulation *time_sim = new TimeSimulation();

  cout <<endl <<"Begin analysis of " <<events->GetEntries()
    <<" events:\n";

  for (int i = 0; i < events->GetEntries(); i++) {

    if(i!=0) break;
    debug::out <<"\n event: " <<i <<std::endl;

    //print and update progress bar
    progress(i, events->GetEntries());

    /* IMPORTANT: cout and printf MUST print strings beginning and
     * ending with a new line (i.e. "\n" or endl) to not overwrite
     * the bar */

		v.resize(nLayer);
		events->GetEntry(i);

    /*
		if (a->GetEntries()>10) {
			printf("\nEvent %d: %d hits\n", i, a->GetEntries());
	   		for (int j = 0; j < a->GetEntries(); j++) {
				TrCluster *cl = (TrCluster *)a->At(j);
				printf("%d) %d %f\n", j, cl->parID, cl->eDep);
				}
		}
    */

		// eDepSegm array initialisation
		stripReset(eDepSegm);
		hitReset(hitPos,Nlayers);

    time_sim->Reset();
    debug::out <<"\nreset\n";

		for (int j = 0; j < a->GetEntries(); j++) {

      debug::out <<"\nentry: " <<j <<std::endl;

			//cout<<endl<<"Entry #"<<i+j<<endl;
			TrCluster *cl = (TrCluster *)a->At(j);

      debug::out <<std::endl;
      debug::print_cl(cl);

      v[cl->layer].push_back(*cl);
			if(cl->parID == 0) hPrimEdep->Fill(cl->eDep); //primary
			if(cl->eDep > 9e-6) hEdep->Fill(cl->eDep); //total

			//Filling the hits array

			hitPos[cl->layer].push_back(cl->pos[cl->segm]);

      // taking hit times
      debug::out <<"\nsplit: ";
      time_sim->SetHit(cl);
      debug::out <<cl->time <<std::endl;

      //get signal example
      if(j==0 && i==0)
      {
        time_sim->GetSignal(vec_current, cl->ladder, cl->strip);
      }

			//Filling the strips with the current energy
			eDepSegm[cl->ladder][cl->strip] += cl->clust[0];

			if(cl->strip == Nstrips-1 && (cl->ladder+1) % Nsquares == 0) //hit on the last strip of the last ladder of the layer row

				continue; //cl->clust[1] energy is lost

			else if(cl->strip==Nstrips-1)
						eDepSegm[cl->ladder+1][0] += cl->clust[1];
					else
						eDepSegm[cl->ladder][cl->strip+1] += cl->clust[1];
		}

		// Sharing of the energy from non-active strips

		if(jump!=1)
      shareEnergy(eDepSegm,jump);

		addNoise(eDepSegm,tr);

		for (int ix = 0; ix < Nladders; ix++) {
			for (int jx = 0; jx < Nstrips; jx+=jump) {


			//Find the boundaries of the clusters

			if(eDepSegm[ix][jx] >= 27e-6) {
				//cout<<"Analysing cluster\n";
				vector_pair<double> strip;

				int i1, j1;
				bool firstPoint = true;
				for(int ii = ix; ii>=0; ii--) {
					for(int jj = Nstrips-1; jj>=jump; jj-= jump) {
						if(firstPoint) {
							jj = jx;
							firstPoint = false;
						}
						i1 = ii;
						j1 = jj;

						/* evaluate if the strip is the first of the first ladder
						* on the layer row
						*/
						bool first_ladder = ii % Nsquares == 0 && jj == 0;

						if(eDepSegm[ii][jj] < 9e-6 || first_ladder)
							goto nextPart1;
					}
				}

				nextPart1:

				firstPoint = true;
				int i2, j2;
				for(int ii = ix; ii<Nladders; ii++) {
					for(int jj = 0; jj<Nstrips; jj+= jump) {
						if(firstPoint) {
							jj = jx;
							firstPoint = false;
						}
						i2 = ii;
						j2 = jj;

						/* evaluate if the strip is the first of the first ladder
						* on the layer row
						*/
						bool first_ladder = ii % Nsquares == 0 && jj == 0;

						if(eDepSegm[ii][jj] < 9e-6 || first_ladder)
							goto nextPart2;
					}
				}

				nextPart2:

				//Filling vector with current cluster

				while((i1*Nstrips)+j1 <= (i2*Nstrips)+j2) {

          if(j1>Nstrips-1 && i1 == Nladders-1)
            break;
          else if(j1>Nstrips-1) {
            i1++;
            j1 = 0;
          }

					double thisPos = ((i1%Nsquares)*squareSide) + (j1*pitch) - (Nsquares*squareSide*0.5);
					strip.push_back(make_pair(thisPos,eDepSegm[i1][j1]));

					j1+=jump;
				}

				for(int k = 0; k < strip.size(); k++) {
					if(strip[k].second >= 27e-6) {

					//Finding the boundaries of the peak

					int k1 = k-1;
					int k2 = k+1;

					while(k1>0 && strip[k1].second >= 27e-6)
						k1--;
					while(k2<(strip.size()-1) && strip[k2].second >= 27e-6)
						k2++;

					//Finding the peak

					int kMax = k1;
					for(int kHold = k1; kHold <= k2; kHold++)
						if(strip[kHold].second>strip[kMax].second)
							kMax = kHold;


					//cout<<"Analysing peak of "<< strip[kMax].second <<"keV\n";

					//Finding neighbour strip

					double kNext;

          if(kMax == 0)
            kNext = 1;
          else if(kMax == strip.size()-1)
            kNext = strip.size()-2;
					else if(strip[kMax+1].second > strip[kMax-1].second)
						kNext = kMax+1;
					else
						kNext = kMax-1;


					//Finding simulated position

					double simPos = ((strip[kMax].first*strip[kMax].second) + (strip[kNext].first*strip[kNext].second)) / (strip[kMax].second + strip[kNext].second);
					int cLayer = ix/(Nsquares*2);

					//Comparing the simulated hit positions with the real ones on the same layer

					for(int m = 0 ; m < hitPos[cLayer].size(); m++)
						segmp->Fill(simPos-hitPos[cLayer][m]);

					//cout<<"Simulated position: "<<simPos<<endl;



					//Advancing within the current cluster
					k = k2;
						}
					}

		//Advancing and looking for other clusters
		ix = i2;
		jx = j2;
				}
			}
	}
}

  //time simulation ended
  delete time_sim;

  cout <<endl <<endl;

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

  for(int i = 0; i < (int) vec_current.size(); ++i)
    outFile->WriteTObject(vec_current[i]);

  outFile->Close();

  debug::end_debug();

	}
