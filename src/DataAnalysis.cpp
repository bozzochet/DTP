
#include "physics.h"
#include "vector2.h"
#include "progress.h"
#include "PosSimulation.h"
#include "TimeSimulation.h"
#include "TrCluster.hh"
#include "Geometry.h"
#include "TimeSegm.h"

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
#include "TMultiGraph.h"

#include "utils/GGSSmartLog.h"

#include <iostream>
#include <vector>

using namespace std;


//global geometric parameters

Geometry GEO;

const int Nlayers = GEO.GetNlayers();
const int Nstrips = GEO.GetNstrips();
const int Nrows = GEO.GetNrows();
const int Nsquares = GEO.GetNsquares();
const int pitch = GEO.GetPitch();
const int squareSide = GEO.GetSquareSide();
const int Nladders = GEO.GetNladders();
const double thickness = GEO.GetThickness();



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

/*
  TGraph *current_ideal;
  TGraph *current_noise;

  TGraph *charge_ideal;
  TGraph *charge_noise;

  TMultiGraph *charge = new TMultiGraph("charge", "charge");
  TMultiGraph *current = new TMultiGraph("current", "current");
*/

  TH1F *htime = new TH1F
    ("htime", "timing; ; entries", 1000, -0.05, 0.05);

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

  PosSimulation *pos_sim = new PosSimulation(&GEO, 2, tr);

  TimeSegm *time_segm = new TimeSegm(&GEO, A, 1);

  //thickness is given in mm; TimeSimulation wants m
  TimeSimulation *time_sim =
    new TimeSimulation(time_segm, thickness*1e-3);

  cout <<endl <<"Begin analysis of " <<events->GetEntries()
    <<" events:\n";

  for (int i = 0; i < events->GetEntries(); i++) {

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

		pos_sim->Reset();

		for (int j = 0; j < a->GetEntries(); j++) {

      //cout<<endl<<"Entry #"<<i+j<<endl;
			TrCluster *cl = (TrCluster *)a->At(j);
			v[cl->layer].push_back(*cl);
			if(cl->parID == 0) hPrimEdep->Fill(cl->eDep); //primary
			if(cl->eDep > 9e-6) hEdep->Fill(cl->eDep); //total

/*
      if(i==0 && j==0)
      {
        charge_ideal = new TGraph();
        time_sim->GetChargeSignal(charge_ideal, cl->eDep*1e+9, false);
        charge_ideal->SetNameTitle("charge_ideal", "ideal charge");

        current_ideal = new TGraph();
        time_sim->GetCurrentSignal(current_ideal, charge_ideal, cl->time*1e-9);
        current_ideal->SetNameTitle("current_ideal", "ideal current");

        charge_noise = new TGraph(*charge_ideal); //copy ideal
        time_sim->AddChargeNoise(charge_noise);
        charge_noise->SetNameTitle("charge_noise", "charge with noise");

        current_noise = new TGraph();
        time_sim->GetCurrentSignal(current_noise, charge_noise, cl->time*1e-9);
        current_noise->SetNameTitle("current_noise", "current with noise");

        charge_noise->SetLineColor(kBlue);
        current_noise->SetLineColor(kBlue);

        charge_ideal->SetLineColor(kRed);
        current_ideal->SetLineColor(kRed);

        charge_noise->SetLineWidth(2);
        current_noise->SetLineWidth(2);

        charge_ideal->SetLineWidth(2);
        current_ideal->SetLineWidth(2);

        charge->Add(charge_noise);
        charge->Add(charge_ideal);

        current->Add(current_noise);
        current->Add(current_ideal);
      }
*/

      /****************************************
       * BUG: strip < 0 from TrCluster object *
       ****************************************/

      if(cl->ladder >= 0 && cl->strip >= 0) //BUG TEMPORARY FIX
        time_segm->SetHit
          (cl->ladder, cl->strip, cl->time * 1e-9, cl->eDep * 1e+9);

      pos_sim->SetHitPos(cl->layer, cl->pos[cl->segm]);
      pos_sim->DepositEnergy(cl->ladder, cl->strip, cl->clust[0]);

			if(cl->strip == Nstrips-1 && (cl->ladder+1) % Nsquares == 0) //hit on the last strip of the last ladder of the layer row

				continue; //cl->clust[1] energy is lost

			else if(cl->strip==Nstrips-1)
        pos_sim->DepositEnergy(cl->ladder+1, 0, cl->clust[1]);
			else
				pos_sim->DepositEnergy(cl->ladder, cl->strip+1, cl->clust[1]);
		}

		// Sharing of the energy from non-active strips

    pos_sim->ShareEnergy();
    pos_sim->AddNoise();
    pos_sim->Segm(segmp);
}

  //time deviations

  for(int i=0; i < time_segm->GetNgroups(); ++i)
  {
    std::vector<mytime_t> true_time;
    time_segm->GetTimes(i, true_time);

    TGraph *current = new TGraph();
    time_sim->GetCurrentSignal(i, current);

    mytime_t meas_time = time_sim->GetTime(current, 0.1);

    delete current;

    for(int j = 0; j < (int) true_time.size(); ++j)
      htime->Fill( (meas_time - true_time[j]) / true_time[j] );
  }

  //simulations ended
  delete pos_sim;
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
/*
  outFile->WriteTObject(current_ideal);
  outFile->WriteTObject(current_noise);
  outFile->WriteTObject(charge_ideal);
  outFile->WriteTObject(charge_noise);
  outFile->WriteTObject(current);
  outFile->WriteTObject(charge);
*/
  outFile->WriteTObject(htime);
  outFile->Close();
	}
