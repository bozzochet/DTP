
#include "physics.h"
#include "vector2.h"
#include "progress.h"
#include "PosSim.h"
#include "TrCluster.hh"
#include "Geometry.h"
#include "measure.h"

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
#include "TString.h"
#include "TBranch.h"

#include "utils/GGSSmartLog.h"
#include "montecarlo/readers/GGSTRootReader.h"

#include <iostream>
#include <vector>

using namespace std;



int main(int argc, char **argv) {

  static const string routineName("analysis");
  GGSSmartLog::verboseLevel = GGSSmartLog::INFO; // Print only INFO messages or more important


  //geometry

  GGSTRootReader reader;
  if (!(reader.Open(argv[1]))) {
    std::cerr << "Cannot open input file " << argv[1] <<std::endl;
    return 1;
  }

  const GGSTGeoParams *GEO = reader.GetGeoParams();
  Geometry *geo = new Geometry;

  geo->Nlayers = GEO->GetIntGeoParam("Nlayers");
  geo->Nstrips = GEO->GetIntGeoParam("Nstrips");
  geo->Nrows = GEO->GetIntGeoParam("Nrows");
  geo->Nsquares = GEO->GetIntGeoParam("Nsquares");
  geo->pitch = 1e-2 * GEO->GetRealGeoParam("pitch");
  geo->thickness = 1e-3 * GEO->GetRealGeoParam("thickness");

  geo->squareSide = geo->pitch * ((double) geo->Nstrips);
  geo->Nladders = geo->Nsquares * geo->Nrows * geo->Nlayers;

  COUT(INFO) <<"Geometric parameters:" <<ENDL;
  COUT(INFO) <<"layers:                 " <<geo->Nlayers <<ENDL;
  COUT(INFO) <<"strips per ladder:      " <<geo->Nstrips <<ENDL;
  COUT(INFO) <<"ladders rows per layer: " <<geo->Nrows  <<ENDL;
  COUT(INFO) <<"squares per side:       " <<geo->Nsquares  <<ENDL;
  COUT(INFO) <<"implant pitch:          " <<geo->pitch  <<ENDL;
  COUT(INFO) <<"layers thickness:       " <<geo->thickness <<ENDL;
  COUT(INFO) <<"squares side:           " <<geo->squareSide <<ENDL;
  COUT(INFO) <<"ladders:                " <<geo->Nladders <<ENDL;

  COUT(INFO) <<ENDL;


  TString inputFileName = argv[2];
  COUT(INFO) <<"Opening input file " <<inputFileName <<"..." <<ENDL;
  auto inFile = TFile::Open(inputFileName);

  TString outFileName = "histos.root";
  COUT(INFO) <<"Recreating output file " <<outFileName <<"..." <<ENDL;
	TFile *outFile = new TFile(outFileName, "recreate");


  COUT(INFO) <<"Creating histos..." <<ENDL;

	TH1F *h = new TH1F("disttemp", "disttemp", 100, -0.5, 0.5);
	TH1F *h1 = new TH1F("bt", "bt", 1000, 0, 2);
	TH1F *h2 = new TH1F("btls", "btls", 1000, 0, 2000);
	TH1F *h3 = new TH1F("btls_log", "btls_log", 1000, 0, 4);
	TH1F *h4 = new TH1F("bt_sim", "bt_sim", 1000, 0, 2);
	TH1F *h5 = new TH1F("btls_sim", "btls_sim", 1500, 0, 15);
	TH1F *h5cut = new TH1F("btls_sim_cut", "btls_sim_cut", 1500, 0, 15);
	TH1F *h5nopri = new TH1F("btls_sim_nopri", "btls_sim_nopri", 1500, 0, 15);
	TH1F *h5nomip = new TH1F("btls_sim_nomip", "btls_sim_nomip", 1500, 0, 15);

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


  //energy

/*
  TH1F *hPrimEdep= new TH1F
  (
    "hPrimEdep", "energy deposited on primary hits;[eV];",1000, -0, 0
  );
  hPrimEdep->SetCanExtend(TH1::kAllAxes);
*/

  TH1F *hEdep = new TH1F
    ("hEdep", "energy deposited on a strip;energy [eV];", 1000, -0, 0);
  hEdep->SetCanExtend(TH1::kAllAxes);

  TH1F *hEdepMeas = new TH1F
  (
    "hEdepMeas", "energy deposited on a strip measures;energy [eV];",
    1000, -0, 0
  );
  hEdepMeas->SetCanExtend(TH1::kAllAxes);

  TH1F *hEdepRes = new TH1F
  (
    "hEdepRes",
    "energy deposited on a strip measurement resolution;[eV];",
    1000, -0, 0
  );
  hEdepRes->SetCanExtend(TH1::kAllAxes);


  //position

	TH1F *hPosRes = new TH1F
    ("hPosRes", "position measurement resolution;[m];", 1000, -0, 0);
  hPosRes->SetCanExtend(TH1::kAllAxes);


  //time

  TH1F *hTimeMeas15= new TH1F
  (
    "hTimeMeas15", "time measures with threshold at 15%;[s];",
    1000, -0, 0
  );
  hTimeMeas15->SetCanExtend(TH1::kAllAxes);

  TH1F *hTimeHit = new TH1F("hTimeHit", "hit times;[s];", 1000, -0, 0);
  hTimeHit->SetCanExtend(TH1::kAllAxes);

  TH1F *hTimeRes15 = new TH1F
  (
    "hTimeRes15",
    "time measurement resolution with threshold at 15%;[s];",
    1000, -0, 0
  );
  hTimeRes15->SetCanExtend(TH1::kAllAxes);


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

  COUT(INFO) <<"Opening TTree object in " <<inputFileName <<"..." <<ENDL;

  TTree *events;
	inFile->GetObject("Data", events);
	events->Print();

	TClonesArray *a = new TClonesArray("TrCluster", 200);
	events->SetBranchAddress("Events", &a);
  TBranch *branch_ev = events->GetBranch("Events");

  measure meas;
  events->SetBranchAddress("Measures", &meas);
  TBranch *branch_meas = events->GetBranch("Measures");


  COUT(INFO) <<"Begin loop over " <<events->GetEntries() <<ENDL;

  int iMeas = 0; //iterator for branch_meas

  for (int i = 0; i < events->GetEntries(); i++) {

		v.resize(geo->Nlayers);

    branch_ev->GetEntry(i);

    /*
		if (a->GetEntries()>10) {
			printf("\nEvent %d: %d hits\n", i, a->GetEntries());
	   		for (int j = 0; j < a->GetEntries(); j++) {
				TrCluster *cl = (TrCluster *)a->At(j);
				printf("%d) %d %f\n", j, cl->parID, cl->eDep);
				}
		}
    */

		for (int j = 0; j < a->GetEntries(); j++) {

      //cout<<endl<<"Entry #"<<i+j<<endl;
			TrCluster *cl = (TrCluster *)a->At(j);

			v[cl->layer].push_back(*cl);

			//if(cl->parID == 0) hPrimEdep->Fill(cl->eDep); //primary

      hTimeHit->Fill(cl->time);


      /* while Events branch work with two indexes (i,j),
       * Measures branch was filled with one index and contains
       * a struct instead of an array. Because of this j element of
       * entry i in Events branch, corresponds to measure on entry
       * iMeas of Measures branch. Look in Digitization.cpp,
       * digitization function */

      branch_meas->GetEntry(iMeas);
      ++iMeas;

      for(int m=0; m<2; ++m)
      {
        if(meas.time[m] >= 0)
        {
          hTimeMeas15->Fill(meas.time[m]);
          hTimeRes15->Fill(meas.time[m] - cl->time);
        }

        hEdep->Fill(cl->clust[m]);
        hEdepMeas->Fill(meas.energy[m]);
        hEdepRes->Fill(meas.energy[m] - cl->clust[m]);
      }

      hPosRes->Fill(meas.position - cl->pos[cl->xy]);

		} //for j
  } //for i


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


  COUT(INFO) <<"Particle identification..." <<ENDL;

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


  COUT(INFO) <<"Writing output..." <<ENDL;

	outFile->WriteTObject(h);
	outFile->WriteTObject(h1);
	outFile->WriteTObject(h2);
	outFile->WriteTObject(h3);
	outFile->WriteTObject(h4);
	outFile->WriteTObject(h5);
	outFile->WriteTObject(h5cut);
	outFile->WriteTObject(h5nopri);
	outFile->WriteTObject(h5nomip);

  outFile->WriteTObject(hprotons);
	outFile->WriteTObject(hneutrons);
	outFile->WriteTObject(hgamma);
	outFile->WriteTObject(hisotopes);
	outFile->WriteTObject(helectron);
	outFile->WriteTObject(hpositron);
	outFile->WriteTObject(helectronmu);
	outFile->WriteTObject(hpi);
	outFile->WriteTObject(hk);

  //outFile->WriteTObject(hPrimEdep);
  outFile->WriteTObject(hEdep);
  outFile->WriteTObject(hEdepMeas);
  outFile->WriteTObject(hEdepRes);

  outFile->WriteTObject(hPosRes);

  outFile->WriteTObject(hTimeHit);
  outFile->WriteTObject(hTimeMeas15);
  outFile->WriteTObject(hTimeRes15);

  outFile->Close();

  COUT(INFO) <<"Output written in " <<outFileName <<ENDL;
	}
