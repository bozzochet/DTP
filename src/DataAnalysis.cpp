
//#include "DEBUG.h"
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
#include "TMath.h"

#include "utils/GGSSmartLog.h"
#include "montecarlo/readers/GGSTRootReader.h"

#include <iostream>
#include <vector>

using namespace std;



int main(int argc, char **argv) {

  //debug::start_debug(); //DEBUG.h

  static const string routineName("DataAnalysis::main");
  GGSSmartLog::verboseLevel = GGSSmartLog::INFO; // Print only INFO messages or more important


  TString inputFileName = argv[1];
  COUT(INFO) <<"Opening input file " <<inputFileName <<"..." <<ENDL;
  auto inFile = TFile::Open(inputFileName);


  TString outFileName;

  if(argc < 3)
    outFileName = "histos.root";
  else
    outFileName = argv[2];

  COUT(INFO) <<"Recreating output file " <<outFileName <<"..." <<ENDL;
  TFile *outFile = new TFile(outFileName, "recreate");

/*
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
*/


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

  TH1F *h_energy = new TH1F
    ("h_energy", "MC energy;energy [eV];", 1000, -150e+6, 150e+6);

  TH1F *h_energy_meas = new TH1F
    ("h_energy_meas", "energy measures;energy [eV];", 1000, 0, 150e+6);

  TH1F *h_energy_res = new TH1F
  (
    "h_energy_res",
    "resolution of energy measurement;E_meas - E_true [eV];",
    1000, -150e+3, 150e+3
  );


  //position

	TH1F *h_pos_res = new TH1F
  (
    "h_pos_res", "resolution position measurement;x_meas - x_true[m];",
    1000, -1, 1
  );


  //time

  TH1F *h_time = new TH1F
    ("h_time", "hit times;log10(t / ns);", 1000, 0, 10);

  TH1F *h_time_meas15= new TH1F
  (
    "h_time_meas15", "time measures (threshold 15%);log10(t / ns);",
    1000, 0, 10
  );

  TH1F *h_time_res15 = new TH1F
  (
    "h_time_res15",
    "resolution of time measurement (threshold 15%);t_meas - t_true [ns];",
    1000, 0, 2
  );


  TH1F *h_time_HIGH_meas15= new TH1F
  (
    "h_time_HIGH_meas15",
    "time measures (E > 0, threshold 15%);log10(t / ns);",
    1000, 0, 10
  );

  TH1F *h_time_HIGH_res15 = new TH1F
  (
    "h_time_HIGH_res15",
    "resolution of time measurement (E > 0, threshold 15%);t_meas - t_true [ns];",
    1000, 0, 2
  );


  //hit time backscattered particles (electrons and protons)

  TH1F *h_time_slow = new TH1F
  (
    "h_time_slow",
    "slowest hit times;log10(t / ns);",
    1000, 0, 10
  );

  TH1F *h_time_meas15_slow = new TH1F
  (
    "h_time_meas15_slow",
    "slowest hit times;log10(t / ns);",
    1000, 0, 10
  );


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


  COUT(INFO) <<"Opening TTree objects in " <<inputFileName <<"..." <<ENDL;

  TTree *events_tree;
  inFile->GetObject("events", events_tree);
  events_tree->Print();

  TClonesArray *a = new TClonesArray("TrCluster", 200);
  events_tree->SetBranchAddress("Events", &a);


  TTree *meas_tree;
  inFile->GetObject("measures", meas_tree);

  measure meas;
  meas_tree->SetBranchAddress("Measures", &meas);


  TTree *geo_tree;
  inFile->GetObject("geometry", geo_tree);

  Geometry geo;
  geo_tree->SetBranchAddress("Geometry", &geo);
  geo_tree->GetEntry(0);

  COUT(INFO) <<ENDL;
  COUT(INFO) <<"=================================" <<ENDL;
  COUT(INFO) <<"Geometric parameters:     " <<ENDL;
  COUT(INFO) <<"  layers:                 " <<geo.Nlayers <<ENDL;
  COUT(INFO) <<"  strips per ladder:      " <<geo.Nstrips <<ENDL;
  COUT(INFO) <<"  ladders rows per layer: " <<geo.Nrows  <<ENDL;
  COUT(INFO) <<"  squares per side:       " <<geo.Nsquares  <<ENDL;
  COUT(INFO) <<"  implant pitch:          " <<geo.pitch  <<ENDL;
  COUT(INFO) <<"  layers thickness:       " <<geo.thickness <<ENDL;
  COUT(INFO) <<"  squares side:           " <<geo.squareSide <<ENDL;
  COUT(INFO) <<"  ladders:                " <<geo.Nladders <<ENDL;
  COUT(INFO) <<"=================================" <<ENDL;


  COUT(INFO) <<ENDL;
  COUT(INFO) <<"Begin loop over " <<events_tree->GetEntries() <<ENDL;


  int iMeas = 0; //iterator for meas_tree


  //lost measures counters

  int energy_lost = 0;
  int time_lost = 0;
  int position_lost = 0;


  for (int i = 0; i < events_tree->GetEntries(); i++) {

		v.resize(geo.Nlayers);

    events_tree->GetEntry(i);

    /*
		if (a->GetEntries()>10) {
			printf("\nEvent %d: %d hits\n", i, a->GetEntries());
	   		for (int j = 0; j < a->GetEntries(); j++) {
				TrCluster *cl = (TrCluster *)a->At(j);
				printf("%d) %d %f\n", j, cl->parID, cl->eDep);
				}
		}
    */


    //measured times for particle i;
    //the slowest is used afterwards to fill h_time_meas15_slow
    //and h_time_slow
    std::vector<mytime_t> v_slow;
    std::vector<mytime_t> v_slow_meas;


		for (int j = 0; j < a->GetEntries(); j++) {

      //cout<<endl<<"Entry #"<<i+j<<endl;
			TrCluster *cl = (TrCluster *)a->At(j);

			v[cl->layer].push_back(*cl);

      h_time->Fill(TMath::Log10(1e+9 * cl->time));

      v_slow.push_back(cl->time);


      /* while Events branch work with two indexes (i,j),
       * Measures branch was filled with one index and contains
       * a struct instead of an array. Because of this j element of
       * entry i in Events branch, corresponds to measure on entry
       * iMeas of Measures branch. Look in Digitization.cpp,
       * digitization function */

      meas_tree->GetEntry(iMeas);
      ++iMeas;

      //scan energies clust and measures

      for(int m=0; m<2; ++m)
      {
        h_energy->Fill(cl->clust[m]);

/* DEBUG.h
        debug::out <<"\ni: " <<i <<" j: " <<j <<" m: " <<m;

        debug::out <<"\n\tE: " <<cl->clust[m];
        debug::out <<"\n\tE + noise: " <<meas.energy[m];
        debug::out <<"\n\tt: " <<cl->time;
        debug::out <<"\n\tt meas: " <<meas.time[m];
        debug::out <<"\n\tpos: " <<cl->pos[cl->xy];
        debug::out <<"\n\tpos meas: " <<meas.position;

        debug::out <<std::endl;
*/

        //analyze valid measures

        if(meas.energy[m] > 0)
        {
          h_energy_meas->Fill(meas.energy[m]);
          h_energy_res->Fill(meas.energy[m] - cl->clust[m]);
        }
        else
          ++energy_lost;


        if(meas.time[m] >= 0)
        {
          h_time_meas15->Fill(TMath::Log10(1e+9 * meas.time[m]));
          h_time_res15->Fill(1e+9 * (meas.time[m] - cl->time));

          v_slow_meas.push_back(meas.time[m]);
        }
        else
          ++time_lost;


        if(meas.time[m] >= 0 && meas.energy[m] > 0)
        {
          h_time_HIGH_meas15->Fill(TMath::Log10(1e+9 * meas.time[m]));
          h_time_HIGH_res15->Fill(1e+9 * (meas.time[m] - cl->time));
        }

      } //for m


      //read only valid position measures without lost ones

      if(TMath::Abs(meas.position) < 1)
      // this if is also a temporary fix for Digitization bug:
      // Digitization.cpp,  line 424
        h_pos_res->Fill(meas.position - cl->pos[cl->xy]);
      else
        ++position_lost;

		} //for j


    //fill slow hit

    h_time_slow->Fill
    (
      TMath::Log10(1e+9 * TMath::MaxElement(v_slow.size(), &v_slow[0]))
    );

    h_time_meas15_slow->Fill
    (
      TMath::Log10
        (1e+9 * TMath::MaxElement(v_slow_meas.size(), &v_slow_meas[0]))
    );

  } //for i


  COUT(INFO) <<ENDL;

  COUT(INFO) <<"Lost energies:  " <<energy_lost <<" on " <<iMeas*2
    <<ENDL;

  COUT(INFO) <<"Lost times:     " <<time_lost <<" on " <<iMeas*2
    <<ENDL;

  COUT(INFO) <<"Lost positions: " <<position_lost <<" on " <<iMeas
    <<ENDL;


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


  COUT(INFO) <<ENDL;
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

/*
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
*/

  outFile->WriteTObject(h_energy);
  outFile->WriteTObject(h_energy_meas);
  outFile->WriteTObject(h_energy_res);

  outFile->WriteTObject(h_pos_res);

  outFile->WriteTObject(h_time);
  outFile->WriteTObject(h_time_meas15);
  //outFile->WriteTObject(h_time_HIGH_meas15);

  outFile->WriteTObject(h_time_slow);
  outFile->WriteTObject(h_time_meas15_slow);

  outFile->WriteTObject(h_time_res15);
  //outFile->WriteTObject(h_time_HIGH_res15);

  outFile->Close();

  COUT(INFO) <<"Output written in " <<outFileName <<ENDL;
  //debug::end_debug();

  return 0;
}
