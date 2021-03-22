//#include "DEBUG.h"
#include "physics.h"
#include "vector2.h"
#include "info.h"
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

  //MC

  TH1F *h_time_MC = new TH1F
    ("h_time_MC", ";log10(t / ns);", 1000, 0, 4);

  TH1F *h_time_MC_slow = new TH1F
    ("h_time_MC_slow", ";log10(t / ns);", 1000, 0, 4);

  TH1F *h_energy_MC = new TH1F
    ("h_energy_MC", ";energy [MeV];", 1000, -150, 150);


  //measures

  TH1F *h_time_meas15= new TH1F
    ("h_time_meas15", ";log10(t / ns);", 1000, 0, 4);

  TH1F *h_time_meas15_slow = new TH1F
    ("h_time_meas15_slow", ";log10(t / ns);", 1000, 0, 4);

  TH1F *h_energy_meas = new TH1F
    ("h_energy_meas", ";energy [MeV];", 1000, 0, 150);


  //resolutions

  TH1F *h_time_res15 = new TH1F
    ("h_time_res15", ";t_meas - t_true [ns];", 1000, 0, 1);

  TH1F *h_energy_res = new TH1F
    ("h_energy_res", ";E_meas - E_true [keV];", 1000, -150, 150);

  TH1F *h_pos_res = new TH1F
    ("h_pos_res", ";x_meas - x_true [cm];", 1000, -20, 20);

  //end of histos


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

  geo.ComputeDerived();
  
  COUT(INFO) <<ENDL;
  COUT(INFO) <<"=================================" <<ENDL;
  COUT(INFO) <<"Geometric parameters:     " <<ENDL;
  COUT(INFO) <<"  Calo side:              " <<geo.CaloSide <<ENDL;
  COUT(INFO) <<"  Calo-Stk gap:           " <<geo.CaloStkGap <<ENDL;
  COUT(INFO) <<"  wafers per side:        " <<geo.Nsquares  <<ENDL;
  COUT(INFO) <<"  ladders per 'column':   " <<geo.Nrows  <<ENDL;
  COUT(INFO) <<"  layers:                 " <<geo.Nlayers <<ENDL;
  COUT(INFO) <<"  gap between layers:     " <<geo.LayerGap <<ENDL;
  COUT(INFO) <<"  gap between planes:     " <<geo.PlaneGap <<ENDL;
  COUT(INFO) <<"  strips per ladder:      " <<geo.Nstrips <<ENDL;
  COUT(INFO) <<"  implant pitch:          " <<geo.pitch  <<ENDL;
  COUT(INFO) <<"  layers thickness:       " <<geo.thickness <<ENDL;
  COUT(INFO) <<"  wafer side:             " <<geo.squareSide <<ENDL;
  COUT(INFO) <<"  total # of ladders:     " <<geo.Nladders <<ENDL;
  COUT(INFO) <<"=================================" <<ENDL;

  COUT(INFO) <<ENDL;
  COUT(INFO) <<"Begin loop over " <<events_tree->GetEntries() <<ENDL;


  int iMeas = 0; //iterator for meas_tree


  //lost measures counters

  int energy_lost = 0;
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
    //and h_time_MC_slow
    std::vector<mytime_t> v_slow;
    std::vector<mytime_t> v_slow_meas;


    for (int j = 0; j < a->GetEntries(); j++) {

      //cout<<endl<<"Entry #"<<i+j<<endl;
      TrCluster *cl = (TrCluster *)a->At(j);

      v[cl->layer].push_back(*cl);

      h_time_MC->Fill(TMath::Log10(1e+9 * cl->time));

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
	  h_energy_MC->Fill(cl->clust[m] * 1e-6);

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
	      h_energy_meas->Fill(meas.energy[m] * 1e-6);
	      h_energy_res->Fill(1e-3 * (meas.energy[m] - cl->clust[m]));
	    }
	  else
	    ++energy_lost;

	  // In the next "if" is used for energy a threshold
	  // proportional to CHARGE_NOISE_DEV_ variable defined in
	  // TimeSim.h.
	  // Digitization executable does not save TimeSim object
	  // parameters used to generate time measures.
	  // Would be better that Digitization saves TimeSim parameters
	  // to read them in analysis.

	  if(meas.time[m] >= 0 && meas.energy[m] > 8*300*ENERGY_COUPLE)
	    {
	      h_time_meas15->Fill(TMath::Log10(1e+9 * meas.time[m]));
	      h_time_res15->Fill(1e+9 * (meas.time[m] - cl->time));

	      v_slow_meas.push_back(meas.time[m]);
	    }

	} //for m


      //read only valid position measures without lost ones

      if(TMath::Abs(meas.position) < 1)
	// this "if" is a temporary fix for Digitization bug:
	// Digitization.cpp,  line 424
        h_pos_res->Fill(1e+2 * (meas.position - cl->pos[cl->xy]));
      else
        ++position_lost;

    } //for j


    //fill slow hit

    h_time_MC_slow->Fill
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

  for (int il = 0; il<(int)(v.size()); il++) {
    for (int hit = 0; hit<(int)(v[il].size()); hit++) {

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

  outFile->WriteTObject(h_time_res15);
  outFile->WriteTObject(h_time_meas15);
  outFile->WriteTObject(h_time_meas15_slow);

  outFile->WriteTObject(h_energy_res);
  outFile->WriteTObject(h_energy_meas);

  //outFile->WriteTObject(h_pos_res);

  outFile->Close();

  COUT(INFO) <<"Output written in " <<outFileName <<ENDL;
  //debug::end_debug();

  return 0;
}
