//#define _DEBUG_

#ifdef _DEBUG_
#include "DEBUG.h"
#endif
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

#define SIGNAL_THRESHOLD 8.0*300.0*ENERGY_COUPLE
#define TIME_RESOLUTION 1.0E-9*0.1

int main(int argc, char **argv) {

#ifdef _DEBUG_
  debug::start_debug(); //DEBUG.h
#endif
  
  static const string routineName("DataAnalysis::main");
  GGSSmartLog::verboseLevel = GGSSmartLog::INFO; // Print only INFO messages or more important

  bool _calo_flag=false;
  TString ExeName = argv[0];
  if (ExeName.EndsWith("_calo")) _calo_flag=true;

  if (!_calo_flag) {
    if (argc<3) {
      COUT(ERROR) << "You need to pass at least an input file name and an output one" << ENDL;
      COUT(ERROR) << argv[0] << " [output file name] <input file name 1> <input file name 2> ..." << ENDL;
      return -1;
    }
  }
  else {
    if(argc < 6) {
      COUT(ERROR) << "You need to pass at least an input file name and an output one, the nominal energy and the selection ones" << ENDL;
      COUT(ERROR) << argv[0] << " [output file name] <input file name 1> <input file name 2> ... [nominal energy] [min deposited energy] [max deposited energy]" << ENDL;
      return -1;
    }
  } 
  
  TString outFileName = argv[1];
  
  TChain *events_tree = new TChain("events");
  TChain *meas_tree = new TChain("measures");
  TChain *calo_tree = new TChain("calorimeter");
  TChain *geo_tree = new TChain("geometry");

  int shift=2;
  for (shift=2; shift<argc; shift++) {
    TString inputFileName = argv[shift];
    //    printf("%d) %s\n", shift, argv[shift]);

    if (inputFileName.IsDigit()) {
      break;//input files are over
    }
    
    COUT(INFO) <<"Opening TTree objects in " <<inputFileName <<"..." <<ENDL;
    
    events_tree->Add(argv[shift]);
    meas_tree->Add(argv[shift]);
    if (_calo_flag) {
      calo_tree->Add(argv[shift]);
    }
    geo_tree->Add(argv[shift]);
  }

  
  COUT(INFO) <<"Recreating output file " <<outFileName <<"..." <<ENDL;
  TFile *outFile = new TFile(outFileName, "recreate");
  outFile->cd();

  COUT(INFO) <<"Creating histos..." <<ENDL;
  
  TH1F *pri_deltat_wrt_mean = new TH1F("pri_deltat_wrt_mean", "#{Delta}t wrt mean (primary)", 100, -0.5, 0.5);
  TH1F *pri_deltat_wrt_start = new TH1F("pri_deltat_wrt_start", "#{Delta}t wrt start (primary)", 1001, -2.0/1000.0, 2);
  TH1F *pri_deltat_wrt_start_longscale = new TH1F("pri_deltat_wrt_start_longscale", "#{Delta}t wrt start (primary)", 1001, -2000.0/1000.0, 2000);
  TH1F *pri_deltat_wrt_start_log = new TH1F("pri_deltat_wrt_start_log", "#{Delta}t wrt start (primary) (log)", 1000, 0, 4);
  
  TH1F *deltat_wrt_start = new TH1F("deltat_wrt_start", "#{Delta}t wrt start", 1001, -2.0/1000.0, 2);
  TH1F *deltat_wrt_start_longscale = new TH1F("deltat_wrt_start_longscale", "#{Delta}t wrt start", 1001, -2000.0/1000.0, 2000);
  TH1F *deltat_wrt_start_log = new TH1F("deltat_wrt_start_log", "#{Delta}t wrt start (log)", 1000, 0, 4);
  
  TH1F *deltat_smeared_wrt_start = new TH1F("deltat_smeared_wrt_start", "#{Delta}t (smeared) wrt start", 1001, -2.0/1000.0, 2);
  TH1F *deltat_smeared_wrt_start_longscale = new TH1F("deltat_smeared_wrt_start_longscale", "#{Delta}t (smeared) wrt start", 1001, -2000.0/1000.0, 2000);
  TH1F *deltat_smeared_wrt_start_log = new TH1F("deltat_smeared_wrt_start_log", "#{Delta}t (smeared) wrt start (log)", 1000, 0, 4);

  TH1F *nopri_deltat_smeared_wrt_start = new TH1F("nopri_deltat_smeared_wrt_start", "#{Delta}t (smeared) wrt start (no primaries)", 1001, -2.0/1000.0, 2);
  TH1F *nopri_deltat_smeared_wrt_start_longscale = new TH1F("nopri_deltat_smeared_wrt_start_longscale", "#{Delta}t (smeared) wrt start (no primaries)", 1001, -2000.0/1000.0, 2000);
  TH1F *nopri_deltat_smeared_wrt_start_log = new TH1F("nopri_deltat_smeared_wrt_start_log", "#{Delta}t (smeared) wrt start (no primaries) (log)", 1000, 0, 4);
  
  TH1F *nomip_deltat_smeared_wrt_start = new TH1F("nomip_deltat_smeared_wrt_start", "#{Delta}t (smeared) wrt start (no MIP)", 1001, -2.0/1000.0, 2);
  TH1F *nomip_deltat_smeared_wrt_start_longscale = new TH1F("nomip_deltat_smeared_wrt_start_longscale", "#{Delta}t (smeared) wrt start (no MIP)", 1001, -2000.0/1000.0, 2000);
  TH1F *nomip_deltat_smeared_wrt_start_log = new TH1F("nomip_deltat_smeared_wrt_start_log", "#{Delta}t (smeared) wrt start (no MIP) (log)", 1000, 0, 4);

  TH1F *slower_deltat_wrt_start = new TH1F("slower_deltat_wrt_start", "#{Delta}t (slower) wrt start", 1001, -2.0/1000.0, 2);
  TH1F *slower_deltat_wrt_start_longscale = new TH1F("slower_deltat_wrt_start_longscale", "#{Delta}t (slower) wrt start", 1001, -2000.0/1000.0, 2000);
  TH1F *slower_deltat_wrt_start_log = new TH1F("slower_deltat_wrt_start_log", "#{Delta}t (slower) wrt start (log)", 1000, 0, 4);

  TH1F *slower_deltat_smeared_wrt_start = new TH1F("slower_deltat_smeared_wrt_start", "#{Delta}t (slower-smeared) wrt start", 1001, -2.0/1000.0, 2);
  TH1F *slower_deltat_smeared_wrt_start_longscale = new TH1F("slower_deltat_smeared_wrt_start_longscale", "#{Delta}t (slower-smeared) wrt start", 1001, -2000.0/1000.0, 2000);
  TH1F *slower_deltat_smeared_wrt_start_log = new TH1F("slower_deltat_smeared_wrt_start_log", "#{Delta}t (slower-smeared) wrt start (log)", 1000, 0, 4);
  
  TH1F *primaries = new TH1F("primaries", "primaries", 1000, 0, 4);
  primaries->SetLineColor(kBlack);
  primaries->SetMarkerColor(kBlack);
  
  TH1F *protons = new TH1F("protons", "protons", 1000, 0, 4);
  protons->SetLineColor(kGreen+2);
  protons->SetMarkerColor(kGreen+2);
  
  TH1F *antip = new TH1F("antip", "antiprotons", 1000, 0, 4);
  antip->SetLineColor(kYellow+2);
  antip->SetMarkerColor(kYellow+2);
  
  TH1F *neutrons = new TH1F("neutrons", "neutrons", 1000, 0, 4);
  neutrons->SetLineColor(kOrange+3);
  neutrons->SetMarkerColor(kOrange+3);
  
  TH1F *gamma = new TH1F("gamma", "gamma", 1000, 0, 4);
  gamma->SetLineColor(kCyan);
  gamma->SetMarkerColor(kCyan);
  
  TH1F *isotopes = new TH1F("isotopes", "isotopes", 1000, 0, 4);
  isotopes->SetLineColor(kOrange-8);
  isotopes->SetMarkerColor(kOrange-8);
  
  TH1F *electrons = new TH1F("electrons", "electrons", 1000, 0, 4);
  electrons->SetLineColor(kBlue);
  electrons->SetMarkerColor(kBlue);
  
  TH1F *positrons = new TH1F("positrons", "positrons", 1000, 0, 4);
  positrons->SetLineColor(kRed+2);
  positrons->SetMarkerColor(kRed+2);
  
  TH1F *muons = new TH1F("muons", "muons", 1000, 0, 4);
  muons->SetLineColor(kAzure+1);
  muons->SetMarkerColor(kAzure+1);
  
  TH1F *pions = new TH1F("pions", "pions", 1000, 0, 4);
  pions->SetLineColor(kOrange+7);
  pions->SetMarkerColor(kOrange+7);
  
  TH1F *kaons = new TH1F("kaons", "kaons", 1000, 0, 4);
  kaons->SetLineColor(kMagenta);
  kaons->SetMarkerColor(kMagenta);
  
  //in the following definitions:
  // - argv[3] must contain energy of the beam (in macros/run.mac)
  // - slow histos store slowest hits for each event
  // - meas histos store quantities measured by Time and Pos libraries
  // - energy_calo hist stores energy deposited in the calorimeter,
  //     while time_calo stores quantities referred to events where
  //     energy deposited in calorimeter is in the range [argv(4),argv(5)]
  
  //argv are expressed in GeV
  energy_t beam_energy = -999999;
  energy_t Ecalo_min = -999999;
  energy_t Ecalo_max = -999999;
  
  if (_calo_flag){
    beam_energy = std::atof(argv[shift]);
    Ecalo_min = std::atof(argv[shift+1]);
    Ecalo_max = std::atof(argv[shift+2]);
    printf("%f %f %f\n", beam_energy, Ecalo_min, Ecalo_max);
    beam_energy *= 1e+9;
    Ecalo_min *= 1e+9;
    Ecalo_max *= 1e+9;
  }
  
  //histos
  
  //MC

  TH1F *h_time_MC = new TH1F
    ("h_time_MC", ";log10(t / ns);", 1000, 0, 4);

  TH1F *h_time_MC_slow = new TH1F
    ("h_time_MC_slow", ";log10(t / ns);", 1000, 0, 4);

  TH1F *h_energy_MC = new TH1F
    ("h_energy_MC", ";energy [MeV];", 1000, -150, 150);


  //measures

  TH1F *h_time_meas15= new TH1F
    ("h_time_meas", ";log10(t / ns);", 1000, 0, 4);

  TH1F *h_time_meas15_slow = new TH1F
    ("h_time_meas_slow", ";log10(t / ns);", 1000, 0, 4);

  TH1F *h_energy_meas = new TH1F
    ("h_energy_meas", ";energy [MeV];", 1000, 0, 150);


  //resolutions

  TH1F *h_time_res15 = new TH1F
    ("h_time_res", ";t_meas - t_true [ns];", 1000, 0, 1);

  TH1F *h_energy_res = new TH1F
    ("h_energy_res", ";E_meas - E_true [keV];", 1000, -150, 150);

  TH1F *h_pos_res = new TH1F
    ("h_pos_res", ";x_meas - x_true [cm];", 1000, -20, 20);

  //calo

  TH1F *h_energy_calo = NULL;
  
  if (_calo_flag) {    
    h_energy_calo = new TH1F("h_energy_calo", ";energy [GeV];", 1000, 0, beam_energy * 1e-9);
  }
  
  //end of histos


  TRandom3 *tr = new TRandom3();
  tr->SetSeed(time(NULL));
  
  events_tree->Print();

  TClonesArray *a = new TClonesArray("TrCluster", 200);
  events_tree->SetBranchAddress("Events", &a);

  measure meas;
  meas_tree->SetBranchAddress("Measures", &meas);

  energy_t Ecalo;
  if (_calo_flag) {
    calo_tree->SetBranchAddress("Events", &Ecalo);
  }
  
  Geometry geo;
  Geometry* p_geo = &geo;
  geo_tree->SetBranchAddress("Geometry", &p_geo);
  /* in Apr 2021 the object was saved wrongly inm Digitazion: skipping for now. To be restored after the MDPI paper
  geo_tree->GetEntry(0);
  */

  geo.ComputeDerived();
  geo.Dump();  

  COUT(INFO) <<ENDL;
  COUT(INFO) <<"Begin loop over " <<events_tree->GetEntries() <<ENDL;

  int iMeas = 0; //iterator for meas_tree

  //lost measures counters

  int energy_lost = 0;//energy measured < 0
  int position_lost = 0;//|position measured| < 1 (to be re-thought)

  for (int i = 0; i < events_tree->GetEntries(); i++) {
    
    vector2<TrCluster> v;
    //    v.resize(geo.Nlayers);
    v.resize(10);//there's a bug and the values of geo are not retrieved correctly
    
    events_tree->GetEntry(i);

    if (_calo_flag) {
      calo_tree->GetEntry(i);
      
      h_energy_calo->Fill(Ecalo * 1e-9);
 
      if (Ecalo < Ecalo_min && Ecalo > Ecalo_max) continue;
    }
    
    /*
      if (a->GetEntries()>10) {
      printf("\nEvent %d: %d hits\n", i, a->GetEntries());
      for (int j = 0; j < a->GetEntries(); j++) {
      TrCluster *cl = (TrCluster *)a->At(j);
      printf("%d) %d %f\n", j, cl->parID, cl->eDep);
      }
      }
    */

    //measured times for primary/event i;
    //the slowest is used afterwards to fill h_time_meas15_slow
    //and h_time_MC_slow
    std::vector<mytime_t> v_times;
    std::vector<mytime_t> v_times_meas;

    for (int j = 0; j < a->GetEntries(); j++) {

      //cout<<endl<<"Entry #"<<i+j<<endl;
      TrCluster *cl = (TrCluster *)a->At(j);

      v[cl->layer].push_back(*cl);

      h_time_MC->Fill(TMath::Log10(1e+9 * cl->time));

      v_times.push_back(cl->time);

      /* while Events branch work with two indexes (i,j),
       * Measures branch was filled with one index and contains
       * a struct instead of an array. Because of this j element of
       * entry i in Events branch, corresponds to measure on entry
       * iMeas of Measures branch. Look in Digitization.cpp,
       * digitization function */

      meas_tree->GetEntry(iMeas);
      ++iMeas;

      //scan energies clust and measures

      for (int m=0; m<2; ++m) { //the clusters are two since there're the two strips around the hit position
	//most likely for the timing this is even correct (at leat when there's no "grouping")
	//but for the energy this is WRONG

	static double ene_true = 0;
	static double ene_meas = 0;
	if (m==0) {
	  ene_true = cl->clust[m];
	  ene_meas = meas.energy[m];
	}
	else {
	  ene_true += cl->clust[m];
	  ene_meas += meas.energy[m];
	}

	if (m==1) {
	  h_energy_MC->Fill(ene_true * 1e-6);
	}
	
#ifdef _DEBUG_
	debug::out <<"\ni: " <<i <<" j: " <<j <<" m: " <<m;
	
	debug::out <<"\n\tE: " <<cl->clust[m];
	debug::out <<"\n\tE + noise: " <<meas.energy[m];
	debug::out <<"\n\tt: " <<cl->time;
	debug::out <<"\n\tt meas: " <<meas.time[m];
	debug::out <<"\n\tpos: " <<cl->pos[cl->xy];
	debug::out <<"\n\tpos meas: " <<meas.position;
	
	debug::out <<std::endl;
#endif	 
	
	//analyze valid measures

	if (m==1) {
	  if(ene_meas > 0) {
	    h_energy_meas->Fill(ene_meas * 1e-6);
	    h_energy_res->Fill(1e-3 * (ene_meas - ene_true));
	  }
	  else {
	    ++energy_lost;
	  }
	}
	
	// In the next "if" is used for energy a threshold
	// proportional to CHARGE_NOISE_DEV_ variable defined in
	// TimeSim.h.
	// Digitization executable does not save TimeSim object
	// parameters used to generate time measures.
	// Would be better that Digitization saves TimeSim parameters
	// to read them in analysis.
	
	if(meas.time[m] >= 0 && meas.energy[m] > SIGNAL_THRESHOLD) {//the cut is on the single energy measurement since the zero suppression is applied on this
	  h_time_meas15->Fill(TMath::Log10(1e+9 * meas.time[m]));
	  h_time_res15->Fill(1e+9 * (meas.time[m] - cl->time));
	  
	  v_times_meas.push_back(meas.time[m]);
	}
	
      } //for m
      
      //read only valid position measures without lost ones
      
      if (TMath::Abs(meas.position) < 1)
	// this "if" is a temporary fix for Digitization bug:
	// Digitization.cpp,  line 424
        h_pos_res->Fill(1e+2 * (meas.position - cl->pos[cl->xy]));
      else
        ++position_lost;
      
    } //for j

    //fill slow hit

    mytime_t slower = TMath::MaxElement(v_times.size(), &v_times[0]);
    mytime_t slower_meas = TMath::MaxElement(v_times_meas.size(), &v_times_meas[0]);
    
    h_time_MC_slow->Fill(TMath::Log10(1e+9 * slower));
    h_time_meas15_slow->Fill(TMath::Log10(1e+9 * slower_meas));
    
    // Find tStart and tMean
    
    float tStart = v[0][0].time;
    float tMean = 0;
    int _n = 0;
    for (auto il : v) {
      for (auto hit : il) {
	if (tStart>hit.time) tStart = hit.time;
	if (hit.parID != 0) continue;
	tMean += hit.time;
	_n++;
      }
    }
    tMean /= _n;
        
    //Particle identification
    
    for (int il = 0; il<(int)(v.size()); il++) { //layer
      for (int hit = 0; hit<(int)(v[il].size()); hit++) { //hit	
	if (v[il][hit].clust[0]>SIGNAL_THRESHOLD || v[il][hit].clust[1]>SIGNAL_THRESHOLD) {
	  
	  if (v[il][hit].parID == 0) {
	    primaries->Fill(log10(1.0 + 1e+9 * (v[il][hit].time - tStart )));
	  }
	  else {
	    if (
		v[il][hit].parPdg == 2212 //proton
		)
	      protons->Fill(log10(1.0 + 1e+9 * (v[il][hit].time - tStart )));
	    
	    if (
		v[il][hit].parPdg == -2212 //antiproton
		)
	      antip->Fill(log10(1.0 + 1e+9 * (v[il][hit].time - tStart)));
	    
	    if (
		v[il][hit].parPdg == 2112 //neeutron
		)
	      neutrons->Fill(log10(1.0 + 1e+9 * (v[il][hit].time - tStart)));
	    
	    if (
		v[il][hit].parPdg == 11 //electron
		)
	      electrons->Fill(log10(1.0 + 1e+9 * (v[il][hit].time - tStart)));
	    
	    if (
		v[il][hit].parPdg == -11 //positron
		)
	      positrons->Fill(log10(1.0 + 1e+9 * (v[il][hit].time - tStart)));
	    
	    if (
		v[il][hit].parPdg == 211 //pi+
		||
		v[il][hit].parPdg == -211 //pi-
		)
	      pions->Fill(log10(1.0 + 1e+9 * (v[il][hit].time - tStart)));
	    
	    if (
		v[il][hit].parPdg == 130 //K0L
		||
		v[il][hit].parPdg == 310 //K0S
		||
		v[il][hit].parPdg == 311 //K0
		||
		v[il][hit].parPdg == 321 //K+ 
		||
		v[il][hit].parPdg == -321 //K-
		)
	      kaons->Fill(log10(1.0 + 1e+9 * (v[il][hit].time - tStart)));
	    
	    if (
		v[il][hit].parPdg == 13 //mu-
		||
		v[il][hit].parPdg == -13 //mu+
		)
	      muons->Fill(log10(1.0 + 1e+9 * (v[il][hit].time - tStart)));
	    
	    if (
		v[il][hit].parPdg == 22 //gamma
		)
	      gamma->Fill(log10(1.0 + 1e+9 * (v[il][hit].time - tStart)));
	    
	    if (
		v[il][hit].parPdg > 1000000000 //isotopes
		)
	      isotopes->Fill(log10(1.0 + 1e+9 * (v[il][hit].time - tStart)));
	  }
	  
	  //time
	  double deltat = v[il][hit].time - tStart;
	  double smeared_deltat = tr->Gaus(deltat, TIME_RESOLUTION);
	  double deltat_mean = v[il][hit].time - tMean;
	  
	  deltat_wrt_start->Fill(1e+9 * deltat);
	  deltat_wrt_start_longscale->Fill(1e+9 * deltat);
	  deltat_wrt_start_log->Fill(log10(1.0 + 1e+9 * deltat));
	  
	  deltat_smeared_wrt_start->Fill(1e+9 * smeared_deltat);
	  deltat_smeared_wrt_start_longscale->Fill(1e+9 * smeared_deltat);
	  deltat_smeared_wrt_start_log->Fill(log10(1.0 + 1e+9 * smeared_deltat));

	  if (v[il][hit].parID == 0) {
	    pri_deltat_wrt_mean->Fill(1e+9 * deltat_mean);
	    pri_deltat_wrt_start->Fill(1e+9 * deltat);
	    pri_deltat_wrt_start_longscale->Fill(1e+9 * deltat);
	    pri_deltat_wrt_start_log->Fill(log10(1.0 + 1e+9 * deltat));
	  }
	  else {
	    nopri_deltat_smeared_wrt_start->Fill(1e+9 * smeared_deltat);
	    nopri_deltat_smeared_wrt_start_longscale->Fill(1e+9 * smeared_deltat);
	    nopri_deltat_smeared_wrt_start_log->Fill(log10(1.0 + 1e+9 * smeared_deltat));
	  }
	  
	  if (v[9].size()>5) {//MD: v[9] is the last layer. Maybe the point is "iff the activity is low"
	    nomip_deltat_smeared_wrt_start->Fill(1e+9 * smeared_deltat);
	    nomip_deltat_smeared_wrt_start_longscale->Fill(1e+9 * smeared_deltat);
	    nomip_deltat_smeared_wrt_start_log->Fill(log10(1.0 + 1e+9 * smeared_deltat));
	  }

	}
      }
    }

    double deltat_slower = slower - tStart;
    double smeared_deltat_slower = tr->Gaus(deltat_slower, TIME_RESOLUTION);

    slower_deltat_wrt_start->Fill(1e+9 * deltat_slower);
    slower_deltat_wrt_start_longscale->Fill(1e+9 * deltat_slower);
    slower_deltat_wrt_start_log->Fill(log10(1.0 + 1e+9 * deltat_slower));
    
    slower_deltat_smeared_wrt_start->Fill(1e+9 * smeared_deltat_slower);
    slower_deltat_smeared_wrt_start_longscale->Fill(1e+9 * smeared_deltat_slower);
    slower_deltat_smeared_wrt_start_log->Fill(log10(1.0 + 1e+9 * smeared_deltat_slower));
    
  } //for i
  
  COUT(INFO) <<ENDL;
  
  COUT(INFO) << "Lost energies:  " <<energy_lost << " on " << iMeas << ENDL;
  COUT(INFO) << "Lost positions: " <<position_lost << " on " << iMeas << ENDL;

  COUT(INFO) <<"Writing output..." <<ENDL;
  
  outFile->Write();

  outFile->Close();

  COUT(INFO) <<"Output written in " <<outFileName <<ENDL;

#ifdef _DEBUG_
  debug::end_debug();
#endif

  return 0;
}
