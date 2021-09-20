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
#include "TH2F.h"
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

#define SIGNAL_THRESHOLD 10000 //eV 
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

  TH1F *primaries = new TH1F ("primaries", "", 40, 0, 6.5);
  primaries->GetXaxis()->SetTitle("Z (cm)");
  primaries->GetYaxis()->SetTitle("Number of hits");
   
  TH1D *electron_hit = new TH1D("electron_hit", "", 40, 0, 40);
  electron_hit->GetXaxis()->SetTitle("Layer");
  electron_hit->GetYaxis()->SetTitle("Number of hits");

  TH1D *positron_hit = new TH1D("positron_hit", "", 40, 0, 40);
  electron_hit->GetXaxis()->SetTitle("Layer");;
  electron_hit->GetYaxis()->SetTitle("Number of hits");


  TH2F *energy_distribution = new TH2F ("energy_distribution", "", 40, 0, 40, 10000, 0, 1000000);
  energy_distribution->GetXaxis()->SetTitle("Layer");
  energy_distribution->GetYaxis()->SetTitle("Energy (eV)");

  TH2F *energy_distribution_electrons = new TH2F ("energy_distribution_electrons", "", 40, 0, 40, 10000, 0, 1000000);
  energy_distribution_electrons->GetXaxis()->SetTitle("Layer");
  energy_distribution_electrons->GetYaxis()->SetTitle("Energy (eV)");


  TH2F *energy_distribution_positrons = new TH2F ("energy_distribution_positrons", "", 40, 0, 40, 10000, 0, 1000000);
  energy_distribution_positrons->GetXaxis()->SetTitle("Layer");
  energy_distribution_positrons->GetYaxis()->SetTitle("Energy (eV)");

  TH1I *n_products = new TH1I("n_products", "number of products", 50, 0, 50);


  
  //MC


  //measures


  //resolutions

  //calo

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
    //std::cout << "ITERATION #: " << i << std::endl;
    
    vector2<TrCluster> v;
    //v.resize(geo.Nlayers);
    v.resize(40);//there's a bug and the values of geo are not retrieved correctly

    events_tree->GetEntry(i);
  
  
    if (_calo_flag) {
      calo_tree->GetEntry(i);
      
      //h_energy_calo->Fill(Ecalo * 1e-9);
 
      if (Ecalo < Ecalo_min && Ecalo > Ecalo_max) continue;
    }
    
    


    /*
      if (a->GetEntries()>10) { // Che qui ci voglia 40?
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
      //std::cout<<"Entry # "<< j <<std::endl;
      TrCluster *cl = (TrCluster *)a->At(j);
      if(cl->firstInteraction == 1 && cl -> parPdg == 22 && cl->primIntPoint[2] > 0 && cl->primIntPoint[2]<6.5) { // Primaries
        primaries -> Fill(cl->primIntPoint[2]);
      }

      if(cl->firstInteraction == 1 && cl -> parPdg == 22) {
        if(cl->numberOfProducts > 2) {
          n_products->Fill(cl->numberOfProducts);
          std::cout << cl-> numberOfProducts << " ";
        }
      }

      /*
      if(cl->firstInteraction == 1 && cl->parPdg == 22 && cl->primIntPoint[2] > 0) {         
        double distance = 0;
        for (int q = 0; q < 40; q++) { // q  < geo.Nlayers
          if(distance < cl->primIntPoint[2] && cl->primIntPoint[2] < distance + 0.015) { //< distance + geo.thickness 
            
            break;
          }
          distance += 0.015 + 0.1512820513;
        }
      }
      */

      v[cl->layer].push_back(*cl);

      //h_time_MC->Fill(TMath::Log10(1e+9 * cl->time));

      v_times.push_back(cl->time);

       //while Events branch work with two indexes (i,j),
       //* Measures branch was filled with one index and contains
       //* a struct instead of an array. Because of this j element of
       //* entry i in Events branch, corresponds to measure on entry
       //* iMeas of Measures branch. Look in Digitization.cpp,
       //* digitization function 


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
	        //h_energy_MC->Fill(ene_true * 1e-6);
	      }
	

	
	//analyze valid measures


	      if (m==1) {
	        if(ene_meas > 0) {
	          //h_energy_meas->Fill(ene_meas * 1e-6);
	          //h_energy_res->Fill(1e-3 * (ene_meas - ene_true));
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
	        //h_time_meas15->Fill(TMath::Log10(1e+9 * meas.time[m]));
	        //h_time_res15->Fill(1e+9 * (meas.time[m] - cl->time));
	  
	        v_times_meas.push_back(meas.time[m]);
	      }
      } //for m
      

      //read only valid position measures without lost ones
      


      if (TMath::Abs(meas.position) < 1) {
	// this "if" is a temporary fix for Digitization bug:
	// Digitization.cpp,  line 424
        //h_pos_res->Fill(1e+2 * (meas.position - cl->pos[cl->xy]));
      }
      else {
        ++position_lost;
      }
      
    } //for j

    //fill slow hit
    //mytime_t slower = TMath::MaxElement(v_times.size(), &v_times[0]);
    //mytime_t slower_meas = TMath::MaxElement(v_times_meas.size(), &v_times_meas[0]); 

    //h_time_MC_slow->Fill(TMath::Log10(1e+9 * slower));
    //h_time_meas15_slow->Fill(TMath::Log10(1e+9 * slower_meas));
    
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
      //std::cout << il << std::endl;
      for (int hit = 0; hit<(int)(v[il].size()); hit++) { //hit
	      if (v[il][hit].eDep > SIGNAL_THRESHOLD) {
          energy_distribution -> Fill(v[il][hit].layer, v[il][hit].eDep);
          if (v[il][hit].parPdg == 11 && v[il][hit].parID == 1) { // Electrons
            electron_hit -> Fill(v[il][hit].layer);
            energy_distribution_electrons -> Fill(v[il][hit].layer, v[il][hit].eDep);
          }
          if (v[il][hit].parPdg == -11 && v[il][hit].parID == 1) { // Positrons
            positron_hit -> Fill(v[il][hit].layer);
            energy_distribution_positrons -> Fill(v[il][hit].layer, v[il][hit].eDep);
          }

	      }
      }
    }
  } 

  COUT(INFO) <<ENDL;

  COUT(INFO) <<"Writing output..." <<ENDL;




  
  outFile->Write();


  outFile->Close();

  COUT(INFO) <<"Output written in " <<outFileName <<ENDL;

#ifdef _DEBUG_
  debug::end_debug();
#endif

  return 0;
}
