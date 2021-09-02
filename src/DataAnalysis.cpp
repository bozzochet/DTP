// Ispirazione da POX, da un punto di vista tecnico tocca consultare però DTP
// Tocca pescare le informazioni dalla digitizzazione, la digitizzazione potrebbe non avere tutte le informazioni necessarie
// per la DataAnalysis
// Cosa produrre:
//
// - Vedere quanti fotoni hanno convertito nei silici, specificando il layer, e nel calorimetro, specificando il cubetto
//
// - Contare, dopo la conversione del fotone, quante volte elettrone e positrone hanno hit nei piani di silicio successivo
//
// - Contare quante volte gli elettroni e i positroni vanno a finire nel calorimetro, e in quale cubo, e quante volte escono
//
// - Contare quanti degli elettroni e positroni che vanno a finire nel calorimetro lo attraversano senza rilasciare energia (?)
//
// - Quanta energia viene rilasciata nel calorimetro da positrone ed elettrone?
//
// - Distinzione del fondo -> Protone che simula le caratteristiche del fotone
//
// - Quante volte l'evento converte nella posizione giusta (che poi definiremo)
// in qualche modo questa rappresenta l'accettanza del rivelatore (calcoli sul taccuino)
// PLOT FINALE: quanti fotoni potrebbe vedere SLA se sta in orbita un certo periodo di tempo guardando Crab Nebula a diverse energie
//
// - (se ci si riesce) Valutazione di elettrone e positrone come singolo oggetto: quando escono troppo vicini non siamo
// in grado di distinguere le particelle -> Quante volte sono talmente vicini da non riuscire a distinguerli?
//
// 
//
// Informazioni ricostruite non servono (ie: posizione ricostruita con la digitizzazione)
//
// vedere POX, fa le stesse cose

// 1) Plot di: Numero di conversioni, fotoni che convertono e in quale silicio che stanno in una certa z
// Istogramma con sulle x il numero di piani e sui bin il numero delle conversioni
// 2) Plot di: Quante hit nei silici hanno generato elettrone e positrone e quante fanno più di 10keV
// 3) Per elettrone e positrone distribuzione dell'energia di deposito
// 4) Quante hit fanno elettrone e positrone (Praticamente identico al numero di hit del plot precedente)
// 5) Istogramma delle energia depositata separata tra elettrone e positrone



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
#include "THStack.h"

#include "utils/GGSSmartLog.h"
#include "montecarlo/readers/GGSTRootReader.h"

#include <iostream>
#include <vector>

using namespace std;

#define SIGNAL_THRESHOLD 0 //8.0*300.0*ENERGY_COUPLE
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
  
  //histos
  // 1) Plot di: Numero di conversioni, fotoni che convertono e in quale silicio che stanno in una certa z
  // Istogramma con sulle x il numero di piani e sui bin il numero delle conversioni
  // 2) Plot di: Quante hit nei silici hanno generato elettrone e positrone e quante fanno più di 10keV
  // 3) Per elettrone e positrone distribuzione dell'energia di deposito
  // 4) Quante hit fanno elettrone e positrone (Praticamente identico al numero di hit del plot precedente)
  // 5) Istogramma delle energia depositata separata tra elettrone e positrone

  TH1D *photon_primaries = new TH1D ("photon_primaries", "Photon primaries", 40, 0, 40);
  TH1D *electron_positron_generation = new TH1D ("electron_positron_generation", "Photon that generated electrons and positrons", 40, 0, 40);
  TH1D *energy_distribution = new TH1D ("energy_distribution", "Energy distribution", 40, 0, 40);
  TH1D *electron_hit = new TH1D("electron_hit", "Electron hit", 40, 0, 40);
  TH1D *positron_hit = new TH1D("positron_hit", "Positron hit", 40, 0, 40);
  TH1D *energy_distribution_electrons = new TH1D("energy_distribution_electrons", "Energy distributions of electrons", 40, 0, 40);
  TH1D *energy_distribution_positrons = new TH1D("energy_distribution_positrons", "Energy distributions of positrons", 40, 0, 40);

  
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


  //-------------------------------------------------------------------------------------------------------
  for (int i = 0; i < events_tree->GetEntries(); i++) {
    //std::cout << "ITERATION #: " << i << std::endl;
    
    vector2<TrCluster> v;
    //v.resize(geo.Nlayers);
    v.resize(40);//there's a bug and the values of geo are not retrieved correctly
    //QUI IL VALORE ERA 10, MA TOCCA METTERCI I LAYER DISPONIBILI CHE SU SLA SONO 40

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
    for (int il = 0; il<(int)(v.size()); il++) { //layer, la dimensione del vector va bene
      for (int hit = 0; hit<(int)(v[il].size()); hit++) { //hit
	      if (v[il][hit].clust[0]>SIGNAL_THRESHOLD || v[il][hit].clust[1]>SIGNAL_THRESHOLD) {
          energy_distribution -> Fill(v[il][hit].layer, v[il][hit].eDep);
          if (v[il][hit].parPdg == 22 && v[il][hit].firstInteraction == 1 && v[il][hit].primIntPoint[2] != -99999) {
            photon_primaries -> Fill(v[il][hit].layer);
          }

          if (v[il][hit].parPdg == 11) {
            electron_hit -> Fill(v[il][hit].layer);
            energy_distribution_electrons -> Fill(v[il][hit].layer, v[il][hit].eDep);
          }
          if (v[il][hit].parPdg == -11) {
            positron_hit -> Fill(v[il][hit].layer);
            energy_distribution_positrons -> Fill(v[il][hit].layer, v[il][hit].eDep);
          }

	      }
      }
    }
  } 
  //-------------------------------------------------------------------------------------

  COUT(INFO) <<ENDL;
  
  //COUT(INFO) << "Lost energies:  " <<energy_lost << " on " << iMeas << ENDL;
  //COUT(INFO) << "Lost positions: " <<position_lost << " on " << iMeas << ENDL;

  COUT(INFO) <<"Writing output..." <<ENDL;




  
  outFile->Write();


  outFile->Close();

  COUT(INFO) <<"Output written in " <<outFileName <<ENDL;

#ifdef _DEBUG_
  debug::end_debug();
#endif

  return 0;
}
