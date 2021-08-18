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

  
  static const string routineName("DataAnalysis::main");
  GGSSmartLog::verboseLevel = GGSSmartLog::INFO; // Print only INFO messages or more important

  // ExeName contiene il nome dell'eseguibile includendo il path delle informazioni
  // in parole povere credo sia questo che l'output della digitization sia letto da queste parti

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
  
  //Definito il nome in uscita del file
  TString outFileName = argv[1];
  

  //TChain rappresenta una collezione di file contenenti oggetti TTree
  //Gli TTree servono a organizzare grosse collezioni di oggetti 

  // Per ora ne togliamo un paio per vedere se funziona
  TChain *events_tree = new TChain("events");
  TChain *meas_tree = new TChain("measures");
  TChain *calo_tree = new TChain("calorimeter");
  TChain *geo_tree = new TChain("geometry");

  // shift viene utilizzata come variabile per prendere i nomi di file di input e assegnare i nomi ai file di output
  int shift=2;
  for (shift=2; shift<argc; shift++) {
    TString inputFileName = argv[shift];
    //    printf("%d) %s\n", shift, argv[shift]);

    if (inputFileName.IsDigit()) {
      break;//input files are over
    }
    
    COUT(INFO) <<"Opening TTree objects in " <<inputFileName <<"..." <<ENDL; // <-------
    
    // QUI VENGONO RACCOLTI I DATI DAL FILE DI INPUT
    events_tree->Add(argv[shift]);
    meas_tree->Add(argv[shift]);
    if (_calo_flag) {
      calo_tree->Add(argv[shift]);
    }
    geo_tree->Add(argv[shift]);
  }

  COUT(INFO) <<"Recreating output file " <<outFileName <<"..." <<ENDL; // <--------
  TFile *outFile = new TFile(outFileName, "recreate");
  outFile->cd();

  COUT(INFO) <<"Creating histos..." <<ENDL;

// 1) Plot di: Numero di conversioni, fotoni che convertono e in quale silicio che stanno in una certa z
// Istogramma con sulle x il numero di piani e sui bin il numero delle conversioni
// 2) Plot di: Quante hit nei silici hanno generato elettrone e positrone e quante fanno più di 10keV
// 3) Per elettrone e positrone distribuzione dell'energia di deposito
// 4) Quante hit fanno elettrone e positrone (Praticamente identico al numero di hit del plot precedente)
// 5) Istogramma delle energia depositata separata tra elettrone e positrone

  TH1D *photon_conversion = new TH1D("photon_conversion", "Number of converted photons", 40, 0, 40);
  TH1D *electron_hit = new TH1D("electron_hit", "Number of electron hit", 40, 0, 40);
  TH1D *positron_hit = new TH1D("positron_hit", "Number of positron hit", 40, 0, 40);
  TH1D *energy_distribution = new TH1D("energy_distribution", "Energy distribution ", 40, 0, 40);
  TH1D *energy_distribution_electron = new TH1D("energy_distribution_electron", "Energy distribution of electrons", 40, 0, 40);
  TH1D *energy_distribution_positron = new TH1D("energy_distribution_positron", "Energy distribution of positrons", 40, 0, 40);

  events_tree->Print();

  // SetBranchAddress -> Change branch address, dealing with clone trees properly
  TClonesArray *a = new TClonesArray("TrCluster", 200); 
  events_tree->SetBranchAddress("Events", &a);

  // Questi dovrebbero essere oggetti di GGS/Geant4
  Geometry geo;
  Geometry* p_geo = &geo;
  geo_tree->SetBranchAddress("Geometry", &p_geo);
  /* in Apr 2021 the object was saved wrongly in Digitazion: skipping for now. To be restored after the MDPI paper
  geo_tree->GetEntry(0);
  */

  geo.ComputeDerived();
  geo.Dump();  

  COUT(INFO) <<ENDL;
  COUT(INFO) <<"Begin loop over " <<events_tree->GetEntries() <<ENDL;

  for (int i = 0; i < events_tree->GetEntries(); i++) {
    
    vector2<TrCluster> v;
    //    v.resize(geo.Nlayers);
    v.resize(40);//there's a bug and the values of geo are not retrieved correctly
    
    for (int j=0; j < a->GetEntries(); j++) {
      TrCluster *cl = (TrCluster *)a->At(j); // Così facendo dovrebbe contare tutte le hit dei fotoni e non devono aver
      if (cl -> parID == 0) {                  // necessariamente convertito ma intanto proviamo
        photon_conversion -> Fill(cl -> layer);
      }
      if(cl -> parPdg == 11) {
        electron_hit -> Fill(cl -> layer);
        energy_distribution_electron -> Fill(cl -> layer, cl -> eDep);
      }
      if(cl -> parPdg == -11) {
        positron_hit -> Fill(cl -> layer); 
        energy_distribution_positron -> Fill(cl -> layer, cl -> eDep);
      }
      energy_distribution -> Fill(cl -> layer, cl -> eDep);
    }
  }
outFile->Write();
outFile->Close();
COUT(INFO) << "Output written in " <<outFileName << ENDL;

  return 0;
}