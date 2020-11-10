
// Executable for time measures analysis

#include "physics.h"
#include "vector2.h"
#include "progress.h"
#include "Geometry.h"
#include "TrCluster.hh"
#include "TimeSim.h"
#include "TimeSegm.h"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TRandom3.h"


#include <iostream>


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


int main(int argc, char **argv)
{

  //input file

  if(argc < 2)
  {
    std::cerr <<"\nfatal error: input file not specified\n";
    return 1;
  }

  TFile *inFile = TFile::Open(argv[1]);

  if(inFile->IsZombie())
  {
    std::cerr <<"\nfatal error: unable to open " <<argv[1] <<std::endl;
    return 1;
  }



  //output file

  TFile *outFile = new TFile("time.root", "recreate");

  if(outFile->IsZombie())
  {
    std::cerr <<"\nfatal error: unable to open output file\n";
    return 1;
  }



  TTree *events;
	inFile->GetObject("Tree", events);
	events->Print();

  TClonesArray *branch = new TClonesArray("TrCluster", 200);
	events->SetBranchAddress("Events", &branch);



  TH1F *h_time = new TH1F
  (
    "h_time", "time measures; t_meas - t_true [s]; entries",
    10000, -4e-10, 0
  );

  TRandom3 *random = new TRandom3(9298);

  TimeSegm *time_segm = new TimeSegm(&GEO, A, 1);

  TimeSim *time_sim = new TimeSim
    (time_segm, random, thickness * 1e-3 /*convert mm to m*/ );



  std::cout <<"\nAnalysis of " <<events->GetEntries() <<" events...\n";

  int negative_lad_strip = 0; //for debug purpose

  for (int i = 0; i < events->GetEntries(); i++)
  {

    //print and update progress bar
    progress(i, events->GetEntries());

    events->GetEntry(i); //get one object (particle or other) trace


    //scan a single particle hits

    for(int j = 0; j < branch->GetEntries(); ++j)

      /****************************************
      * BUG: strip < 0 from TrCluster object *
      ****************************************/

      if //TEMPORARY FIX
      (
        ((TrCluster*) branch->At(j)) ->ladder >= 0
        && ((TrCluster*) branch->At(j)) ->strip >= 0
      )
      {
        time_segm->SetHit
        (
          ((TrCluster*) branch->At(j)) ->ladder,
          ((TrCluster*) branch->At(j)) ->strip,

          //convert ns to s
          ((TrCluster*) branch->At(j)) ->time * 1e-9,

          //convert GeV to eV
          ((TrCluster*) branch->At(j)) ->eDep * 1e+9
        );
      }

      else
        ++negative_lad_strip; //for debug purpose


    //constant fraction of current signal

    for(int k=0; k < time_segm->GetNgroups(); ++k)
    {
      //print and update progress bar
      //progress(analyzed_hits, total_hits);

      std::vector<mytime_t> true_time;
      time_segm->GetTimes(k, true_time);

      if(true_time.size() == 0)
        continue;
        //HERE COULD BE GENERATED NOISE SIGNAL WITHOUT HIT

      TGraph *current = new TGraph();
      time_sim->GetCurrentSignal(k, current);

      mytime_t meas_time = time_sim->GetTime(current, 0.1);

      delete current;

      for(int l = 0; l < (int) true_time.size(); ++l)
        h_time->Fill( meas_time - true_time[l] );

      //analyzed_hits += true_time.size();
    }

    time_segm->Clear(); //clear particle traces
  }

  delete time_sim;
  delete time_segm;

  std::cout <<std::endl <<std::endl;


  //debug

  if(negative_lad_strip > 0)
  {
    std::cerr <<"[DEBUG] negative strips or ladders:\t";
    std::cerr <<negative_lad_strip <<std::endl <<std::endl;
  }


  // write output

  outFile->WriteTObject(h_time);
  outFile->Close();

  std::cout <<"Results written in:\ttime.root\n\n";

  return 0;
}
