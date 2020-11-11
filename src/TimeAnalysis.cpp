
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


int analyze
(
  const char *input, const char *output,
  const time_segm &S, const int &J
);


int main(int argc, char **argv)
{
  if(argc < 2)
  {
    std::cerr <<"\nfatal error: input file not specified\n";
    return 1;
  }

  time_segm S[4] = {A, A, B, B};
  int J[4] = {1, 10, 1, 10};
  const char* output[4] = {"time--A1.root", "time--A10.root", "time--B1.root", "time--B10.root"};

  for(int i=0; i < 4; ++i)
  {
    std::cout <<std::endl <<(char) S[i] <<std::to_string(J[i])
      <<" ANALYSIS:\n\n";

    if(analyze(argv[1], output[i], S[i], J[i]) != 0)
      return 1;
  }

  return 0;
}

int analyze
(
  const char *input, const char *output,
  const time_segm &S, const int &J
)
{

  TFile *inFile = TFile::Open(input);

  if(inFile->IsZombie())
  {
    std::cerr <<"\nfatal error: unable to open " <<input <<std::endl;
    return 1;
  }



  //output file

  TFile *outFile = new TFile(output, "recreate");

  if(outFile->IsZombie())
  {
    std::cerr <<"\nfatal error: unable to open output file\n";
    return 1;
  }



  TTree *events;
	inFile->GetObject("Tree", events);
	//events->Print();

  TClonesArray *branch = new TClonesArray("TrCluster", 200);
	events->SetBranchAddress("Events", &branch);



  //time measures with segm and noise
  TH1F *h_time = new TH1F
  (
    "h_time", "time measures; t_meas - t_true [s]; entries",
    10000, -1e-10, 1e-10
  );

  //time measures without segmentation and without noise
  TH1F *h_ideal = new TH1F
  (
    "h_ideal", "time measures; t_meas - t_true [s]; entries",
    10000, -1e-10, 1e-10
  );

  //confront signal t0 with t_hit
  TH1F *h_signal_zero = new TH1F
  (
    "h_signal_zero",
    "hit time for signal; t_signal - t_true [s]; entries",
    10000, -1e-10, 1e-10
  );



  TRandom3 *random = new TRandom3(9298);

  TimeSegm *time_segm = new TimeSegm(&GEO, S, J);

  TimeSim *time_sim = new TimeSim(time_segm, random);



  std::cout <<"\nAnalysis of " <<events->GetEntries() <<" events...\n";

  int negative_lad_strip = 0; //for debug purpose

  for (int i = 0; i < events->GetEntries(); i++)
  {

    //print and update progress bar
    progress(i, events->GetEntries());

    events->GetEntry(i); //get one object (particle or other) trace


    //scan a single particle hits

    for(int j = 0; j < branch->GetEntries(); ++j)
    {
      TrCluster *cl = (TrCluster*) branch->At(j);

      /****************************************
      * BUG: strip < 0 from TrCluster object *
      ****************************************/

      if(cl->ladder >= 0 && cl->strip >= 0) //TEMPORARY FIX
      {
        time_segm->SetHit
        (
          cl ->ladder, cl ->strip,
          cl ->time * 1e-9, //convert ns to s
          cl ->eDep * 1e+9 //convert GeV to eV
        );

        TGraph *current = new TGraph();
        TGraph *charge = new TGraph();

        time_sim->GetChargeSignal(cl->eDep * 1e+9, charge, false);

        time_sim->GetCurrentSignal(cl->time * 1e-9, current, charge);

        h_ideal->Fill
          (time_sim->GetTime(current, 0.1) - cl->time * 1e-9);

        mytime_t t0;
        current_t tmp;

        current->GetPoint(0, t0, tmp);

        h_signal_zero->Fill(t0 - cl->time * 1e-9);

        delete current;
        delete charge;
      }

      else
        ++negative_lad_strip; //for debug purpose

    } //scan single particle hits


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
  outFile->WriteTObject(h_ideal);
  outFile->WriteTObject(h_signal_zero);

  outFile->Close();

  std::cout <<"Results written in:\t" <<output <<"\n\n";

  return 0;
}
