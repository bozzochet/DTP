
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
#include <ctime>



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
  (const char*, const char*, const time_segm&, const int&, const bool);
void analyze_noise(TTree*, TClonesArray*, TimeSim*, TH1F*);
void analyze_segm(TTree*, TClonesArray*, TimeSegm*, TimeSim*, TH1F*);



int main(int argc, char **argv)
{
  if(argc < 2)
  {
    std::cerr <<"\nfatal error: input file not specified\n";
    return 1;
  }

  time_segm S[2] = {A,B};
  int J[2] = {10,10};

  const char* output[2] = {"time--A10.root", "time--B10.root"};

  for(int i = -1; i < 0; ++i)
  {

    std::clock_t start = std::clock();


    if(i == -1)
    {
      std::cout <<"\nANALYSIS WITHOUT SEGM:\n\n";
      analyze(argv[1], "time--nosegm.root", A, 1, false);
    }

    else
    {
      std::cout <<"\nANALYSIS WITH SEGM "
        <<(char) S[i] <<std::to_string(J[i]) <<":\n\n";

      analyze(argv[1], output[i], S[i], J[i], true);
    }


    int sec = (std::clock() - start) / CLOCKS_PER_SEC;

    if(sec < 1) continue;

    std::cout <<"task completed in ";

    if(sec / 60 >= 60)
      std::cout << sec / 3600 << ":" << sec % 3600 / 60 << ":"
        << (sec % 3600) % 60 ;

    else if(sec >= 60)
      std::cout << sec / 60 << ":" << sec % 60 <<" min\n";

    else
      std::cout << sec <<" secondi\n";

  }

  return 0;
}



int analyze
(
  const char *input, const char *output,
  const time_segm &S, const int &J, const bool segm
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
  TH1F *h_segm = new TH1F
  (
    "h_segm",
    "time measures with noise and segm; t_meas - t_true [s]; entries",
    10000, -1e-10, 1e-10
  );

  h_segm->SetBit(TH1::kAutoBinPTwo);


  //time measures with noise only
  TH1F *h_noise = new TH1F
  (
    "h_noise",
    "time measures with noise; t_meas - t_true [s]; entries",
    10000, -3e-12, 3e-12
  );

  h_noise->SetBit(TH1::kAutoBinPTwo);



  TRandom3 *random = new TRandom3(9298);

  TimeSegm *time_segm = new TimeSegm(&GEO, S, J);

  TimeSim *time_sim = new TimeSim(time_segm, random);



  if( !segm)
    analyze_noise(events, branch, time_sim, h_noise);
  else
  {
    analyze_segm(events, branch, time_segm, time_sim, h_segm);
    time_segm->Clear(); //clear particle traces
  }


  delete time_sim;
  delete time_segm;
  delete random;

  delete events;
  delete branch;

  std::cout <<std::endl <<std::endl;


  // write output

  if(segm)
    outFile->WriteTObject(h_segm);
  else
    outFile->WriteTObject(h_noise);

  outFile->Close();

  delete outFile;

  //delete h_segm;  //cause program crash
  //delete h_noise;  //cause program crash

  std::cout <<"Results written in:\t" <<output <<"\n\n";

  return 0;
}


//constant fraction with noise
void analyze_noise
(
  TTree *events, TClonesArray *branch,
  TimeSim *time_sim, TH1F *h_noise
)
{

  std::cout <<"\nAnalysis of " <<events->GetEntries() <<" events...\n";


  for (int i = 0; i < events->GetEntries(); i++)
  {

    //print and update progress bar
    progress(i, events->GetEntries());

    events->GetEntry(i); //get one particle trace

    for(int j = 0; j < branch->GetEntries(); ++j)
    {
      TrCluster *cl = (TrCluster*) branch->At(j);

      /****************************************
      * BUG: strip < 0 from TrCluster object *
      ****************************************/

      if(cl->ladder >= 0 && cl->strip >= 0) //TEMPORARY FIX
      {
        TGraph *current = new TGraph();
        TGraph *charge = new TGraph();

        time_sim->GetChargeSignal(cl->eDep * 1e+9, charge);

        time_sim->GetCurrentSignal(cl->time * 1e-9, current, charge);

        h_noise->Fill
          (time_sim->GetTime(current, 0.1) - cl->time * 1e-9);

        delete current;
        delete charge;
      }
    }
  }
}


//constant fraction with segm and noise
void analyze_segm
(
  TTree *events, TClonesArray *branch,
  TimeSegm *time_segm, TimeSim *time_sim, TH1F *h_segm
)
{

  std::cout <<"\nAnalysis of " <<events->GetEntries() <<" events...\n";


  for (int i = 0; i < events->GetEntries(); i++)
  {

    //print and update progress bar
    progress(i, events->GetEntries());

    events->GetEntry(i); //get one particle trace

    for(int j = 0; j < branch->GetEntries(); ++j)
    {
      TrCluster *cl = (TrCluster*) branch->At(j);

      if(cl->ladder >= 0 && cl->strip >= 0) //TEMPORARY FIX
        time_segm->SetHit
        (
          cl ->ladder, cl ->strip,
          cl ->time * 1e-9, //convert ns to s
          cl ->eDep * 1e+9 //convert GeV to eV
        );
    }


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

      for(int m = 0; m < (int) true_time.size(); ++m)
        h_segm->Fill( meas_time - true_time[m] );
    }
  }
}
