
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
#include "TMultiGraph.h"
#include "TGaxis.h"

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



void charge_meas
(
  const char *input, const char *output,
  const double &threshold
);

void get_example(const char *input, const char *output);



int main(int argc, char **argv)
{
  if(argc < 3)
  {
    std::cerr <<"\nfatal error: missing arguments\n";
    return 1;
  }


  std::clock_t start = std::clock();


  if(std::strcmp(argv[2], "example") == 0)
  {
    get_example(argv[1], "time--example.root");
  }


  else if(std::strcmp(argv[2], "charge") == 0)
  {
    for(int i=1; i <= std::stod(argv[4]) / std::stod(argv[3]); ++i)
    {
      double threshold = std::stod(argv[3]) * i;

      std::string output = "time--charge--thresh.root";

      output.insert
        (output.length()-5, std::to_string((int) (threshold*100)) );

      std::cout <<"\n\nCONSTANT FRACTION OF CHARGE WITH THRESHOLD "
        << threshold <<"\n";

      charge_meas(argv[1], output.c_str(), threshold);
    }
  }


  std::cout <<"\nexecution completed in ";

  int sec = (std::clock() - start) / CLOCKS_PER_SEC;

  if(sec / 60 >= 60)
    std::cout << sec / 3600 << "h " << sec % 3600 / 60 << "min "
      << (sec % 3600) % 60 <<"s\n";

  else if(sec >= 60)
    std::cout << sec / 60 << "min " << sec % 60 <<"s\n";

  else
    std::cout << sec <<"s\n";


  return 0;
}



void charge_meas
(
  const char *input, const char *output,
  const double &threshold
)
{

  TFile *inFile = TFile::Open(input);

  if(inFile->IsZombie())
  {
    std::cerr <<"\nfatal error: unable to open " <<input <<std::endl;
    exit(1);
  }


  //output file

  TFile *outFile = new TFile(output, "recreate");

  if(outFile->IsZombie())
  {
    std::cerr <<"\nfatal error: unable to open output file\n";
    exit(1);
  }


  TTree *events;
	inFile->GetObject("Tree", events);
	//events->Print();

  TClonesArray *branch = new TClonesArray("TrCluster", 200);
	events->SetBranchAddress("Events", &branch);


  //time measures with noise
  TH1F *h_meas = new TH1F("h_meas", " ", 1e+3, -0, 0);

  //time measures without noise
  TH1F *h_ideal = new TH1F("h_ideal", " ", 1e+3, -0, 0);

  //ideal charge collected
  TH1F *h_charge_ideal = new TH1F("h_charge_ideal", " ", 1e+3, -0, 0);

  //charge collected without noise
  TH1F *h_charge = new TH1F("h_charge", " ", 1e+3, -0, 0);

  h_ideal->SetTitle
    ("time measures without noise; t_meas - t_true [s]; entries");

  h_meas->SetTitle
    ("time measures with noise; t_meas - t_true [s]; entries");

  h_charge_ideal->SetTitle
    ("ideal charge collected;charge [C];entries");

  h_charge->SetTitle
    ("charge collected with noise;charge [C];entries");

  h_ideal->SetCanExtend(TH1::kAllAxes);
  h_meas->SetCanExtend(TH1::kAllAxes);
  h_charge_ideal->SetCanExtend(TH1::kAllAxes);
  h_charge->SetCanExtend(TH1::kAllAxes);

  h_ideal->SetLineColor(kRed);
  h_charge_ideal->SetLineColor(kRed);


  TRandom3 *random = new TRandom3(9298);

  TimeSim *time_sim = new TimeSim(NULL, random);


  std::cout <<"Analysis of " <<events->GetEntries() <<" events...\n";

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
        for(int k = 0; k < 2; ++k)
        {
          TGraph *charge = new TGraph();
          TGraph *charge_ideal = new TGraph();

          charge_t q;

          h_charge_ideal->Fill
          (
            q = time_sim->GetChargeSignal
              (cl->clust[k] * 1e+9, charge_ideal, false)
          );

          for(int m=0; m < charge_ideal->GetN(); ++m)
            charge->SetPoint
              (m, charge_ideal->GetX()[m], charge_ideal->GetY()[m]);

          h_charge->Fill(q + time_sim->AddChargeNoise(charge));


          mytime_t t_true = time_sim->GetMeas(charge_ideal, threshold);

          mytime_t t_meas = time_sim->GetMeas(charge, threshold);

          if(t_true != -9999)
            h_ideal->Fill(t_true);

          if(t_meas != -9999)
            h_meas->Fill(t_meas);


          delete charge;
          delete charge_ideal;

        }//for k
      } //if lad >= 0 strip >= 0

    } //for j
  } //for i


  delete time_sim;
  delete random;

  delete events;
  delete branch;

  std::cout <<std::endl <<std::endl;

  std::cout <<"Mean: " <<h_meas->GetMean();
  std::cout <<"\nStd Dev: " <<h_meas->GetStdDev() <<std::endl;


  // write output

  outFile->WriteTObject(h_meas);
  outFile->WriteTObject(h_ideal);

  outFile->WriteTObject(h_charge);
  outFile->WriteTObject(h_charge_ideal);

  outFile->Close();

  delete outFile;

  std::cout <<"\nResults written in " <<output;
}



void get_example(const char *input, const char *output)
{
  TFile *inFile = TFile::Open(input);

  if(inFile->IsZombie())
  {
    std::cerr <<"\nfatal error: unable to open " <<input <<std::endl;
    exit(1);
  }


  //output file

  TFile *outFile = new TFile(output, "recreate");

  if(outFile->IsZombie())
  {
    std::cerr <<"\nfatal error: unable to open output file\n";
    exit(1);
  }


  TTree *events;
  inFile->GetObject("Tree", events);
  //events->Print();

  TClonesArray *branch = new TClonesArray("TrCluster", 200);
  events->SetBranchAddress("Events", &branch);


  TRandom3 *random = new TRandom3(9298);

  TimeSim *time_sim = new TimeSim(NULL, random);

  TGraph *charge = new TGraph();
  TGraph *charge_ideal = new TGraph();
  TGraph *current = new TGraph();
  TGraph *current_ideal = new TGraph();


  events->GetEntry(0); //get one particle trace

  TrCluster *cl = (TrCluster*) branch->At(1);

      /****************************************
      * BUG: strip < 0 from TrCluster object *
      ****************************************/

  if(cl->ladder >= 0 && cl->strip >= 0) //TEMPORARY FIX
  {
    time_sim->GetChargeSignal
      (cl->clust[0] * 1e+9, charge_ideal, false);

    for(int h = 0; h < charge_ideal->GetN(); ++h)
    {
      mytime_t t;
      charge_t q;

      charge_ideal->GetPoint(h, t, q);
      charge->SetPoint(charge->GetN(), t, q);
    }

    time_sim->AddChargeNoise(charge);

    time_sim->GetCurrentSignal
      (cl->time * 1e-9, current_ideal, charge_ideal);
    time_sim->GetCurrentSignal(cl->time * 1e-9, current, charge);
  }


  delete time_sim;
  delete random;

  delete events;
  delete branch;

  std::cout <<std::endl <<std::endl;


  charge->SetLineColor(kBlue);
  charge->SetLineWidth(2);

  charge_ideal->SetLineColor(kRed);
  charge_ideal->SetLineWidth(2);

  current->SetLineColor(kBlue);
  current->SetLineWidth(2);

  current_ideal->SetLineColor(kRed);
  current_ideal->SetLineWidth(2);

  TMultiGraph *mg_current = new TMultiGraph();
  TMultiGraph *mg_charge = new TMultiGraph();

  mg_current->Add(current);
  mg_current->Add(current_ideal);

  mg_charge->Add(charge);
  mg_charge->Add(charge_ideal);


  // write output

  outFile->WriteTObject(mg_current);
  outFile->WriteTObject(mg_charge);

  outFile->Close();

  delete outFile;

  std::cout <<"examples written in:\t" <<output <<"\n";
}
