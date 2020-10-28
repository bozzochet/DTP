
#include "TimeSimulation.h"


TimeSimulation::TimeSimulation()
{
  /* open weightfield2 output ROOT file:
   * originally in <weightfield2_folder>/sensors/graph/parameters.root
   * moved in <DTP_source>/data/weightfield2.root
   */

  //user must execute executables files in build or installation folder
  TFile *f = new TFile("data/weightfield2.root");

  //check if file is correctly open
  if(f->IsZombie())
  {
    std::cerr
      <<"\n[TIME SIMULATION] fatal error: unable to open ROOT file\n";
    return;
  }

  //current signals are contained in a TCanvas object called "currents"
  TList *list =
    ( (TCanvas*) f->Get("currents") )->GetListOfPrimitives();


  //get total current TGraph
  TGraph *current;

  for(int i=0; i < (int) list->GetEntries(); ++i)
  {
    TGraph *obj = (TGraph*) list->At(i);

    if(
      std::strcmp(obj->ClassName(), "TGraph") == 0 //obj is a TGraph
      && obj->GetLineColor() == 3 //TGraph line color is green (total current is green)
    )
    {
      current = obj;
      break;
    }
  }

  x_ = new TVectorD();
  y_ = new TVectorD();

  //copy current TGraph points in TVectors x_ , y_
  x_->Use(current->GetN(), current->GetX());
  y_->Use(current->GetN(), current->GetY());

  f->Close();
  delete f;
}


TimeSimulation::~TimeSimulation()
{
  delete x_;
  delete y_;
}


/*
TGraph* TimeSimulation::GetSignal(const int &lad, const int &s)
{
  TGraph *gr = new TGraph();

  //hit times for strip passed
  std::vector<mytime_t> vec = Stopwatch::GetTimes(lad,s);

  //WORK IN PROGRESS
}
*/
