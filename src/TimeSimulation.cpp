
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

  X_ = new TVectorD();
  Y_ = new TVectorD();

  //copy current TGraph points in TVectors X_ , Y_
  X_->Use(current->GetN(), current->GetX());
  Y_->Use(current->GetN(), current->GetY());

  f->Close();
  delete f;
}


TimeSimulation::~TimeSimulation()
{
  delete X_;
  delete Y_;
}
