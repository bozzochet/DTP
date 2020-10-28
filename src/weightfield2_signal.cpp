
#include "TFile.h"
#include "TCanvas.h"
#include "TList.h"
#include "TObject.h"
#include "TGraph.h"

#include <iostream>
#include <string>
#include <cstring>

int main(int argc, char *argv[])
{
  if(argc < 1)
  {
    std::cerr <<"\nfatal error: ROOT filename must be passed as argument\n";
    return 1;
  }
  TFile *f = new TFile(argv[1]);

  if(f->IsZombie())
  {
    std::cerr <<"\nfatal error: unable to open ROOT file\n";
    return 1;
  }

  std::cout <<"\nobjects in file:\n";
  f->Map();

  TCanvas *data = (TCanvas*) f->Get("currents");

  TList *list;

  std::cout
    << "\ncurrents canvas contains "
    << ( list = data->GetListOfPrimitives() )->GetEntries()
    << " objects:\n";

  std::cout << "\tClassName\tTitle\n";

  for(int i=0; i < (int) list->GetEntries(); ++i)
    std::cout
      << i << ":\t"
      << list->At(i)->ClassName() << "\t"
      << list->At(i)->GetTitle() << std::endl;

  std::cout << std::endl;

  TGraph *current;

  std::cout << "\nidentifying total current\n";

  for(int i=0; i < (int) list->GetEntries(); ++i) {

    TGraph *obj = (TGraph*) list->At(i);

    std::cout
      << i << ":\t" << obj->ClassName() <<"\t";

    if( std::strcmp(obj->ClassName(), "TGraph") == 0)
      std::cout << obj->GetLineColor();

    std::cout << std::endl;

    if(
      std::strcmp(obj->ClassName(), "TGraph") == 0
      && obj->GetLineColor() == 3
    )
    {
      current = obj;
      std::cout << "\ntotal current is object at " << i << std::endl;
      break;
    }
  }

  TCanvas *c = new TCanvas("c","weightfield2 simulation",600,600);
  c->cd();

  current->Draw();

  c->SaveAs("total_current.png");

  f->Close();
  delete f;
}
