
#include "TimeSimulation.h"


TimeSimulation::TimeSimulation()
{
  charge_ = new TF1(
    "charge",
    "[0] * (1 - TMath::Exp(-[1]*x) )",
    T_START_, T_END_
  );

  charge_->SetParName(0, "Q");
  charge_->SetParName(1, "k");
  charge_->SetParameter("k", 1.0 / T_CAPACITOR_);

  up_ = new TF1("up", "[0]*x", T_START_, T_END_);

  up_->SetParName(0, "Slew rate");
  up_->SetParameter("Slew rate", SLEW_RATE_);
}

/*
TF1* TimeSimulation::SetExp(const char* filename)
{
  TFile *file = new TFile(filename);

  //check if file is correctly open
  if(file->IsZombie())
  {
    std::cerr
      <<"\n[TIME SIMULATION] fatal error: unable to open ROOT file\n";
    return NULL;
  }

  //current signals are contained in a TCanvas object called "currents"
  TList *list =
    ( (TCanvas*) file->Get("currents") )->GetListOfPrimitives();


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


  //get exp parameters

  TF1 *exp = new TF1("exp", "expo(0)", 0, 3e-8);  //weightfield2 simulate a signal from 0 to 30 ns = 3e-8 s

  current->Fit(exp);


  file->Close();
  delete file;

  return exp;
}
*/

/*
void TimeSimulation::SetEnergy(const vector2<double> &vec)
{
  for(int lad = 0; lad < Nladders; ++lad)
    for(int strip = 0; strip < Nstrips; ++strip)
    {
      energy_t energy = vec[lad][strip] * 1e-9; //convert GeV to eV
      //GGS energies are expressed in GeV

      if(energy != 0)
        energies_ [absStrip(lad,strip)] = energy;
    }
}
*/

void TimeSimulation::AddSignal
  (TH1F *hist, const energy_t &hitEnergy, const mytime_t &hitTime)
{
  charge_t Q = hitEnergy / ENERGY_COUPLE_
    * FOND_CHARGE_;

  charge_->SetParameter("Q", Q);

  TGraph *down = new TGraph(charge_, "d"); //derivative of charge_

  // fill hist with up_

  //t_peak is delta_t between hit and peak of current
  mytime_t t_peak = down->Eval(0) / up_->GetParameter("Slew rate");

  for( mytime_t t = 0; t < t_peak; t += T_SAMPLING_ )
    hist->Fill(t + hitTime, up_->Eval(t));

  // fill with down_

  for(int i = 0; i < (int) down->GetN(); ++i)
  {
    mytime_t t;
    current_t I;

    down->GetPoint(i, t, I);
    hist->Fill(t + hitTime + t_peak, I);
  }

}



std::vector<TH1F*>* TimeSimulation::GetSignal
  (const int &lad, const int &strip)
{
  std::vector<TH1F*> *vec = new std::vector<TH1F*>;

  std::string name = "current(:)";
  name.insert(8, std::to_string(lad));
  name.insert(name.length()-1, std::to_string(strip));

  std::string title = "current signal [ladder:  strip: ]";
  title.insert(24, std::to_string(lad));
  title.insert(title.length()-1, std::to_string(strip));

  /* name and title variables are used again near the end of
   * this method for another hist; they will be replaced for the needs
   */

  TH1F *hist = new TH1F
    ( name.c_str(), title.c_str(), N_BINS_, T_START_, T_END_ );

  //hist is stored in vec[0]
  vec->push_back(hist);

  hist->GetXaxis()->SetTitle("run time [s]");
  hist->GetYaxis()->SetTitle("current [A]");

  //draw hist with points instead of using bars
  hist->SetLineColor(kWhite);
  hist->SetMarkerColor(kBlue);
  hist->SetMarkerStyle(kFullDotMedium);

  //GENERATE NOISE HERE

  if(energies_.find(absStrip(lad,strip)) == energies_.end())
    return vec; //return hist with noise only

  // fill hist with signal and push_back charge collection graphs

  for(int i = 0; i < (int) energies_[absStrip(lad,strip)].size(); ++i)
  {
    AddSignal
    (
      hist,
      energies_[absStrip(lad,strip)][i],
      times_[absStrip(lad,strip)][i]
    );

    //hist for charge collected

    name = "charge(:)";
    name.insert(name.length()-3, std::to_string(i));
    name.insert(name.length()-2, std::to_string(lad));
    name.insert(name.length()-1, std::to_string(strip));

    title = "charge collected [hit:  ladder:  strip: ]";
    title.insert(title.length()-18, std::to_string(i));
    title.insert(title.length()-9, std::to_string(lad));
    title.insert(title.length()-1, std::to_string(strip));

    TH1F *charge = new TH1F
      (name.c_str(), title.c_str(), N_BINS_, T_START_, T_END_);

    charge->GetXaxis()->SetTitle("time from hit [s]");
    charge->GetYaxis()->SetTitle("charge collected [C]");

    //draw charge with points instead of using bars
    charge->SetLineColor(kWhite);
    charge->SetMarkerColor(kBlue);
    charge->SetMarkerStyle(kFullDotMedium);

    for( mytime_t t = 0; t < T_END_; t += T_SAMPLING_ )
      charge->Fill(t, charge_->Eval(t));

    vec->push_back(charge);
  }

  return vec;
}


TimeSimulation::~TimeSimulation()
{
  delete up_;
  delete charge_;
}
