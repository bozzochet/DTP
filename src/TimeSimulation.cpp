
#include "TimeSimulation.h"


TimeSimulation::TimeSimulation()
{
  /* open weightfield2 output ROOT file:
   * originally in <weightfield2_folder>/sensors/graph/parameters.root
   * moved in <DTP_source>/data/weightfield2.root
   */

  exp_ = SetExp("data/weightfield2.root");

  line_ = new TF1("line", "[0]*x", 0, T_PEAK_);
  line_->SetParName(0, "Slope");
}


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

  TF1 *exp = new TF1("exp", "expo(0)", 0, 3e-5);

  current->Fit(exp);

  k_exp_ = exp->GetParameter(1);


  file->Close();
  delete file;

  return exp;
}


TimeSimulation::~TimeSimulation()
{
  delete line_;
  delete exp_;
}


void TimeSimulation::SetEnergy(const vector2<double> &vec)
{
  for(int lad = 0; lad < Nladders; ++lad)
    for(int strip = 0; strip < Nstrips; ++strip)
    {
      energy_t energy = vec[lad][strip];

      if(energy != 0)
        energy_ [absStrip(lad,strip)] = energy;
    }
}

/*
void TimeSimulation::AddSignal
  (
    TH1F *hist,
    const energy_t &energy,
    const std::vector<mytime_t> &times
  )
{

  charge_t Q = energy / ENERGY_COUPLE_
    * FOND_CHARGE_;

  // set integral equal to Q collected by strip

  double peak = 2.0 * TMath::Abs(exp_->GetParameter("Slope")) * Q
    / (T_PEAK_ * TMath::Abs(exp_->GetParameter("Slope")) + 2);

    std::cout <<"\npeak: " <<peak <<"\n";

  line_->SetParameter("Slope", peak / T_PEAK_);

  exp_->SetParameter
    ( "Constant", TMath::Log(peak) - exp_->GetParameter(1) * T_PEAK_ );

  std::cout <<"\nconstant: " << exp_->GetParameter(0) <<"\n"
    <<"k_exp: " <<exp_->GetParameter(1) <<"\n";

  line_->SetRange(T_START_, T_END_);
  exp_->SetRange(T_START_, T_END_);

  //sample signal and fill hist

  for(int i = 0; i < (int) times.size(); ++i)
    for(mytime_t t = times[i]; t < T_END_; t += T_SAMPLING_)
      if(t - times[i] < T_PEAK_)
        hist->Fill(t, line_->Eval(t - times[i]));
      else
        hist->Fill(t, exp_->Eval(t - times[i]));
}
*/


TH1F* TimeSimulation::GetSignal(const int &lad, const int &strip)
{
  std::string name = "current(:)";
  name.insert(8, std::to_string(lad));
  name.insert(name.length()-1, std::to_string(strip));

  std::string title = "ladder:  strip: ";
  title.insert(8, std::to_string(lad));
  title += std::to_string(strip);

  TH1F *hist = new TH1F(name.c_str(), title.c_str(),
    N_BINS_, T_START_, T_END_);

  hist->GetXaxis()->SetTitle("time [s]");
  hist->GetYaxis()->SetTitle("current [???]");

  //draw hist as a function instead of using bars
  hist->SetLineColor(kWhite);
  hist->SetMarkerColor(kBlue);
  hist->SetMarkerStyle(kFullDotMedium);

  //GENERATE NOISE HERE

  if(energy_.find(absStrip(lad,strip)) == energy_.end())
    return hist; //return hist with noise only

/*
  AddSignal
    (hist, energy_[absStrip(lad,strip)], times_[absStrip(lad,strip)]);
*/

  return hist;
}
