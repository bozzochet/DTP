
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

  random_ = new TRandom3(9298);

  deviation_ = new TH1F("charge_dev", "charge deviation", 1000, 0, 1);

  deviation_->GetXaxis()->SetTitle("relative deviation");
  deviation_->GetYaxis()->SetTitle("Entries");
}


TimeSimulation::~TimeSimulation()
{
  delete up_;
  delete charge_;
  delete deviation_;
  delete random_;
}


TGraph* TimeSimulation::GetChargeSignal(const energy_t &energy)
{
  charge_t Q = energy / ENERGY_COUPLE_ * FOND_CHARGE_;

  std::cout <<"\nQ = " << Q <<" C\n";

  charge_->SetParameter("Q", Q);

//get ideal charge signal

  TGraph* real_charge = new TGraph(charge_);

//add noise

  for(int i=0; i < real_charge->GetN(); ++i)
  {
    mytime_t t;
    charge_t q;

    real_charge->GetPoint(i,t,q);

    charge_t q0 = q;
    q += (q / 10.0) * random_->Uniform();

    deviation_->Fill( (q - q0) / q0);

    real_charge->SetPoint(i,t,q);
  }

  real_charge->Sort();

  return real_charge;
}


void TimeSimulation::AddSignal
  (TGraph *gr, const TGraph *charge, const mytime_t &hitTime)
{

  std::cout <<"\nhit time = " << hitTime <<" s\n";

// fill with charge derivative

  current_t peak;
  time_t t_peak; //delta_t between hitTime and current = peak

  for(int i=0; i < charge->GetN()-1; ++i)
  {
    mytime_t t1, t2;
    charge_t q1, q2;

    //charge is sorted
    charge->GetPoint(i, t1, q1);
    charge->GetPoint(i+1, t2, q2);

    if(i==0)
    {
      peak = (q2 - q1) / (t2 - t1);
      t_peak = peak / up_->GetParameter("Slew rate");
    }

    gr->SetPoint(
      gr->GetN(),
      hitTime + t_peak + (t2 + t1)*0.5,
      (q2 - q1) / (t2 - t1)
    );

    /* SetPoint could be replaced by AddPoint but it's beeing used
     * 6.18 ROOT version, where AddPoint was not defined yet */
  }

// fill gr with up_

  for( mytime_t t = 0; t < t_peak; t += T_SAMPLING_ )
    gr->SetPoint(gr->GetN(), t + hitTime, up_->Eval(t));

    /* SetPoint could be replaced by AddPoint but it's beeing used
     * 6.18 ROOT version, where AddPoint was not defined yet */

  gr->Sort();
}


void TimeSimulation::GetSignal
  (std::vector<TGraph*> &vec, const int &lad, const int &strip)
{
  std::string name = "current(:)";
  name.insert(8, std::to_string(lad));
  name.insert(name.length()-1, std::to_string(strip));

  std::string title = "current signal [ladder:  strip: ]";
  title.insert(24, std::to_string(lad));
  title.insert(title.length()-1, std::to_string(strip));

  /* name and title variables are used again near the end of
   * this method for another graph; they will be replaced for the needs
   */

  TGraph *gr = new TGraph();

  gr->SetNameTitle(name.c_str(), title.c_str());

  //gr is stored in vec[0]
  vec.push_back(gr);

  gr->GetXaxis()->SetTitle("run time [s]");
  gr->GetYaxis()->SetTitle("current [A]");

  std::vector<energy_t> *energies ;
  bool dark = false;

  if(energies_.find(absStrip(lad,strip)) == energies_.end())
  {
    dark = true; //no hit => strip is dark: fill with noise only
    energies = new std::vector<energy_t> (1,0);
  }
  else
    energies = &energies_[absStrip(lad,strip)];

// fill gr with signal and push_back charge collection graphs

  for(int i = 0; i < (int) energies->size(); ++i)
  {
    //charge collected with noise
    TGraph *charge = GetChargeSignal(energies->at(i));

    name = "charge(:)";
    name.insert(name.length()-3, std::to_string(i));
    name.insert(name.length()-2, std::to_string(lad));
    name.insert(name.length()-1, std::to_string(strip));

    title = "charge collected [hit:  ladder:  strip: ]";
    title.insert(title.length()-18, std::to_string(i));
    title.insert(title.length()-9, std::to_string(lad));
    title.insert(title.length()-1, std::to_string(strip));

    charge->SetNameTitle(name.c_str(), title.c_str());

    charge->GetXaxis()->SetTitle("time from hit [s]");
    charge->GetYaxis()->SetTitle("charge collected [C]");

    vec.push_back(charge);

    AddSignal(gr, charge, times_[absStrip(lad,strip)][i] );
  }

  if(dark)
    delete energies;
}
