
#include "TimeSimulation.h"


TimeSimulation::TimeSimulation()
{
  up_ = new TF1("up", "[0]*x", T_START_, T_END_);

  up_->SetParName(0, "Slew rate");
  up_->SetParameter("Slew rate", SLEW_RATE_);

  random_ = new TRandom3(9298);
}


TimeSimulation::~TimeSimulation()
{
  delete up_;
  delete random_;
}


void TimeSimulation::AddChargeSignal(TGraph *signal, const TF1 *ideal)
{
  for(mytime_t t = 0; t < ideal->GetXmax(); t += T_SAMPLING_ )
    signal->SetPoint(signal->GetN(), t, ideal->Eval(t));

  signal->Sort();
}

charge_t TimeSimulation::AddChargeNoise(TGraph *signal, const charge_t &Q)
{
  TH1F *uniform = new TH1F
    ("h_uniform", "uniform", signal->GetN()-1, 0, 1e+6);

  for(int i = 0; i < 1e+4; ++i)
    uniform->Fill(random_->Uniform(0, 1e+6));

  int integral = 0; //integral of uniform

  //charge collected is 0 at t=0 => no noise at t=0 => start from i=1
  for(int i = 1; i < signal->GetN(); ++i)
  {
    charge_t q;
    mytime_t t;

    signal->GetPoint(i, t, q);

    integral += uniform->GetBinContent(i);

    charge_t dq = Q * CHARGE_NOISE_ *
      integral / uniform->GetSumOfWeights();

    //make dq a multiple of fondamental charge
    dq = TMath::Floor(dq / FOND_CHARGE) * FOND_CHARGE;

    signal->SetPoint(i, t, q + dq);
  }

  signal->Sort();

  return Q * CHARGE_NOISE_ ;
}


charge_t TimeSimulation::GetChargeSignal
  (TGraph *signal, const energy_t &energy, const bool noise)
{
  signal->SetNameTitle("charge", "charge collected");
  signal->GetXaxis()->SetTitle("time from hit [s]");
  signal->GetYaxis()->SetTitle("charge collected [C]");

  //cumulative function of ideal charge collected

  charge_t Q = GetChargeFromEnergy(energy);

  TF1 *charge = new TF1(
    "charge",
    "[0] * (1 - TMath::Exp(-[1]*x) )",
    0, - T_CAPACITOR_ * TMath::Log(1-STOP_CHARGE_FRACTION_)
  );

  charge->SetParName(0, "Q");
  charge->SetParName(1, "k");

  charge->SetParameter("k", 1.0 / T_CAPACITOR_);
  charge->SetParameter("Q", Q);

  AddChargeSignal(signal, charge);
  delete charge;

  if(noise) Q += AddChargeNoise(signal, Q);

  return Q;
}


void TimeSimulation::AddCurrentSignal
  (TGraph *signal, const TGraph *charge, const mytime_t &hitTime)
{

//fill with charge derivative

  current_t peak;
  mytime_t t_peak; //delta_t between hitTime and current = peak

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

    signal->SetPoint(
      signal->GetN(),
      hitTime + t_peak + (t2 + t1)*0.5,
      (q2 - q1) / (t2 - t1)
    );

  }

//fill signal with up_

  for( mytime_t t = 0; t < t_peak; t += T_SAMPLING_ )
    signal->SetPoint(signal->GetN(), t + hitTime, up_->Eval(t));

  signal->Sort();
}


void TimeSimulation::GetCurrentSignal
  (TGraph *signal, const TGraph *charge, const mytime_t &hitTime)
{
  signal->SetNameTitle("current", "current signal");
  signal->GetXaxis()->SetTitle("run time [s]");
  signal->GetYaxis()->SetTitle("current [A]");

  AddCurrentSignal(signal, charge, hitTime);
}
