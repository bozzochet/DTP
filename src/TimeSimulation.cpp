
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
}


TimeSimulation::~TimeSimulation()
{
  delete up_;
  delete charge_;
  delete random_;
}


charge_t TimeSimulation::GetChargeSignal
  (TGraph *real_charge, const energy_t &energy, const bool fill_graph)
{

//set graph

  if(fill_graph)
  {
    real_charge->SetNameTitle("charge", "charge collected");
    real_charge->GetXaxis()->SetTitle("time from hit [s]");
    real_charge->GetYaxis()->SetTitle("charge collected [C]");
  }

//compute signal and deviations

  charge_t Q = energy / ENERGY_COUPLE_ * FOND_CHARGE_;
  charge_->SetParameter("Q", Q);

//signal
  charge_t q = 0;

  for( mytime_t t = T_SAMPLING_; t < T_END_; t += T_SAMPLING_ )
  {
    //ideal charge collected
    charge_t dq = charge_->Eval(t) - charge_->Eval(t - T_SAMPLING_);

    charge_t q_noise = dq / 10.0 * random_->Uniform();

    //charge collected at step t
    q += dq + q_noise;

    if(fill_graph)
      real_charge->SetPoint(real_charge->GetN(), t, q);
  }

  if(fill_graph)
    real_charge->Sort();

  return Q;
}


void TimeSimulation::AddSignal
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


void TimeSimulation::GetSignal
  (TGraph *signal, const TGraph *charge, const mytime_t &hitTime)
{

//set graph

  signal->SetNameTitle("current", "current signal");
  signal->GetXaxis()->SetTitle("run time [s]");
  signal->GetYaxis()->SetTitle("current [A]");

//fill signal

  AddSignal(signal, charge, hitTime);
}
