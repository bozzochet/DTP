
#include "TimeSimulation.h"


void TimeSimulation::AddChargeSignal(signal_t *signal, const signal_fun_t *ideal)
{
  for(mytime_t t = 0; t < ideal->GetXmax(); t += T_SAMPLING_ )
    signal->SetPoint(signal->GetN(), t, ideal->Eval(t));
}


charge_t TimeSimulation::AddChargeNoise(signal_t *signal)
{
  hist_t *uniform = new hist_t
    ("h_uniform", "uniform", signal->GetN()-1, 0, 1e+6);

  for(int i = 0; i < 1e+4; ++i)
    uniform->Fill(random_->Uniform(0, 1e+6));

  int integral = 0; //integral of uniform
  charge_t q_noise = 0; //total noise added

  //charge collected is 0 at t=0 => no noise at t=0 => start from i=1
  for(int i = 1; i < signal->GetN(); ++i)
  {
    charge_t q;
    mytime_t t;

    signal->GetPoint(i, t, q);

    integral += uniform->GetBinContent(i);

    charge_t dq = noise_->GetChargeNoise() *
      integral / uniform->GetSumOfWeights();

    //make dq a multiple of fondamental charge
    dq = TMath::Floor(dq / FOND_CHARGE) * FOND_CHARGE;

    q_noise += dq;

    signal->SetPoint(i, t, q + dq);
  }

  delete uniform;

  return q_noise;
}


charge_t TimeSimulation::GetChargeSignal
  (signal_t *signal, const energy_t &energy, const bool noise)
{
  signal->SetNameTitle("charge", "charge collected");
  signal->GetXaxis()->SetTitle("time from hit [s]");
  signal->GetYaxis()->SetTitle("charge collected [C]");

  charge_t Q = GetChargeFromEnergy(energy);

  //cumulative function of ideal charge collected

  signal_fun_t *charge = new signal_fun_t(
    "charge",
    "[0] * (1 - TMath::Exp(-[1]*x) )",
    0, - T_CAPACITOR_ * TMath::Log(1 - STOP_CHARGE_FRACTION_)
  );

  charge->SetParameter(0, Q);
  charge->SetParameter(1, 1.0 / T_CAPACITOR_);

  AddChargeSignal(signal, charge);

  delete charge;

  if(noise) Q += AddChargeNoise(signal);

  return Q;
}


void TimeSimulation::AddCurrentSignal
  (signal_t *signal, const signal_t *charge, const mytime_t &hitTime)
{

//fill with charge derivative

  //position of current graph first point
  current_t I_first;
  mytime_t t_first;

  for(int i=0; i < charge->GetN()-1; ++i)
  {
    mytime_t t1, t2;
    charge_t q1, q2;

    //charge is sorted
    charge->GetPoint(i, t1, q1);
    charge->GetPoint(i+1, t2, q2);

    if(i==0)
    {
      I_first = (q2 - q1) / (t2 - t1);
      t_first = I_first / SLEW_RATE_;
    }

    signal->SetPoint(
      signal->GetN(),
      hitTime + t_first + (t2 + t1)*0.5,
      (q2 - q1) / (t2 - t1)
    );

  }

//fill first part of signal

  signal_fun_t *line = new signal_fun_t("line", "[0]*x", 0, t_first);
  line->SetParameter(0, SLEW_RATE_);

  for( mytime_t t = 0; t < t_first; t += T_SAMPLING_ )
    signal->InsertPointBefore
      ( t / T_SAMPLING_ , t + hitTime, line->Eval(t));

  delete line;
}


void TimeSimulation::GetCurrentSignal
  (signal_t *signal, const signal_t *charge, const mytime_t &hitTime)
{
  signal->SetNameTitle("current", "current signal");
  signal->GetXaxis()->SetTitle("run time [s]");
  signal->GetYaxis()->SetTitle("current [A]");

  AddCurrentSignal(signal, charge, hitTime);
}
