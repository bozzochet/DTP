
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
  (const energy_t &energy, signal_t *signal, const bool noise)
{
  if(signal->GetN() != 0)
    return 0;

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
  (mytime_t hitTime, signal_t *signal, const signal_t *charge)
{
  //align hitTime to samples
  hitTime = TMath::Floor(hitTime / T_SAMPLING_) * T_SAMPLING_;

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
  (mytime_t hitTime, signal_t *signal, const signal_t *charge)
{
  if(signal->GetN() != 0)
    return;

  signal->SetNameTitle("current", "current signal");
  signal->GetXaxis()->SetTitle("run time [s]");
  signal->GetYaxis()->SetTitle("current [A]");

  AddCurrentSignal(hitTime, signal, charge);
}


void TimeSimulation::GetCurrentSignal(signal_t *signal, const int &i)
{
  std::map <mytime_t, energy_t> time_energy_;
  segm_->GetHits(time_energy_, i);

  for(auto it = time_energy_.begin(); it != time_energy_.end(); ++it)
  {
    signal_t *charge = new signal_t();
    GetChargeSignal(it->second, charge);

    signal_t *current = new signal_t();
    GetCurrentSignal(it->first, current, charge);

    SumCurrentSignal(signal, current);

    delete charge;
    delete current;
  }
}


void TimeSimulation::SumCurrentSignal
  (signal_t *sum, const signal_t *add)
{

  // get max and min and bin length (t_sample)

  mytime_t t_sum_sample, t_sum_min, t_sum_max;
  mytime_t t_add_sample, t_add_min, t_add_max;
  current_t tmp;

  sum->GetPoint(0, t_sum_min, tmp);
  sum->GetPoint(1, t_sum_sample, tmp);
  sum->GetPoint(sum->GetN()-1, t_sum_max, tmp);

  add->GetPoint(0, t_add_min, tmp);
  add->GetPoint(1, t_add_sample, tmp);
  add->GetPoint(add->GetN()-1, t_add_max, tmp);

  if(t_add_sample != t_sum_sample)
    return;

  mytime_t t_sample = t_sum_sample;

  //sum

  for(int i=0; i < add->GetN(); ++i)
  {
    mytime_t t_add = 0, t_sum = 0;
    current_t I_add = 0, I_sum = 0;

    add->GetPoint(i, t_add, I_add);

    int i_sum = (t_add - t_sum_min) / t_sample;

    /* "add" signal point is < of every "sum" signal point: add zeros
     * at the begin of "sum" signal */

    for(mytime_t t = t_sum_min - t_sample; i_sum < 0; t -= t_sample)
    {
      sum->InsertPointBefore(1, t, 0);
      t_sum_min = t;
      i_sum = (t_add - t_sum_min) / t_sample;
    }

    // evaluate if "add" signal point is inside "sum" range or above

    if(i_sum < sum->GetN())
      sum->GetPoint(i_sum, t_sum, I_sum);

    /* append 0 to "sum" signal until "sum" range reaches "add" signal
     * point */

    else
      for
      (
        mytime_t t = t_sum_max + t_sample;
        i_sum > sum->GetN();
        t += t_sample
      )
      {
        sum->SetPoint(sum->GetN(), t, 0);
        t_sum_max = t; //this does not affect t in the next iterations
      }

    //if "add" point was above, i_sum should be equal to sum->GetN()

    sum->SetPoint(i_sum, t_add + t_sum, I_add + I_sum);
  }

}


mytime_t TimeSimulation::GetTime
  (const signal_t *current, const double threshold)
{
  const current_t I_max = current->GetMaximum();

  for(int i = 0; i < current->GetN(); ++i)
  {
    mytime_t t;
    current_t I;

    current->GetPoint(i, t, I);

    if(I > threshold * I_max)
      return t;
  }

  return -9999;
}
