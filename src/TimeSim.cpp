
#include "TimeSim.h"


charge_t TimeSim::AddChargeNoise(signal_t *signal)
{
  //total charge noise
  charge_t Q_NOISE =
    random_->Gaus(0, 2*CHARGE_NOISE_) -random_->Gaus(0, CHARGE_NOISE_);

  //make Q_NOISE a multiple of fondamental charge
  Q_NOISE = TMath::Floor(Q_NOISE / FOND_CHARGE) * FOND_CHARGE;

  mytime_t T = TMath::MaxElement(signal->GetN(), signal->GetX());

  //cumulative noise already added
  charge_t q_noise = 0;

  //charge to add every t_sample
  charge_t dq =  Q_NOISE / T * T_SAMPLING_;

  //charge collected is 0 at t=0 => no noise at t=0 => start from i=1
  for(int i = 1; i < signal->GetN(); ++i)
  {
    charge_t q;
    mytime_t t;

    signal->GetPoint(i, t, q);
    signal->SetPoint(i, t, q + q_noise + dq);

    q_noise += dq;
  }

  return q_noise;
}


charge_t TimeSim::GetChargeSignal
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


void TimeSim::AddCurrentSignal
  (const mytime_t &hitTime, signal_t *signal, const signal_t *charge)
{
  //fill with charge derivative

  for(int i=0; i < charge->GetN()-1; ++i)
  {
    mytime_t t1, t2;
    charge_t q1, q2;

    //charge is sorted
    charge->GetPoint(i, t1, q1);
    charge->GetPoint(i+1, t2, q2);

    signal->SetPoint(
      signal->GetN(),
      hitTime + (t2 + t1)*0.5,
      (q2 - q1) / (t2 - t1)
    );
  }

}


void TimeSim::GetCurrentSignal
  (const mytime_t &hitTime, signal_t *signal, const signal_t *charge)
{
  if(signal->GetN() != 0)
    return;

  signal->SetNameTitle("current", "current signal");
  signal->GetXaxis()->SetTitle("run time [s]");
  signal->GetYaxis()->SetTitle("current [A]");

  AddCurrentSignal(hitTime, signal, charge);
}


void TimeSim::GetCurrentSignal(const int &gr, signal_t *signal)
{
  std::map <mytime_t, energy_t> time_energy;
  segm_->GetHits(gr, time_energy);

  for(auto it = time_energy.begin(); it != time_energy.end(); ++it)
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


void TimeSim::SumCurrentSignal
  (signal_t *sum, const signal_t *add)
{

  // get max and min and bin length (t_sample)

  mytime_t t_sample=0;

  mytime_t t_sum_min=0, t_sum_max=0;
  mytime_t t_add_min=0, t_add_max=0;

  current_t tmp;

  sum->GetPoint(0, t_sum_min, tmp);
  sum->GetPoint(sum->GetN()-1, t_sum_max, tmp);

  add->GetPoint(0, t_add_min, tmp);
  add->GetPoint(add->GetN()-1, t_add_max, tmp);


  if(sum->GetN() == 0) //copy add in sum
  {
    for(int i = 0; i < add->GetN(); ++i)
    {
      mytime_t t;
      current_t I;

      add->GetPoint(i, t, I);
      sum->GetPoint(i, t, I);
    }

    return;
  }

  else //get sampling time
    sum->GetPoint(1, t_sample, tmp);


  //"sum" is not void => add points

  for(int i=0; i < add->GetN(); ++i)
  {
    mytime_t t_add = 0, t_sum = 0;
    current_t I_add = 0, I_sum = 0;

    add->GetPoint(i, t_add, I_add);

    double i_sum = (t_add - t_sum_min) / t_sample;

    /* "add" points are not in general on the same t of "sum" ;
     * they need to be aligned to "sum" samples */

     i_sum =
      i_sum - TMath::FloorNint(i_sum) < 0.5 ?
      TMath::FloorNint(i_sum) : TMath::CeilNint(i_sum);



    /* if "add" signal point is < of every "sum" signal point
     * (i.e. i_sum < 0): add zeros at the begin of "sum" signal */

    for(mytime_t t = t_sum_min - t_sample; i_sum < 0; t -= t_sample)
    {
      sum->InsertPointBefore(0, t, 0);
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

        t_sum_max = t;
        //this does not affect next iterations of this cycle
      }



    /* if "add" point was above "sum" range, i_sum should be equal
     * to sum->GetN() now, otherwise it is < sum->GetN() */

    sum->SetPoint(i_sum, t_add + t_sum, I_add + I_sum);
  }

}


mytime_t TimeSim::GetMeas
  (const signal_t *signal, const double threshold)
{
  /* signal could be charge of current; anyway charge_t and current_t
   * are both double => it is used logical type double directly */

  const double peak =
    TMath::MaxElement(signal->GetN(), signal->GetY());

  for(int i = 0; i < signal->GetN(); ++i)
  {
    mytime_t t;
    double y;

    signal->GetPoint(i, t, y);

    if(y > threshold * peak)
      return t;
  }

  return -9999;
}
