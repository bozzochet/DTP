
#include "TimeSim.h"


TimeSim::signal_fun_t* TimeSim::GetChargeNoiseFun
  (const mytime_t &t_min, const mytime_t &t_max, const charge_t &Q)
{
  signal_fun_t *noise = new signal_fun_t
  (
    "cumulative noise charge collected;time [s];charge [C]",
    "[0] * (x - [1])", t_min, t_max
  );

  noise->SetParameter(0, Q / (t_max - t_min));
  noise->SetParameter(1, t_min);

  return noise;
}


charge_t TimeSim::AddChargeNoise(signal_t *signal)
{
  //total charge noise
  charge_t Q_NOISE = GetChargeNoise();

  signal_fun_t *noise = GetChargeNoiseFun
  (
    TMath::MinElement(signal->GetN(), signal->GetX()),
    TMath::MaxElement(signal->GetN(), signal->GetX()),
    Q_NOISE
  );

  for(int i = 0; i < signal->GetN(); ++i)
  {
    charge_t q;
    mytime_t t;

    signal->GetPoint(i, t, q);
    signal->SetPoint(i, t, q + noise->Eval(t));
  }

  return Q_NOISE;
}


TimeSim::signal_fun_t* TimeSim::GetChargeIdealFun
  (const mytime_t &hitTime, const charge_t &Q)
{
  signal_fun_t *charge = new signal_fun_t
  (
    "charge",
    "[0] * (1 - TMath::Exp( -[1]*(x-[2]) ) )",
    hitTime, hitTime - T_CAPACITOR_*TMath::Log(1-STOP_CHARGE_FRACTION_)
  );

  charge->SetParameter(0, Q);
  charge->SetParameter(1, 1.0 / T_CAPACITOR_);
  charge->SetParameter(2, hitTime);

  return charge;
}


charge_t TimeSim::GetChargeSignal
(
  const mytime_t &hitTime, const energy_t &energy,
  signal_t *signal, const bool noise
)
{
  if(signal->GetN() != 0)
    return 0;

  signal->SetNameTitle("charge", "charge collected");
  signal->GetXaxis()->SetTitle("time [s]");
  signal->GetYaxis()->SetTitle("charge collected [C]");

  charge_t Q = GetChargeFromEnergy(energy);

  signal_fun_t *charge = GetChargeIdealFun(hitTime, Q);

  AddChargeSignal(signal, charge);

  delete charge;

  if(noise) Q += AddChargeNoise(signal);

  return Q;
}


void TimeSim::AddCurrentSignal
  (signal_t *signal, const signal_t *charge)
{
  //fill with charge derivative

  for(int i=0; i < charge->GetN()-1; ++i)
  {
    mytime_t t1, t2;
    charge_t q1, q2;

    //charge is sorted
    charge->GetPoint(i, t1, q1);
    charge->GetPoint(i+1, t2, q2);

    signal->SetPoint
      (signal->GetN(), (t2 + t1)*0.5, (q2 - q1) / (t2 - t1) );
  }

}


void TimeSim::GetCurrentSignal
  (signal_t *signal, const signal_t *charge)
{
  if(signal->GetN() != 0)
    return;

  signal->SetNameTitle("current", "current signal");
  signal->GetXaxis()->SetTitle("run time [s]");
  signal->GetYaxis()->SetTitle("current [A]");

  AddCurrentSignal(signal, charge);
}

/* THIS METHOD NEEDS TO BE REVISITED: signal generated has non noise;
 * how to collect noise ?
 *
int TimeSim::GetChargeSignal(const int &gr, signal_t *signal)
{
  std::map <mytime_t, energy_t> time_energy;
  segm_->GetHits(gr, time_energy);

  if( (int) time_energy.size() == 0)
    return 0; //no hit on this group


  std::vector<signal_fun_t*> charge;
  std::vector<mytime_t> hit_time;

  //get ideal charge functions and store them and hit times in vecs

  for(auto it = time_energy.begin(); it != time_energy.end(); ++it)
  {
    charge.push_back
    (
      GetChargeIdealFun(it->first, GetChargeFromEnergy(it->second))
    );

    hit_time.push_back(it->first);
  }


  //sample group charge signal and fill signal parameter

  for
  (
    mytime_t t = hit_time[0];
    t <= charge.back()->GetXmax();
    t += T_SAMPLING_
  )
  {
    charge_t q = 0; //q(t)


    //sum single hit contributions to q(t)

    for(int i=0; i < (int) charge.size(); ++i)
      if(t >= charge[i]->GetXmin() && t <= charge[i]->GetXmax())
        //t is in charge[i] range
        q += charge[i]->Eval(t);


    signal->SetPoint(signal->GetN(), t, q);

  } //for t

  return (int) time_energy.size();
}
*/


/* DEPRECATED : look at GetChargeSignal(gr,signal) above; it is
 * commented but it has a better approach
 *
 *
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

    //"add" points are not in general on the same t of "sum" ;
    // they need to be aligned to "sum" samples

     i_sum =
      i_sum - TMath::FloorNint(i_sum) < 0.5 ?
      TMath::FloorNint(i_sum) : TMath::CeilNint(i_sum);



    // if "add" signal point is < of every "sum" signal point
    // (i.e. i_sum < 0): add zeros at the begin of "sum" signal

    for(mytime_t t = t_sum_min - t_sample; i_sum < 0; t -= t_sample)
    {
      sum->InsertPointBefore(0, t, 0);
      t_sum_min = t;
      i_sum = (t_add - t_sum_min) / t_sample;
    }



    // evaluate if "add" signal point is inside "sum" range or above

    if(i_sum < sum->GetN())
      sum->GetPoint(i_sum, t_sum, I_sum);

    // append 0 to "sum" signal until "sum" range reaches "add" signal
    // point

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



    // if "add" point was above "sum" range, i_sum should be equal
    // to sum->GetN() now, otherwise it is < sum->GetN()

    sum->SetPoint(i_sum, t_add + t_sum, I_add + I_sum);
  }

}
*/


mytime_t TimeSim::GetMeas
  (const signal_t *charge, const double threshold)
{
  const charge_t peak =
    TMath::MaxElement(charge->GetN(), charge->GetY());

  for(int i = 0; i < charge->GetN(); ++i)
  {
    mytime_t t;
    double y;

    charge->GetPoint(i, t, y);

    if(y > threshold * peak)
      return t;
  }

  return -9999;
}
