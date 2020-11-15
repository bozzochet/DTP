
/**********************************************************************

  TimeSim class simulates signal of current in time
  generated in the strips by particles hits.

  The signal simulated, line from 0 to peak and an exponential
  descent, based on collected charge.

  Charge collecting process follows the capacitor charging law and
  noise is added to charge signal.

 *********************************************************************/

#ifndef TIME_SIM_INCLUDE
#define TIME_SIM_INCLUDE


#include "physics.h"
#include "absStrip.h"
#include "TrCluster.hh"
#include "vector2.h"
#include "TimeSegm.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TList.h"
#include "TObject.h"
#include "TGraph.h"
#include "TF1.h"
#include "TAxis.h"
#include "TRandom3.h"
#include "TH1F.h"

#include <vector>
#include <string>
#include <cstring>
#include <iostream>


class TimeSim
{

  typedef TGraph signal_t;
  typedef TH1F hist_t;
  typedef TF1 signal_fun_t;
  typedef TRandom3 random_gen_t;


// const variables

  //tau of strip as a capacitor
  const mytime_t T_CAPACITOR_ = 1e-9;

  //charge function is stopped when reaching this fraction of peak
  const double STOP_CHARGE_FRACTION_ = 0.999;

  //signal parameters
  const mytime_t T_SAMPLING_ = 1e-12;

  /* noise deviation:
   *   pair/um * thickness * fond_charge */
  const charge_t CHARGE_NOISE_ = 8 * 300 * FOND_CHARGE;



// variables

  TimeSegm *segm_ = NULL;

  //random generator
  random_gen_t *random_ = NULL;



// methods

  /* add to signal a single hit current generated using up_ and
   * charge signal obtained by GetChargeSignal; signal points are
   * sorted after this method execution*/
  void AddCurrentSignal(signal_t *signal, const signal_t *charge);


  // signal points are sorted after this method execution
  inline void AddChargeSignal
    (signal_t *signal, const signal_fun_t *ideal)
  {
    for
    (
      mytime_t t = ideal->GetXmin();
      t <= ideal->GetXmax();
      t += T_SAMPLING_
    )
      signal->SetPoint(signal->GetN(), t, ideal->Eval(t));
  }


  //cumulative function of ideal charge collected
  signal_fun_t* GetIdealChargeFun(const mytime_t&, const charge_t&);


  /* add second signal to first one passed;
   * signals passed MUST BE sorted and sampled with same time;
   * first signal is sorted after this method execution */
  //void SumCurrentSignal(signal_t*, const signal_t*);



public:

  TimeSim
    (TimeSegm *segm, random_gen_t *random)
  { segm_ = segm; random_ = random; }


  inline charge_t GetChargeFromEnergy(const energy_t &E)
  { return TMath::Floor(E / ENERGY_COUPLE) * FOND_CHARGE; }


  inline charge_t GetChargeNoise()
  {
    return
      TMath::Floor
      (
        (random_->Gaus(0, 2*CHARGE_NOISE_)-random_->Gaus(0, CHARGE_NOISE_)) / FOND_CHARGE
      ) * FOND_CHARGE;
  }


  /* add noise to signal passed and return total charge noise
   * collected;
   * IMPORTANT: signals passed MUST be sorted */
  charge_t AddChargeNoise(signal_t *signal);


  /* generate charge signal in time with noise;
   * return total charge collected;
   * signal passed MUST BE VOID;
   * signal points are sorted after this method execution */
  charge_t GetChargeSignal
  (
    const mytime_t &hitTime, const energy_t&,
    signal_t *signal, const bool noise = true
  );


  /* get current signal on strip ( #ladder, #strip) based on charge
   * collected;
   * signal passed MUST BE VOID;
   * signal points are sorted after this method execution */
  void GetCurrentSignal(signal_t *signal, const signal_t *charge);


  /* get charge signal of i-th group of strips;
   * signal passed MUST BE VOID;
   * signal points are sorted after this method execution;
   * return number of hits for this group */
   int GetChargeSignal(const int &i, signal_t *signal);


   /* get time when charge becomes > threshold_fraction * peak
    * IMPORTANT: charge MUST be sorted */
   mytime_t GetMeas
    (const signal_t *charge, const double threshold_fraction = 0.15);

};


#endif //include guard
