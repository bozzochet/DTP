
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

  //slew rate of line_
  const double SLEW_RATE_ = 1e+5;

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
  void AddCurrentSignal
  (
    const mytime_t &hitTime, signal_t *signal, const signal_t *charge
  );


  // signal points are sorted after this method execution
  inline void AddChargeSignal
    (signal_t *signal, const signal_fun_t *ideal)
  {
    for(mytime_t t = 0; t < ideal->GetXmax(); t += T_SAMPLING_ )
      signal->SetPoint(signal->GetN(), t, ideal->Eval(t));
  }


  inline charge_t GetChargeFromEnergy(const energy_t &E)
  { return E / ENERGY_COUPLE * FOND_CHARGE; }


  bool Trigger(double&, const signal_t*, const int&, const double&);


  /* add second signal to first one passed;
   * signals passed MUST BE sorted and sampled with same time;
   * first signal is sorted after this method execution */
  void SumCurrentSignal(signal_t*, const signal_t*);



public:

  TimeSim
    (TimeSegm *segm, random_gen_t *random)
  { segm_ = segm; random_ = random; }


  /* add noise to signal passed and return total charge noise
   * collected;
   * IMPORTANT: signals passed MUST be sorted */
  charge_t AddChargeNoise(signal_t *signal);


  /* generate charge signal in time with noise;
   * return charge created by hit;
   * signal passed MUST BE VOID;
   * signal points are sorted after this method execution */
  charge_t GetChargeSignal
  (
    const energy_t&, signal_t *signal, const bool noise = true
  );


  /* get current signal on strip ( #ladder, #strip) based on charge
   * collected;
   * signal passed MUST BE VOID;
   * signal points are sorted after this method execution */
  void GetCurrentSignal
  (
    const mytime_t &hitTime, signal_t *signal, const signal_t *charge
  );


  /* get current signal i-th group of strips;
   * signal passed MUST BE VOID;
   * signal points are sorted after this method execution */
   void GetCurrentSignal(const int &i, signal_t *signal);


   /* get time when current becomes > threshold_fraction * peak
    * IMPORTANT: current MUST be sorted */
   mytime_t GetTime
    (const signal_t *current, const double threshold_fraction = 0.1);

};


#endif //include guard
