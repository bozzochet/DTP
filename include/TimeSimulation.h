
/**********************************************************************

  TimeSimulation class simulates signal of current in time
  generated in the strips by particles hits.

  The signal simulated, line from 0 to peak and an exponential
  descent, based on collected charge.

  Charge collecting process follows the capacitor charging law and
  noise is added to charge signal.

 *********************************************************************/

#ifndef TIME_SIMULATION_INCLUDE
#define TIME_SIMULATION_INCLUDE


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





class TimeSimulation
{

  typedef TGraph signal_t;
  typedef TH1F hist_t;
  typedef TF1 signal_fun_t;
  typedef TRandom3 random_gen_t;


  //make CHARGE_NOISE_ not editable after initialization
  class Noise
  {
    charge_t CHARGE_NOISE_ = 0;

  public:

    Noise(const length_t &thickness)
    // Q_noise ~ 8 um^(-1) * thickness * e
    { CHARGE_NOISE_ = thickness * (8 * 1e+6) * FOND_CHARGE; }

    inline charge_t GetChargeNoise()
    { return CHARGE_NOISE_; }
  };


// const variables

  //slew rate of line_
  const double SLEW_RATE_ = 1e+5;

  //tau of strip as a capacitor
  const mytime_t T_CAPACITOR_ = 1e-9;

  //charge function is stopped when reaching this fraction of peak
  const double STOP_CHARGE_FRACTION_ = 0.999;

  //signal parameters
  const mytime_t T_SAMPLING_ = 1e-11;


// variables

  TimeSegm *segm_ = NULL;

  Noise *noise_ = NULL;

  //random generator
  random_gen_t *random_ = NULL;


// methods

  /* add to signal a single hit current generated using up_ and
   * charge signal obtained by GetChargeSignal; signal points are
   * sorted after this method execution*/
  void AddCurrentSignal
  (
    signal_t *signal, const signal_t *charge,
    const mytime_t &hitTime
  );

  // signal points are sorted after this method execution
  void AddChargeSignal(signal_t *signal, const signal_fun_t *ideal);

  inline charge_t GetChargeFromEnergy(const energy_t &E)
  { return E / ENERGY_COUPLE * FOND_CHARGE; }

  bool Trigger(double&, const signal_t*, const int&, const double&);


public:

  TimeSimulation(TimeSegm *segm, const double &thickness)
  {
    segm_ = segm;
    noise_ = new Noise(thickness);
    random_ = new random_gen_t();
  };

  virtual ~TimeSimulation()
  {
    delete noise_;
    delete random_;
  }

  /* add noise to signal passed and return total charge noise
   * collected;
   * IMPORTANT: signals passed MUST be sorted */
  charge_t AddChargeNoise(signal_t *signal);

  inline charge_t GetTotalChargeNoise()
  { return noise_->GetChargeNoise(); }

  /* generate charge signal in time with noise;
   * return charge created by hit; signal points are sorted after this
   * method execution */
  charge_t GetChargeSignal
  (
    signal_t *signal, const energy_t&,
    const bool noise = true
  );

  /* get current signal on strip ( #ladder, #strip) based on charge
   * collected; signal points are sorted after this method execution */
  void GetCurrentSignal
  (
    signal_t *signal, const signal_t *charge, const mytime_t &hitTime
  );

  /*
   * IMPLEMENTATION NEEDED: in case of multiple hits on the same strip
   * GetSignal is able to generate the multiple signals and append them
   * to the total signal returnes; at this moment the signals are NOT
   * summed (as has to be). Multiple currents are just drawn on the
   * same graph
   *
   */

   /* get time when current becomes > threshold_fraction * peak
    * IMPORTANT: current MUST be sorted */
   mytime_t GetTime
    (const signal_t *current, const double threshold_fraction = 0.1);

};


#endif //include guard
