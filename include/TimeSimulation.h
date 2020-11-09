
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
  const double SLEW_RATE_ = 1e+4;

  //tau of strip as a capacitor
  const mytime_t T_CAPACITOR_ = 1e-9;

  //charge function is stopped when reaching this fraction of peak
  const double STOP_CHARGE_FRACTION_ = 0.999;

  //signal parameters
  const mytime_t T_START_ = 0;
  const mytime_t T_END_ = 5e-9;
  const int N_BINS_ = 1000;
  const mytime_t T_SAMPLING_ = (T_END_ - T_START_) / (double) N_BINS_;


// variables

  Noise *noise_ = NULL;

  signal_fun_t *up_ = NULL; //slope of current signal to the peak

  //random generator
  random_gen_t *random_ = NULL;


// methods

  /* add to signal a single hit current generated using up_ and
   * charge signal obtained by GetChargeSignal */
  void AddCurrentSignal
    (signal_t *signal, const signal_t *charge, const mytime_t &hitTime);

  void AddChargeSignal(signal_t *signal, const signal_fun_t *ideal);

  //return total charge noise
  charge_t AddChargeNoise(signal_t *signal);

  inline charge_t GetChargeFromEnergy(const energy_t &E)
  { return E / ENERGY_COUPLE * FOND_CHARGE; }


public:

  TimeSimulation(const double&);
  virtual ~TimeSimulation();

  //add noise to signal passed
  inline charge_t GetChargeNoise(signal_t *signal)
  { return AddChargeNoise(signal); }

  inline charge_t GetTotalChargeNoise()
  { return noise_->GetChargeNoise(); }

  /* generate charge signal in time with noise;
   * return charge created by hit */
  charge_t GetChargeSignal
    (signal_t *signal, const energy_t&, const bool noise = true);

  /* get current signal on strip ( #ladder, #strip) based on charge
   * collected */
  void GetCurrentSignal
    (signal_t *signal, const signal_t *charge, const mytime_t &hitTime);

  /*
   * NEED IMPLEMENTATION: in case of multiple hits on the same strip
   * GetSignal is able to generate the multiple signals and append them
   * to the total signal returnes; at this moment the signals are NOT
   * summed (as has to be). Multiple currents are just drawn on the
   * same graph
   *
   */
};


#endif //include guard
