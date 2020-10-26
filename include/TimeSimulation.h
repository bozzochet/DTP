
/**********************************************************************

  TimeSimulation class simulates the time development of the current
  signal generated in the strip by the passage of a particle.

  The signal simulated, line from 0 to peak and an exponential
  descent, based on T_PEAK_ and T_RELAX_ variables. The current is
  normalized to the maximum current, so the peak of a single signal
  has a value of 1; anyway, the sum of more than one signal could
  give a value > 1.

 *********************************************************************/

#ifndef TIME_SIMULATION_INCLUDE
#define TIME_SIMULATION_INCLUDE


#include "global.h"
#include "types.h"
#include "Stopwatch.h"

#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TAxis.h"

#include <vector>
#include <string>
#include <iostream>


class TimeSimulation
{
  //signal histogram parameters
  const int N_BINS_ = 1000;
  const mytime_t T_MIN_ = 0;
  const mytime_t T_MAX_ = 4;
  const double BIN_LENGTH_ = (T_MAX_ - T_MIN_) / (double) N_BINS_;

  //signal generation parameters
  const double PEAK_VALUE_ = 1; //signal normalize to peak current
  const mytime_t T_PEAK_ = 0.1; // = 100 ps ; time to reach peak
  const double SLEW_RATE_ = PEAK_VALUE_ / T_PEAK_ ;
  const double K_EXP_ = 1; // k in formula exp(-kx) for signal descent
  const double ZERO_THRESH_ = 0.05; //value under wich exp is considered = 0


  //ideal signal
  TF1 *signal_up_ = NULL; //from 0 to peak
  TF1 *signal_down_ = NULL; //from peak to 0


/*
  //fill hist with noise
  void add_noise(TH1D *hist);
*/

  //fill hist with simulated signal
  void add_signal(TH1D *hist, const std::vector<mytime_t> &times);


public:

  //create TF1 functions
  TimeSimulation();

  ~TimeSimulation();

  //return hist with simulated signal of strip passed
  TH1D* get_signal(const int &ladder, const int &strip, Stopwatch*);
};


#endif //include guard
