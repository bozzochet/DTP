
/*****************************************************************
  Stopwatch class simulates the time development of the current
  signal generated in the strip by the passage of a particle.

  GGS times are stored in a vector associated to a strip through a
  map.

  The signal simulated, for a first implementation, is a sawtooth
  based on T_PEAK_ and T_RELAX_ variables. The current is normalized
  to the maximum current, so the peak of a single signal has a value
  of 1; anyway, the sum of more than one signal could give a
  value > 1.

 ****************************************************************/

#ifndef STOPWATCH_INCLUDE
#define STOPWATCH_INCLUDE

#include "global.h"

#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TAxis.h"

#include <map>
#include <vector>
#include <iostream>
#include <string>

typedef int abs_strip_t;

typedef double mytime_t; //times are expressed in ns

typedef std::map <abs_strip_t, std::vector <mytime_t> > times_map_t;

class Stopwatch
{
  //signal histogram parameters
  const int N_BINS_ = 1000;
  const mytime_t T_MIN_ = 0;
  const mytime_t T_MAX_ = 4;
  const double BIN_LENGTH_ = (T_MAX_ - T_MIN_) / (double) N_BINS_;

  //time to reach the peak after a hit
  const mytime_t T_PEAK_ = 0.5;

  //relaxing time between hit and when current returns to 0
  const mytime_t T_RELAX_ = 2;


  //active strips are one every jump, beginning from a row
  int jump_ = 1;

  /* a hit happened between strip s and strip s+1 (s = absolute ID)
   * is stored in original_ with key s and value t */
  times_map_t original_;

  //active_ register times of active strips
  times_map_t active_;

  //ideal signal
  TF1 *signal_up_ = NULL; //from 0 to peak
  TF1 *signal_down_ = NULL; //from peak to 0


  //absolute strip ID
  inline abs_strip_t abs_strip(const int &ladder, const int &strip)
  { return ladder*Nstrips + strip; }

  inline bool is_active(const abs_strip_t &strip)
  {
    if(strip % Nrows*Nstrips % jump_ == 0) return true;
    return false;
  }

/*
  //fill hist with noise
  void add_noise(TH1D *hist);
*/

  //fill hist with simulated signal
  void add_signal(TH1D *hist, const std::vector<mytime_t> &times);


public:

  //get jump and set TF1 private variables
  Stopwatch(const int &jump);

  ~Stopwatch();

  //take time on strip
  inline void split
    (const int &ladder, const int &strip, const mytime_t &t)
  { original_[abs_strip(ladder,strip)].push_back(t); }

  //get hit times on (ladder,strip); if not a hit return void vector
  inline std::vector<mytime_t> time
    (const int &ladder, const int &strip)
  {
    if(active_.find(abs_strip(ladder,strip)) != active_.end())
      return active_[abs_strip(ladder,strip)];
    return {}; //void vector
  }

  inline void reset()
  { original_.clear(); active_.clear(); }

  //all times_map taken
  void stop();

  //return hist with simulated signal of strip passed
  TH1D* get_signal(const int &ladder, const int &strip);

};

#endif //include guard
