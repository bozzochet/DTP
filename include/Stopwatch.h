
/*****************************************************************
  Stopwatch class simulate detector time profiling

  BACKSCATTERING NOT SUPPORTED:
    one time for strip at maximum is stored

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

typedef double mytime_t;

typedef std::map <abs_strip_t, mytime_t> times_map_t;

class Stopwatch
{
  //active strips are one every jump, beginning from a row
  int jump_ = 1;

  /* a hit happened between strip s and strip s+1 (s = absolute ID)
   * is stored in original_ with key s and value t */
  times_map_t original_;

  //active_ register times of active strips
  times_map_t active_;

  //signal histogram parameters
  const int N_BINS_ = 1000;
  const mytime_t T_MIN_ = 0;
  const mytime_t T_MAX_ = 4;
  const double BIN_LENGTH_ = (T_MAX_ - T_MIN_) / (double) N_BINS_;

  //time to reach the peak after a hit
  const mytime_t T_PEAK_ = 0.5;

  //relaxing time between hit and when current returns to 0
  const mytime_t T_RELAX_ = 2;

  //ideal signal
  TF1 *signal_up_ = NULL; //from 0 to peak
  TF1 *signal_down_ = NULL; //from peak to 0

/*********************************************************************/

  //absolute strip ID
  inline abs_strip_t abs_strip(const int &ladder, const int &strip)
  { return ladder*Nstrips + strip; }

  inline bool is_active(const abs_strip_t &strip)
  {
    if(strip % Nrows*Nstrips % jump_ == 0) return true;
    return false;
  }

  //fill hist with simulated signal
  void add_signal(TH1D *hist, const mytime_t &hitTime);

/*********************************************************************/

public:

  //get jump and set TF1 private variables
  Stopwatch(const int &jump);
  
  ~Stopwatch();

  //take time on strip
  inline void split
    (const int &ladder, const int &strip, const mytime_t &t)
  { original_[abs_strip(ladder,strip)] = t; }

  //all times_map taken
  void stop();

  //get hit time on (ladder,strip); if not a hit return negative value
  inline mytime_t time(const int &ladder, const int &strip)
  {
    if(active_.find(abs_strip(ladder,strip)) != active_.end())
      return active_[abs_strip(ladder,strip)];
    return -1;
  }

  inline void reset()
  { original_.clear(); active_.clear(); }

  TH1D* get_signal(const int &ladder, const int &strip);

};

#endif //include guard
