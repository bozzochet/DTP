
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

/*********************************************************************/

  //absolute strip ID
  inline abs_strip_t abs_strip(const int &ladder, const int &strip)
  { return ladder*Nstrips + strip; }

  inline bool is_active(const abs_strip_t &strip)
  {
    if(strip % Nrows*Nstrips % jump_ == 0) return true;
    return false;
  }

/*********************************************************************/

public:

  Stopwatch(const int &jump)
  { jump_ = jump; }

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

};

#endif //include guard
