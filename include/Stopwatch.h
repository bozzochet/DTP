
// Stopwatch handles GGS hit times

#ifndef STOPWATCH_INCLUDE
#define STOPWATCH_INCLUDE


#include "global.h"
#include "abs_strip.h"
#include "mytime.h"

#include <map>
#include <vector>


typedef std::map <abs_strip_t, std::vector <mytime_t> > times_map_t;

class Stopwatch
{
  //active strips are one every jump, beginning from a row
  int jump_ = 1;

  /* a hit happened between strip s and strip s+1 (s = absolute ID)
   * is stored in original_ with key s and value t */
  times_map_t original_;

  //active_ register times of active strips
  times_map_t active_;

  inline bool is_active(const abs_strip_t &strip)
  {
    if(strip % Nrows*Nstrips % jump_ == 0) return true;
    return false;
  }


public:

  //take time on strip
  inline void split
    (const int &ladder, const int &strip, const mytime_t &t)
  { original_[abs_strip(ladder,strip)].push_back(t); }

  //get hit times on (ladder,strip); if not a hit return void vector
  inline std::vector<mytime_t> time(const int &ladder, const int &strip)
  {
    if(active_.find(abs_strip(ladder,strip)) != active_.end())
      return active_[abs_strip(ladder,strip)];
    return {}; //void vector
  }

  inline void reset()
  { original_.clear(); active_.clear(); }

};


#endif //include guard
