
// Stopwatch handles GGS hit times

#ifndef STOPWATCH_INCLUDE
#define STOPWATCH_INCLUDE


#include "global.h"
#include "types.h"
#include "absStrip.h"

#include <map>
#include <vector>


class Stopwatch
{

protected:

  //associate time to strip
  typedef std::map <abs_strip_t, std::vector<mytime_t> > times_map_t;


  /* a hit happened between strip s and strip s+1 (s = absolute ID)
   * is stored in times_ with key s and value t */
  times_map_t times_;


public:

  //take time on strip
  inline void Split
    (const int &ladder, const int &strip, const mytime_t &t)
  { times_[absStrip(ladder,strip)].push_back(t); }

  //get hit times on (ladder,strip); if not a hit return void vector
  inline std::vector<mytime_t> GetTimes
    (const int &ladder, const int &strip)
  {
    if(times_.find(absStrip(ladder,strip)) != times_.end())
      return times_[absStrip(ladder,strip)];
    return {}; //void vector
  }

  inline void Reset()
  { times_.clear(); }

};


#endif //include guard
