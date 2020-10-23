
#include "Stopwatch.h"


void Stopwatch::stop()
{

  for(auto it = original_.begin(); it != original_.end(); ++it)
  {
    if(is_active(it->first))
      active_[it->first] = original_[it->first];
    else
    {
      //get index of active strips before and after strip it->first
      abs_strip_t active_before =
        it->first - it->first % (Nsquares*Nstrips) % jump_ ;

      abs_strip_t active_after = active_before + jump_ ;

      active_[active_before] = original_[it->first];
      active_[active_after] = original_[it->first];
    }
  }

}
