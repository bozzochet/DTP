
//types for time value

#ifndef MYTIME_INCLUDE
#define MYTIME_INCLUDE


#include "abs_strip.h"

#include <map>
#include <vector>


//time
typedef double mytime_t;

//time associated to strip
typedef std::map <abs_strip_t, std::vector<mytime_t> > times_map_t;


#endif
