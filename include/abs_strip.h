
// handling absolute strip IDs

#ifndef ABS_STRIP_INCLUDE
#define ABS_STRIP_INCLUDE


#include "global.h"


//absolute strip ID
typedef int abs_strip_t;

//convert relative strip ID to absolute
inline abs_strip_t abs_strip(const int &ladder, const int &strip)
{ return ladder*Nstrips + strip; }


#endif
