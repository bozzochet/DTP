
#ifndef ABS_STRIP_INCLUDE
#define ABS_STRIP_INCLUDE


extern const int Nstrips; //defined in global.h

//absolute strip ID; look to abs_strip function for definition
typedef int abs_strip_t;

inline abs_strip_t abs_strip(const int &ladder, const int &strip)
{ return ladder*Nstrips + strip; }


#endif
