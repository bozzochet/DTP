
#ifndef PROGRESS_INCLUDE
#define PROGRESS_INCLUDE


#include <iostream>
#include <string>
#include <ctime>


//print a progress bar
void progress
(
  const std::clock_t &time_from_start,
  const int &partial, const int &tot
);


#endif
