
/*********************************************************************
  info namespace contains function progress to print a progress bar
  with extimated time for the process in execution and elapsed_time
  function to print total time took by process to complete the task
*********************************************************************/

#ifndef INFO_INCLUDE
#define INFO_INCLUDE


#include <iostream>
#include <string>
#include <ctime>


namespace info
{

  //print a progress bar
  void progress
    (const std::clock_t &start, const int &partial, const int &tot);

  //print total elapsed time
  void elapsed_time(const std::clock_t &start);

};


#endif
