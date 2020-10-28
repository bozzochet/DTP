
/**********************************************************************

  TimeSimulation class simulates the time development of the current
  signal generated in the strip by the passage of a particle.

  The signal simulated, line from 0 to peak and an exponential
  descent, based on T_PEAK_ and T_RELAX_ variables. The current is
  normalized to the maximum current, so the peak of a single signal
  has a value of 1; anyway, the sum of more than one signal could
  give a value > 1.

 *********************************************************************/

#ifndef TIME_SIMULATION_INCLUDE
#define TIME_SIMULATION_INCLUDE


#include "global.h"
#include "types.h"
#include "absStrip.h"
#include "Stopwatch.h"
#include "TrCluster.hh"

#include "TFile.h"
#include "TCanvas.h"
#include "TList.h"
#include "TObject.h"
#include "TGraph.h"
#include "TVectorD.h"

#include <vector>
#include <string>
#include <cstring>
#include <iostream>


class TimeSimulation : protected Stopwatch
{
  //current signal
  TVectorD *x_ = NULL;
  TVectorD *y_ = NULL;


public:

  //read weightfield2 current signal from its output root file
  TimeSimulation();

  ~TimeSimulation();

  //store hit times
  inline void SetHit(const TrCluster *cl)
  { Stopwatch::Split(cl->ladder, cl->strip, cl->time); }

  //Stopwatch::Reset is protected in this class
  virtual inline void Reset()
  { Stopwatch::Reset(); }

  //return graph with current signal on strip ( #ladder, #strip)
  //TGraph* GetSignal(const int &ladder, const int &strip);

};


#endif //include guard
