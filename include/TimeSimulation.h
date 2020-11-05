
/**********************************************************************

  TimeSimulation class simulates the time development of the current
  signal generated in the strip by the passage of a particle.

  The signal simulated, line from 0 to peak and an exponential
  descent, based on collected charge.

 *********************************************************************/

#ifndef TIME_SIMULATION_INCLUDE
#define TIME_SIMULATION_INCLUDE


#include "globals_and_types.h"
#include "absStrip.h"
#include "Stopwatch.h"
#include "TrCluster.hh"
#include "vector2.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TList.h"
#include "TObject.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1F.h"
#include "TAxis.h"

#include <vector>
#include <string>
#include <cstring>
#include <iostream>


class TimeSimulation : protected Stopwatch
{

// typedefs

  typedef std::map <abs_strip_t, std::vector<energy_t> >
    energies_map_t;


// const variables

  const double SLEW_RATE_ = 1e+9; //slew rate of line_
  const mytime_t T_CAPACITOR_ = 1e-9; //tau of strip as a capacitor
  const mytime_t T_START_ = 0;
  const mytime_t T_END_ = 5e-9;
  const int N_BINS_ = 1000; //signal hist bins
  const mytime_t T_SAMPLING_ = T_END_ / (double) N_BINS_; //hist sampling time


// variables

  //current signal components
  TF1 *up_ = NULL;
  TGraph *down_ = NULL;
  TF1 *charge_ = NULL; //charge collected in time

  //energy deposit
  energies_map_t energies_;


// methods

  //fit exponential to weightfield2 output signal
  TF1* SetExp(const char*);

  //add signal simulated with line_ and exp_ , based on on energy_
  void AddSignal(TH1F*, const energy_t&, const mytime_t&);


public:

  //set line_ and exp_
  TimeSimulation();

  ~TimeSimulation();

  //store hit times; TrCluster store times in ns
  inline void SetHit(const TrCluster *cl)
  {
    Stopwatch::Split(cl->ladder, cl->strip, cl->time * 1e-9);
    energies_[absStrip(cl->ladder, cl->strip)].push_back(cl->eDep);
  }

  //Stopwatch::Reset is protected in this class
  virtual inline void Reset()
  { Stopwatch::Reset(); }

  //return graph with current signal on strip ( #ladder, #strip)
  TH1F* GetSignal(const int &ladder, const int &strip);

};


#endif //include guard
