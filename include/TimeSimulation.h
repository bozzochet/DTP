
/**********************************************************************

  TimeSimulation class simulates signal of current in time
  generated in the strips by particles hits.

  The signal simulated, line from 0 to peak and an exponential
  descent, based on collected charge.

  Charge collecting process follows the capacitor charging law and
  noise is added to charge signal.

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
#include "TAxis.h"
#include "TRandom3.h"
#include "TH1F.h"

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

  //slew rate of line_
  const double SLEW_RATE_ = 1e+4;

  //tau of strip as a capacitor
  const mytime_t T_CAPACITOR_ = 1e-9;

  //signal hist parameters
  const mytime_t T_START_ = 0;
  const mytime_t T_END_ = 5e-9;
  const int N_BINS_ = 1000;

  //hist sampling time
  const mytime_t T_SAMPLING_ = (T_END_ - T_START_) / (double) N_BINS_;


// variables

  TF1 *up_ = NULL; //slope of current signal to the peak
  TF1 *charge_ = NULL; //charge collected in time without noise

  //energy deposit by particles
  energies_map_t energies_;

  //random generator
  TRandom3 *random_ = NULL;

  //hist for analysis
  TH1F *deviation_ = NULL;


// methods

  //generate charge signal in time with noise
  TGraph* GetChargeSignal(const energy_t&);

  /* add to signal a single hit current generated using up_ and
   * charge signal obtained by GetChargeSignal */
  void AddSignal
    (TGraph *signal, const TGraph *charge, const mytime_t &hitTime);


public:

  TimeSimulation();
  virtual ~TimeSimulation();

  //store hit times; TrCluster store times in ns
  inline void SetHit(const TrCluster *cl)
  {
    Stopwatch::Split(cl->ladder, cl->strip, cl->time * 1e-9);
    energies_[absStrip(cl->ladder, cl->strip)].push_back(cl->eDep);
  }

  //Stopwatch::Reset is protected in this class
  virtual inline void Reset()
  { Stopwatch::Reset(); }

  /* get current signal on strip ( #ladder, #strip) in vec begin
   * and collected charge graphs for every hit*/
  void GetSignal
    (std::vector<TGraph*> &vec, const int &ladder, const int &strip);

  /*
   * NEED IMPLEMENTATION: in case of multiple hits on the same strip
   * GetSignal is able to generate the multiple signals and append them
   * to the total signal returnes; at this moment the signals are NOT
   * summed (as has to be). Multiple currents are just drawn on the
   * same graph
   *
   */

};


#endif //include guard
