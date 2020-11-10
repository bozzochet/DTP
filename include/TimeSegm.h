
/*********************************************************************

  time segmentation functions and variables;

  there are implemented three types (A,B,C) of segm:

    A - measures are given considering groups of N strips as one

    B - measures are given considering i-th strip, the (i+N)-th strip
      etc... as one.

    C - measures are given considering Si pads placed on the backplate
      of the layer and oriented perpendicular to the strip direction

 ********************************************************************/

#ifndef TIME_SEGM_INCLUDE
#define TIME_SEGM_INCLUDE


#include "physics.h"
#include "Geometry.h"

#include "TMath.h"
#include "TGraph.h"

#include <string>
#include <vector>
#include <map>


enum time_segm : const char
{
  NONE = '_',
  A = 'A',
  B = 'B',
  C = 'C'
};


class TimeSegm
{

// nested classes

  class Group //group of strips depending on segm choosen
  {
    const char *name_ = "_NONE_";

    TGraph *time_energy_; //store points (hit time, energy deposited)

    bool sorted_ = false;

  public:

    Group(const char *name)
    { name_ = name; time_energy_ = new TGraph(); }

    ~Group()
    { delete time_energy_; }

    inline void SetHit(const mytime_t &t, const energy_t &E)
    {
      sorted_ = false;
      time_energy_->SetPoint(time_energy_->GetN(), t, E);
    }

    inline const char* GetName()
    { return name_; }

    /* get time_energy_ points;
     * after execution, map passed contains points sorted */
    void GetHits(std::map <mytime_t, energy_t>&);

    void GetTimes(std::vector<mytime_t>&);

  };


// variables

  Geometry *geo_ = NULL;

  time_segm S_ = NONE;
  int jump_ = 0;
  double side_ = 0;

  int Ngroups_row_ = 0; //groups per row

  std::vector<Group*> group_;


//methods

public:

  TimeSegm
  (
    Geometry *geo, const time_segm S, const int jump,
    const double side = 0
  );

  ~TimeSegm();

  inline time_segm GetSegm()
  { return S_; }

  inline int GetJump()
  { return jump_; }

  inline double GetSide()
  { return side_; }

  inline int GetNgroups()
  { return group_.size(); }

  inline void GetHits(const int &i, std::map<mytime_t, energy_t> &m)
  { group_[i]->GetHits(m); }

  inline void GetTimes(const int &i, std::vector<mytime_t> &v)
  { group_[i]->GetTimes(v); }

  void SetHit
    (
      const int &lad, const int &strip,
      const mytime_t &t, const energy_t &E
    );

};


#endif //include guard
