
//
// Andrea Serpolla
//   a.serpolla@icloud.com
//

/*********************************************************************

  Time segmentation functions and variables;

  there are implemented three types (A,B,C) of segm:

    A - measures are given considering groups of N consecutive strips
      as one

    B - measures are given considering i-th strip, the (i+N)-th strip
      etc... as one

    C - measures are given considering Si pads placed on the backplate
      of the layer and oriented perpendicular to the strip direction


  ID definitions:

    - group A0 includes strips {0, ..., N} on ladder 0

    - group B0 includes strips {0, N, 2N, ...} on ladder 0

    - if layer 0 contains strips developing along y direction and
      occupies region [x,x+dx] * [y,y+dy], group C0 is placed in region
      [x,x+dx] * [y,y+side] (in this case dx = square side)
      where side is the parameter of C segm.

  parameter for A and B segm is jump: strips are grouped every N=jump

  in every segm, after filling one ladder, the numbering follows
  strips and/or ladder IDs

 ********************************************************************/

#ifndef TIME_SEGM_INCLUDE
#define TIME_SEGM_INCLUDE


#include "physics.h"
#include "Geometry.h"
#include "TrCluster.hh"

#include "TMath.h"
#include "TGraph.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>


enum time_segm : const char
{
  NONE = '_',
  A = 'A',
  B = 'B',
  C = 'C'
};


class TimeSegm
{

  class Group //group of strips depending on segm choosen
  {

    const char *name_ = "_NONE_";

    TGraph *time_energy_; //store points (hit time, energy deposited)
    TGraph *time_noise_; //store points (hit time, noise collected)

    bool sorted_ = false;


  public:

    Group(const char *name)
    {
      name_ = name;
      time_energy_ = new TGraph();
      time_noise_ = new TGraph();
    }

    ~Group()
    { delete time_energy_; delete time_noise_; }


    inline void Sort()
    {
      if( !sorted_)
      {
        time_energy_->Sort();
        time_noise_->Sort();
        sorted_ = true;
      }
    }


    inline void SetHit
      (const mytime_t &t, const energy_t &E, const energy_t &noise)
    {
      sorted_ = false;
      time_energy_->SetPoint(time_energy_->GetN(), t, E);
      time_noise_->SetPoint(time_noise_->GetN(), t, noise);
    }


    inline const char* GetName()
    { return name_; }


    inline void Reset()
    {
      delete time_energy_; time_energy_ = new TGraph();
      delete time_noise_; time_noise_ = new TGraph();
    }


    inline void GetTimes(std::vector<mytime_t> &v)
    {
      Sort();

      for(int i = 0; i < time_energy_->GetN(); ++i)
        v.push_back(time_energy_->GetX()[i]);
    }


    /* get time_energy_ points;
     * after execution, map passed contains points sorted */
    void GetHits
      (std::map <mytime_t, energy_t>&, std::map <mytime_t, energy_t>&);


  }; //class Group


// variables

  Geometry *geo_   = NULL;

  time_segm S_     = NONE;
  int jump_        = -1;           //A or B segm
  length_t side_   = -1;           //C segm

  int Ngroups_lad_ = -1;           //groups per ladder

  std::vector<Group*> group_;


//methods

  //get void space between C segm padds
  length_t GetVoidSpace();


public:

  //if S == A or B, only jump is read; otherwise only side
  TimeSegm
  (
    Geometry *geo, const time_segm S, const int jump,
    const double side = 0
  );

  ~TimeSegm()
  {
    for(int i = 0; i < (int) group_.size(); ++i)
      delete group_[i];
  }


  inline time_segm GetSegm()
  { return S_; }


  inline int GetJump()
  { return jump_; }


  inline double GetSide()
  { return side_; }


  inline int GetNgroups()
  { return group_.size(); }


  inline void GetHits
  (
    const int &i,
    std::map<mytime_t, energy_t> &E,
    std::map<mytime_t, energy_t> &noise
  )
  { group_[i]->GetHits(E, noise); }


  inline void GetTimes(const int &i, std::vector<mytime_t> &v)
  { group_[i]->GetTimes(v); }


  inline void Reset()
  {
    for(int i = 0; i < (int) group_.size(); ++i)
      group_[i]->Reset();
  }


  //individuates right group and call Group::SetHit; return group ID
  int SetHit(TrCluster *cl, const energy_t &noise);

};


#endif //include guard
