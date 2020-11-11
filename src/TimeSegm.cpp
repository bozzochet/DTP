
#include "TimeSegm.h"

TimeSegm::TimeSegm
  (
    Geometry *geo,
    const time_segm S, const int jump, const double side
  )
{
  geo_ = geo;

  if( S_ == NONE && (S == A || S == B) && jump > 0 )
  {
    S_ = S;
    jump_ = jump;

    //compute groups on a single ladder
    Ngroups_lad_ = S == A ? geo_->GetNstrips() / jump : jump;

    if(geo_->GetNstrips() % jump != 0 && S == A)
      ++Ngroups_lad_; //add one group for module division


    /* if Nstrips * Nsquares is not a multiple of jump:
     *
     *  - A segm:
     *      last group contains a number of strips < jump;
     *
     *  - B segm:
     *      B groups from 0 to (Nstrips * Nsquares) % jump
     *      contain
     *
     *            ceil( (Nstrips * Nsquares) / jump )
     *
     *      strips while the others contain
     *
     *            floor( (Nstirps * Nsquares) / jump )
     */

    //number of total groups
    const int N = Ngroups_lad_ * geo_->GetNladders();

    for(int i=0; i < N; ++i)
    {
      char c = S_;

      std::string name = std::to_string(i);
      name.insert(0, &c);

      group_.push_back( new Group(name.c_str()) );
    }
  }

  else if(S == C && side > 0)
  {
    S_ = S;
    side_ = side;

    //SEGM C
  }
}


TimeSegm::~TimeSegm()
{
  for(int i = 0; i < (int) group_.size(); ++i)
    delete group_[i];
}


void TimeSegm::SetHit
  (
    const int &lad, const int &strip,
    const mytime_t &t, const energy_t &E
  )
{
  int i = -1;

  if(S_ == A)
    i = (lad - 1) * Ngroups_lad_ + strip / jump_ ;

  else if(S_ == B)
    i = (lad - 1) * Ngroups_lad_ + strip % jump_;

  //else if(S_ == C)
    //SEGM C

  else //error
  {
    std::cout <<"\n[TIMESEGM] fatal error on SetHit: ";
    std::cout <<"S_ set badly\n";
    std::cout <<"segm: " <<S_ <<std::endl;

    exit(1);
  }


  //error
  if( i >= (int) group_.size() || i < 0)
  {
    std::cout <<"\n[TIMESEGM] fatal error on SetHit: ";
    std::cout <<"call out of vector size\n";

    std::cout <<"i: " <<i <<" segm: " << (char) S_ <<" size: "
      <<group_.size();
    std::cout <<"\nlad: " <<lad <<" Nsquares: " <<geo_->GetNsquares();
    std::cout <<"\nNgroups per row: " <<Ngroups_lad_ <<" strip: "
      <<strip;
    std::cout <<" jump: " <<jump_ <<std::endl;

    exit(1);
  }


  group_[i]->SetHit(t,E);
}


void TimeSegm::Group::GetHits(std::map <mytime_t, energy_t> &m)
{
  if( !sorted_)
  {
    time_energy_->Sort();
    sorted_ = true;
  }

  for(int i = 0; i < time_energy_->GetN(); ++i)
  {
    mytime_t t;
    energy_t E;

    time_energy_->GetPoint(i, t, E);

    m[t] = E;
  }
}


void TimeSegm::Group::GetTimes(std::vector<mytime_t> &v)
{
  if( !sorted_)
  {
    time_energy_->Sort();
    sorted_ = true;
  }

  for(int i = 0; i < time_energy_->GetN(); ++i)
  {
    mytime_t t;
    energy_t E;

    time_energy_->GetPoint(i, t, E);

    v.push_back(t);
  }
}
