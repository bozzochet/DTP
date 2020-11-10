
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

    //compute groups on a single row
    Ngroups_ =

      S == A ?
        TMath::Ceil
          (
            (double) (geo_->GetNstrips() * geo_->GetNsquares())
            / (double) jump
          )
      :
        jump
    ;

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

    const int N = Ngroups_ * geo_->GetNrows() * geo_->GetNlayers();

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
    i = (lad % geo_->GetNsquares()) * Ngroups_ + strip / jump_;

  else if(S_ == B)
    i = (lad % geo_->GetNsquares()) * Ngroups_ + strip % jump_;

  //else if(S_ == C)
    //SEGM C

  else
  {
    std::cout <<"\n[TIMESEGM] fatal error on SetHit\n";
    exit(1);
  }

  group_[i]->SetHit(t,E);
}
