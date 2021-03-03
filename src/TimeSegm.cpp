
#include "TimeSegm.h"


TimeSegm::TimeSegm
  (
    Geometry *geo,
    const time_segm S, const int jump, const length_t side
  )
{
  geo_ = geo;

  if( S_ == NONE && (S == A || S == B) && jump > 0 )
  {
    S_ = S;
    jump_ = jump;

    //compute groups on a single ladder
    Ngroups_lad_ = S == A ? geo_->Nstrips / jump : jump;

    if(geo_->Nstrips % jump != 0 && S == A)
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
  }

  else if(S == C && side > 0)
  {
    S_ = S;
    side_ = side;

    //compute groups on a single ladder
    Ngroups_lad_ = geo_->squareSide / side_ ;
  }

  else
  {
    std::cerr <<"[TIMESEGM] segm passed not acceptable: " <<(char)S
      <<" jump: " <<jump <<" side: " <<side <<std::endl;
    exit(1);
  }


  //create groups

  //number of total groups
  const int N = Ngroups_lad_ * geo_->Nladders;

  for(int i=0; i < N; ++i)
  {
    char c = S_;

    std::string name = std::to_string(i);
    name.insert(0, &c);

    group_.push_back( new Group(name.c_str()) );
  }
}



length_t TimeSegm::GetVoidSpace()
{
  return
  (
    //total ladder side along strip direction
    ((double) geo_->Nsquares) * geo_->squareSide
      / ((double) geo_->Nrows)

    //pads sides
    - Ngroups_lad_ * side_
  )
  / (Ngroups_lad_ - 1); //number of spaces between pads
}



int TimeSegm::SetHit(TrCluster *cl)
{
  int i = -1;

  if(S_ == A)
    i = cl->ladder * Ngroups_lad_ + cl->strip / jump_ ;

  else if(S_ == B)
    i = cl->ladder * Ngroups_lad_ + cl->strip % jump_;

  else if(S_ == C)
  {
    length_t side_void = side_ + GetVoidSpace();

    length_t pos; //WHICH UNIT USED FOR LENGTH BY TRCLUSTER ?

    if(cl->xy == 0)
      pos = cl->pos[1];
    else
      pos = cl->pos[0];

    //in GGS coordinates, 0 is in the middle of the layer
    pos += geo_->Nsquares * geo_->squareSide * 0.5;
    //now pos belongs to [0, Nsquares * squareSide]


    //possible bugs: check if pos is in [0, Nsquares * squareSide]
    if(pos > geo_->Nsquares * geo_->squareSide)
      pos = geo_->Nsquares * geo_->squareSide;


    i = TMath::FloorNint(pos / side_void);

    //hit in void space between i-th and (i+1)-th pads
    if(pos - side_void*i > side_)
    {
      //??????
    }
  }

  else //error
  {
    std::cerr <<"\n[TIMESEGM] fatal error on SetHit: ";
    std::cerr <<"S_ set badly\n";
    std::cerr <<"segm: " <<S_ <<std::endl;

    exit(1);
  }


  //error
  if( i >= (int) group_.size() || i < 0)
  {
    std::cerr <<"\n[TIMESEGM] fatal error on SetHit: ";
    std::cerr <<"call out of vector size\n";

    std::cerr <<"i: " <<i <<" segm: " << (char) S_ <<" size: "
      <<group_.size();

    std::cerr <<"\nlad: " <<cl->ladder <<" Nsquares: "
      <<geo_->Nsquares;

    std::cerr <<"\nNgroups per row: " <<Ngroups_lad_ <<" strip: "
      <<cl->strip;

    std::cerr <<" jump: " <<jump_ <<std::endl;

    exit(1);
  }


  group_[i]->SetHit(cl->time, cl->eDep);
}



void TimeSegm::Group::GetHits(std::map <mytime_t, energy_t> &m)
{
  Sort();

  for(int i = 0; i < time_energy_->GetN(); ++i)
  {
    mytime_t t;
    energy_t E;

    time_energy_->GetPoint(i, t, E);

    m[t] = E;
  }
}
