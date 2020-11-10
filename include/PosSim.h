
#ifndef POS_SIM
#define POS_SIM


#include "physics.h"
#include "vector2.h"
#include "Geometry.h"

#include "TRandom3.h"
#include "TH1F.h"

#include <vector>

using namespace std;


class PosSim
{
  Geometry *geo_ = NULL;
  TRandom3 *random_ = NULL;
  vector2<double> *eDepSegm = NULL;
	vector2<double> *hitPos = NULL;

  int jump = 0;


  inline void SetVectors()
  {
    eDepSegm = new vector2<double>
      (
        geo_->GetNladders(), vector<double>(geo_->GetNstrips())
      );

  	hitPos = new vector2<double> (geo_->GetNlayers());
  }

  void GetCluster(int&, int&, int&, int&);

  void FillCluster(vector_pair<double>&, int, int, int, int);

  void FillHist(TH1F*, const vector_pair<double>&, const int);


public:

  PosSim(Geometry *geo, int j, TRandom3 *r)
  { geo_ = geo; SetVectors(); jump = j; random_ = r; }

  ~PosSim()
  { delete eDepSegm; delete hitPos; }

  inline void Reset()
  { delete eDepSegm; delete hitPos; SetVectors(); }

  inline void SetHitPos(const int &layer, const double &pos)
  { (*hitPos)[layer].push_back(pos); }

  inline void DepositEnergy
    (const int &ladder, const int &strip, const double &energy)
  { (*eDepSegm)[ladder][strip] += energy; }

  //share energy between active strips: one every jump
  void ShareEnergy();

  void AddNoise();

  void Segm(TH1F*);
};


#endif //include guard
