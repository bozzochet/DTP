

/* Andrea S. : this lib was created copying functions and code parts
 *    written in DataAnalysis.cpp, (not originally meant for beeing
 *    organized as here) to simulate position measures.
 *    I did not optimize this lib for the new organization of the
 *    project, leaving it where possible as the original version,
 *    modifying only the necessary to make it work with the rest of the
 *    project.
 *
 * The original behaviour is in legacy branch
 */

#ifndef POS_SIM
#define POS_SIM


#include "physics.h"
#include "vector2.h"
#include "Geometry.h"

//#include "TRandom3.h"
//#include "TH1F.h"

#include <vector>
#include <iostream>

using namespace std;


class PosSim
{
  Geometry *geo_ = NULL;
  //TRandom3 *random_ = NULL;
  vector2<double> eDepSegm;
	vector2<double> hitPos;

  int jump = 0;


  inline void ResizeVectors()
  {
    eDepSegm.resize(geo_->Nladders, vector<double>(geo_->Nstrips));
  	hitPos.resize(geo_->Nlayers, vector<double>(0));
  }

  void GetCluster(int&, int&, int&, int&);

  void FillCluster(vector_pair<double>&, int, int, int, int);

  double GetSimPos(const vector_pair<double>&, const int);


public:

  PosSim(Geometry *geo, int j /*, TRandom3 *r*/)
  { geo_ = geo; ResizeVectors(); jump = j; /*random_ = r;*/ }

  inline void Clear()
  { eDepSegm.clear(); hitPos.clear(); ResizeVectors(); }

  inline void SetHitPos(const int &layer, const double &pos)
  { hitPos[layer].push_back(pos); }

  inline void DepositEnergy
    (const int &ladder, const int &strip, const double &energy)
  { eDepSegm[ladder][strip] += energy; }

  //share energy between active strips: one every jump
  void ShareEnergy();

  double GetMeas();
};


#endif //include guard
