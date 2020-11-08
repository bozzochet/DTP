
#ifndef POS_SIMULATION
#define POS_SIMULATION


#include "globals_and_types.h"
#include "vector2.h"

#include "TRandom3.h"
#include "TH1F.h"

#include <vector>

using namespace std;


class PosSimulation
{

  TRandom3 *random_ = NULL;
  vector2<double> *eDepSegm = NULL;
	vector2<double> *hitPos = NULL;

  int jump = 0;


  inline void SetVectors()
  {
    eDepSegm = new vector2<double> (Nladders, vector<double>(Nstrips));
  	hitPos = new vector2<double> (Nlayers);
  }

  void GetCluster(int&, int&, int&, int&);

  void FillCluster(vector_pair<double>&, int, int, int, int);

  void FillHist(TH1F*, const vector_pair<double>&, const int);


public:

  PosSimulation(int j, TRandom3 *r)
  { SetVectors(); jump = j; random_ = r; }

  ~PosSimulation()
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
