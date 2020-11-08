
#ifndef TR_CLUSTER_INCLUDE
#define TR_CLUSTER_INCLUDE

#include <iostream>
#include "TObject.h"

class TrCluster:public TObject{

public:
	int segm;
	double pos[3];
	double time;
	double eDep;
	double spRes;
	double tRes;
	int parID;
	int parPdg;
	int layer;
	int ID;
	int strip;
	int ladder;
	double clust[2];

ClassDef(TrCluster, 2)

	void Dump(){
	std::cout << "segm = " << segm << std::endl;
	std::cout << "pos = (" << pos[0] << " , " << pos[1] << "," << pos[2] << ")" << std::endl;
	std::cout << "time = " << time << std::endl;
	std::cout << "eDep = " << eDep << std::endl;
	std::cout << "spRes = " << spRes << std::endl;
	std::cout << "tRes = " << tRes << std::endl;
	std::cout << "parID = " << parID << std::endl;
	std::cout << "parPdg = " << parPdg << std::endl;
	std::cout << "layer = " << layer << std::endl;
	std::cout << "ID = " << ID << std::endl;
	std::cout << "strip = " << strip << std::endl;
	std::cout << "ladder = " << ladder << std::endl;
	std::cout << "cluster = (" << clust[0] << " , " << clust[1] << std::endl;
	}
};

#endif
