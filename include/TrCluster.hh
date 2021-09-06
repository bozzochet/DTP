#ifndef TR_CLUSTER_INCLUDE
#define TR_CLUSTER_INCLUDE

#include <iostream>
#include "TObject.h"

class TrCluster:public TObject{

public:
	int xy; //0 if x, 1 if y
	double pos[3]; //m
	double time; //s
	double eDep; //eV
	//double spRes;
	//double tRes;
	int parID;
	int parPdg;
	int layer;
	int ID;
	int strip;
	int ladder;
	double clust[2]; //eV


	float primIntPoint[3];
	int firstInteraction; //not primary particle == 0, primary particle == 1
	
	int isDead; //short



	void Dump() const {
	std::cout << "xy = " << xy << std::endl;
	std::cout << "pos = (" << pos[0] << " , " << pos[1] << "," << pos[2] << ")" << std::endl;
	std::cout << "time = " << time << std::endl;
	std::cout << "eDep = " << eDep << std::endl;
	//std::cout << "spRes = " << spRes << std::endl;
	//std::cout << "tRes = " << tRes << std::endl;
	std::cout << "parID = " << parID << std::endl;
	std::cout << "parPdg = " << parPdg << std::endl;
	std::cout << "layer = " << layer << std::endl;
	std::cout << "ID = " << ID << std::endl;
	std::cout << "strip = " << strip << std::endl;
	std::cout << "ladder = " << ladder << std::endl;
	std::cout << "cluster = (" << clust[0] << " , " << clust[1] << std::endl;
	}

  ClassDef(TrCluster, 2)
};

#endif