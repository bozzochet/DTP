class TrCluster:public TObject{
public:
	double segm;
	double pos[2];
	double time;
	double eDep;
	double spRes;
	double tRes;
	int parID;
	int layer;

ClassDef(TrCluster, 1)
};
