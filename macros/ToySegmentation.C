std::vector< std::pair<double, bool> > Segmentation(int readout_step=1, double implant_pitch=100) {

  double readout_pitch = implant_pitch*readout_step;

  std::vector< std::pair<double, bool> > strip_pos;

  strip_pos.push_back(std::make_pair(-readout_pitch, true));

  int ii=0;
  while (1) {
    ii++;
    std::pair<double, bool> _pair; 
    _pair.first=strip_pos[ii-1].first + implant_pitch;
    if (ii%readout_step) _pair.second=false;
    else _pair.second=true;
    strip_pos.push_back(_pair);
    if ((strip_pos[ii].first+implant_pitch)>readout_pitch) break;
  }

  // for (int ii=0; ii<(int)(strip_pos.size()); ii++) {
  //   printf("%d) %f (%d)\n", ii, strip_pos[ii].first, strip_pos[ii].second);
  // }
  
  return strip_pos;
}

std::vector< std::pair<double, bool> > ChargeSharing(double pos, std::vector< std::pair<double, bool> > strip_pos, double Ene=1.0, int kSharingType=0){
  // kSharingType = 0; linear
  // kSharingType = 1; only closer
  // kSharingType = 2; 50% - 50% between th two closest

  std::vector< std::pair<double, bool> > ene_dep(strip_pos.size());

  double Eleft = 0;
  double Eright = 0;

  for (int ii=1; ii<(int)(strip_pos.size()); ii++) {
    double left = strip_pos[ii-1].first;
    double right = strip_pos[ii].first;
    if (pos>left && pos<=right) {
      if (kSharingType==1) {
	if (fabs(pos-left)<fabs(pos-right)) {
	  Eleft = Ene;
	  Eright = 0;
	}
	else {
	  Eleft = 0;
	  Eright = Ene;
	}
	ene_dep[ii-1].first = Eleft;
	ene_dep[ii].first = Eright;
      }
      else if (kSharingType==2) {
	Eright = 0.5*Ene;
	Eleft = 0.5*Ene;
	ene_dep[ii-1].first = Eleft;
	ene_dep[ii].first = Eright;
      }
      else {
	//      if (kSharingType==0) {
	Eright = Ene*(pos-left)/(right-left);
	Eleft = Ene*(right-pos)/(right-left);
	ene_dep[ii-1].first = Eleft;
	ene_dep[ii].first = Eright;	
      }
    }
    ene_dep[ii-1].second = strip_pos[ii-1].second;
    ene_dep[ii].second = strip_pos[ii].second;
  }

  //  printf("pos=%f -> Eleft=%f, Eright=%f\n", pos, Eleft, Eright);
  
  return ene_dep;
}

std::vector< std::pair<double, bool> > ChargeCoupling(std::vector< std::pair<double, bool> > ene_dep){

  std::vector< std::pair<double, bool> > ene_readout(ene_dep);

  // 4% sharing, up to the 3rd neighbour...
  for (int ii=0; ii<(int)(ene_readout.size()); ii++) {
    double ene = ene_readout[ii].first;
    if (ii>0) {
      ene_dep[ii-1].first += 0.04*ene;
      ene_dep[ii].first   -= 0.04*ene;
    }
    if (ii<(int)(ene_readout.size())-1) {
      ene_dep[ii+1].first += 0.04*ene;
      ene_dep[ii].first   -= 0.04*ene;
    }
    if (ii>1) {
      ene_dep[ii-2].first += 0.04*0.04*ene;
      ene_dep[ii].first   -= 0.04*0.04*ene;
    }
    if (ii<(int)(ene_readout.size())-2) {
      ene_dep[ii+2].first += 0.04*0.04*ene;
      ene_dep[ii].first   -= 0.04*0.04*ene;
    }
    if (ii>2) {
      ene_dep[ii-3].first += 0.04*0.04*0.04*ene;
      ene_dep[ii].first   -= 0.04*0.04*0.04*ene;
    }
    if (ii<(int)(ene_readout.size())-3) {
      ene_dep[ii+3].first += 0.04*0.04*0.04*ene;
      ene_dep[ii].first   -= 0.04*0.04*0.04*ene;
    }
  }

  // "fluence" up to the readout... 
  for (int ll=0; ll<30; ll++) {
    ene_readout = ene_dep;
    for (int ii=0; ii<(int)(ene_readout.size()); ii++) {
      double ene = ene_readout[ii].first;
      if (!ene_readout[ii].second && ii!=0) {//non readout: it works if first and last are readout
	ene_dep[ii-1].first += 0.5*ene;
	ene_dep[ii+1].first += 0.5*ene;
	ene_dep[ii].first   -= ene;
      }
    }
  }
  
  // for (int ii=0; ii<(int)(ene_readout.size()); ii++) {
  //   printf("%d) %f (%d)\n", ii, ene_readout[ii].first, ene_readout[ii].second);
  // }
  
  return ene_readout;
}

std::vector<double> RemoveNotRead(std::vector< std::pair<double, bool> > vec){

  std::vector<double> vec_stripped;
  
  for (int ii=0; ii<(int)(vec.size()); ii++) {
    if ((bool)(vec[ii].second)) {
      vec_stripped.push_back(vec[ii].first);
    }
  }
  
  // for (int ii=0; ii<(int)(vec_stripped.size()); ii++) {
  //   printf("%d) %f\n", ii, vec_stripped[ii]);
  // }
  
  return vec_stripped;
}

std::vector<double> AddNoise(std::vector<double> ene_readout, double ene) {

  std::vector<double> ene_readout_withnoise(ene_readout);
  
  double noise = ene/10.0;
  //  printf("%f\n", noise);

  for (int ii=0; ii<(int)(ene_readout.size()); ii++) {
    double ene = gRandom->Gaus(ene_readout[ii], noise);
    //    printf("ene = %f\n", ene);
    if (ene<0) ene=0;
    ene_readout_withnoise[ii] = ene;
  }
  
  // for (int ii=0; ii<(int)(ene_readout.size()); ii++) {
  //   printf("%d) %f -> %f\n", ii, ene_readout[ii], ene_readout_withnoise[ii]);
  // }

  return ene_readout_withnoise;
}

double Baricenter(std::vector<double> strip_pos, std::vector<double> ene_readout){

  int seed = -999;
  int neigh = -999;

  double highest = -999;
  for (int ii=0; ii<(int)(ene_readout.size()); ii++) {
    if (ene_readout[ii]>highest) {
      highest = ene_readout[ii];
      seed = ii;
    }
  }

  double second = -999;
  for (int ii=0; ii<(int)(ene_readout.size()); ii++) {
    if (ene_readout[ii]>second && ii!=seed) {
      second = ene_readout[ii];
      neigh = ii;
    }
  }
  
  // printf("%f * %f\n", strip_pos[seed], ene_readout[seed]);
  // printf("%f * %f\n", strip_pos[neigh], ene_readout[neigh]);
  double reco_pos = strip_pos[seed]*ene_readout[seed] + strip_pos[neigh]*ene_readout[neigh];
  reco_pos/=(ene_readout[seed] + ene_readout[neigh]);
  
  //  printf("reco_pos = %f\n", reco_pos);
 
  return reco_pos;
}
  
void ToySegmentation() {
  
  printf("**** pitch: *****\n");
  double implant_pitch=50.0;
  double readout_step=3;
  printf("pitch = %f, %f\n", implant_pitch, readout_step*implant_pitch);
  printf("*****************\n");

  printf("**** ene: *****\n");
  double tot_ene=30.0;
  printf("ene = %f\n", tot_ene);
  printf("*****************\n");
  
  int nentries = 1000000;
  TH1F* h = new TH1F("h", "h", 1000, -readout_step*implant_pitch, readout_step*implant_pitch);

  for (int nn=0; nn<nentries; nn++) {
    
    if (nentries==1) printf("**** pos: *****\n");
    double pos=15.0;
    if (nentries>1) pos=gRandom->Uniform(-readout_step*implant_pitch/2.0, readout_step*implant_pitch/2.0);
    if (nentries==1) printf("pos = %f\n", pos);
    if (nentries==1) printf("*****************\n");
    
    if (nentries==1) printf("**** strip: *****\n");
    std::vector< std::pair<double, bool> > implant_strip_pos = Segmentation(readout_step, implant_pitch);
    for (int ii=0; ii<(int)(implant_strip_pos.size()); ii++) {
      if (nentries==1) printf("%d) %f (%d)\n", ii, implant_strip_pos[ii].first, implant_strip_pos[ii].second);
    }
    if (nentries==1) printf("*****************\n");

    if (nentries==1) printf("**** depo: *****\n");
    std::vector< std::pair<double, bool> > ene_dep = ChargeSharing(pos, implant_strip_pos, tot_ene);
    for (int ii=0; ii<(int)(ene_dep.size()); ii++) {
      if (nentries==1) printf("%d) %f (%d)\n", ii, ene_dep[ii].first, ene_dep[ii].second);
    }
    if (nentries==1) printf("*****************\n");

    if (nentries==1) printf("**** collected: *****\n");
    std::vector< std::pair<double, bool> > ene_collected = ChargeCoupling(ene_dep);
    for (int ii=0; ii<(int)(ene_collected.size()); ii++) {
      if (nentries==1) printf("%d) %f (%d) \n", ii, ene_collected[ii].first, ene_collected[ii].second);
    }
    if (nentries==1) printf("*****************\n");

    if (nentries==1) printf("**** remove not-read: *****\n");
    std::vector<double> strip_pos = RemoveNotRead(implant_strip_pos);
    std::vector<double> ene_readout = RemoveNotRead(ene_collected);
    for (int ii=0; ii<(int)(ene_readout.size()); ii++) {
      if (nentries==1) printf("%d) %f %f\n", ii, strip_pos[ii], ene_readout[ii]);
    }
    if (nentries==1) printf("*****************\n");

    if (nentries==1) printf("**** add noise: *****\n");
    std::vector<double> ene_readout_withnoise = AddNoise(ene_readout, tot_ene);
    for (int ii=0; ii<(int)(ene_readout_withnoise.size()); ii++) {
      if (nentries==1) printf("%d) %f\n", ii, ene_readout_withnoise[ii]);
    }
    if (nentries==1) printf("*****************\n");
    
    if (nentries==1) printf("**** pos reco: *****\n");
    double reco_pos=Baricenter(strip_pos, ene_readout_withnoise);
    if (nentries==1) printf("reco_pos = %f\n", reco_pos);
    if (nentries==1) printf("*****************\n");

    h->Fill(reco_pos-pos);
  }

  h->Draw();
  
  return;
}
