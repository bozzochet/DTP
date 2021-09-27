void ComputeEntryExitPos(double pos, double theta, double thickness, double& entrypos, double& exitpos);

std::vector< std::pair<double, bool> > Segmentation(int readout_step=1, double implant_pitch=100);
std::vector< std::pair<double, bool> > ChargeSharing(double pos, std::vector< std::pair<double, bool> > strip_pos, double Ene=86400.0, int kSharingType=3, double theta=0.0, double thickness=300.0);
std::vector< std::pair<double, bool> > ChargeCoupling(std::vector< std::pair<double, bool> > ene_dep, bool kFourPercentCoupling=true);
std::vector<double> RemoveNotRead(std::vector< std::pair<double, bool> > vec);
std::vector<double> AddNoise(std::vector<double> ene_readout);
double Baricenter(std::vector<double> strip_pos, std::vector<double> ene_readout);
  
void ToySegmentation() {
  
  printf("**** pitch: *****\n");
  // // POX
  // double implant_pitch=50.0;
  // double readout_step=3;
  // double thickness=150;
  // AMS-S
  double implant_pitch=27.5;
  double readout_step=4;
  double thickness=300;
  //
  printf("pitch = %f, %f\n", implant_pitch, readout_step*implant_pitch);
  printf("*****************\n");
  printf("thickness = %f\n", thickness);
  printf("*****************\n");
  
  printf("**** ene: *****\n");
  double vert_ene=thickness*80*3.6;//eV
  printf("ene = %f\n", vert_ene);
  printf("*****************\n");
  
  int nentries = 1000000;
  //  int nentries = 1;

  TH1F* hreso = new TH1F("hreso", "reso", 1000, -readout_step*implant_pitch, readout_step*implant_pitch);

  TH1F* hpos = new TH1F("hpos", "pos", 1000, -readout_step*implant_pitch, readout_step*implant_pitch);
  TH2F* hresovspos = new TH2F("hresovspos", "resovspos",  1000, -readout_step*implant_pitch, readout_step*implant_pitch, 1000, -readout_step*implant_pitch, readout_step*implant_pitch);
  
  TH1F* hentrypos = new TH1F("hentrypos", "entrypos", 1000, -2.0*readout_step*implant_pitch, 2.0*readout_step*implant_pitch);
  TH2F* hresovsentrypos = new TH2F("hresovsentrypos", "resovsentrypos",  1000, -2.0*readout_step*implant_pitch, 2.0*readout_step*implant_pitch, 1000, -readout_step*implant_pitch, readout_step*implant_pitch);
  
  TH1F* hexitpos = new TH1F("hexitpos", "exitpos", 1000, -2.0*readout_step*implant_pitch, 2.0*readout_step*implant_pitch);
  TH2F* hresovsexitpos = new TH2F("hresovsexitpos", "resovsexitpos",  1000, -2.0*readout_step*implant_pitch, 2.0*readout_step*implant_pitch, 1000, -readout_step*implant_pitch, readout_step*implant_pitch);

  TH1F* hangle = new TH1F("hangle", "angle", 1000, -360, 360);
  TH2F* hresovsangle = new TH2F("hresovsangle", "resovsangle", 1000, -360, 360, 1000, -readout_step*implant_pitch, readout_step*implant_pitch);
  
  TH1F* hene = new TH1F("hene", "ene", 1000, 4, 7);
  TH2F* hresovsene = new TH2F("hresovsene", "resovsene", 1000, 4, 7, 1000, -readout_step*implant_pitch, readout_step*implant_pitch);
  
  for (int nn=0; nn<nentries; nn++) {
    
    if (nentries==1) printf("**** pos: *****\n");
    double pos=15.0;
    //    double pos=100.0;
    if (nentries>1) pos=gRandom->Uniform(-readout_step*implant_pitch/2.0, readout_step*implant_pitch/2.0);
    if (nentries==1) {
      printf("pos = %f\n", pos);
      printf("*****************\n");
    }
    hpos->Fill(pos);
    
    if (nentries==1) printf("**** theta: *****\n");
    // generation features | theta is wrt positive z-axis
    const double PI= TMath::Pi();
    Double_t THETA1= PI/2.0;
    Double_t THETA2= PI;
    Double_t PHI1= 0;
    Double_t PHI2= 2.0*PI;
    double cossqtheta=1.0;
    if (nentries>1) cossqtheta=gRandom->Uniform(cos(THETA1)*cos(THETA1), cos(THETA2)*cos(THETA2));
    double phi = 0.0;
    phi = gRandom->Uniform(PHI1, PHI2);
    double theta = 0.0;
    theta = TMath::ACos(TMath::Sqrt(cossqtheta));
    double theta_x = TMath::ATan2(TMath::Sin(phi)*TMath::Sin(theta), TMath::Cos(theta));
    if (nentries==1) {
      printf("theta_x = %f\n", TMath::RadToDeg()*theta_x);
      printf("*****************\n");
    }
    hangle->Fill(TMath::RadToDeg()*theta_x);
    
    double entrypos = 0.0;
    double exitpos = 0.0;
    ComputeEntryExitPos(pos, theta_x, thickness, entrypos, exitpos);
    hentrypos->Fill(entrypos);
    hexitpos->Fill(exitpos);
    
    double tot_ene=vert_ene/fabs(TMath::Cos(theta));//here better to put theta (i.e. treat the deposit in 3D)!
    hene->Fill(log10(tot_ene));
    
    if (nentries==1) printf("**** strip: *****\n");
    std::vector< std::pair<double, bool> > implant_strip_pos = Segmentation(readout_step, implant_pitch);
    for (int ii=0; ii<(int)(implant_strip_pos.size()); ii++) {
      if (nentries==1) printf("%d) %f (%d)\n", ii, implant_strip_pos[ii].first, implant_strip_pos[ii].second);
    }
    if (nentries==1) printf("*****************\n");

    if (nentries==1) printf("**** depo: *****\n");
    //    std::vector< std::pair<double, bool> > ene_dep = ChargeSharing(pos, implant_strip_pos, tot_ene, theta_x, thickness);
    std::vector< std::pair<double, bool> > ene_dep = ChargeSharing(pos, implant_strip_pos, tot_ene, 3, theta_x, thickness);
    for (int ii=0; ii<(int)(ene_dep.size()); ii++) {
      if (nentries==1) printf("%d) %f (%d)\n", ii, ene_dep[ii].first, ene_dep[ii].second);
    }
    if (nentries==1) printf("*****************\n");

    if (nentries==1) printf("**** collected: *****\n");
    std::vector< std::pair<double, bool> > ene_collected = ChargeCoupling(ene_dep);
    //    std::vector< std::pair<double, bool> > ene_collected = ChargeCoupling(ene_dep, false);
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
    std::vector<double> ene_readout_withnoise = AddNoise(ene_readout);
    for (int ii=0; ii<(int)(ene_readout_withnoise.size()); ii++) {
      if (nentries==1) printf("%d) %f\n", ii, ene_readout_withnoise[ii]);
    }
    if (nentries==1) printf("*****************\n");
    
    if (nentries==1) printf("**** pos reco: *****\n");
    double reco_pos=Baricenter(strip_pos, ene_readout_withnoise);
    if (nentries==1) printf("reco_pos = %f\n", reco_pos);
    if (nentries==1) printf("*****************\n");

    hreso->Fill(reco_pos-pos);

    hresovsene->Fill(log10(tot_ene), reco_pos-pos);
    hresovsentrypos->Fill(entrypos, reco_pos-pos);
    hresovsexitpos->Fill(exitpos, reco_pos-pos);
    hresovsangle->Fill(TMath::RadToDeg()*theta_x, reco_pos-pos);
    hresovspos->Fill(pos, reco_pos-pos);
  }

  TCanvas* creso = new TCanvas("creso", "reso");
  hreso->Draw();

  TCanvas* cpos = new TCanvas("cpos", "pos");
  hpos->Draw();
  TCanvas* cresovspos = new TCanvas("cresovspos", "reso vs pos");
  hresovspos->Draw("colz");
  
  TCanvas* cangle = new TCanvas("cangle", "theta_x");
  hangle->Draw();
  TCanvas* cresovsangle = new TCanvas("cresovsangle", "reso vs theta_x");
  hresovsangle->Draw("colz");
  
  TCanvas* cene = new TCanvas("cene", "ene");
  hene->Draw();
  TCanvas* cresovsene = new TCanvas("cresovsene", "reso vs ene");
  hresovsene->Draw("colz");
  
  TCanvas* centrypos = new TCanvas("centrypos", "entrypos");
  hentrypos->Draw();
  TCanvas* cresovsentrypos = new TCanvas("cresovsentrypos", "reso vs entrypos");
  hresovsentrypos->Draw("colz");
  
  TCanvas* cexitpos = new TCanvas("cexitpos", "exitpos");
  hexitpos->Draw();
  TCanvas* cresovsexitpos = new TCanvas("cresovsexitpos", "reso vs exitpos");
  hresovsexitpos->Draw("colz");
  
  return;
}

std::vector< std::pair<double, bool> > Segmentation(int readout_step, double implant_pitch) {

  double howmanypitchs=3.0;
  
  double readout_pitch = implant_pitch*readout_step;

  std::vector< std::pair<double, bool> > strip_pos;

  strip_pos.push_back(std::make_pair(-howmanypitchs*readout_pitch, true));

  int ii=0;
  while (1) {
    ii++;
    std::pair<double, bool> _pair; 
    _pair.first=strip_pos[ii-1].first + implant_pitch;
    if (ii%readout_step) _pair.second=false;
    else _pair.second=true;
    strip_pos.push_back(_pair);
    if ((strip_pos[ii].first+implant_pitch)>howmanypitchs*readout_pitch) break;
  }

  // for (int ii=0; ii<(int)(strip_pos.size()); ii++) {
  //   printf("%d) %f (%d)\n", ii, strip_pos[ii].first, strip_pos[ii].second);
  // }
  
  return strip_pos;
}

void ComputeEntryExitPos(double pos, double theta, double thickness, double& entrypos, double& exitpos){

  entrypos = pos;
  exitpos = pos;
  
  double deltapos = 0.5*thickness*TMath::Sin(theta);
  entrypos-=deltapos;
  exitpos+=deltapos;
    
  //  printf("%f (%f) -> {%f, %f, %f}\n", TMath::RadToDeg()*theta-180, TMath::Sin(theta), entrypos, pos, exitpos);

  return;
}
  
std::vector< std::pair<double, bool> > ChargeSharing(double pos, std::vector< std::pair<double, bool> > strip_pos, double Ene, int kSharingType, double theta, double thickness){
  // kSharingType = 0; linear
  // kSharingType = 1; only closer
  // kSharingType = 2; 50% - 50% between th two closest
  // kSharingType = 3; split based on entry and exit position wrt to strip one
  
  std::vector< std::pair<double, bool> > ene_dep(strip_pos.size());

  if (kSharingType==0 || kSharingType==1 || kSharingType==2) {
    double Eleft = 0;
    double Eright = 0;
    
    for (int ii=1; ii<(int)(strip_pos.size()); ii++) {
      double left = strip_pos[ii-1].first;
      double right = strip_pos[ii].first;
      if (pos>left && pos<=right) {//in this cases this can happen just once
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
  }
  else {//kSharingType==3
    double entrypos = 0.0;
    double exitpos = 0.0;
    ComputeEntryExitPos(pos, theta, thickness, entrypos, exitpos);
        
    for (int ii=0; ii<(int)(strip_pos.size()); ii++) {
      ene_dep[ii].first = 0.0;
      ene_dep[ii].second = strip_pos[ii].second;
    }
    
    for (int ii=0; ii<(int)(strip_pos.size()); ii++) {
      //      printf("%d (%d)\n", ii, (int)(strip_pos.size()));
      
      double E = 0;
      double strip = strip_pos[ii].first;
      double pitch = 0.0;
      if (ii==((int)(strip_pos.size())-1)) {//last
	double last = strip_pos[ii-1].first;
	pitch = strip-last;
      }
      else {
	double next = strip_pos[ii+1].first;
	pitch = next-strip;
      }

      double path = 0.0;
      if (entrypos<exitpos) {//entering from top left toward bottom right
	if (exitpos<strip-0.5*pitch || entrypos>strip+0.5*pitch) continue;
	double leftborder = std::max(entrypos, strip-0.5*pitch);
	double rightborder = std::min(exitpos, strip+0.5*pitch);
	path = rightborder-leftborder;
      }
      else {//entering from top right toward bottom left
	if (exitpos>strip+0.5*pitch || entrypos<strip-0.5*pitch) continue;
	double leftborder = std::max(exitpos, strip-0.5*pitch);
	double rightborder = std::min(entrypos, strip+0.5*pitch);
	path = rightborder-leftborder;
      }

      E = Ene*(path)/(exitpos-entrypos);
      ene_dep[ii].first = E;
    }
  }
    
  return ene_dep;
}

std::vector< std::pair<double, bool> > ChargeCoupling(std::vector< std::pair<double, bool> > ene_dep, bool kFourPercentCoupling){

  std::vector< std::pair<double, bool> > ene_readout(ene_dep);

  if (kFourPercentCoupling) {
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
  }
    
  // "fluence" up to the readout... 
  for (int ll=0; ll<50; ll++) {
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

std::vector<double> AddNoise(std::vector<double> ene_readout) {

  std::vector<double> ene_readout_withnoise(ene_readout);
  
  double noise = 8600.0;//eV
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

