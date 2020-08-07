std::vector<double> Segmentation(int readout_step=1, double implant_pitch=100) {

  double readout_pitch = implant_pitch*readout_step;

  std::vector<double> strip_pos(readout_step+1); 

  strip_pos[0] = -readout_pitch/2.0;
  
  for (int ii=1; ii<=readout_step; ii++) {
    strip_pos[ii] = strip_pos[ii-1] + implant_pitch; 
  }

  // for (int ii=0; ii<(int)(strip_pos.size()); ii++) {
  //   printf("%d) %f\n", ii, strip_pos[ii]);
  // }
  
  return strip_pos;
}

std::vector<double> ChargeDeposit(double pos, std::vector<double> strip_pos, double Ene=1.0){

  std::vector<double> ene_dep(strip_pos.size());

  double Eleft = 0;
  double Eright = 0;

  for (int ii=1; ii<(int)(strip_pos.size()); ii++) {
    double left = strip_pos[ii-1];
    double right = strip_pos[ii];
    if (pos>left && pos<=right) {
      Eright = Ene*(pos-left)/(right-left);
      Eleft = Ene*(right-pos)/(right-left);
      ene_dep[ii-1] = Eleft;
      ene_dep[ii] = Eright;
    }
  }

  //  printf("pos=%f -> Eleft=%f, Eright=%f\n", pos, Eleft, Eright);
  
  return ene_dep;
}

std::vector<double> ChargeSharing(std::vector<double> ene_dep){

  std::vector<double> ene_readout(ene_dep);

  for (int ll=0; ll<30; ll++) {
    ene_readout = ene_dep;
    
    for (int ii=1; ii<(int)(ene_readout.size())-1; ii++) {
      double ene = ene_readout[ii];
      ene_dep[ii-1] += 0.5*ene;
      ene_dep[ii+1] += 0.5*ene;
      ene_dep[ii]   -= ene;
    }
  }
  
  // for (int ii=0; ii<(int)(ene_readout.size()); ii++) {
  //   printf("%d) %f\n", ii, ene_readout[ii]);
  // }
  
  return ene_readout;
}

std::vector<double> RemoveNotRead(std::vector<double> vec, int readout_step){

  int n = (int)(((int)(vec.size())-1)/readout_step)+1;

  std::vector<double> vec_stripped(n);

  for (int ii=0; ii<(int)(vec_stripped.size()); ii++) {
    vec_stripped[ii] = vec[ii*readout_step];
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

  int n = (int)(strip_pos.size())-1;

  // printf("%f * %f\n", strip_pos[0], ene_readout[0]);
  // printf("%f * %f\n", strip_pos[n], ene_readout[n]);
  double reco_pos = strip_pos[0]*ene_readout[0] + strip_pos[n]*ene_readout[n];
  reco_pos/=(ene_readout[0] + ene_readout[n]);
  
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
  
  int nentries = 10000;
  TH1F* h = new TH1F("h", "h", 1000, -readout_step*implant_pitch, readout_step*implant_pitch);

  for (int nn=0; nn<nentries; nn++) {
    
    if (nentries==1) printf("**** pos: *****\n");
    double pos=15.0;
    if (nentries>1) pos=gRandom->Uniform(-readout_step*implant_pitch/2.0, readout_step*implant_pitch/2.0);
    if (nentries==1) printf("pos = %f\n", pos);
    if (nentries==1) printf("*****************\n");
    
    if (nentries==1) printf("**** strip: *****\n");
    std::vector<double> implant_strip_pos = Segmentation(readout_step, implant_pitch); 
    for (int ii=0; ii<(int)(implant_strip_pos.size()); ii++) {
      if (nentries==1) printf("%d) %f\n", ii, implant_strip_pos[ii]);
    }
    if (nentries==1) printf("*****************\n");
    
    if (nentries==1) printf("**** depo: *****\n");
    std::vector<double> ene_dep = ChargeDeposit(pos, implant_strip_pos, tot_ene);
    for (int ii=0; ii<(int)(ene_dep.size()); ii++) {
      if (nentries==1) printf("%d) %f\n", ii, ene_dep[ii]);
    }
    if (nentries==1) printf("*****************\n");
    
    if (nentries==1) printf("**** collected: *****\n");
    std::vector<double> ene_collected = ChargeSharing(ene_dep);
    for (int ii=0; ii<(int)(ene_collected.size()); ii++) {
      if (nentries==1) printf("%d) %f\n", ii, ene_collected[ii]);
    }
    if (nentries==1) printf("*****************\n");
    
    if (nentries==1) printf("**** remove not-read: *****\n");
    std::vector<double> strip_pos = RemoveNotRead(implant_strip_pos, readout_step);
    std::vector<double> ene_readout = RemoveNotRead(ene_collected, readout_step);
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
