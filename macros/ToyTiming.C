void ToyNoiseVsTime(){

  TH1F* h_vs_t = new TH1F("h_vs_t", "h_vs_t; t; Events", 10000000, 0, 1000000);

  for (int ii=0; ii<10000; ii++) {
    h_vs_t->Fill(gRandom->Uniform(0, 1000000));
  }

  TCanvas* c_vs_t = new TCanvas();
  h_vs_t->Draw();

  TH1F* h_int_vs_t = new TH1F("h_int_vs_t", "h_int_vs_t; t; Events", 10000000, 0, 1000000);
  TH1F* h_vs_dt = new TH1F("h_vs_dt", "h_vs_dt; dt; Event", 1000, 0, 10000);
  TH1F* h_vs_df = new TH1F("h_vs_df", "h_vs_df; df; Event", 500, 0, 0.01);

  double entriessofar;
  double deltat=0;
  for (int ii=1; ii<=h_vs_t->GetNbinsX()+1; ii++) {
    double bc = h_vs_t->GetBinContent(ii);
    entriessofar+=bc;
    if (ii<h_vs_t->GetNbinsX()+1) h_int_vs_t->SetBinContent(ii, entriessofar);
    if (bc>0) {
      h_vs_dt->Fill(deltat+1);
      h_vs_df->Fill(1.0/(deltat+1));
      deltat=0;
    }
    else {
      deltat++;
    }
  }

  TCanvas* c_int_vs_t = new TCanvas();
  c_int_vs_t->SetLogy();
  h_int_vs_t->Draw();
  
  TCanvas* c_vs_dt = new TCanvas();
  c_vs_dt->SetLogy();
  h_vs_dt->Draw();
  h_vs_dt->Fit("expo", "", "", 0, 10000);

  TCanvas* c_vs_df = new TCanvas();
  c_vs_df->SetLogy();
  h_vs_df->Draw();
  h_vs_df->Fit("expo", "", "", 0.0005, 0.1);

  //#define _FFT

#ifdef _FFT
  //Compute the transform and look at the magnitude of the output
  TH1* h_vs_f = NULL;
  TVirtualFFT::SetTransform(0);
  //  h_vs_t->Rebin(1000000);
  h_vs_f = h_vs_t->FFT(h_vs_f, "MAG R2C M");
  
  h_vs_f->SetBit(TH1::kNoStats);
  h_vs_f->SetBit(TH1::kNoTitle);
#else
  TH1F* h_vs_f = new TH1F("h_vs_f", "h_vs_f; f; Event", 100000, 1.0E-6, 0.001);

  for (int ii=1; ii<=h_int_vs_t->GetNbinsX(); ii++) {
    double bc = h_int_vs_t->GetBinContent(ii);
    double xx = h_int_vs_t->GetBinCenter(ii);
    double f = 1.0/xx;
    h_vs_f->Fill(f, bc);
  }  
#endif
    
  TCanvas* c_vs_f = new TCanvas();
  c_vs_f->SetLogx();
  c_vs_f->SetLogy();
  h_vs_f->Draw("hist");
  TF1* f = new TF1("f", "[0]*TMath::Power(x, [1])", 1.0E-6, 0.001);
  f->SetNpx(10000);
  f->SetLineColor(kRed+2);
  f->SetLineWidth(3);
  f->SetParameter(0, 1.0E-5);
  f->SetParameter(1, -2.0);
  h_vs_f->Fit(f, "", "", 1.0E-6, 1.0E-4);
  f->Draw("same");
  // h_vs_f->Draw("same");
  
}

void ToyTiming() {

  ToyNoiseVsTime();

  return;
}
