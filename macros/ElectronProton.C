
void ElectronProton()
{
  int E[3] = {10, 100, 1000}; //beam energies in GeV

  gStyle->SetOptTitle(0);

  std::string label[2] = {"meas15","slow"};

  for(int i=0; i<2; ++i)
    for(int j=0; j<3; ++j)
    {
      std::string p_file = "histos_proton_";
      p_file += std::to_string(E[j]);
      p_file += ".root";

      std::string e_file = "histos_electron_";
      e_file += std::to_string(E[j]);
      e_file += ".root";

      std::string object = "h_time_";
      object += label[i];

      TH1F *h_p = (TH1F*) TFile::Open(p_file.c_str())->Get(object.c_str());
      TH1F *h_e = (TH1F*) TFile::Open(e_file.c_str())->Get(object.c_str());

      h_p->SetTitle("proton");
      h_e->SetTitle("electron");

      h_p->SetBit(TH1::kNoStats);
      h_e->SetBit(TH1::kNoStats);

      h_e->SetLineColor(kRed);

      std::string title = label[i];
      title += "_";
      title += std::to_string(E[j]);

      TCanvas *c = new TCanvas(title.c_str(), title.c_str(), 600, 400);

      c->SetLogy();

      h_p->Draw();

      if(j!=0) //at 10 GeV same graphs for electrons and protons
        h_e->Draw("same");

    } //for i,j
}
