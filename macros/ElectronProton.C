
void ElectronProton(int max)
{
  int E[4] = {10, 100, 1000, 10000}; //beam energies in GeV

  gStyle->SetOptTitle(0);

  std::string label[2] = {"meas15","slow"};

  for(int i=0; i<2; ++i)
    for(int j=0; j<max; ++j)
    {
      std::string p_file = "p/";
      p_file += std::to_string(E[j]);
      p_file += "/histos.root";

      std::string e_file = "e/";
      e_file += std::to_string(E[j]);
      e_file += "/histos.root";

      std::string object = "h_time_";
      object += label[i];

      TH1F *h_p =
        (TH1F*) TFile::Open(p_file.c_str())->Get(object.c_str());

      TH1F *h_e =
        (TH1F*) TFile::Open(e_file.c_str())->Get(object.c_str());

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
      h_e->Draw("same");

    } //for i,j
}
