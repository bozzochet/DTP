
void PrintGraph(const char *object, const int &E_e, const int &E_p)
{
  std::string p_file = "p/";
  p_file += std::to_string(E_p);
  p_file += "/histos.root";

  std::string e_file = "e/";
  e_file += std::to_string(E_e);
  e_file += "/histos.root";

  TH1F *h_p =
    (TH1F*) TFile::Open(p_file.c_str())->Get(object);

  TH1F *h_e =
    (TH1F*) TFile::Open(e_file.c_str())->Get(object);

  h_p->SetTitle("proton");
  h_e->SetTitle("electron");

  h_p->SetBit(TH1::kNoStats);
  h_e->SetBit(TH1::kNoStats);

  h_e->SetLineColor(kRed);

  std::string title = object;
  title += "_";
  title += std::to_string(E_e);
  title += "_";
  title += std::to_string(E_p);

  TCanvas *c = new TCanvas(title.c_str(), title.c_str(), 600, 400);

  c->SetLogy();

  h_e->DrawNormalized();
  h_p->DrawNormalized("same");
}


void ElectronProton(bool calo, int E_e, int E_p /*in GeV*/)
{
  gStyle->SetOptTitle(0);

  std::string label[4] = {"meas15", "meas15_slow", "calo", "calo_slow"};

  int max = calo ? 4 : 2;
  
  for(int i = 0; i < max; ++i)
  {
    std::string title = "h_time_";
    title += label[i];

    PrintGraph(title.c_str(), E_e, E_p);
  }
}
