
#include "TimeSimulation.h"


TimeSimulation::TimeSimulation()
{
  signal_up_ = new TF1("up", "[0]*x", 0, T_PEAK_);
  signal_up_->SetParameter(0, SLEW_RATE_ );

  //formula for signal_down_
  std::string formula = "TMath::Exp(-*x)";
  formula.insert(12, std::to_string(K_EXP_));

  double x_max = -1.0 / K_EXP_ * TMath::Log(ZERO_THRESH_);
  signal_down_ = new TF1("down", formula.c_str(), 0, x_max);
}


TimeSimulation::~TimeSimulation()
{
  delete signal_up_;
  delete signal_down_;
}


/*
void Stopwatch::add_noise(TH1D *hist)
{
    for(mytime_t t = T_MIN_; t < T_MAX_; t += BIN_LENGTH_)
      hist->Fill( NOISE );
}
*/


void TimeSimulation::add_signal
  (TH1D *hist, const std::vector<mytime_t> &times)
{
  for(int i = 0; i < (int) times.size(); ++i)
    for(
      mytime_t t = 0;
      t < times[i] + signal_down_->GetXmax();
      t += BIN_LENGTH_
    )

      if(t < T_PEAK_)
        hist->Fill(t+times[i], signal_up_->Eval(t));
      else
        hist->Fill(t+times[i], signal_down_->Eval(t-T_PEAK_) );
}


TH1D* TimeSimulation::get_signal
  (const int &lad, const int &s, Stopwatch *watch)
{
  //hist title
  std::string title = "ladder:  strip: ";
  title.insert(8, std::to_string(lad) );
  title.insert(title.length(), std::to_string(s) );

  TH1D *hist = new TH1D(
    "current", title.c_str(), N_BINS_, T_MIN_, T_MAX_
  );

  hist->GetXaxis()->SetTitle("time [ns]");
  hist->GetYaxis()->SetTitle("current / peak_current");

  //draw hist as a function instead of using bars
  hist->SetLineColor(kWhite);
  hist->SetMarkerColor(kBlue);
  hist->SetMarkerStyle(kFullDotMedium);

  //add_noise(hist);

  if(watch->time(lad,s).size() == 0)
  {
    //if no hit on strip ...
    std::cout <<"\nVECTOR VOID\n";
    return hist; // ... return hist with noise only
  }

  add_signal(hist, watch->time(lad,s));
  return hist;
}
