
#include "Stopwatch.h"


Stopwatch::Stopwatch(const int &jump)
{
  jump_ = jump;

  //set TF1 functions to describe simulated signal

  signal_up_ = new TF1("up", "[0]*x", 0, T_PEAK_);
  signal_up_->SetParameter(0, SLEW_RATE_UP_ );

  //formula for signal_down_
  std::string formula = "TMath::Exp(-*x)";
  formula.insert(12, std::to_string(K_EXP_));

  double x_max = -1.0 / K_EXP_ * TMath::Log(ZERO_THRESH_);
  signal_down_ = new TF1("down", formula.c_str(), 0, x_max);
}


Stopwatch::~Stopwatch()
{
  delete signal_up_;
  delete signal_down_;
}


void Stopwatch::stop()
{

  for(auto it = original_.begin(); it != original_.end(); ++it)
  {
    if(is_active(it->first))
      active_[it->first] = original_[it->first];
    else
    {
      //get index of active strips before and after strip it->first
      abs_strip_t active_before =
        it->first - it->first % (Nsquares*Nstrips) % jump_ ;

      abs_strip_t active_after = active_before + jump_ ;

      active_[active_before] = original_[it->first];
      active_[active_after] = original_[it->first];
    }
  }

}

/*
void Stopwatch::add_noise(TH1D *hist)
{
    for(mytime_t t = T_MIN_; t < T_MAX_; t += BIN_LENGTH_)
      hist->Fill( NOISE );
}
*/

void Stopwatch::add_signal
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


TH1D* Stopwatch::get_signal(const int &lad, const int &s)
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

  abs_strip_t strip = abs_strip(lad,s);

  if(original_.find(strip) == original_.end())
    //if no signal at strip passed ...
    return hist; // ... return hist with noise only

  add_signal(hist, original_[strip]);
  return hist;
}
