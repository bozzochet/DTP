
#include "info.h"

void info::progress
  (const std::clock_t &start, const int &n, const int &N)
{

    double frac = (double) n / (double) N ;
    if(n == N-1) frac = 1; //fix ending with 99%

    //don't update if same percentage //NOT WORKING
    if( (int)(frac*100) == ((n-1) / N * 100) )
      return;

    int length = 40; //number of characters to use to print bar
    int completed = frac*length;

    std::string bar = "  [";

    for(int i=0; i < completed; ++i)
      bar += "=";

    if(completed < length)
      bar += ">";

    for(int j=completed+1; j < length; ++j)
      bar += "-";

    //percentage

    bar += "]  " + std::to_string((int)(frac*100)) + "%";

    std::cout <<"\r" <<bar <<"  ";

    if(frac*100 < 100) std::cout <<" ";
    if(frac*100 < 10) std::cout <<" ";


    //extimated time

    double v =
      (double) n
      / (double) ( (std::clock() - start) / CLOCKS_PER_SEC );

    int extimated = (N-n) / v;

    if(extimated / 3600 > 0)
    {
      if(extimated / 3600 < 10) std::cout <<" ";
      std::cout <<extimated / 3600 <<"h ";
      extimated %= 3600;
    }
    else
      std::cout <<"    ";


    if(extimated / 60 > 0)
    {
      if(extimated / 60 < 10) std::cout <<" ";
      std::cout <<extimated / 60 <<"min ";
      extimated %= 60;
    }
    else
      std::cout <<"      ";


    if(extimated > 0 && n+1 != N)
    {
      if(extimated < 10) std::cout <<" ";
      std::cout <<extimated <<"s ";
    }
    else
      std::cout <<"   ";


    if(n+1 == N)
      std::cout <<std::endl;
}


void info::elapsed_time(const std::clock_t &start)
{
  int sec = (std::clock() - start) / CLOCKS_PER_SEC;

  if(sec != 0)
    std::cout <<"elapsed time: ";
  else
    std::cout <<std::endl;

  if(sec >= 3600)
  {
    std::cout <<sec/3600 <<"h ";
    sec %= 3600;
  }

  if(sec >= 60)
  {
    std::cout <<sec/60 <<"min ";
    sec %= 60;
  }

  if(sec > 0)
    std::cout <<sec <<"s " <<std::endl;
}
