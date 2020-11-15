
#include "progress.h"

void progress(const std::clock_t &time, const int &n, const int &N)
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

    if(n+1 == N)
    {
      std::cout <<std::endl;
      return; //don't print extimated time
    }

    //extimated time

    double v = ((double) n) / ((double)(time / CLOCKS_PER_SEC));
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


    if(extimated > 0)
    {
      if(extimated < 10) std::cout <<" ";
      std::cout <<extimated <<"s ";
    }
}
