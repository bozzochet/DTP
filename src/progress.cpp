
#include "progress.h"

void progress(const int &n, const int &N) {

    double frac = (double) n / (double) N ;
    if(frac > 0.99) frac = 1; //fix ending with 99%

    int length = 40; //number of characters to use to print bar
    int completed = frac*length;

    string bar = "  [";

    for(int i=0; i < completed; ++i)
      bar += "=";

    if(completed < length)
      bar += ">";

    for(int j=completed+1; j < length; ++j)
      bar += "-";

    //percentage

    bar += "]  " + to_string((int)(frac*100)) + "%";

    cout <<"\r" <<bar;
}
