
#include "PosSimulation.h"


void PosSimulation::ShareEnergy()
{
  std::vector<double> fill;
  for (int ix = 0; ix < (*eDepSegm).size(); ix++) {
    for (int jx = 0; jx < (*eDepSegm)[0].size(); jx+=jump) {
      fill.clear();
		  fill.shrink_to_fit();
		  int i0 = ix;
		  int j0 = jx;

		  //Saving the energy from non-active strips

			for(int jp = 1; jp<jump; jp++) {
  			if(j0!=(*eDepSegm)[0].size()-1) {
  			  j0++;
		  	}

			  else {
				  i0++;
				  j0 = 0;
			  }

			  if(i0==(*eDepSegm).size())
				  break;

				if((*eDepSegm)[i0][j0]!=0) {
  				fill.push_back((*eDepSegm)[i0][j0]);
	  			(*eDepSegm)[i0][j0] = 0;
				}
		  }

		  //Moving to right-side strip

		  if(j0!=(*eDepSegm).size()-1) {
			  j0++;
			}
  		else if( (i0+1) % Nsquares != 0) {
				i0++;
				j0 = 0;
			}
  		else {
	  	/* no active strips between (ix,jx) and layer
		  * row end: distribute the entire energy to
		  * (ix,jx) strip
		  */
  		  i0 = ix;
	  	  j0 = jx;
		  }

      //Distributing the energy to left and right-side strips

		  for(int s = 0; s < fill.size(); s++) {
			  (*eDepSegm)[ix][jx] += fill[s]/(2*(s+1));
			  (*eDepSegm)[i0][j0] += fill[s]/(2*(fill.size()-s));
		  }

		}
  }

}


void PosSimulation::AddNoise()
{

	for (int ii = 0; ii < (*eDepSegm).size(); ii++) {
		for (int jj = 0; jj < (*eDepSegm)[ii].size(); jj++) {
			double fluct = random_->Gaus(0, 9e-6);

			(*eDepSegm)[ii][jj] += fluct;
			}
		}
}


void PosSimulation::Segm(TH1F *segmp)
{

  for (int ix = 0; ix < Nladders; ix++)
    for (int jx = 0; jx < Nstrips; jx+=jump) {


      //Find the boundaries of the clusters

      if((*eDepSegm)[ix][jx] < 27e-6)
        continue;

      //cout<<"Analysing cluster\n";
      vector_pair<double> strip;

      int i1, j1;
      bool firstPoint = true;
      for(int ii = ix; ii>=0; ii--) {
        for(int jj = Nstrips-1; jj>=jump; jj-= jump) {
          if(firstPoint) {
            jj = jx;
            firstPoint = false;
          }
          i1 = ii;
          j1 = jj;

          /* evaluate if the strip is the first of the first ladder
           * on the layer row
           */
          bool first_ladder = ii % Nsquares == 0 && jj == 0;

          if((*eDepSegm)[ii][jj] < 9e-6 || first_ladder)
            goto nextPart1;
        }
      }

      nextPart1:

      firstPoint = true;
      int i2, j2;
      for(int ii = ix; ii<Nladders; ii++) {
        for(int jj = 0; jj<Nstrips; jj+= jump) {
          if(firstPoint) {
            jj = jx;
            firstPoint = false;
          }
          i2 = ii;
          j2 = jj;

          /* evaluate if the strip is the first of the first ladder
           * on the layer row
           */
          bool first_ladder = ii % Nsquares == 0 && jj == 0;

          if((*eDepSegm)[ii][jj] < 9e-6 || first_ladder)
            goto nextPart2;
        }
      }

      nextPart2:

      //Filling vector with current cluster

      while((i1*Nstrips)+j1 <= (i2*Nstrips)+j2) {

        if(j1>Nstrips-1 && i1 == Nladders-1)
          break;
        else if(j1>Nstrips-1) {
          i1++;
          j1 = 0;
        }

        double thisPos = ((i1%Nsquares)*squareSide) + (j1*pitch) - (Nsquares*squareSide*0.5);
        strip.push_back(make_pair(thisPos,(*eDepSegm)[i1][j1]));

        j1+=jump;
      }

      for(int k = 0; k < strip.size(); k++) {

        if(strip[k].second < 27e-6)
          continue;

        //Finding the boundaries of the peak

        int k1 = k-1;
        int k2 = k+1;

        while(k1>0 && strip[k1].second >= 27e-6)
          k1--;
        while(k2<(strip.size()-1) && strip[k2].second >= 27e-6)
          k2++;

        //Finding the peak

        int kMax = k1;
        for(int kHold = k1; kHold <= k2; kHold++)
          if(strip[kHold].second>strip[kMax].second)
            kMax = kHold;


        //cout<<"Analysing peak of "<< strip[kMax].second <<"keV\n";

        //Finding neighbour strip

        double kNext;

        if(kMax == 0)
          kNext = 1;
        else if(kMax == strip.size()-1)
          kNext = strip.size()-2;
        else if(strip[kMax+1].second > strip[kMax-1].second)
          kNext = kMax+1;
        else
          kNext = kMax-1;


        //Finding simulated position

        double simPos = ((strip[kMax].first*strip[kMax].second) + (strip[kNext].first*strip[kNext].second)) / (strip[kMax].second + strip[kNext].second);
        int cLayer = ix/(Nsquares*Nrows);

        //Comparing the simulated hit positions with the real ones on the same layer

        for(int m = 0 ; m < (*hitPos)[cLayer].size(); m++)
          segmp->Fill(simPos-(*hitPos)[cLayer][m]);

        //cout<<"Simulated position: "<<simPos<<endl;

        //Advancing within the current cluster
        k = k2;
      }

      //Advancing and looking for other clusters
      ix = i2;
      jx = j2;
    }

}
