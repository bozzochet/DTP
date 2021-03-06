#include "PosSim.h"


void PosSim::ShareEnergy()
{
  std::vector<double> fill;
  for (int ix = 0; ix < (int)((*eDepSegm).size()); ix++) {
    for (int jx = 0; jx < (int)((*eDepSegm)[0].size()); jx+=jump) {
      fill.clear();
      fill.shrink_to_fit();
      int i0 = ix;
      int j0 = jx;

      //Saving the energy from non-active strips

      for(int jp = 1; jp<jump; jp++) {
	if(j0!=(int)((*eDepSegm)[0].size())-1) {
	  j0++;
	}

	else {
	  i0++;
	  j0 = 0;
	}

	if(i0==(int)((*eDepSegm).size()))
	  break;

	if((*eDepSegm)[i0][j0]!=0) {
	  fill.push_back((*eDepSegm)[i0][j0]);
	  (*eDepSegm)[i0][j0] = 0;
	}
      }

      //Moving to right-side strip

      if(j0!=(int)((*eDepSegm).size())-1) {
	j0++;
      }
      else if( (i0+1) % geo_->Nsquares != 0) {
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

      for(int s = 0; s < (int)(fill.size()); s++) {
	(*eDepSegm)[ix][jx] += fill[s]/(2*(s+1));
	(*eDepSegm)[i0][j0] += fill[s]/(2*(fill.size()-s));
      }

    }
  }

}


/*
  void PosSim::AddNoise()
  {

  for (int ii = 0; ii < (*eDepSegm).size(); ii++) {
  for (int jj = 0; jj < (*eDepSegm)[ii].size(); jj++) {
  double fluct = random_->Gaus(0, 9e+3); //9kev

  (*eDepSegm)[ii][jj] += fluct;
  }
  }
  }
*/


void PosSim::GetCluster(int &i1, int &j1, int &i2, int &j2)
{
  bool firstPoint = true;

  int row = i1 / geo_->Nsquares; //get number of row which ladder belongs to

  for(int ii = i1; ii / geo_->Nsquares == row && ii >= 0; ii--)
    for(int jj = geo_->Nstrips-1; jj>=jump; jj-= jump) {

      if(firstPoint) {
        jj = j1;
        firstPoint = false;
      }

      if((*eDepSegm)[ii][jj] < 9e+3)
        break;

      i1 = ii;
      j1 = jj;
    }

  firstPoint = true;

  for(int ii = i2; ii / geo_->Nsquares == row && ii >= 0; ii++)
    for(int jj = 0; jj<geo_->Nstrips; jj+= jump) {

      if(firstPoint) {
        jj = j2;
        firstPoint = false;
      }

      if((*eDepSegm)[ii][jj] < 9e+3)
        break;

      i2 = ii;
      j2 = jj;
    }
}


void PosSim::FillCluster
(vector_pair<double> &strip, int i1, int j1, int i2, int j2)
{
  while((i1*geo_->Nstrips)+j1 <= (i2*geo_->Nstrips)+j2) {

    if(j1>geo_->Nstrips-1 && i1 == geo_->Nladders-1)
      break;
    else if(j1>geo_->Nstrips-1) {
      i1++;
      j1 = 0;
    }

    double thisPos = ((i1%geo_->Nsquares)*geo_->squareSide) + (j1*geo_->pitch) - (geo_->Nsquares*geo_->squareSide*0.5);
    strip.push_back(make_pair(thisPos,(*eDepSegm)[i1][j1]));

    j1+=jump;
  }
}


double PosSim::GetSimPos
(const vector_pair<double> &strip, const int layer)
{
  for(int k = 0; k < (int)(strip.size()); k++) {

    if(strip[k].second < 27e+3)
      continue;

    //Finding the boundaries of the peak

    int k1 = k-1;
    int k2 = k+1;

    while(k1>0 && strip[k1].second >= 27e+3)
      k1--;
    while(k2<(int)((strip.size())-1) && strip[k2].second >= 27e+3)
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
    else if(kMax == (int)(strip.size())-1)
      kNext = strip.size()-2;
    else if(strip[kMax+1].second > strip[kMax-1].second)
      kNext = kMax+1;
    else
      kNext = kMax-1;


    //Finding simulated position

    /* Andrea S.: in Digitization this lib is used one hit at a time.
     * Because of this I expect this method will find just one simPos.
     * Below I return for this reason. */

    double simPos = ((strip[kMax].first*strip[kMax].second) + (strip[kNext].first*strip[kNext].second)) / (strip[kMax].second + strip[kNext].second);

    return simPos;

    //Comparing the simulated hit positions with the real ones on the same layer
    /*
      for(int m = 0 ; m < (*hitPos)[layer].size(); m++)
      segmp->Fill(simPos-(*hitPos)[layer][m]);
    */
    //cout<<"Simulated position: "<<simPos<<endl;
    
    //Advancing within the current cluster
    k = k2;
    // FIX ME
  }
}


double PosSim::GetMeas()
{
  for (int ix = 0; ix < geo_->Nladders; ix++)
    for (int jx = 0; jx < geo_->Nstrips; jx+=jump) {

      //Find the boundaries of the clusters

      if((*eDepSegm)[ix][jx] < 27e+3)
        continue;

      //cout<<"Analysing cluster\n";

      int i1 = ix;
      int i2 = ix;

      int j1 = jx;
      int j2 = jx;

      GetCluster(i1,j1,i2,j2);

      //Filling vector with current cluster

      vector_pair<double> strip;
      FillCluster(strip, i1, j1, i2, j2);

      return GetSimPos(strip, ix/(geo_->Nsquares*geo_->Nrows));

      //Advancing and looking for other clusters
      ix = i2;
      jx = j2;
      // FIX ME, se return of GetSimPos
    }
}
