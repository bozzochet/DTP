
#include "segm.h"

//active strips are one every jump
void one_every_N::shareEnergy(vector2<double> &array, const int &jump)
{
  std::vector<double> fill;
  for (int ix = 0; ix < array.size(); ix++) {
    for (int jx = 0; jx < array[0].size(); jx+=jump) {
      fill.clear();
		  fill.shrink_to_fit();
		  int i0 = ix;
		  int j0 = jx;

		  //Saving the energy from non-active strips

			for(int jp = 1; jp<jump; jp++) {
  			if(j0!=array[0].size()-1) {
  			  j0++;
		  	}

			  else {
				  i0++;
				  j0 = 0;
			  }

			  if(i0==array.size())
				  break;

				if(array[i0][j0]!=0) {
  				fill.push_back(array[i0][j0]);
	  			array[i0][j0] = 0;
				}
		  }

		  //Moving to right-side strip

		  if(j0!=array.size()-1) {
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
			  array[ix][jx] += fill[s]/(2*(s+1));
			  array[i0][j0] += fill[s]/(2*(fill.size()-s));
		  }

		}
  }

}
