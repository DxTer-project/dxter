/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2015, The University of Texas and Bryan Marker

    DxTer is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DxTer is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.               

    You should have received a copy of the GNU General Public License
    along with DxTer.  If not, see <http://www.gnu.org/licenses/>.
*/



#include "transform.h"
#include "costs.h"
#include <stdio.h>
#include <string>
#include <iomanip>
#include <sstream>
#include <cmath>

using namespace std;

Cost ALPHA(ALPHAVAL);
Cost BETA(BETAVAL);
Cost GAMMA(GAMMAVAL);
Cost PSIW(PSIWVAL);
Cost PSIR(PSIRVAL);
Cost CACHEMISS(CACHEMISSVAL);
#if TWOD
Size P(PVAL);
Size R(RVAL);
Size C(CVAL);
#elif DOTENSORS
Size GridLens[NUM_GRID_DIMS] = {8, 4, 4, 4};
#endif
#if DOELEM
Size ELEM_BS(ELEM_BSVAL);
#elif DOBLIS
Size BLIS_MC_BS(BLIS_MC_BSVAL);
Size BLIS_KC_BS(BLIS_KC_BSVAL);
Size BLIS_NC_BS(BLIS_NC_BSVAL);
Size BLIS_OUTER_BS(BLIS_OUTER_BSVAL);
#elif DOTENSORS
Size TENSOR_BS(TENSOR_BSVAL);
#endif
Cost TWO(2);
Cost ONE(1);
Cost ZERO(0);

inline double log2(double in)
{
  return ceil(log(in) / log(2.0));
}

Cost AllGather(Size totalSize, Size numProcs)
{
  if (!totalSize)
    return 0;
  if (totalSize < 0) {
    cout << "totalSize in AllGather\n";
    throw;
  }
  return ( log2(numProcs) * ALPHA )
    + (((numProcs - ONE) / numProcs) * totalSize * BETA);
}

Cost AllReduce(Size totalSize, Size numProcs)
{
  if (!totalSize)
    return 0;
  if (totalSize < 0) {
    cout << "totalSize in AllReduce\n";
    throw;
  }
  return ( log2(numProcs) * ALPHA )
    + (TWO * ((numProcs - ONE) / numProcs) * totalSize * BETA)
    + (((numProcs - ONE) / numProcs) * totalSize * GAMMA);
}

Cost ReduceScatter(Size totalSize, Size numProcs)
{
  if (!totalSize)
    return 0;
  if (totalSize < 0) {
    cout << "totalSize in ReduceScatter\n";
    throw;
  }
      
  return (log2(numProcs) * ALPHA)
    + (((numProcs - ONE) / numProcs) * totalSize * (BETA + GAMMA));  
}

Cost SendRecv(Size totalSize)
{
  if (!totalSize)
    return 0;
  if (totalSize < 0) {
    cout << "totalSize in SendRecv\n";
    throw;
  }
  return ALPHA + (totalSize * BETA);
}

Cost AllToAll(Size totalSize, Size numProcs)
{
  if (!totalSize)
    return 0;
  if (totalSize < 0) {
    cout << "totalSize in AllToAll\n";
    throw;
  }
  if (!(totalSize/numProcs)){
    cout << "nope\n";
    throw;
  }
  return (numProcs - 1) * SendRecv(totalSize/numProcs);
}

Cost CopyCost(Size inner, Size outer, Size readLdim, Size writeLdim, bool readCacheMissOnOuterLoop, bool writeCacheMissOnOuterLoop)
{
  if (!outer || !inner)
    return 0;
  if (outer < 0 || inner < 0 || readLdim < 0 || writeLdim < 0) {
    cout << "Copy Cost problem\n";
    throw;
  }
      
  Cost cost = 0;
  if (max(readLdim, writeLdim) < 16)
    cost = inner * outer * PSIR;
  else
    cost = inner * outer * (PSIR + 10 * GAMMAVAL);
  if (readCacheMissOnOuterLoop || writeCacheMissOnOuterLoop)
    cost += (outer * CACHEMISS);
  return cost;
}
