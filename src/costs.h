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



#pragma once

#include <string>
#include "sizes.h"

using namespace std;

#define GAMMAVAL 1
#define BETAVAL 100 * GAMMAVAL
#define ALPHAVAL 120000*GAMMAVAL
#define PSIWVAL .15*GAMMAVAL
#define PSIRVAL .1*GAMMAVAL
#if DOELEM
#define RVAL 40
#define CVAL 40
#define ELEM_BSVAL 128
#elif DOBLIS
#define RVAL 1
#define CVAL 1
#define NUMCORESPERL2 2
#define NUML2PERPROC 3
#define NUMPROCS 4
#define BLIS_MC_BSVAL 368
#define BLIS_KC_BSVAL 256
#define BLIS_NC_BSVAL 4096
#define BLIS_OUTER_BSVAL 256
#elif DOTENSORS
#define TENSOR_BSVAL 32
#define NUM_GRID_DIMS 4
#elif DOLLDLA
#define RVAL 40
#define CVAL 40

#define MU_VAR_NAME "MUVALUE"
#endif 
#if TWOD
#define PVAL RVAL*CVAL
#endif
#define CACHEMISSVAL 50 * GAMMAVAL


extern Cost ALPHA;
extern Cost BETA;
extern Cost GAMMA;
extern Cost PSIW;
extern Cost PSIR;
extern Cost CACHEMISS;
#if TWOD
extern Size P;
extern Size R;
extern Size C;
#elif DOTENSORS
extern Size GridLens[NUM_GRID_DIMS];
#endif
#if DOELEM
extern Size ELEM_BS;
#elif DOBLIS
extern Size BLIS_MC_BS;
extern Size BLIS_KC_BS;
extern Size BLIS_NC_BS;
extern Size BLIS_OUTER_BS;
#elif DOTENSORS
extern Size TENSOR_BS;
#endif
extern Cost ONE;
extern Cost ZERO;
extern Cost TWO;

Cost AllGather(Size totalSize, Size numProcs);
Cost AllReduce(Size totalSize, Size numProcs);
Cost ReduceScatter(Size totalSize, Size numProcs);
Cost SendRecv(Size totalSize);
Cost AllToAll(Size totalSize, Size numProcs);
Cost CopyCost(Size inner, Size outer, Size readLdim, Size writeLdim, bool cacheMissOnOuterLoop = false, bool writeCacheMissOnOuterLoop = false);

typedef vector<Cost> CostVec;
typedef CostVec::iterator CostVecIter;


Size GridModeLens(const DimVec &modes);
