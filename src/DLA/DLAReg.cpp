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



#include "DLAReg.h"
#include "universe.h"
#include "loopSupport.h"
#include "critSect.h"

#define Reg(NODE) Universe::RegCons(NODE::GetClass(), &(NODE::BlankInst))

#if DOTENSORS
void RegAllTensorNodes()
{
  Reg(Contraction);
}
#elif DOBLIS||DOELEM
void RegAllDLANodes()
{
  Reg(TriInv);
  Reg(Axpy);
  Reg(Scal);
  Reg(ConstScal);
#if DOELEM
  Reg(Chol);
  Reg(RedistNode);
  Reg(SumScatterNode);
  Reg(SumScatterFrom);
  Reg(SumOverCommNode);
  Reg(LocalSymmAcc);
  Reg(MakeTrapNode);
  Reg(LocalTrmmAcc);
  Reg(ViewAroundDiag);
  Reg(ViewAroundDiagCombine);
  Reg(ViewPan);
#endif
  Reg(Gemm);
  Reg(InputNode);
  //  Reg(ConstVal);
  Reg(TempVarNode);
  Reg(OutputNode);
  Reg(ScaleNode);
  Reg(ScaleTrapNode);
  Reg(Hemm);
  Reg(Her2k);
  Reg(Tri2k);
  Reg(Herk);
  Reg(TriRK);
#if DOELEM
  Reg(Hetrmm);
#endif
#ifndef SKIPTWOSIDED
  Reg(TwoSidedTrxm);
#endif
  Reg(CombineSingleIter);
  Reg(SplitSingleIter);
  Reg(LoopTunnel);
#if 0
  Reg(CritSectTunnel);
#endif
  Reg(Tunnel);
  Reg(Trxm);
#if DOBLIS
  Reg(Transpose);
  Reg(TrxmBP);
#endif
  Reg(LU);
  Reg(PanelLU);
  Reg(ViewTL);
  Reg(ViewTLCombine);
}
#elif DOLLDLA
void RegAllLLDLANodes()
{

}
#elif DOBOOL
#elif DORQO
#else
asdklfja
#endif
