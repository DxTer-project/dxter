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

#include "madd.h"

#if DOLLDLA

#include "exampleUtils.h"
#include "miscellaneousExamples.h"
#include "mvmul.h"
#include "partition.h"
#include "recombine.h"
#include "setToZero.h"
#include "svmul.h"
#include "LLDLATranspose.h"
#include "vadd.h"
#include "vvdot.h"

RealPSet* BasicMultiAssign(Type dataType, int m, int n) {
  auto Ain = InputTunnel("A",
			 m, n,
			 1, m,
			 dataType);

  auto q = InputTunnel("q",
		       m, 1,
		       1, m,
		       dataType);


  auto p = InputTunnel("p",
		       n, 1,
		       1, n,
		       dataType);

  auto s = InputTunnel("s",
		       n, 1,
		       1, n,
		       dataType);

  auto r = InputTunnel("r",
		       m, 1,
		       1, m,
		       dataType);

  auto qZ = new SetToZero(ABSLAYER);
  qZ->AddInput(q, 0);

  auto sZ = new SetToZero(ABSLAYER);
  sZ->AddInput(s, 0);

  auto Aqz = new MVMul(ABSLAYER);
  Aqz->AddInputs(6,
		 Ain, 0,
		 p, 0,
		 qZ, 0);

  auto AT = new LLDLATranspose(ABSLAYER);
  AT->AddInput(Ain, 0);

  auto ATrs = new MVMul(ABSLAYER);
  ATrs->AddInputs(6,
		  AT, 0,
		  r, 0,
		  sZ, 0);

  Poss* innerPoss = new Poss(2, Aqz, ATrs, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode* out1 = new OutputNode;
  out1->AddInput(innerSet->OutTun(0), 0);

  OutputNode* out2 = new OutputNode;
  out2->AddInput(innerSet->OutTun(1), 0);

  Poss* outerPoss = new Poss(2, out1, out2, true);
  RealPSet* outerSet = new RealPSet(outerPoss);

  return outerSet;
}

RealPSet* LGenCompareL2(Type dataType, int m, int n) {
  auto x = InputTunnel("X",
		       n, 1,
		       1, n,
		       dataType);

  auto a = InputTunnel("A",
		       m, n,
		       1, m,
		       dataType);

  auto b = InputTunnel("B",
		       n, m,
		       1, n,
		       dataType);

  auto y = InputTunnel("Y",
		       m, 1,
		       1, m,
		       dataType);

  auto bTrans = new LLDLATranspose(ABSLAYER);
  bTrans->AddInput(b, 0);

  auto aPlusBTrans = new MAdd(ABSLAYER);
  aPlusBTrans->AddInputs(4,
			 a, 0,
			 bTrans, 0);

  auto zeroY = new SetToZero(ABSLAYER);
  zeroY->AddInput(y, 0);

  auto result = new MVMul(ABSLAYER);
  result->AddInputs(6,
		    aPlusBTrans, 0,
		    x, 0,
		    zeroY, 0);

  return WrapInPSet(result);
}

RealPSet* VMulAddBenchmark(Type dataType, int m) {
  auto alpha = InputTunnel("alpha",
			     1, 1,
			     1, 1,
			     dataType);

  auto beta = InputTunnel("beta",
			    1, 1,
			    1, 1,
			    dataType);

  auto x = InputTunnel("X",
			 m, 1,
			 1, m,
			 dataType);

  auto y = InputTunnel("Y",
			 m, 1,
			 1, m,
			 dataType);

  auto z = InputTunnel("Z",
			 m, 1,
			 1, m,
			 dataType);

  auto zPlusY = new VAdd(ABSLAYER, COLVECTOR);
  zPlusY->AddInputs(4,
		    z, 0,
		    y, 0);

  auto betaZY = new SVMul(ABSLAYER);
  betaZY->AddInputs(4,
		    beta, 0,
		    zPlusY, 0);

  auto alphaX = new SVMul(ABSLAYER);
  alphaX->AddInputs(4,
		    alpha, 0,
		    x, 0);

  auto finalY = new VAdd(ABSLAYER, COLVECTOR);
  finalY->AddInputs(4,
		    alphaX, 0,
		    betaZY, 0);

  Poss *innerPoss = new Poss(finalY, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);

  return outerSet;
  //  return WrapInPSet(finalY);
}

RealPSet* GenSizeColSVMul(Type dataType, int m)
{
  auto Ain = new InputNode("A input", m, 1, "A",
			   m, 1,
			   "ANumRows", "ANumCols",
			   "ARowStride", "AColStride", dataType);

  auto xIn = new InputNode("x input", 1, 1, "X",
			   m, 1,
			   "XNumRows", "XNumCols",
			   "XRowStride", "XColStride", dataType);

  auto tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Size partSplitPoint = m - (m % arch->VecRegWidth(dataType));
  cout << "Part split point = " << partSplitPoint << endl;
  Partition* part =
    new Partition(ABSLAYER, VERTICAL, partSplitPoint);
  part->AddInput(tunA, 0);

  auto tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  SVMul* topSVMul = new SVMul(ABSLAYER);
  topSVMul->AddInputs(4,
		      tunX, 0,
		      part, 0);

  SVMul* bottomSVMul = new SVMul(ABSLAYER);
  bottomSVMul->AddInputs(4,
			 tunX, 0,
			 part, 1);

  Recombine* recombine = new Recombine(ABSLAYER, VERTICAL);
  recombine->AddInputs(6,
		       topSVMul, 0,
		       bottomSVMul, 0,
		       tunA, 0);

  Poss *innerPoss = new Poss(recombine, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);

  return outerSet;
}

RealPSet* MVMul2(Type dataType, int m, int n, int p)
{
  auto xIn = new InputNode("x input", p, 1, "X",
				 1, p,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  auto yIn = new InputNode("y input", n, 1, "Y",
				 1, n,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

  auto zIn = new InputNode("z input", m, 1, "Z",
				 1, m,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride", dataType);

  auto AIn = new InputNode("a input", m, n, "A",
				 1, m,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride", dataType);

  auto BIn = new InputNode("b input", n, p, "B",
				 1, n,
				 "BNumRows", "BNumCols",
				 "BRowStride", "BColStride", dataType);
  
  auto tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  auto tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  auto tunZ = new Tunnel(POSSTUNIN);
  tunZ->AddInput(zIn, 0);

  auto tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(AIn, 0);

  auto tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(BIn, 0);

  MVMul* mvmul1 = new MVMul(ABSLAYER);
  mvmul1->AddInputs(6,
		    tunB, 0,
		    tunX, 0,
		    tunY, 0);

  MVMul* mvmul2 = new MVMul(ABSLAYER);
  mvmul2->AddInputs(6,
		    tunA, 0,
		    mvmul1, 0,
		    tunZ, 0);

  Poss* innerPoss = new Poss(mvmul2, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* MAdd2(Type dataType, int m, int n)
{
  auto xIn = new InputNode("x input", m, n, "X",
				 1, m,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  auto yIn = new InputNode("y input", m, n, "Y",
				 1, m,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

  auto zIn = new InputNode("z input", m, n, "Z",
				 1, m,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride", dataType);
  
  auto tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  auto tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  auto tunZ = new Tunnel(POSSTUNIN);
  tunZ->AddInput(zIn, 0);

  MAdd* madd1 = new MAdd(ABSLAYER);
  madd1->AddInputs(4,
		  tunX, 0,
		  tunY, 0);

  MAdd* madd2 = new MAdd(ABSLAYER);
  madd2->AddInputs(4,
		   tunZ, 0,
		   madd1, 0);

  Poss* innerPoss = new Poss(madd2, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* VAdd2(Type dataType, int m)
{
  auto xIn = new InputNode("x input", m, 1, "X",
				 1, m,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  auto yIn = new InputNode("y input", m, 1, "Y",
				 1, m,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

  auto zIn = new InputNode("z input", m, 1, "Z",
				 1, m,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride", dataType);
  
  auto tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  auto tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  auto tunZ = new Tunnel(POSSTUNIN);
  tunZ->AddInput(zIn, 0);

  VAdd* vadd1 = new VAdd(ABSLAYER, COLVECTOR);
  vadd1->AddInputs(4,
		  tunX, 0,
		  tunY, 0);

  VAdd* vadd2 = new VAdd(ABSLAYER, COLVECTOR);
  vadd2->AddInputs(4,
		   tunZ, 0,
		   vadd1, 0);

  Poss* innerPoss = new Poss(vadd2, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;

}

RealPSet* DoubleGemm(Type dataType, Trans transA, Trans transB, int m, int n, int p, int k)
{
  InputNode *Ain = new InputNode("A input",  m, n, "A",
				 n, 1,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride", dataType);

  InputNode *Bin = new InputNode("B input", n, p, "B",
				 p, 1,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride", dataType);

  InputNode *Cin = new InputNode("C input",  p, k, "C",
				 k, 1,
				 "CNumRows","CNumCols",
				 "CRowStride","CColStride", dataType);

  Tunnel *tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain,0);

  Tunnel *tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(Bin,0);

  Tunnel *tunC = new Tunnel(POSSTUNIN);
  tunC->AddInput(Cin,0);

  Gemm *gemm1 = new Gemm(ABSLAYER, transA, transB, COEFONE, COEFONE, dataType);
  gemm1->AddInputs(6,
		  tunA,0,
		  tunB,0,
		  tunC,0);

  Gemm *gemm2 = new Gemm(ABSLAYER, transA, transB, COEFONE, COEFONE, dataType);
  gemm2->AddInputs(6,
		  tunA,0,
		  tunB,0,
		  gemm1,0);

  Poss *innerPoss = new Poss(gemm2,true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Cout,true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

#endif //DOLLDLA
