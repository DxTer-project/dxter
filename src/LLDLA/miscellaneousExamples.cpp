#include "madd.h"
#include "miscellaneousExamples.h"
#include "mvmul.h"
#include "partition.h"
#include "recombine.h"
#include "svmul.h"
#include "vadd.h"
#include "vvdot.h"

#if DOLLDLA

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

  SVMul* topSVMul = new SVMul(COLVECTOR, ABSLAYER);
  topSVMul->AddInputs(4,
		      tunX, 0,
		      part, 0);

  SVMul* bottomSVMul = new SVMul(COLVECTOR, ABSLAYER);
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

  VAdd* vadd1 = new VAdd(COLVECTOR, ABSLAYER);
  vadd1->AddInputs(4,
		  tunX, 0,
		  tunY, 0);

  VAdd* vadd2 = new VAdd(COLVECTOR, ABSLAYER);
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
