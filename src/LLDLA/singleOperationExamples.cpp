#include "singleOperationExamples.h"

#if DOLLDLA

RealPSet* DotTest(Type dataType, int m)
{
  InputNode* Ain = new InputNode("A", 1, m,
				 m, 1,
				 dataType);

  InputNode* Bin = new InputNode("B", m, 1,
				 m, 1,
				 dataType);

  InputNode* Cin = new InputNode("C", 1, 1,
				 m, 1,
				 dataType);

  Tunnel *tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain,0);

  Tunnel *tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(Bin,0);

  Tunnel *tunC = new Tunnel(POSSTUNIN);
  tunC->AddInput(Cin,0);

  VVDot* dot = new VVDot(ABSLAYER);
  dot->AddInputs(6,
		 tunA, 0,
		 tunB, 0,
		 tunC, 0);

  Poss *innerPoss = new Poss(dot, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* MVMulTest(Type dataType, bool transpose, int m, int n)
{
  InputNode* Ain;
  if (transpose) {
    Ain = new InputNode("A", n, m,
			1, n,
			dataType);
  } else {
    Ain = new InputNode("A", m, n,
			1, m,
			dataType);
  }

  InputNode* xIn = new InputNode("x", n, 1,
				 1, n,
				 dataType);

  InputNode* yIn = new InputNode("y", m, 1,
				 1, m,
				 dataType);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  MVMul* mvmul = new MVMul(ABSLAYER);
  LLDLATranspose* trans = new LLDLATranspose(ABSLAYER);
  if (transpose) {
    trans->AddInputs(2,
		     tunA, 0);
    mvmul->AddInputs(6,
		     trans, 0,
		     tunX, 0,
		     tunY, 0);
  } else {
    mvmul->AddInputs(6,
		     tunA, 0,
		     tunX, 0,
		     tunY, 0);
  }
  Poss *innerPoss = new Poss(mvmul, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* MAddTest(Type dataType, int m, int n)
{
  InputNode* Ain = new InputNode("A", m, n,
				 1, m,
				 dataType);

  InputNode* Bin = new InputNode("B", m, n,
				 1, m,
				 dataType);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(Bin, 0);

  MAdd* madd = new MAdd(ABSLAYER);
  madd->AddInputs(4,
		 tunA, 0,
		 tunB, 0);

  Poss *innerPoss = new Poss(madd, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* SVMulTest(Type dataType, VecType vecType, int m)
{
  int numRows, numCols, rowStride, colStride;
  if (vecType == ROWVECTOR) {
    numRows = 1;
    numCols = m;
    rowStride = m;
    colStride = 1;
  } else {
    numRows = m;
    numCols = 1;
    rowStride = 1;
    colStride = m;
  }

  InputNode* Ain = new InputNode("A", numRows, numCols,
				 rowStride, colStride,
				 dataType);

  InputNode* xIn = new InputNode("x", 1, 1,
				 1, 1,
				 dataType);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  SVMul* svmul = new SVMul(vecType, ABSLAYER);
  svmul->AddInputs(4,
		   tunX, 0,
		   tunA, 0);

  Poss *innerPoss = new Poss(svmul, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* SMMulTest(Type dataType, int m, int n)
{
  InputNode* Ain = new InputNode("A", m, n,
				 n, 1,
				 dataType);

  InputNode* xIn = new InputNode("x", 1, 1,
				 1, 1,
				 dataType);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  SMMul* smmul = new SMMul(ABSLAYER);
  smmul->AddInputs(4,
		   tunX, 0,
		   tunA, 0);

  Poss *innerPoss = new Poss(smmul, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* VMMulTest(Type dataType, int m, int n)
{
  InputNode* Ain = new InputNode("A", m, n,
				 n, 1,
				 dataType);

  InputNode* xIn = new InputNode("x", 1, m,
				 m, 1,
				 dataType);

  InputNode* yIn = new InputNode("y", 1, n,
				 n, 1,
				 dataType);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  VMMul* vmmul = new VMMul(ABSLAYER);
  vmmul->AddInputs(6,
		   tunX, 0,
		   tunA, 0,
		   tunY, 0);

  Poss *innerPoss = new Poss(vmmul, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* VAddTest(Type dataType, VecType vecType, int m)
{

  int nRows, nCols;
  if (vecType == ROWVECTOR) {
    nRows = 1;
    nCols = m;
  } else {
    nRows = m;
    nCols = 1;
  }

  InputNode* xIn = new InputNode("x", nRows, nCols,
				 nCols, nRows,
				 dataType);

  InputNode* yIn = new InputNode("y", nRows, nCols,
				 nCols, nRows,
				 dataType);
  
  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  VAdd* vadd = new VAdd(vecType, ABSLAYER);
  vadd->AddInputs(4,
		  tunX, 0,
		  tunY, 0);

  Poss* innerPoss = new Poss(vadd, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

#endif //DOLLDLA
