#include "singleOperationExamples.h"

#if DOLLDLA

RealPSet* GemmExample(Type dataType, Trans transA, Trans transB, int m, int n, int p)
{
  InputNode *Ain= new InputNode("A input", m, p, "A",
				 1, m,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride", dataType);

  InputNode *Bin = new InputNode("B input", p, n, "B",
				 1, p,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride", dataType);

  InputNode *Cin = new InputNode("C input", m, n, "C",
				 1, m,
				 "CNumRows","CNumCols",
				 "CRowStride","CColStride", dataType);

  Tunnel *tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel *tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(Bin, 0);

  Tunnel *tunC = new Tunnel(POSSTUNIN);
  tunC->AddInput(Cin, 0);

  Gemm *gemm = new Gemm(ABSLAYER, transA, transB, COEFONE, COEFONE, dataType);
  gemm->AddInputs(6,
		  tunA, 0,
		  tunB, 0,
		  tunC, 0);

  Poss *innerPoss = new Poss(gemm,true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Cout,true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* DotExample(Type dataType, int m)
{
  InputNode* Ain = new InputNode("A input", 1, m, "A", 
				 m, 1,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride", dataType);

  InputNode* Bin = new InputNode("B input", m, 1, "B", 
				 m, 1,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride", dataType);

  InputNode* Cin = new InputNode("C input", 1, 1, "C", 
				 m, 1,
				 "CNumRows","CNumCols",
				 "CRowStride","CColStride", dataType);

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

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* MVMulExample(Type dataType, bool transpose, int m, int n)
{
  InputNode* Ain;
  if (transpose) {
    Ain = new InputNode("A input", n, m, "A",
			1, n,
			"ANumRows", "ANumCols",
			"ARowStride", "AColStride", dataType);
  } else {
    Ain = new InputNode("A input", m, n, "A",
			1, m,
			"ANumRows", "ANumCols",
			"ARowStride", "AColStride", dataType);
  }

  InputNode* xIn = new InputNode("x input", n, 1, "X",
				 1, n,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);
  InputNode* yIn = new InputNode("y input", m, 1, "Y",
				 1, m,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

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

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* MAddExample(Type dataType, int m, int n)
{
  InputNode* Ain = new InputNode("A input", m, n, "A", 
				 1, m,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride", dataType);

  InputNode* Bin = new InputNode("B input", m, n, "B", 
				 1, m,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride", dataType);

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

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* SVMulRowExample(Type dataType, int m)
{
  InputNode* Ain = new InputNode("A input", 1, m, "A",
				 m, 1,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride", dataType);

  InputNode* xIn = new InputNode("x input", 1, 1, "X",
				 m, 1,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  SVMul* svmul = new SVMul(ROWVECTOR, ABSLAYER);
  svmul->AddInputs(4,
		   tunX, 0,
		   tunA, 0);

  Poss *innerPoss = new Poss(svmul, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* SVMulColExample(Type dataType, int m)
{
  InputNode* Ain = new InputNode("A input", m, 1, "A",
				 m, 1,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride", dataType);

  InputNode* xIn = new InputNode("x input", 1, 1, "X",
				 m, 1,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  SVMul* svmul = new SVMul(COLVECTOR, ABSLAYER);
  svmul->AddInputs(4,
		   tunX, 0,
		   tunA, 0);

  Poss *innerPoss = new Poss(svmul, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* SMMulExample(Type dataType, int m, int n)
{
  InputNode* Ain = new InputNode("A input", m, n, "A",
				 n, 1,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride", dataType);
  InputNode* xIn = new InputNode("x input", 1, 1, "X",
				 1, 1,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

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

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* VMMulExample(Type dataType, int m, int n)
{
  InputNode* Ain = new InputNode("A input", m, n, "A",
				 n, 1,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride", dataType);
  InputNode* xIn = new InputNode("x input", 1, m, "X",
				 m, 1,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  InputNode* yIn = new InputNode("y input", 1, n, "Y",
				 n, 1,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

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

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* VAddExample(Type dataType, int m)
{
  InputNode* xIn = new InputNode("x input", m, 1, "X",
				 1, m,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  InputNode* yIn = new InputNode("y input", m, 1, "Y",
				 1, m,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);
  
  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  VAdd* vadd = new VAdd(COLVECTOR, ABSLAYER);
  vadd->AddInputs(4,
		  tunX, 0,
		  tunY, 0);

  Poss* innerPoss = new Poss(vadd, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;

}

#endif //DOLLDLA
