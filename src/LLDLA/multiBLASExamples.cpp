#include "multiBLASExamples.h"

#if DOLLDLA

RealPSet* MAddGemm(Type dataType, int m, int n, int p)
{
  InputNode* alphaIn = new InputNode("alpha input", 1, 1, "Alpha",
				     1, m,
				     "AlphaNumRows", "AlphaNumCols",
				     "AlphaRowStride", "AlphaColStride", dataType);

  InputNode* betaIn = new InputNode("beta input", 1, 1, "Beta",
				    1, m,
				    "BetaNumRows", "BetaNumCols",
				    "BetaRowStride", "BetaColStride", dataType);

  InputNode* A0In = new InputNode("A0 input", m, n, "A0",
				 n, 1,
				 "A0NumRows", "A0NumCols",
				 "A0RowStride", "A0ColStride", dataType);

  InputNode* A1In = new InputNode("A1 input", m, n, "A1",
				 n, 1,
				 "A1NumRows", "A1NumCols",
				 "A1RowStride", "A1ColStride", dataType);

  InputNode* BIn = new InputNode("B input", m, p, "B",
				 1, m,
				 "BNumRows", "BNumCols",
				 "BRowStride", "BColStride", dataType);

  InputNode* CIn = new InputNode("C input", m, p, "C",
				 1, m,
				 "CNumRows", "CNumCols",
				 "CRowStride", "CColStride", dataType);

  Tunnel* tunA0 = new Tunnel(POSSTUNIN);
  tunA0->AddInput(A0In, 0);

  Tunnel* tunA1 = new Tunnel(POSSTUNIN);
  tunA1->AddInput(A1In, 0);

  Tunnel* tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(BIn, 0);

  Tunnel* tunC = new Tunnel(POSSTUNIN);
  tunC->AddInput(CIn, 0);

  Tunnel* tunAlpha = new Tunnel(POSSTUNIN);
  tunAlpha->AddInput(alphaIn, 0);

  Tunnel* tunBeta = new Tunnel(POSSTUNIN);
  tunBeta->AddInput(betaIn, 0);

  MAdd* madd = new MAdd(ABSLAYER);
  madd->AddInputs(4,
		  tunA0, 0,
		  tunA1, 0);

  /*  LLDLATranspose* trans = new LLDLATranspose(ABSLAYER);
  trans->AddInput(madd, 0);

  SMMul* alphaA = new SMMul(ABSLAYER);
  alphaA->AddInputs(4,
		    tunAlpha, 0,
		    trans, 0);

  SMMul* betaC = new SMMul(ABSLAYER);
  betaC->AddInputs(4,
		   tunBeta, 0,
		   tunC, 0);

  Gemm* gemm = new Gemm(ABSLAYER, NORMAL, NORMAL, COEFONE, COEFONE, dataType);
  gemm->AddInputs(6,
		  alphaA, 0,
		  tunB, 0,
		  tunC, 0);*/

  Poss* innerPoss = new Poss(madd, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* Gemv2(Type dataType, int m, int n, int k)
{
  InputNode* alphaIn = new InputNode("alpha input", 1, 1, "Alpha",
				     1, m,
				     "AlphaNumRows", "AlphaNumCols",
				     "AlphaRowStride", "AlphaColStride", dataType);

  InputNode* betaIn = new InputNode("beta input", 1, 1, "Beta",
				    1, m,
				    "BetaNumRows", "BetaNumCols",
				    "BetaRowStride", "BetaColStride", dataType);

  InputNode* AIn = new InputNode("A input", m, n, "A",
				 n, 1,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride", dataType);

  InputNode* BIn = new InputNode("B input", m, k, "B",
				 1, m,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

  InputNode* xIn = new InputNode("y input", n, 1, "Y",
				 1, n,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride", dataType);

  InputNode* yIn = new InputNode("w input", k, 1, "Y",
				 1, k,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

  InputNode* zIn = new InputNode("z input", m, 1, "Z",
				 1, m,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride", dataType);

  InputNode* wIn = new InputNode("w input", m, 1, "W",
				 1, m,
				 "WNumRows", "WNumCols",
				 "WRowStride", "WColStride", dataType);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  Tunnel* tunW = new Tunnel(POSSTUNIN);
  tunW->AddInput(wIn, 0);

  Tunnel* tunZ = new Tunnel(POSSTUNIN);
  tunZ->AddInput(zIn, 0);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(AIn, 0);

  Tunnel* tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(BIn, 0);

  Tunnel* tunAlpha = new Tunnel(POSSTUNIN);
  tunAlpha->AddInput(alphaIn, 0);

  Tunnel* tunBeta = new Tunnel(POSSTUNIN);
  tunBeta->AddInput(betaIn, 0);

  MVMul* aX = new MVMul(ABSLAYER);
  aX->AddInputs(6,
		   tunA, 0,
		   tunX, 0,
		   tunZ, 0);

  MVMul* bY = new MVMul(ABSLAYER);
  bY->AddInputs(6,
		tunB, 0,
		tunY, 0,
		tunW, 0);

  SVMul* alphaAX = new SVMul(COLVECTOR, ABSLAYER);
  alphaAX->AddInputs(4,
		    tunAlpha, 0,
		    aX, 0);

  SVMul* betaBY = new SVMul(COLVECTOR, ABSLAYER);
  betaBY->AddInputs(4,
		    tunBeta, 0,
		    bY, 0);
  
  VAdd* sum = new VAdd(COLVECTOR, ABSLAYER);
  sum->AddInputs(4,
		 alphaAX, 0,
		 betaBY, 0);

  Poss* innerPoss = new Poss(sum, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* VMVMulExample(Type dataType, int m, int n)
{
  InputNode* Ain = new InputNode("A input", m, n, "A",
				 1, m,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride", dataType);

  InputNode* xIn = new InputNode("x input", n, 1, "X",
				 1, n,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  InputNode* zIn = new InputNode("z input", m, 1, "Z",
				 1, m,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride", dataType);

  InputNode* yIn = new InputNode("y input", 1, m, "Y",
				 1, 1,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

  InputNode* wIn = new InputNode("w input", 1, 1, "W",
				 1, 1,
				 "WNumRows", "WNumCols",
				 "WRowStride", "WColStride", dataType);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunZ = new Tunnel(POSSTUNIN);
  tunZ->AddInput(zIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  Tunnel* tunW = new Tunnel(POSSTUNIN);
  tunW->AddInput(wIn, 0);

  MVMul* mvmul = new MVMul(ABSLAYER);
  mvmul->AddInputs(6,
		   tunA, 0,
		   tunX, 0,
		   tunZ, 0);

  VVDot* vvdot = new VVDot(ABSLAYER);
  vvdot->AddInputs(6,
		   tunY, 0,
		   mvmul, 0,
		   tunW, 0);

  Poss* innerPoss = new Poss(vvdot, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

#endif // DOLLDLA
