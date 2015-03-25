#include "LLDLATranspose.h"
#include "madd.h"
#include "multiBLASExamples.h"
#include "mvmul.h"
#include "setToZero.h"
#include "smmul.h"
#include "svmul.h"
#include "vadd.h"
#include "vvdot.h"

#if DOLLDLA

RealPSet* Gesummv(Type dataType, int m, int n)
{
  InputNode* alphaIn = new InputNode("Alpha",
				     1, 1,
				     1, 1,
				     dataType);

  InputNode* betaIn = new InputNode("Beta",
				    1, 1,
				    1, 1,
				    dataType);

  InputNode* AIn = new InputNode("A",
				 m, n,
				 1, m,
				 dataType);

  InputNode* BIn = new InputNode("B",
				 m, n,
				 1, m,
				 dataType);

  InputNode* YIn = new InputNode("Y",
				 m, 1,
				 1, m,
				 dataType);

  InputNode* XIn = new InputNode("X",
				 n, 1,
				 1, n,
				 dataType);

  Tunnel* tunAlpha = new Tunnel(POSSTUNIN);
  tunAlpha->AddInput(alphaIn, 0);

  Tunnel* tunBeta = new Tunnel(POSSTUNIN);
  tunBeta->AddInput(betaIn, 0);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(AIn, 0);

  Tunnel* tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(BIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(YIn, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(XIn, 0);

  auto zeroY = new SetToZero(ABSLAYER);
  zeroY->AddInput(tunY, 0);

  SMMul* alphaA = new SMMul(ABSLAYER);
  alphaA->AddInputs(4,
		    tunAlpha, 0,
		    tunA, 0);

  SMMul* betaB = new SMMul(ABSLAYER);
  betaB->AddInputs(4,
		   tunBeta, 0,
		   tunB, 0);

  MVMul* alphaY = new MVMul(ABSLAYER);
  alphaY->AddInputs(6,
		    alphaA, 0,
		    tunX, 0,
		    zeroY, 0);

  MVMul* betaBY = new MVMul(ABSLAYER);
  betaBY->AddInputs(6,
		    betaB, 0,
		    tunX, 0,
		    alphaY, 0);

  Poss* innerPoss = new Poss(betaBY, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* Gemam(Type dataType, int m, int n, int p)
{
  InputNode* alphaIn = new InputNode("Alpha",
				     1, 1,
				     1, m,
				     dataType);

  InputNode* betaIn = new InputNode("Beta",
				    1, 1,
				    1, m,
				    dataType);

  InputNode* A0In = new InputNode("A0",
				  m, n,
				  1, m,
				  dataType);

  InputNode* A1In = new InputNode("A1",
				  m, n,
				  1, m,
				  dataType);

  InputNode* BIn = new InputNode("B",
				 m, p,
				 1, m,
				 dataType);

  InputNode* CIn = new InputNode("C",
				 n, p,
				 1, p,
				 dataType);

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

  LLDLATranspose* trans = new LLDLATranspose(ABSLAYER);
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
		  betaC, 0);

  Poss* innerPoss = new Poss(gemm, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* Gemv2(Type dataType, int m, int n, int k)
{
  InputNode* alphaIn = new InputNode("Alpha",
				     1, 1,
				     1, m,
				     dataType);

  InputNode* betaIn = new InputNode("Beta",
				    1, 1,
				    1, m,
				    dataType);

  InputNode* AIn = new InputNode("A",
				 m, n,
				 n, 1,
				 dataType);

  InputNode* BIn = new InputNode("B",
				 m, k,
				 1, m,
				 dataType);

  InputNode* xIn = new InputNode("Y",
				 n, 1,
				 1, n,
				 dataType);

  InputNode* yIn = new InputNode("W",
				 k, 1,
				 1, k,
				 dataType);

  InputNode* zIn = new InputNode("Z",
				 m, 1,
				 1, m,
				 dataType);

  InputNode* wIn = new InputNode("W",
				 m, 1,
				 1, m,
				 dataType);

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

  SVMul* alphaAX = new SVMul(ABSLAYER);
  alphaAX->AddInputs(4,
		    tunAlpha, 0,
		    aX, 0);

  SVMul* betaBY = new SVMul(ABSLAYER);
  betaBY->AddInputs(4,
		    tunBeta, 0,
		    bY, 0);
  
  VAdd* sum = new VAdd(ABSLAYER, COLVECTOR);
  sum->AddInputs(4,
		 alphaAX, 0,
		 betaBY, 0);

  Poss* innerPoss = new Poss(sum, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* VMVMul(Type dataType, int m, int n)
{
  InputNode* Ain = new InputNode("A",
				 m, n,
				 1, m,
				 dataType);

  InputNode* xIn = new InputNode("X",
				 n, 1,
				 1, n,
				 dataType);

  InputNode* zIn = new InputNode("Z",
				 m, 1,
				 1, m,
				 dataType);

  InputNode* yIn = new InputNode("Y",
				 1, m,
				 1, 1,
				 dataType);

  InputNode* wIn = new InputNode("W",
				 1, 1,
				 1, 1,
				 dataType);

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

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

#endif // DOLLDLA
