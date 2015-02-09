#include "blasExamples.h"
#include "DLAReg.h"
#include "LLDLATranspose.h"
#include "localInput.h"
#include "mvmul.h"
#include "smmul.h"
#include "svmul.h"
#include "vadd.h"


#if DOLLDLA

RealPSet* GemmTest(Type dataType, Trans transA, Trans transB, int m, int n, int p)
{
  InputNode* Ain;
  InputNode* Bin;
  InputNode* Cin;

  if (transA == NORMAL) {
    Ain = new InputNode("A", m, p,
			1, m,
			dataType);
  } else {
    Ain = new InputNode("A", p, m,
			1, p,
			dataType);
  }

  if (transB == NORMAL) {
    Bin = new InputNode("B", p, n,
			1, p,
			dataType);
  } else {
    Bin = new InputNode("B", n, p,
			1, n,
			dataType);
  }

  Cin = new InputNode("C", m, n,
		      1, m,
		      dataType);

  auto alpha = new InputNode("alpha", 1, 1,
			     1, 1,
			     dataType);

  auto beta = new InputNode("beta", 1, 1,
			    1, 1,
			    dataType);

  auto tunAlpha = new Tunnel(POSSTUNIN);
  tunAlpha->AddInput(alpha, 0);

  auto tunBeta = new Tunnel(POSSTUNIN);
  tunBeta->AddInput(beta, 0);

  Tunnel *tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel *tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(Bin, 0);

  Tunnel *tunC = new Tunnel(POSSTUNIN);
  tunC->AddInput(Cin, 0);

  auto betaC = new SMMul(ABSLAYER);
  betaC->AddInputs(4,
		   tunBeta, 0,
		   tunC, 0);

  auto transposeA = new LLDLATranspose(ABSLAYER);
  if (transA == TRANS) {
    transposeA->AddInput(tunA, 0);
  }

  auto transposeB = new LLDLATranspose(ABSLAYER);
  if (transB == TRANS) {
    transposeB->AddInput(tunB, 0);
  }

  SMMul* alphaA = new SMMul(ABSLAYER);
  if (transA == TRANS) {
    alphaA->AddInputs(4,
		      tunAlpha, 0,
		      transposeA, 0);
  } else {
    alphaA->AddInputs(4,
		     tunAlpha, 0,
		     tunA, 0);
  }

  Gemm *gemm = new Gemm(ABSLAYER, NORMAL, NORMAL, COEFONE, COEFONE, dataType);
  gemm->AddInput(alphaA, 0);
  if (transB == NORMAL) {
    gemm->AddInput(tunB, 0);
  } else {
    gemm->AddInput(transposeB, 0);
  }
  gemm->AddInput(betaC, 0);

  Poss *innerPoss = new Poss(gemm,true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Cout,true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* Axpy(Type dataType, VecType vType, int m)
{
  auto* alphaIn = new InputNode("alpha", 1, 1,
				1, m,
				dataType);

  int nCols, nRows;
  if (vType == COLVECTOR) {
    nCols = 1;
    nRows = m;
  } else {
    nCols = m;
    nRows = 1;
  }

  auto* xIn = new InputNode("x", nRows, nCols,
			    nCols, nRows,
			    dataType);

  auto* yIn = new InputNode("y", nRows, nCols,
			    nCols, nRows,
			    dataType);

  auto* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  auto* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  auto* tunAlpha = new Tunnel(POSSTUNIN);
  tunAlpha->AddInput(alphaIn, 0);

  auto* axMul = new SVMul(vType, ABSLAYER);
  axMul->AddInputs(4,
		   tunAlpha, 0,
		   tunX, 0);

  auto* axPlusY = new VAdd(vType, ABSLAYER);
  axPlusY->AddInputs(4,
		     axMul, 0,
		     tunY, 0);

  auto* innerPoss = new Poss(axPlusY, true);
  auto* innerSet = new RealPSet(innerPoss);

  auto* Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  auto* outerPoss = new Poss(Cout, true);
  auto* outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* Gemv(Type dataType, bool transpose, int m, int n)
{
  auto xIn = new InputNode("X", n, 1,
			    1, n,
			    dataType);

  auto yIn = new InputNode("Y", m, 1,
			    1, m,
			    dataType);

  auto vIn = new LocalInput("V",
			    m, 1,
			    1, m,
			    dataType);

  InputNode* AIn;
  if (transpose) {
    AIn = new InputNode("A", n, m,
			1, n,
			dataType);
  } else {
    AIn = new InputNode("A", m, n,
			1, m,
			dataType);
  }

  auto alphaIn = new InputNode("Alpha",
				1, 1,
				1, 1,
				dataType);

  auto betaIn = new InputNode("Beta", 1, 1,
			       1, 1,
			       dataType);

  auto tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  auto tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  auto tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(AIn, 0);

  auto tunV = new Tunnel(POSSTUNIN);
  tunV->AddInput(vIn, 0);

  auto tunAlpha = new Tunnel(POSSTUNIN);
  tunAlpha->AddInput(alphaIn, 0);

  auto tunBeta = new Tunnel(POSSTUNIN);
  tunBeta->AddInput(betaIn, 0);

  auto by = new SVMul(COLVECTOR, ABSLAYER);
  by->AddInputs(4,
		tunBeta, 0,
		tunY, 0);

  auto trans = new LLDLATranspose(ABSLAYER);
  auto axMul = new MVMul(ABSLAYER);
  if (transpose) {
    trans->AddInputs(2,
		     tunA, 0);
    axMul->AddInputs(6,
		     trans, 0,
		     tunX, 0,
		     tunV, 0);
  } else {
    axMul->AddInputs(6,
		     tunA, 0,
		     tunX, 0,
		     tunV, 0);
    delete trans;
  }

  auto alphaAXMul = new SVMul(COLVECTOR, ABSLAYER);
  alphaAXMul->AddInputs(4,
			tunAlpha, 0,
			axMul, 0);

  auto sumVecs = new VAdd(COLVECTOR, ABSLAYER);
  sumVecs->AddInputs(4,
		     alphaAXMul, 0,
		     by, 0);

  auto innerPoss = new Poss(sumVecs, true);
  auto innerSet = new RealPSet(innerPoss);

  auto Cout = new OutputNode;
  Cout->AddInput(innerSet->OutTun(0), 0);

  auto outerPoss = new Poss(Cout, true);
  auto outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

#endif //DOLLDLA
