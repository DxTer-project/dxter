#include "blasExamples.h"

#if DOLLDLA

RealPSet* Axpy(Type dataType, VecType vType, int m)
{
  auto* alphaIn = new InputNode("alpha input", 1, 1, "alpha",
				1, m,
				"alphaNumRows", "alphaNumCols",
				"alphaRowStride", "alphaColStride",
				dataType);

  int nCols, nRows;
  if (vType == COLVECTOR) {
    nCols = 1;
    nRows = m;
  } else {
    nCols = m;
    nRows = 1;
  }

  auto* xIn = new InputNode("x input", nRows, nCols, "x",
			    nCols, nRows,
			    "xNumRows", "xNumCols",
			    "xRowStride", "xColStride",
			    dataType);

  auto* yIn = new InputNode("y input", nRows, nCols, "y",
			    nCols, nRows,
			    "yNumRows", "yNumCols",
			    "yRowStride", "yColStride",
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

  auto* Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  auto* outerPoss = new Poss(Cout, true);
  auto* outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* Gemv(Type dataType, bool transpose, int m, int n)
{
  auto* xIn = new InputNode("x input", n, 1, "X",
				 1, n,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  auto* yIn = new InputNode("y input", m, 1, "Y",
				 1, m,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

  auto* zIn = new InputNode("z input", m, 1, "Z",
				 1, m,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride", dataType);


  InputNode* AIn;
  if (transpose) {
    AIn = new InputNode("a input", n, m, "A",
			1, n,
			"ANumRows", "ANumCols",
			"ARowStride", "AColStride", dataType);
  } else {
    AIn = new InputNode("a input", m, n, "A",
			1, m,
			"ANumRows", "ANumCols",
			"ARowStride", "AColStride", dataType);
  }

  auto* alphaIn = new InputNode("alpha input", 1, 1, "Alpha",
				     1, m,
				     "AlphaNumRows", "AlphaNumCols",
				     "AlphaRowStride", "AlphaColStride", dataType);

  auto* betaIn = new InputNode("beta input", 1, 1, "Beta",
				    1, m,
				    "BetaNumRows", "BetaNumCols",
				    "BetaRowStride", "BetaColStride", dataType);

  auto* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  auto* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  auto* tunZ = new Tunnel(POSSTUNIN);
  tunZ->AddInput(zIn, 0);

  auto* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(AIn, 0);

  auto* tunAlpha = new Tunnel(POSSTUNIN);
  tunAlpha->AddInput(alphaIn, 0);

  auto* tunBeta = new Tunnel(POSSTUNIN);
  tunBeta->AddInput(betaIn, 0);

  auto* by = new SVMul(COLVECTOR, ABSLAYER);
  by->AddInputs(4,
		tunBeta, 0,
		tunY, 0);

  auto* trans = new LLDLATranspose(ABSLAYER);
  auto* axMul = new MVMul(ABSLAYER);
  if (transpose) {
    trans->AddInputs(2,
		     tunA, 0);
    axMul->AddInputs(6,
		     trans, 0,
		     tunX, 0,
		     tunZ, 0);
  } else {
    axMul->AddInputs(6,
		     tunA, 0,
		     tunX, 0,
		     tunZ, 0);
  }

  auto* alphaAXMul = new SVMul(COLVECTOR, ABSLAYER);
  alphaAXMul->AddInputs(4,
			tunAlpha, 0,
			axMul, 0);

  auto* sumVecs = new VAdd(COLVECTOR, ABSLAYER);
  sumVecs->AddInputs(4,
		     alphaAXMul, 0,
		     by, 0);

  auto* innerPoss = new Poss(sumVecs, true);
  auto* innerSet = new RealPSet(innerPoss);

  auto* Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  auto* outerPoss = new Poss(Cout, true);
  auto* outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

#endif //DOLLDLA
