#include "blasExamples.h"

#if DOLLDLA

RealPSet* GemvExample(Type dataType, bool transpose, int m, int n)
{
  InputNode* xIn = new InputNode("x input", n, 1, "X",
				 1, n,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  InputNode* yIn = new InputNode("y input", m, 1, "Y",
				 1, m,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

  InputNode* zIn = new InputNode("z input", m, 1, "Z",
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

  InputNode* alphaIn = new InputNode("alpha input", 1, 1, "Alpha",
				     1, m,
				     "AlphaNumRows", "AlphaNumCols",
				     "AlphaRowStride", "AlphaColStride", dataType);

  InputNode* betaIn = new InputNode("beta input", 1, 1, "Beta",
				    1, m,
				    "BetaNumRows", "BetaNumCols",
				    "BetaRowStride", "BetaColStride", dataType);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  Tunnel* tunZ = new Tunnel(POSSTUNIN);
  tunZ->AddInput(zIn, 0);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(AIn, 0);

  Tunnel* tunAlpha = new Tunnel(POSSTUNIN);
  tunAlpha->AddInput(alphaIn, 0);

  Tunnel* tunBeta = new Tunnel(POSSTUNIN);
  tunBeta->AddInput(betaIn, 0);

  SVMul* by = new SVMul(COLVECTOR, ABSLAYER);
  by->AddInputs(4,
		tunBeta, 0,
		tunY, 0);

  LLDLATranspose* trans = new LLDLATranspose(ABSLAYER);
  MVMul* axMul = new MVMul(ABSLAYER);
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

  SVMul* alphaAXMul = new SVMul(COLVECTOR, ABSLAYER);
  alphaAXMul->AddInputs(4,
			tunAlpha, 0,
			axMul, 0);

  VAdd* sumVecs = new VAdd(COLVECTOR, ABSLAYER);
  sumVecs->AddInputs(4,
		     alphaAXMul, 0,
		     by, 0);

  Poss* innerPoss = new Poss(sumVecs, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

#endif //DOLLDLA
