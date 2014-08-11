/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2014, The University of Texas and Bryan Marker

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

#include "base.h"
#include "costs.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#include "transform.h"
#include "loopSupport.h"
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <climits>

#if DOLLDLA

#include "LLDLAGemm.h"
#include "DLAReg.h"
#include "madd.h"
#include "vadd.h"
#include "vmmul.h"
#include "mvmul.h"
#include "vvdot.h"

#define DOEMPIRICALEVAL 1
#define PRINTCOSTS 1

#define DOLOOPUNROLLING 1
#define DO2MUTRANSFORMATIONS 1
#define DO3MUTRANSFORMATIONS 0

#include <sstream>

#include "driverUtils.h"
#include "debug.h"
#include "LLDLAGemmTransformations.h"
#include "smmul.h"
#include "svmul.h"
#include "runtimeEvaluation.h"
#include "loopUnrolling.h"

Size one = 1;
Size smallSize = 12;
Size medSize = 36;
Size bigSize = 1000;
//Size bs = ELEM_BS;

RealPSet* VAddExample();
RealPSet* VMVMulExample();
RealPSet* SMMulExample();
RealPSet* VMMulExample();
RealPSet* SVMulRowExample();
RealPSet* SVMulColExample();
RealPSet* MVMulExample();
RealPSet* MAddExample();
RealPSet* DotExample();
RealPSet* GemmExample();
RealPSet* DoubleGemmExample();

Trans transA, transB;

ImplementationMap ImpStrMap(Universe *uni)
{
  ImplementationMap impMap;
  GraphNum i;
  for (i = 1; i <= uni->TotalCount(); i++) {
    std::stringbuf sbuf;
    std::ostream out(&sbuf);
    IndStream istream = IndStream(&out, LLDLASTREAM);
    uni->Print(istream, i);
    impMap.insert(NumImplementationPair(i, sbuf.str()));
  }
  return impMap;
}

void PrintImpMap(ImplementationRuntimeMap &impTimes)
{
  ImplementationRuntimeMapIter mit;
  for (mit = impTimes.begin(); mit != impTimes.end(); ++mit) {
    TimeVecIter vit;
    cout << "IMPLEMENTATION # " << std::to_string(mit->first) << endl;
    for (vit = mit->second.begin(); vit != mit->second.end(); ++vit) {
      cout << std::to_string(*vit) << endl;
    }
    cout << endl;
  }
}

GraphNum PrintImpMapInFlops(ImplementationRuntimeMap &impTimes, double flopCost, int chunkSize) {
  /***************************************************************************
   * WARNING: These numbers are processor specific to Dillon's machine in GDC
   ***************************************************************************/
  double ticksPerSec = 1.0e6;
  double peakFLOPS = 2.7e9 * 8;//30e9;
  GraphNum bestImpNum = 0;
  double bestFLOPS = 0;
  ImplementationRuntimeMapIter mit;
  for (mit = impTimes.begin(); mit != impTimes.end(); ++mit) {
    TimeVecIter vit;
    cout << "IMPLEMENTATION # " << std::to_string(mit->first) << endl;
    for (vit = mit->second.begin(); vit != mit->second.end(); ++vit) {
      double totalFlops = flopCost * chunkSize;
      double totalTimeInSecs = *vit / ticksPerSec;
      double actualFLOPS = totalFlops / totalTimeInSecs;
      double pctPeak = (actualFLOPS / peakFLOPS) * 100;
      if (actualFLOPS > bestFLOPS) {
	bestFLOPS = actualFLOPS;
	bestImpNum = mit->first;
      }
      cout << "FLOPS = " << std::to_string(actualFLOPS) << "\t%Peak = " << std::to_string(pctPeak) << endl;
      /*      if (pctPeak > 100) {
	cout << "pctPeak > 100\n";
	throw;
	}*/
    }
    cout << endl;
  }
  cout << "Best flops achieved: " << std::to_string(bestFLOPS) << endl;
  cout << "Best percent of peak: " << std::to_string((bestFLOPS / peakFLOPS) * 100) << endl;
  return bestImpNum;      
}

void AddGemmTrans()
{
    // Convert gemm into loop over mvmul
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmToMVMul(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  // Transform gemm into loop over vmmuls
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmToVMMul(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  //Introduces loops in the m, n, and k dimensions, respectively
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLAMu), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLAMu), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLAMu), LLDLALOOPPHASE);

#if DO2MUTRANSFORMATIONS
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLA2Mu), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLA2Mu), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLA2Mu), LLDLALOOPPHASE);
#endif 

#if DO3MUTRANSFORMATIONS
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLA3Mu), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLA3Mu), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLA3Mu), LLDLALOOPPHASE);
#endif
  
  //Lowers the layer tag of a Gemm node that is USELLDLAMU in all three dimensions
  Universe::AddTrans(Gemm::GetClass(), new LLDAGemmLowerLayer(ABSLAYER, LLDLAMIDLAYER, LLDLAMu.GetSize()), LLDLALOOPPHASE);
  return;
}

void AddVVDotTrans()
{
  // Vector dot product transforms
  Universe::AddTrans(VVDot::GetClass(), new VVDotLowerLayer(ABSLAYER, LLDLAMIDLAYER, LLDLAMu.GetSize()), LLDLALOOPPHASE);
  
  Universe::AddTrans(VVDot::GetClass(), new VVDotLoopRef(ABSLAYER, ABSLAYER, LLDLAMu), LLDLALOOPPHASE);

  Universe::AddTrans(VVDot::GetClass(), new VVDotToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  return;
}

void AddMAddTrans()
{
  // Transformers for Matrix Matrix add
  Universe::AddTrans(MAdd::GetClass(), new MAddLowerLayer(ABSLAYER, LLDLAMIDLAYER, LLDLAMu.GetSize()), LLDLALOOPPHASE);

  // Introduce loop in M dimension
  Universe::AddTrans(MAdd::GetClass(), new MAddLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMu), LLDLALOOPPHASE);

  // Introduce loop in N dimension
  Universe::AddTrans(MAdd::GetClass(), new MAddLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMu), LLDLALOOPPHASE);

  // Convert to register arithmetic
  Universe::AddTrans(MAdd::GetClass(), new MAddToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  return;
}

void AddMVMulTrans()
{
  // Transformers for matrix vector multiply
  // Introduce loop in M dimension
  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMu), LLDLALOOPPHASE);

  // Introduce loop in N dimension
  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMu), LLDLALOOPPHASE);

  // Convert mvmul to vector arithmetic
  Universe::AddTrans(MVMul::GetClass(), new MVMulToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  // Lower layer tag
  Universe::AddTrans(MVMul::GetClass(), new MVMulLowerLayer(ABSLAYER, LLDLAMIDLAYER, LLDLAMu.GetSize()), LLDLALOOPPHASE);

  return;
}

void AddSMMulTrans()
{
  //Introduces loops in the m and n dimension for SMMul
  Universe::AddTrans(SMMul::GetClass(), new SMulLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMu), LLDLALOOPPHASE);

  Universe::AddTrans(SMMul::GetClass(), new SMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMu), LLDLALOOPPHASE);

  //Lowers the layer tag of a SMMul node that is USELLDLAMU in both dimensions
  Universe::AddTrans(SMMul::GetClass(), new SMulLowerLayer(ABSLAYER, LLDLAMIDLAYER, LLDLAMu.GetSize()), LLDLALOOPPHASE);

  return;
}

void AddUnrollingTrans()
{

#if DOLOOPUNROLLING
#if DO2MUTRANSFORMATIONS
  Universe::AddTrans(SplitSingleIter::GetClass(),
		     new FullyUnrollLoop(2), LLDLALOOPUNROLLPHASE);
#endif // DO2MUTRANSFORMATIONS

#if DO3MUTRANSFORMATIONS
  Universe::AddTrans(SplitSingleIter::GetClass(), 
		     new FullyUnrollLoop(3), LLDLALOOPUNROLLPHASE);
#endif // DO3MUTRANSFORMATIONS
#endif // DOLOOPUNROLLING

  return;
}

void AddSVMulTrans()
{
  // Transformers for scalar vector multiply
  Universe::AddTrans(SVMul::GetClass(), new SVMulLoopRef(ABSLAYER, ABSLAYER, COLVECTOR, LLDLAMu), LLDLALOOPPHASE);

  Universe::AddTrans(SVMul::GetClass(), new SVMulLoopRef(ABSLAYER, ABSLAYER, ROWVECTOR, LLDLAMu), LLDLALOOPPHASE);

  Universe::AddTrans(SVMul::GetClass(), new SVMulToRegArith(ABSLAYER, ABSLAYER, ROWVECTOR), LLDLALOOPPHASE);

  Universe::AddTrans(SVMul::GetClass(), new SVMulToRegArith(ABSLAYER, ABSLAYER, COLVECTOR), LLDLALOOPPHASE);

  Universe::AddTrans(SVMul::GetClass(), new SVMulLowerLayer(ABSLAYER, LLDLAMIDLAYER, LLDLAMu.GetSize()), LLDLALOOPPHASE);

  return;
}

void AddVMMulTrans()
{
  // Transformers for vector matrix multiply
  Universe::AddTrans(VMMul::GetClass(), new VMMulToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(VMMul::GetClass(), new VMMulLoopRef(ABSLAYER, ABSLAYER, DIMK, LLDLAMu), LLDLALOOPPHASE);

  Universe::AddTrans(VMMul::GetClass(), new VMMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMu), LLDLALOOPPHASE);

  Universe::AddTrans(VMMul::GetClass(), new VMMulLowerLayer(ABSLAYER, LLDLAMIDLAYER, LLDLAMu.GetSize()), LLDLALOOPPHASE);

  return;
}

void AddVAddTrans()
{
  Universe::AddTrans(VAdd::GetClass(), new VAddLoopRef(ABSLAYER, ABSLAYER, COLVECTOR, LLDLAMu), LLDLALOOPPHASE);

  Universe::AddTrans(VAdd::GetClass(), new VAddLoopRef(ABSLAYER, ABSLAYER, ROWVECTOR, LLDLAMu), LLDLALOOPPHASE);

  Universe::AddTrans(VAdd::GetClass(), new VAddLowerLayer(ABSLAYER, LLDLAMIDLAYER, LLDLAMu.GetSize()), LLDLALOOPPHASE);

  return;
}

void AddTrans()
{
  AddGemmTrans();
  AddVVDotTrans();
  AddMAddTrans();
  AddMVMulTrans();
  AddSMMulTrans();
  AddSVMulTrans();
  AddVMMulTrans();
  AddVAddTrans();

  AddUnrollingTrans();
  
}

void AddSimplifiers()
{ 
  //Replaces a Gemm node with a LLDLAGemm node
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmToPrim(LLDLAMIDLAYER, LLDLAPRIMITIVELAYER), SIMP);

  //Lowers the layer tag of a SMMul node that is USELLDLAMU in both dimensions
  Universe::AddTrans(SMMul::GetClass(), new SMulLowerLayer(LLDLAMIDLAYER, LLDLAPRIMITIVELAYER, LLDLAMu.GetSize()), SIMP);

  //Changes Gemm with transposition to non-transposed version use Transpose nodes
  Universe::AddTrans(Gemm::GetClass(), new GemmTransToNotTrans(LLDLAMIDLAYER), SIMP);

  // Lowers the layer tag of a VVDot node
  Universe::AddTrans(VVDot::GetClass(), new VVDotLowerLayer(LLDLAMIDLAYER, LLDLAPRIMITIVELAYER, LLDLAMu.GetSize()), SIMP);

  // Lowers the layer tag of a MAdd node
  Universe::AddTrans(MAdd::GetClass(), new MAddLowerLayer(LLDLAMIDLAYER, LLDLAPRIMITIVELAYER, LLDLAMu.GetSize()), SIMP);

  // Lowers the layer tag of a MVMul node
  Universe::AddTrans(MVMul::GetClass(), new MVMulLowerLayer(LLDLAMIDLAYER, LLDLAPRIMITIVELAYER, LLDLAMu.GetSize()), SIMP);

  // Lower layer tag on VMMul node
  Universe::AddTrans(VMMul::GetClass(), new VMMulLowerLayer(LLDLAMIDLAYER, LLDLAPRIMITIVELAYER, LLDLAMu.GetSize()), SIMP);

  // Lowers the tag of an SVMul node
  Universe::AddTrans(SVMul::GetClass(), new SVMulLowerLayer(LLDLAMIDLAYER, LLDLAPRIMITIVELAYER, LLDLAMu.GetSize()), SIMP);

  // Lowers the tag of a vadd node
  Universe::AddTrans(VAdd::GetClass(), new VAddLowerLayer(LLDLAMIDLAYER, LLDLAPRIMITIVELAYER, LLDLAMu.GetSize()), SIMP);
}

void Usage()
{
  cout << "./driver arg1 arg2 ...\n";
  cout <<" arg1 == 0  -> Load from file arg1\n";
  cout <<"         1  -> Gemm  N/T N/T\n";
  cout <<"         2  -> Double Gemm  N/T N/T\n";
  cout <<"         3  -> Dot prod\n";
  cout <<"         4  -> Matrix add\n";
  cout <<"         5  -> Matrix vector multiply\n";
  cout <<"         6  -> Scalar column vector multiply\n";
  cout <<"         7  -> Scalar row vector multiply\n";
  cout <<"         8  -> Vector matrix multiply\n";
  cout <<"         9  -> Scalar matrix multiply\n";
  cout <<"        10  -> Vector add\n";
  cout <<"        11  -> Vector matrix vector multiply\n";
}

int main(int argc, const char* argv[])
{
#ifdef _OPENMP
  omp_set_num_threads(1);
  omp_set_nested(true);
#endif

  //  PrintType printType = CODE;
  int numIters = -1;
  RealPSet* (*algFunc)();
  //  GraphNum whichGraph = 0;
  int algNum;
  string fileName;
  string opName;
  Cost flopCost = 0;

  if(argc < 2) {
    Usage();
    return 0;
  }
  else {
    algNum = atoi(argv[1]);
    switch(algNum) {
    case(1):
      if (argc != 4) {
	Usage();
	return 0;
      }
      opName = "dxt_gemm";
      algFunc = GemmExample;
      transA = CharToTrans(*argv[2]);
      transB = CharToTrans(*argv[3]);
      break;
    case(2):
      if (argc != 4) {
	Usage();
	return 0;
      }
      opName = "dxt_double_gemm";
      algFunc = DoubleGemmExample;
      transA = CharToTrans(*argv[2]);
      transB = CharToTrans(*argv[3]);
      break;
    case(3):
      if (argc != 2) {
	Usage();
	return 0;
      }
      opName = "dxt_dot";
      algFunc = DotExample;
      break;
    case(4):
      if (argc != 2) {
	Usage();
	return 0;
      }
      opName = "dxt_madd";
      algFunc = MAddExample;
      break;
    case(5):
      if (argc != 2) {
	Usage();
	return 0;
      }
      opName = "dxt_mvmul";
      algFunc = MVMulExample;
      break;
    case(6):
      if (argc != 2) {
	Usage();
	return 0;
      }
      opName = "dxt_sv_col_mul";
      algFunc = SVMulColExample;
      break;
    case(7):
      if (argc != 2) {
	Usage();
	return 0;
      }
      opName = "dxt_sv_row_mul";
      algFunc = SVMulRowExample;
      break;
    case(8):
      if (argc != 2) {
	Usage();
	return 0;
      }
      opName = "dxt_vmmul";
      algFunc = VMMulExample;
      break;
    case(9):
      if (argc != 2) {
	Usage();
	return 0;
      }
      opName = "dxt_smmul";
      algFunc = SMMulExample;
      break;
    case(10):
      if (argc != 2) {
	Usage();
	return 0;
      }
      opName = "dxt_vadd";
      algFunc = VAddExample;
      break;
    case(11):
      if (argc != 2) {
	Usage();
	return 0;
      }
      opName = "dxt_vmvmul";
      algFunc = VMVMulExample;
      break;
    default:
      Usage();
      return 0;
    }
  }

  RegAllLLDLANodes();
  AddTrans();
  AddSimplifiers();

  Universe uni;
  time_t start, start2, end;
  uni.PrintStats();
  string absImpStr;
  if (algNum==0) {
    time(&start);
    uni.Init(fileName);
    time(&end);
    cout << "Unflatten took " << difftime(end,start) << " seconds\n";
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
  else {
    cout << "Creating startSet\n";
    RealPSet *startSet = algFunc();
    cout << "Created startSet\n";
    uni.Init(startSet);
    cout << "Initialized universe\n";
    uni.Prop();
    GraphIter graphIter(startSet->m_posses.begin()->second);
    cout << "Printing evaluation code\n";
    flopCost = graphIter.EvalAndSetBest();
    // Print abstract implementation to string for use in testing
    // EXTREMELY HACKY, I could not figure out how to redirect an
    // ostream to a string
    std::stringstream ss;
    IndStream optOut(&ss, LLDLASTREAM);
    graphIter.PrintRoot(optOut, 0, true, startSet);
    absImpStr = ss.str();
    cout << "IMPLEMENTATION FOR CORRECTNESS CHECK:\n" << absImpStr;
    cout << "Flops for operation = " << std::to_string(flopCost) << endl;
    time(&start);
  }


#if DOLLDLALOOPPHASE
  if (CurrPhase == LLDLALOOPPHASE) {
    cout << "Expanding LLDLA loop phase\n";
    uni.Expand(-1, LLDLALOOPPHASE, LLDLACull);
    time(&end);
    cout << "LLDLALOOP phase took " << difftime(end,start) << " seconds\n";
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif


#if DOLLDLALOOPUNROLLPHASE
  if (CurrPhase == LLDLALOOPUNROLLPHASE) {
    cout << "LLDLALOOPUNROLL phase\n";
    uni.Expand(-1, LLDLALOOPUNROLLPHASE, LLDLACull);
    time(&end);
    cout << "LLDLALOOPUNROLL phase took " << difftime(end,start) << " seconds\n";
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

#if DOLLDLAPRIMPHASE
  if (CurrPhase == LLDLAPRIMPHASE) {
    cout << "Expanding LL DLA prim phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, LLDLAPRIMPHASE, LLDLACull);
    time(&end);
    cout << "LLDLAPRIM phase took " << difftime(end,start2) << " seconds\n";

    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

  cout << "Full expansion took " << difftime(end,start) << " seconds\n";
  cout.flush();

#if DOEMPIRICALEVAL  
  cout << "Writing all implementations to runtime eval files\n";

  int chunkSize = 3000;
  int numIterations = 1;
  RuntimeTest rtest(opName, uni.m_argNames, uni.m_declarationVectors, uni.m_constantDefines, numIterations, chunkSize);
  string evalDirName = "runtimeEvaluation";
  RuntimeEvaluator evaler = RuntimeEvaluator(evalDirName);
  cout << "About to evaluate\n";
  ImplementationRuntimeMap impMap = evaler.EvaluateImplementationsWithCorrectnessCheck(rtest, ImpStrMap(&uni), absImpStr);
  cout << "Done evaluating\n";
  GraphNum best = PrintImpMapInFlops(impMap, flopCost, chunkSize);
  cout << "All implementations printed\n";
  cout << "Best times";

#endif //DOEMPIRICALEVAL

#if 1
  uni.PrintAll(algNum, best);
#else
  uni.PrintBest();
#endif

#if PRINTCOSTS  
  uni.PrintCosts(impMap);
#endif

  /*  if (whichGraph <= 0)
      uni.PrintAll();
      else
      uni.Print(cout, CODE, whichGraph); */

  return 0;
}

RealPSet* VAddExample()
{
  InputNode* xIn = new InputNode("x input", bigSize, 1, "X",
				 1, bigSize,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride");

  InputNode* yIn = new InputNode("y input", bigSize, 1, "Y",
				 1, bigSize,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride");
  
  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  VAdd* vadd = new VAdd(COLVECTOR, ABSLAYER, REAL);
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

RealPSet* VMVMulExample()
{
  InputNode* Ain = new InputNode("A input", bigSize, 4, "A",
				 1, 4,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride");

  InputNode* xIn = new InputNode("x input", 4, 1, "X",
				 1, 4,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride");

  InputNode* zIn = new InputNode("z input", bigSize, 1, "Z",
				 1, bigSize,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride");

  InputNode* yIn = new InputNode("y input", 1, bigSize, "Y",
				 1, 1,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride");

  InputNode* wIn = new InputNode("w input", 1, 1, "W",
				 1, 1,
				 "WNumRows", "WNumCols",
				 "WRowStride", "WColStride");

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

  MVMul* mvmul = new MVMul(ABSLAYER, REAL);
  mvmul->AddInputs(6,
		   tunA, 0,
		   tunX, 0,
		   tunZ, 0);

  VVDot* vvdot = new VVDot(ABSLAYER, REAL);
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

RealPSet* SMMulExample()
{
  InputNode* Ain = new InputNode("A input", medSize, medSize, "A",
				 medSize, 1,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride");
  InputNode* xIn = new InputNode("x input", 1, 1, "X",
				 medSize, 1,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride");

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  SMMul* smmul = new SMMul(ABSLAYER, REAL);
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

RealPSet* VMMulExample()
{
  InputNode* Ain = new InputNode("A input", medSize, medSize, "A",
				 medSize, 1,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride");
  InputNode* xIn = new InputNode("x input", 1, medSize, "X",
				 medSize, 1,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride");
  InputNode* yIn = new InputNode("y input", 1, medSize, "Y",
				 medSize, 1,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride");

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  VMMul* vmmul = new VMMul(ABSLAYER, REAL);
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

RealPSet* SVMulRowExample()
{
  InputNode* Ain = new InputNode("A input", 1, medSize, "A",
				 medSize, 1,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride");
  InputNode* xIn = new InputNode("x input", 1, 1, "X",
				 medSize, 1,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride");

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  SVMul* svmul = new SVMul(ROWVECTOR, ABSLAYER, REAL);
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

RealPSet* SVMulColExample()
{
  InputNode* Ain = new InputNode("A input", medSize, 1, "A",
				 medSize, 1,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride");
  InputNode* xIn = new InputNode("x input", 1, 1, "X",
				 medSize, 1,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride");

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  SVMul* svmul = new SVMul(COLVECTOR, ABSLAYER, REAL);
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

RealPSet* MVMulExample()
{
  InputNode* Ain = new InputNode("A input", 4, bigSize, "A",
				 1, 4,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride");
  InputNode* xIn = new InputNode("x input", bigSize, 1, "X",
				 1, bigSize,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride");
  InputNode* yIn = new InputNode("y input", 4, 1, "Y",
				 1, 4,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride");

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  MVMul* mvmul = new MVMul(ABSLAYER, REAL);
  mvmul->AddInputs(6,
		   tunA, 0,
		   tunX, 0,
		   tunY, 0);

  Poss *innerPoss = new Poss(mvmul, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* MAddExample()
{
  InputNode* Ain = new InputNode("A input", medSize, medSize, "A", 
				 medSize, 1,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride");
  InputNode* Bin = new InputNode("B input", medSize, medSize, "B", 
				 medSize, 1,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride");
  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(Bin, 0);

  MAdd* madd = new MAdd(ABSLAYER, REAL);
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

RealPSet* DotExample()
{
  InputNode* Ain = new InputNode("A input", 1, medSize, "A", 
				 medSize, 1,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride");
  InputNode* Bin = new InputNode("B input", medSize, 1, "B", 
				 medSize, 1,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride");
  InputNode* Cin = new InputNode("C input", 1, 1, "C", 
				 medSize, 1,
				 "CNumRows","CNumCols",
				 "CRowStride","CColStride");

  Tunnel *tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain,0);

  Tunnel *tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(Bin,0);

  Tunnel *tunC = new Tunnel(POSSTUNIN);
  tunC->AddInput(Cin,0);

  VVDot* dot = new VVDot(ABSLAYER, REAL);
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

RealPSet* GemmExample()
{
  InputNode *Ain = new InputNode("A input", 4, 4, "A",
				 4, 1,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride");
  InputNode *Bin = new InputNode("B input", 4, 352, "B",
				 352, 1,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride");
  InputNode *Cin = new InputNode("C input", 4, 352, "C",
				 352, 1,
				 "CNumRows","CNumCols",
				 "CRowStride","CColStride");

  Tunnel *tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel *tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(Bin, 0);

  Tunnel *tunC = new Tunnel(POSSTUNIN);
  tunC->AddInput(Cin, 0);

  Gemm *gemm = new Gemm(ABSLAYER, transA, transB, COEFONE, COEFONE, REAL);
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

RealPSet* DoubleGemmExample()
{
  InputNode *Ain = new InputNode("A input",  4, bigSize, "A",
				 bigSize, 1,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride");
  InputNode *Bin = new InputNode("B input", bigSize, 4, "B",
				 4, 1,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride");
  InputNode *Cin = new InputNode("C input",  4, 4, "C",
				 4, 1,
				 "CNumRows","CNumCols",
				 "CRowStride","CColStride");

  Tunnel *tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain,0);

  Tunnel *tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(Bin,0);

  Tunnel *tunC = new Tunnel(POSSTUNIN);
  tunC->AddInput(Cin,0);

  Gemm *gemm1 = new Gemm(ABSLAYER, transA, transB, COEFONE, COEFONE, REAL);
  gemm1->AddInputs(6,
		  tunA,0,
		  tunB,0,
		  tunC,0);

  Gemm *gemm2 = new Gemm(ABSLAYER, transA, transB, COEFONE, COEFONE, REAL);
  gemm2->AddInputs(6,
		  tunA,0,
		  tunB,0,
		  gemm1,0);

  Poss *innerPoss = new Poss(gemm2,true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Cout,true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}


#endif //DOLLDLA
