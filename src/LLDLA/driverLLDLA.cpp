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
#include "driverUtils.h"
#include "debug.h"
#include "LLDLAGemmTransformations.h"
#include "smmul.h"
#include "svmul.h"
#include "runtimeEvaluation.h"
#include "loopUnrolling.h"

#define DOEMPIRICALEVAL 1
#define PRINTCOSTS 1

#define DOLOOPUNROLLING 1
#define DO2MUTRANSFORMATIONS 1
#define DO3MUTRANSFORMATIONS 1
#define DO16MUTRANSFORMATIONS 1
#define DOLARGEMUTRANSFORMATIONS 0

#include <sstream>

Size one = 1;
Size smallSize = 4;
Size medSize = 32;
Size bigSize = 1024;
//Size bs = ELEM_BS;

RealPSet* GemvExample();
RealPSet* MVMul2Example();
RealPSet* MAdd2Example();
RealPSet* VAdd2Example();
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

Architecture* arch;
Type dataType = REAL_DOUBLE;

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
    cout << "IMPLEMENTATION # " << std::to_string((long long int) mit->first) << endl;
    for (vit = mit->second.begin(); vit != mit->second.end(); ++vit) {
      cout << std::to_string((long double) *vit) << endl;
    }
    cout << endl;
  }
}

GraphNum PrintImpMapInFlops(Type type, ImplementationRuntimeMap &impTimes, double flopCost, int chunkSize) {
  double peakFLOPS = arch->FlopsPerCycle(type) * arch->CyclesPerSecond();
  GraphNum bestImpNum = 0;
  double bestFLOPS = 0;
  ImplementationRuntimeMapIter mit;
  for (mit = impTimes.begin(); mit != impTimes.end(); ++mit) {
    TimeVecIter vit;
    cout << "IMPLEMENTATION # " << std::to_string((long long int) mit->first) << endl;
    for (vit = mit->second.begin(); vit != mit->second.end(); ++vit) {
      double totalFlops = flopCost * chunkSize;
      double totalTimeInSecs = *vit;
      double actualFLOPS = totalFlops / totalTimeInSecs;
      double pctPeak = (actualFLOPS / peakFLOPS) * 100;
      if (actualFLOPS > bestFLOPS) {
	bestFLOPS = actualFLOPS;
	bestImpNum = mit->first;
      }
      cout << "GFLOPS = " << std::to_string((long double) actualFLOPS / 1.0e9) << "\t%Peak = " << std::to_string((long double) pctPeak) << endl;
      /*      if (pctPeak > 100) {
	cout << "pctPeak > 100\n";
	throw;
	}*/
    }
    cout << endl;
  }
  cout << "Best GFLOPS achieved: " << std::to_string((long double) bestFLOPS / 1.0e9) << endl;
  cout << "Best percent of peak: " << std::to_string((long double) (bestFLOPS / peakFLOPS) * 100) << endl;
  return bestImpNum;      
}

void AddGemmTrans()
{
    // Convert gemm into loop over mvmul
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmToMVMul(ABSLAYER, ABSLAYER, REAL_SINGLE), LLDLALOOPPHASE);

  // Transform gemm into loop over vmmuls
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmToVMMul(ABSLAYER, ABSLAYER, REAL_SINGLE), LLDLALOOPPHASE);

    // Convert gemm into loop over mvmul
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmToMVMul(ABSLAYER, ABSLAYER, REAL_DOUBLE), LLDLALOOPPHASE);

  // Transform gemm into loop over vmmuls
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmToVMMul(ABSLAYER, ABSLAYER, REAL_DOUBLE), LLDLALOOPPHASE);

  //Introduces loops in the m, n, and k dimensions, respectively
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);

#if DO2MUTRANSFORMATIONS
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLA2MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLA2MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLA2MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
#endif

#if DO3MUTRANSFORMATIONS
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLA3MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLA3MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLA3MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
#endif

  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLAMuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLAMuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLAMuSingle, REAL_DOUBLE), LLDLALOOPPHASE);

#if DO2MUTRANSFORMATIONS
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLA2MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLA2MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLA2MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
#endif

#if DO3MUTRANSFORMATIONS
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLA3MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLA3MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLA3MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
#endif

  return;
}

void AddVVDotTrans()
{
  //  Universe::AddTrans(VVDot::GetClass(), new VVDotToRegArith(ABSLAYER, ABSLAYER, REAL_SINGLE), LLDLALOOPPHASE);

  Universe::AddTrans(VVDot::GetClass(), new VVDotToRegArith(ABSLAYER, ABSLAYER, REAL_DOUBLE), LLDLALOOPPHASE);
  return;
}

void AddMAddTrans()
{
    Universe::AddTrans(MAdd::GetClass(), new MAddLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);

    Universe::AddTrans(MAdd::GetClass(), new MAddLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);

    Universe::AddTrans(MAdd::GetClass(), new MAddToRegArith(ABSLAYER, ABSLAYER, REAL_SINGLE), LLDLALOOPPHASE);

    Universe::AddTrans(MAdd::GetClass(), new MAddLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMuDouble, REAL_DOUBLE), LLDLALOOPPHASE);

    Universe::AddTrans(MAdd::GetClass(), new MAddLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuDouble, REAL_DOUBLE), LLDLALOOPPHASE);

    Universe::AddTrans(MAdd::GetClass(), new MAddToRegArith(ABSLAYER, ABSLAYER, REAL_DOUBLE), LLDLALOOPPHASE);

  return;
}

void AddMVMulTrans()
{
  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);

  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);

  Universe::AddTrans(MVMul::GetClass(), new MVMulToRegArith(ABSLAYER, ABSLAYER, REAL_SINGLE), LLDLALOOPPHASE);

  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMuDouble, REAL_DOUBLE), LLDLALOOPPHASE);

  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuDouble, REAL_DOUBLE), LLDLALOOPPHASE);

  Universe::AddTrans(MVMul::GetClass(), new MVMulToRegArith(ABSLAYER, ABSLAYER, REAL_DOUBLE), LLDLALOOPPHASE);

  return;
}

void AddSMMulTrans()
{
  //Introduces loops in the m and n dimension for SMMul
  Universe::AddTrans(SMMul::GetClass(), new SMulLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMuSingle), LLDLALOOPPHASE);

  Universe::AddTrans(SMMul::GetClass(), new SMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuDouble), LLDLALOOPPHASE);

  return;
}

void AddUnrollingTrans()
{

#if DOLOOPUNROLLING
#if DO2MUTRANSFORMATIONS
  Universe::AddTrans(SplitSingleIter::GetClass(),
		     new CompactlyUnrollLoop(2), LLDLALOOPUNROLLPHASE);
#endif // DO2MUTRANSFORMATIONS

#if DO3MUTRANSFORMATIONS
  Universe::AddTrans(SplitSingleIter::GetClass(), 
		     new CompactlyUnrollLoop(3), LLDLALOOPUNROLLPHASE);
#endif // DO3MUTRANSFORMATIONS

#if DO16MUTRANSFORMATIONS
  Universe::AddTrans(SplitSingleIter::GetClass(), 
		     new CompactlyUnrollLoop(16), LLDLALOOPUNROLLPHASE);
#endif // DO3MUTRANSFORMATIONS

#if DOLARGEMUTRANSFORMATIONS
  Universe::AddTrans(SplitSingleIter::GetClass(), 
		     new CompactlyUnrollLoop(bigSize / LLDLA_MU), LLDLALOOPUNROLLPHASE);
#endif // DOLARGEMUTRANSFORMATIONS

#endif // DOLOOPUNROLLING

  return;
}

void AddSVMulTrans()
{
  Universe::AddTrans(SVMul::GetClass(), new SVMulToRegArith(ABSLAYER, ABSLAYER, ROWVECTOR, REAL_SINGLE), LLDLALOOPPHASE);

  Universe::AddTrans(SVMul::GetClass(), new SVMulToRegArith(ABSLAYER, ABSLAYER, COLVECTOR, REAL_SINGLE), LLDLALOOPPHASE);

  Universe::AddTrans(SVMul::GetClass(), new SVMulToRegArith(ABSLAYER, ABSLAYER, ROWVECTOR, REAL_DOUBLE), LLDLALOOPPHASE);

  Universe::AddTrans(SVMul::GetClass(), new SVMulToRegArith(ABSLAYER, ABSLAYER, COLVECTOR, REAL_DOUBLE), LLDLALOOPPHASE);

  return;
}

void AddVMMulTrans()
{
  // Transformers for vector matrix multiply
  Universe::AddTrans(VMMul::GetClass(), new VMMulToRegArith(ABSLAYER, ABSLAYER, REAL_SINGLE), LLDLALOOPPHASE);

  Universe::AddTrans(VMMul::GetClass(), new VMMulLoopRef(ABSLAYER, ABSLAYER, DIMK, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);

 Universe::AddTrans(VMMul::GetClass(), new VMMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);

  Universe::AddTrans(VMMul::GetClass(), new VMMulToRegArith(ABSLAYER, ABSLAYER, REAL_DOUBLE), LLDLALOOPPHASE);

  Universe::AddTrans(VMMul::GetClass(), new VMMulLoopRef(ABSLAYER, ABSLAYER, DIMK, LLDLAMuDouble, REAL_DOUBLE), LLDLALOOPPHASE);

  Universe::AddTrans(VMMul::GetClass(), new VMMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuDouble, REAL_DOUBLE), LLDLALOOPPHASE);

  return;
}

void AddVAddTrans()
{
  Universe::AddTrans(VAdd::GetClass(), new VAddToRegArith(ABSLAYER, ABSLAYER, REAL_SINGLE), LLDLALOOPPHASE);

  Universe::AddTrans(VAdd::GetClass(), new VAddToRegArith(ABSLAYER, ABSLAYER, REAL_DOUBLE), LLDLALOOPPHASE);

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
  cout <<"        11  -> Vector add twice\n";
  cout <<"        12  -> Vector matrix vector multiply\n";
  cout <<"        13  -> Matrix add twice\n";
  cout <<"        14  -> Matrix vector multiply twice\n";
  cout <<"        15  -> Gemv\n";
}

int main(int argc, const char* argv[])
{
#ifdef _OPENMP
  omp_set_num_threads(1);
  omp_set_nested(true);
#endif

  arch = new AMDEngSample();
  //  PrintType printType = CODE;
  int numIters = -1;
  RealPSet* (*algFunc)();
  //  GraphNum whichGraph = 0;
  int algNum;
  string fileName;
  string opName;
  Cost flopCost = 0;
  printf("dataType == REAL_DOUBLE ? %d\n", dataType == REAL_DOUBLE);
  printf("dataType == REAL_SINGLE ? %d\n", dataType == REAL_SINGLE);

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
      opName = "dxt_vadd2";
      algFunc = VAdd2Example;
      break;
    case(12):
      if (argc != 2) {
	Usage();
	return 0;
      }
      opName = "dxt_vmvmul";
      algFunc = VMVMulExample;
      break;
    case(13):
      if (argc != 2) {
	Usage();
	return 0;
      }
      opName = "dxt_madd2";
      algFunc = MAdd2Example;
      break;
    case(14):
      if (argc != 2) {
	Usage();
	return 0;
      }
      opName = "dxt_mvmul2";
      algFunc = MVMul2Example;
      break;
    case(15):
      if (argc != 2) {
	Usage();
	return 0;
      }
      opName = "dxt_gemv";
      algFunc = GemvExample;
      break;
    default:
      Usage();
      return 0;
    }
  }

  RegAllLLDLANodes();
  AddTrans();

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
    cout << "Flops for operation = " << std::to_string((long double) flopCost) << endl;
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

  int chunkSize = 300;
  int numIterations = 1;
  RuntimeTest rtest(dataType, opName, uni.m_argNames, uni.m_declarationVectors, uni.m_constantDefines, numIterations, chunkSize);
  string evalDirName = "runtimeEvaluation";
  RuntimeEvaluator evaler = RuntimeEvaluator(evalDirName);
  cout << "About to evaluate\n";
  ImplementationRuntimeMap impMap = evaler.EvaluateImplementationsWithCorrectnessCheck(rtest, ImpStrMap(&uni), absImpStr);
  cout << "Done evaluating\n";
  GraphNum best = PrintImpMapInFlops(dataType, impMap, flopCost, chunkSize);
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

RealPSet* GemvExample()
{
  InputNode* xIn = new InputNode("x input", 4, 1, "X",
				 1, 4,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride");

  InputNode* yIn = new InputNode("y input", medSize, 1, "Y",
				 1, medSize,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride");

  InputNode* zIn = new InputNode("z input", medSize, 1, "Z",
				 1, medSize,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride");

  InputNode* AIn = new InputNode("a input", medSize, 4, "A",
				 1, medSize,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride");

  InputNode* alphaIn = new InputNode("alpha input", 1, 1, "Alpha",
				     1, medSize,
				     "AlphaNumRows", "AlphaNumCols",
				     "AlphaRowStride", "AlphaColStride");

  InputNode* betaIn = new InputNode("beta input", 1, 1, "Beta",
				    1, medSize,
				    "BetaNumRows", "BetaNumCols",
				    "BetaRowStride", "BetaColStride");

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

  SVMul* by = new SVMul(COLVECTOR, ABSLAYER, dataType);
  by->AddInputs(4,
		tunBeta, 0,
		tunY, 0);

  MVMul* axMul = new MVMul(ABSLAYER, dataType);
  axMul->AddInputs(6,
		   tunA, 0,
		   tunX, 0,
		   tunZ, 0);

  SVMul* alphaAXMul = new SVMul(COLVECTOR, ABSLAYER, dataType);
  alphaAXMul->AddInputs(4,
			tunAlpha, 0,
			axMul, 0);

  VAdd* sumVecs = new VAdd(COLVECTOR, ABSLAYER, dataType);
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

RealPSet* MVMul2Example()
{
  InputNode* xIn = new InputNode("x input", medSize, 1, "X",
				 1, medSize,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride");

  InputNode* yIn = new InputNode("y input", medSize, 1, "Y",
				 1, medSize,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride");

  InputNode* zIn = new InputNode("z input", medSize, 1, "Z",
				 1, medSize,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride");

  InputNode* AIn = new InputNode("a input", medSize, medSize, "A",
				 1, medSize,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride");

  InputNode* BIn = new InputNode("b input", medSize, medSize, "B",
				 1, medSize,
				 "BNumRows", "BNumCols",
				 "BRowStride", "BColStride");
  
  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  Tunnel* tunZ = new Tunnel(POSSTUNIN);
  tunZ->AddInput(zIn, 0);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(AIn, 0);

  Tunnel* tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(BIn, 0);

  MVMul* mvmul1 = new MVMul(ABSLAYER, dataType);
  mvmul1->AddInputs(6,
		    tunB, 0,
		    tunX, 0,
		    tunY, 0);

  MVMul* mvmul2 = new MVMul(ABSLAYER, dataType);
  mvmul2->AddInputs(6,
		    tunA, 0,
		    mvmul1, 0,
		    tunZ, 0);

  Poss* innerPoss = new Poss(mvmul2, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* MAdd2Example()
{
  InputNode* xIn = new InputNode("x input", medSize, medSize, "X",
				 1, medSize,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride");

  InputNode* yIn = new InputNode("y input", medSize, medSize, "Y",
				 1, medSize,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride");

  InputNode* zIn = new InputNode("z input", medSize, medSize, "Z",
				 1, medSize,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride");
  
  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  Tunnel* tunZ = new Tunnel(POSSTUNIN);
  tunZ->AddInput(zIn, 0);

  MAdd* madd1 = new MAdd(ABSLAYER, dataType);
  madd1->AddInputs(4,
		  tunX, 0,
		  tunY, 0);

  MAdd* madd2 = new MAdd(ABSLAYER, dataType);
  madd2->AddInputs(4,
		   tunZ, 0,
		   madd1, 0);

  Poss* innerPoss = new Poss(madd2, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* VAdd2Example()
{
  InputNode* xIn = new InputNode("x input", bigSize, 1, "X",
				 1, bigSize,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride");

  InputNode* yIn = new InputNode("y input", bigSize, 1, "Y",
				 1, bigSize,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride");

  InputNode* zIn = new InputNode("z input", bigSize, 1, "Z",
				 1, bigSize,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride");
  
  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  Tunnel* tunZ = new Tunnel(POSSTUNIN);
  tunZ->AddInput(zIn, 0);

  VAdd* vadd1 = new VAdd(COLVECTOR, ABSLAYER, dataType);
  vadd1->AddInputs(4,
		  tunX, 0,
		  tunY, 0);

  VAdd* vadd2 = new VAdd(COLVECTOR, ABSLAYER, dataType);
  vadd2->AddInputs(4,
		   tunZ, 0,
		   vadd1, 0);

  Poss* innerPoss = new Poss(vadd2, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;

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

  VAdd* vadd = new VAdd(COLVECTOR, ABSLAYER, dataType);
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
  InputNode* Ain = new InputNode("A input", medSize, 8, "A",
				 1, 8,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride");

  InputNode* xIn = new InputNode("x input", 8, 1, "X",
				 1, 8,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride");

  InputNode* zIn = new InputNode("z input", medSize, 1, "Z",
				 1, medSize,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride");

  InputNode* yIn = new InputNode("y input", 1, medSize, "Y",
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

  MVMul* mvmul = new MVMul(ABSLAYER, dataType);
  mvmul->AddInputs(6,
		   tunA, 0,
		   tunX, 0,
		   tunZ, 0);

  VVDot* vvdot = new VVDot(ABSLAYER, dataType);
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

  SMMul* smmul = new SMMul(ABSLAYER, dataType);
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

  VMMul* vmmul = new VMMul(ABSLAYER, dataType);
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

  SVMul* svmul = new SVMul(ROWVECTOR, ABSLAYER, dataType);
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

  SVMul* svmul = new SVMul(COLVECTOR, ABSLAYER, dataType);
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
  InputNode* Ain = new InputNode("A input", 16, medSize, "A",
				 1, 16,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride");
  InputNode* xIn = new InputNode("x input", medSize, 1, "X",
				 1, medSize,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride");
  InputNode* yIn = new InputNode("y input", 16, 1, "Y",
				 1, 16,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride");

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  MVMul* mvmul = new MVMul(ABSLAYER, dataType);
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
				 1, medSize,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride");
  InputNode* Bin = new InputNode("B input", medSize, medSize, "B", 
				 1, medSize,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride");
  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(Bin, 0);

  MAdd* madd = new MAdd(ABSLAYER, dataType);
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

  VVDot* dot = new VVDot(ABSLAYER, dataType);
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
  InputNode *Ain = new InputNode("A input", smallSize, bigSize, "A",
				 bigSize, 1,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride");
  InputNode *Bin = new InputNode("B input", bigSize, smallSize, "B",
				 smallSize, 1,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride");
  InputNode *Cin = new InputNode("C input", smallSize, smallSize, "C",
				 smallSize, 1,
				 "CNumRows","CNumCols",
				 "CRowStride","CColStride");

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

RealPSet* DoubleGemmExample()
{
  InputNode *Ain = new InputNode("A input",  8, bigSize, "A",
				 bigSize, 1,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride");
  InputNode *Bin = new InputNode("B input", bigSize, 8, "B",
				 8, 1,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride");
  InputNode *Cin = new InputNode("C input",  8, 8, "C",
				 8, 1,
				 "CNumRows","CNumCols",
				 "CRowStride","CColStride");

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

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Cout,true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}


#endif //DOLLDLA
