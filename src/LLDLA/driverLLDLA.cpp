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

#include "allTransformations.h"
#include "base.h"
#include "blasExamples.h"
#include "costs.h"
#include "driverUtils.h"
#include "debug.h"
#include "DLAReg.h"
#include "loopUnrolling.h"
#include "LLDLAGemm.h"
#include "LLDLAGemmTransformations.h"
#include "miscellaneousExamples.h"
#include "multiBLASExamples.h"
#include "runtimeEvaluation.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#include "singleOperationExamples.h"
#include "transform.h"
#include "LLDLATranspose.h"
#include "transpose.h"
#include "loopSupport.h"
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <climits>

#if DOLLDLA

#define DOEMPIRICALEVAL 1
#define PRINTCOSTS 1

#define DOCOMPACTLOOPUNROLLING 0
#define DO2MUTRANSFORMATIONS 1
#define DO3MUTRANSFORMATIONS 1
#define DO16MUTRANSFORMATIONS 1
#define DOLARGEMUTRANSFORMATIONS 0

#define DOPARTIALLOOPUNROLLING 1
#define PARTIALUNROLLINGSTARTCOEF 2
#define PARTIALUNROLLINGENDCOEF 16

#if DOCOMPACTLOOPUNROLLING + DOPARTIALLOOPUNROLLING > 1
do you really want to do compact unrolling and partial unrolling?
#endif

#include <sstream>

void MuNNMuGemmResults(Type precision);
double RunExample(int algNum, RealPSet* algPSet, Type precision, string opName);

Trans transA, transB;

Architecture* arch;

ImplementationMap* ImpStrMap(Universe *uni)
{
  ImplementationMap* impMap = new ImplementationMap();
  GraphNum i;
  for (i = 1; i <= uni->TotalCount(); i++) {
    std::stringbuf sbuf;
    std::ostream out(&sbuf);
    IndStream istream = IndStream(&out, LLDLASTREAM);
    uni->Print(istream, i);
    impMap->insert(NumImplementationPair(i, sbuf.str()));
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

double BestFlopsPerCycle(Type type, ImplementationRuntimeMap &impTimes, double flopCost) {
  double peakFlopsPerCycle = arch->FlopsPerCycle(type);
  double bestFlopsPerCycle = 0;
  ImplementationRuntimeMapIter mit;
  for (mit = impTimes.begin(); mit != impTimes.end(); ++mit) {
    TimeVecIter vit;
    for (vit = mit->second.begin(); vit != mit->second.end(); ++vit) {
      double totalTimeInCycles = *vit;
      double actualFlopsPerCycle = flopCost / totalTimeInCycles;
      double pctPeak = (actualFlopsPerCycle / peakFlopsPerCycle) * 100;
      if (actualFlopsPerCycle > bestFlopsPerCycle) {
	bestFlopsPerCycle = actualFlopsPerCycle;
      }
      if (pctPeak > 100) {
	cout << "pctPeak > 100\n";
	throw;
      }
    }
  }
  return bestFlopsPerCycle;
}

GraphNum PrintImpMapInFlops(Type type, ImplementationRuntimeMap &impTimes, double flopCost) {
  double peakFlopsPerCycle = arch->FlopsPerCycle(type);
  GraphNum bestImpNum = 0;
  double bestFlopsPerCycle = 0;
  ImplementationRuntimeMapIter mit;
  for (mit = impTimes.begin(); mit != impTimes.end(); ++mit) {
    TimeVecIter vit;
    cout << "IMPLEMENTATION # " << std::to_string((long long int) mit->first) << endl;
    for (vit = mit->second.begin(); vit != mit->second.end(); ++vit) {
      double totalTimeInCycles = *vit;
      double actualFlopsPerCycle = flopCost / totalTimeInCycles;
      double pctPeak = (actualFlopsPerCycle / peakFlopsPerCycle) * 100;
      if (actualFlopsPerCycle > bestFlopsPerCycle) {
	bestFlopsPerCycle = actualFlopsPerCycle;
	bestImpNum = mit->first;
      }
      cout << "Flops per cycle = " << std::to_string((long double) actualFlopsPerCycle);
      cout << "\t%Peak = " << std::to_string((long double) pctPeak) << endl;
      if (pctPeak > 100) {
	cout << "pctPeak > 100\n";
	throw;
      }
    }
    cout << endl;
  }
  cout << "Best flops/cycle achieved: " << std::to_string((long double) bestFlopsPerCycle) << endl;
  cout << "Best percent of peak: " << std::to_string((long double) (bestFlopsPerCycle / peakFlopsPerCycle) * 100) << endl;
  return bestImpNum;
}

void Usage()
{
  cout <<"\n";
  cout << "./driver arg1 arg2 ...\n";
  cout <<"\n";
  cout <<"arg1 == 0  -> Load from file arg1\n";
  cout <<"\n";
  cout <<"Single Operation Examples\n";
  cout <<"         3  -> Dot prod F/D M\n";
  cout <<"         4  -> Matrix add F/D M N\n";
  cout <<"         5  -> Matrix vector multiply N/T F/D M N\n";
  cout <<"         6  -> Scalar column vector multiply F/D M\n";
  cout <<"         7  -> Scalar row vector multiply F/D M\n";
  cout <<"         8  -> Vector matrix multiply F/D M N\n";
  cout <<"         9  -> Scalar matrix multiply F/D M N\n";
  cout <<"        10  -> Vector add C/R F/D M\n";
  cout <<"        16  -> Gen Size Col Vector SVMul F/D M\n";
  cout <<"\n";
  cout <<"BLAS Examples\n";
  cout <<"         1  -> Gemm  N/T N/T F/D M N P\n";
  cout <<"        15  -> Gemv N/T F/D M N\n";
  cout <<"        17  -> Axpy C/R F/D M\n";
  cout <<"\n";
  cout <<"Miscellaneous Examples\n";
  cout <<"         2  -> Double Gemm  N/T N/T F/D M N P K\n";
  cout <<"        11  -> Vector add twice F/D M\n";
  cout <<"        12  -> Vector matrix vector multiply F/D M N\n";
  cout <<"        13  -> Matrix add twice F/D M N\n";
  cout <<"        14  -> Matrix vector multiply twice F/D M N P\n";
  cout <<"        18  -> alpha*(A0 + A1)^T*B + beta*C F/D M N P\n";
  cout <<"        19  -> alpha*A*x + beta*B*x + y F/D M N\n";
  cout <<"\n";
}

int main(int argc, const char* argv[])
{
#ifdef _OPENMP
  omp_set_num_threads(1);
  omp_set_nested(true);
#endif

  int m, n, p, k;
  Type precision;
  VecType vecType;

  arch = new HaswellMacbook();

  //  PrintType printType = CODE;

  RealPSet* algPSet;
  //  GraphNum whichGraph = 0;
  int algNum;
  string opName;

  if(argc < 2) {
    Usage();
    return 0;
  }
  else {
    algNum = atoi(argv[1]);
    switch(algNum) {
    case(1):
      if (argc != 8) {
	Usage();
	return 0;
      }
      opName = "dxt_gemm";
      transA = CharToTrans(*argv[2]);
      transB = CharToTrans(*argv[3]);
      precision = CharToType(*argv[4]);
      m = atoi(argv[5]);
      n = atoi(argv[6]);
      p = atoi(argv[7]);
      algPSet = GemmExample(precision, transA, transB, m, n, p);
      break;
    case(2):
      if (argc != 9) {
	Usage();
	return 0;
      }
      opName = "dxt_double_gemm";
      precision = CharToType(*argv[4]);
      m = atoi(argv[5]);
      n = atoi(argv[6]);
      p = atoi(argv[7]);
      k = atoi(argv[8]);
      transA = CharToTrans(*argv[2]);
      transB = CharToTrans(*argv[3]);
      algPSet = DoubleGemmExample(precision, transA, transB, m, n, p, k);
      break;
    case(3):
      if (argc != 4) {
	Usage();
	return 0;
      }
      opName = "dxt_dot";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      algPSet = DotExample(precision, m);
      break;
    case(4):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_madd";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      algPSet = MAddExample(precision, m, n);
      break;
    case(5):
      if (argc != 6) {
	Usage();
	return 0;
      }
      opName = "dxt_mvmul";
      precision = CharToType(*argv[3]);
      m = atoi(argv[4]);
      n = atoi(argv[5]);
      if (TRANS == CharToTrans(*argv[2])) {
	algPSet = MVMulExample(precision, true, m, n);
      } else {
	algPSet = MVMulExample(precision, false, m, n);
      }
      break;
    case(6):
      if (argc != 4) {
	Usage();
	return 0;
      }
      opName = "dxt_sv_col_mul";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      algPSet = SVMulColExample(precision, m);
      break;
    case(7):
      if (argc != 4) {
	Usage();
	return 0;
      }
      precision = CharToType(*argv[2]);
      m = atoi(argv[4]);
      opName = "dxt_sv_row_mul";
      algPSet = SVMulRowExample(precision, m);
      break;
    case(8):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_vmmul";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      algPSet = VMMulExample(precision, m, n);
      break;
    case(9):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_smmul";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      algPSet = SMMulExample(precision, m, n);
      break;
    case(10):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_vadd";
      vecType = CharToVecType(*argv[2]);
      precision = CharToType(*argv[3]);
      m = atoi(argv[4]);
      algPSet = VAddExample(precision, vecType, m);
      break;
    case(11):
      if (argc != 4) {
	Usage();
	return 0;
      }
      opName = "dxt_vadd2";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      algPSet = VAdd2Example(precision, m);
      break;
    case(12):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_vmvmul";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      algPSet = VMVMulExample(precision, m, n);
      break;
    case(13):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_madd2";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      algPSet = MAdd2Example(precision, m, n);
      break;
    case(14):
      if (argc != 6) {
	Usage();
	return 0;
      }
      opName = "dxt_mvmul2";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      p = atoi(argv[5]);
      algPSet = MVMul2Example(precision, m, n, p);
      break;
    case(15):
      if (argc != 6) {
	Usage();
	return 0;
      }
      opName = "dxt_gemv";
      precision = CharToType(*argv[3]);
      m = atoi(argv[4]);
      n = atoi(argv[5]);
      if (TRANS == CharToTrans(*argv[2])) {
	algPSet = Gemv(precision, true, m, n);
      } else {
	algPSet = Gemv(precision, false, m, n);
      }
      break;
    case(16):
      if (argc != 4) {
	Usage();
	return 0;
      }
      opName = "dxt_sv_col_mul_gen";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      algPSet = GenSizeColSVMul(precision, m);
      break;
    case(17):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_saxpy";
      vecType = CharToVecType(*argv[2]);
      precision = CharToType(*argv[3]);
      m = atoi(argv[4]);
      algPSet = Axpy(precision, vecType, m);
      break;
    case(18):
      if (argc != 6) {
	Usage();
	return 0;
      }
      opName = "dxt_sgemam";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      p = atoi(argv[5]);
      algPSet = Gemam(precision, m, n, p);
      break;
    case(19):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_sgemam";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      algPSet = Gesummv(precision, m, n);
      break;
    default:
      Usage();
      return 0;
    }
  }

  RunExample(algNum, algPSet, precision, opName);
  return 0;
}

double RunExample(int algNum, RealPSet* algPSet, Type precision, string opName)
{
  RegAllLLDLANodes();
  AddTransformations();

  int numIters = -1;
  Cost flopCost = 0;
  Universe uni;
  time_t start, start2, end;
  string absImpStr;

  uni.PrintStats();

  cout << "Creating startSet\n";

  RealPSet *startSet = algPSet;
  
  cout << "Created startSet\n";

  uni.Init(startSet);
  
  cout << "Initialized universe\n";
  
  uni.Prop();
  GraphIter* graphIter = new GraphIter(startSet->m_posses.begin()->second);
  cout << "Printing evaluation code\n";
  flopCost = graphIter->EvalAndSetBest();
  // Print abstract implementation to string for use in testing
  // EXTREMELY HACKY, I could not figure out how to redirect an
  // ostream to a string
  std::stringstream ss;
  IndStream optOut(&ss, LLDLASTREAM);
  graphIter->PrintRoot(optOut, 0, true, startSet);
  absImpStr = ss.str();

  cout << "IMPLEMENTATION FOR CORRECTNESS CHECK:\n" << absImpStr;
  cout << "Flops for operation = " << std::to_string((long double) flopCost) << endl;

  time(&start);

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

  uni.CullWorstPerformers(0.50, 0);

#if DOEMPIRICALEVAL  
  cout << "Writing all implementations to runtime eval files\n";

  int numIterations = 10;
  RuntimeTest rtest(precision, opName, uni.m_argNames, uni.m_declarationVectors, uni.m_constantDefines, numIterations);
  string evalDirName = "runtimeEvaluation";
  RuntimeEvaluator evaler = RuntimeEvaluator(evalDirName);
  cout << "About to evaluate\n";
  ImplementationRuntimeMap impMap = evaler.EvaluateImplementationsWithCorrectnessCheck(rtest, ImpStrMap(&uni), absImpStr);

  cout << "Done evaluating\n";
  GraphNum best = PrintImpMapInFlops(precision, impMap, flopCost);
#endif //DOEMPIRICALEVAL

#if 1
  uni.PrintAll(algNum, best);
#else
  uni.PrintBest();
#endif

#if PRINTCOSTS
  uni.PrintCosts(impMap);
#endif

  double bestFPS = BestFlopsPerCycle(precision, impMap, flopCost);
  return bestFPS;
}

#endif //DOLLDLA
