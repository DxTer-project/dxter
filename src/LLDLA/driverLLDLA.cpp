/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2015, The University of Texas and Bryan Marker

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
#include "benchmarkMenu.h"
#include "blasExamples.h"
#include "costs.h"
#include "driverUtils.h"
#include "debug.h"
#include "DLAReg.h"
#include "exampleRunner.h"
#include "loopUnrolling.h"
#include "LLDLAGemm.h"
#include "LLDLAGemmTransformations.h"
#include "miscellaneousExamples.h"
#include "multiBLASExamples.h"
#include "problemInstanceStats.h"
#include "runtimeEvaluation.h"

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

Trans transA, transB;

Architecture* arch;

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

GraphNum PrintImpMapStats(Type type, ImplementationRuntimeMap &impTimes, double flopCost) {

  double peakFlopsPerCycle = arch->FlopsPerCycle(type);
  GraphNum bestImpNum = 0;
  double overallBestAvgFlopsPerCycle = 0;

  for (auto mit : impTimes) {
    double avgFlopsPerCycle = 0.0;
    double numRuns = 0.0;

    cout << "IMPLEMENTATION # " << std::to_string((long long int) mit.first) << endl;
    for (auto vit : mit.second) {
      double timeInCycles = vit;
      double actualFlopsPerCycle = flopCost / timeInCycles;
      avgFlopsPerCycle += actualFlopsPerCycle;
      numRuns += 1.0;
    }

    avgFlopsPerCycle = avgFlopsPerCycle / numRuns;

    if (avgFlopsPerCycle > overallBestAvgFlopsPerCycle) {
      overallBestAvgFlopsPerCycle = avgFlopsPerCycle;
      bestImpNum = mit.first;
    }

    cout << "Avg. Flops per cycle = " << std::to_string((long double) avgFlopsPerCycle);
    double avgPctPeak = (avgFlopsPerCycle / peakFlopsPerCycle) * 100;
    cout << "\t%Peak = " << std::to_string((long double) avgPctPeak) << endl;
  }
  cout << "Best avg. flops/cycle achieved: " << std::to_string((long double) overallBestAvgFlopsPerCycle) << endl;
  cout << "Best avg. percent of peak: " << std::to_string((long double) (overallBestAvgFlopsPerCycle / peakFlopsPerCycle) * 100) << endl;
  return bestImpNum;
}

void Usage()
{
  cout <<"\n";
  cout << "./driver arg1 arg2 ...\n";
  cout <<"\n";
  cout <<"arg1 == 0   -> View benchmarks\n";
  cout <<"\n";
  cout <<"Single Operation Examples\n";
  cout <<"         3  -> Dot prod F/D M\n";
  cout <<"         4  -> Matrix add F/D M N\n";
  cout <<"         5  -> Matrix vector multiply N/T F/D M N\n";
  cout <<"         6  -> Scalar vector multiply C/R F/D M\n";
  cout <<"         7  -> Vector matrix multiply F/D M N\n";
  cout <<"         8  -> Scalar matrix multiply F/D M N\n";
  cout <<"         9  -> Vector add C/R F/D M\n";
  cout <<"        15  -> Gen Size Col Vector SVMul F/D M\n";
  cout <<"\n";
  cout <<"BLAS Examples\n";
  cout <<"         1  -> Gemm  N/T N/T F/D M N P\n";
  cout <<"        14  -> Gemv N/T F/D M N\n";
  cout <<"        16  -> Axpy C/R F/D M\n";
  cout <<"\n";
  cout <<"Miscellaneous Examples\n";
  cout <<"         2  -> Double Gemm  N/T N/T F/D M N P K\n";
  cout <<"        10  -> Vector add twice F/D M\n";
  cout <<"        11  -> Vector matrix vector multiply F/D M N\n";
  cout <<"        12  -> Matrix add twice F/D M N\n";
  cout <<"        13  -> Matrix vector multiply twice F/D M N P\n";
  cout <<"        17  -> alpha*(A0 + A1)^T*B + beta*C F/D M N P\n";
  cout <<"        18  -> alpha*A*x + beta*B*x F/D M N\n";
  cout <<"        19  -> y <- Ax F/D M N\n";
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
  ProblemInstance problemInstance;

  arch = new HaswellMacbook();

  RealPSet* algPSet;
  int algNum;
  string opName;

  if (argc == 2 && *argv[1] == '0') {
    BenchmarkMenu();
  } else if(argc < 2) {
    Usage();
    return 0;
  } else {
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
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      problemInstance.AddDimension(p, "p");
      algPSet = GemmTest(precision, transA, transB, m, n, p);
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
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      problemInstance.AddDimension(p, "p");
      problemInstance.AddDimension(k, "k");
      algPSet = DoubleGemm(precision, transA, transB, m, n, p, k);
      break;
    case(3):
      if (argc != 4) {
	Usage();
	return 0;
      }
      opName = "dxt_dot";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      problemInstance.AddDimension(m, "m");
      algPSet = DotTest(precision, m);
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
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = MAddTest(precision, m, n);
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
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      if (TRANS == CharToTrans(*argv[2])) {
	algPSet = MVMulTest(precision, true, m, n);
      } else {
	algPSet = MVMulTest(precision, false, m, n);
      }
      break;
    case(6):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_sv_mul";
      vecType = CharToVecType(*argv[2]);
      precision = CharToType(*argv[3]);
      m = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      algPSet = SVMulTest(precision, vecType, m);
      break;
    case(7):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_vmmul";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = VMMulTest(precision, m, n);
      break;
    case(8):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_smmul";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = SMMulTest(precision, m, n);
      break;
    case(9):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_vadd";
      vecType = CharToVecType(*argv[2]);
      precision = CharToType(*argv[3]);
      m = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      algPSet = VAddTest(precision, vecType, m);
      break;
    case(10):
      if (argc != 4) {
	Usage();
	return 0;
      }
      opName = "dxt_vadd2";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      problemInstance.AddDimension(m, "m");
      algPSet = VAdd2(precision, m);
      break;
    case(11):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_vmvmul";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = VMVMul(precision, m, n);
      break;
    case(12):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_madd2";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = MAdd2(precision, m, n);
      break;
    case(13):
      if (argc != 6) {
	Usage();
	return 0;
      }
      opName = "dxt_mvmul2";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      p = atoi(argv[5]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      problemInstance.AddDimension(p, "p");
      algPSet = MVMul2(precision, m, n, p);
      break;
    case(14):
      if (argc != 6) {
	Usage();
	return 0;
      }
      opName = "dxt_gemv";
      precision = CharToType(*argv[3]);
      m = atoi(argv[4]);
      n = atoi(argv[5]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      if (TRANS == CharToTrans(*argv[2])) {
	algPSet = Gemv(precision, true, m, n);
      } else {
	algPSet = Gemv(precision, false, m, n);
      }
      break;
    case(15):
      if (argc != 4) {
	Usage();
	return 0;
      }
      opName = "dxt_sv_col_mul_gen";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      problemInstance.AddDimension(m, "m");
      algPSet = GenSizeColSVMul(precision, m);
      break;
    case(16):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_saxpy";
      vecType = CharToVecType(*argv[2]);
      precision = CharToType(*argv[3]);
      m = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      algPSet = Axpy(precision, vecType, m);
      break;
    case(17):
      if (argc != 6) {
	Usage();
	return 0;
      }
      opName = "dxt_sgemam";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      p = atoi(argv[5]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      problemInstance.AddDimension(p, "p");
      algPSet = Gemam(precision, m, n, p);
      break;
    case(18):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_sgemam";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = Gesummv(precision, m, n);
      break;
    case(19):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_zeroMVMul";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = SetToZeroTest(precision, m, n);
      break;
    default:
      Usage();
      return 0;
    }

    problemInstance.SetType(precision);
    problemInstance.SetName(opName);
    RunExample(algNum, algPSet, &problemInstance);
  }
  return 0;
}

#endif //DOLLDLA
