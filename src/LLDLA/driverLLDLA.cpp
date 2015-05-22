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

#include "benchmarkMenu.h"
#include "blasExamples.h"
#include "driverMenu.h"
#include "driverSettings.h"
#include "driverUtils.h"
#include "miscellaneousExamples.h"
#include "multiBLASExamples.h"
#include "nodeTestExamples.h"
#include "problemRunner.h"
#include "singleOperationExamples.h"
#include "testSuites.h"

#if DOLLDLA

Trans transA, transB;

int main(int argc, const char* argv[])
{
#ifdef _OPENMP
  omp_set_num_threads(1);
  omp_set_nested(true);
#endif

  SetUpGlobalState();

  int m, n, p, k;
  Type precision;
  VecType vecType;
  ProblemInstance problemInstance;

  RealPSet* algPSet;
  int algNum;
  string opName;

  if (argc == 2 && *argv[1] == '0') {
    BenchmarkMenu();
  } else if(argc < 2) {
    PrintMainMenu();
    TearDownGlobalState(); return 0;
  } else {
    algNum = atoi(argv[1]);
    switch(algNum) {
    case(1):
      if (argc != 8) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
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
	PrintMainMenu();
	TearDownGlobalState(); return 0;
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
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_dot";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      problemInstance.AddDimension(m, "m");
      algPSet = DotTest(precision, m);
      break;
    case(4):
      if (argc != 5) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
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
	PrintMainMenu();
	TearDownGlobalState(); return 0;
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
	PrintMainMenu();
	TearDownGlobalState(); return 0;
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
	PrintMainMenu();
	TearDownGlobalState(); return 0;
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
	PrintMainMenu();
	TearDownGlobalState(); return 0;
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
	PrintMainMenu();
	TearDownGlobalState(); return 0;
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
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_vadd2";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      problemInstance.AddDimension(m, "m");
      algPSet = VAdd2(precision, m);
      break;
    case(11):
      if (argc != 5) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
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
	PrintMainMenu();
	TearDownGlobalState(); return 0;
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
	PrintMainMenu();
	TearDownGlobalState(); return 0;
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
	PrintMainMenu();
	TearDownGlobalState(); return 0;
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
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_sv_col_mul_gen";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      problemInstance.AddDimension(m, "m");
      algPSet = GenSizeColSVMul(precision, m);
      break;
    case(16):
      if (argc != 5) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
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
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_sgemam";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      p = atoi(argv[5]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      problemInstance.AddDimension(p, "p");
      algPSet = LGenCompareL3(precision, m, n, p);
      break;
    case(18):
      if (argc != 5) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_Gesummv";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = Gesummv(precision, m, n);
      break;
    case(19):
      if (argc != 5) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_zeroMVMul";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = SetToZeroTest(precision, m, n);
      break;
    case(20):
      if (argc != 4) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_pack_test";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      problemInstance.AddDimension(m, "m");
      algPSet = PackTest(precision, m);
      break;
    case(21):
      if (argc != 5) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_copy_test";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = CopyTest(precision, m, n);
      break;
    case(22):
      if (argc != 4) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_vertical_partition_recombine_test";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      problemInstance.AddDimension(m, "m");
      algPSet = VerticalPartitionRecombineTest(precision, m);
      break;
    case(23):
      if (argc != 4) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_horizontal_partition_recombine_test";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      problemInstance.AddDimension(m, "m");
      algPSet = HorizontalPartitionRecombineTest(precision, m);
      break;
    case(24):
      if (argc != 4) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_vertical_refined_pack_test";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      problemInstance.AddDimension(m, "m");
      algPSet = VerticalRefinedPackTest(precision, m);
      break;
    case(25):
      if (argc != 4) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_vertical_pack_unpack_test";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      problemInstance.AddDimension(m, "m");
      algPSet = VerticalPackUnpackTest(precision, m);
      break;
    case(26):
      if (argc != 5) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_vertical_2D_pack_unpack_test";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = TwoDVerticalPackUnpackTest(precision, m, n);
      break;
    case(27):
      if (argc != 5) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_vertical_2D_unpack_test";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = TwoDVerticalUnpackTest(precision, m, n);
      break;
    case(28):
      if (argc != 5) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_horizontal_2D_unpack_test";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = TwoDHorizontalUnpackTest(precision, m, n);
      break;
    case(29):
      if (argc != 5) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_horizontal_2D_copy_test";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = HorizontalCopyTest(precision, m, n);
      break;
    case(30):
      if (argc != 2) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      BasicNoRuntimeEvalTests();
      TearDownGlobalState(); return 0;
    case(31):
      if (argc != 4) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_vmuladd_benchmark";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      problemInstance.AddDimension(m, "m");
      algPSet = LGenCompareL1(precision, m);
      break;
    case(32):
      if (argc != 5) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "lgen_comparison_l2";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = LGenCompareL2(precision, m, n);
      break;
    case(33):
      if (argc != 6) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_mmul";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      p = atoi(argv[5]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      problemInstance.AddDimension(p, "p");
      algPSet = MMMulTest(precision, m, n, p);
      break;
    case(34):
      if (argc != 5) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_basic_multi_assign";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = BasicMultiAssign(precision, m, n);
      break;
    case(35):
      if (argc != 5) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "dxt_trsml";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      problemInstance.AddDimension(m, "m");
      problemInstance.AddDimension(n, "n");
      algPSet = TRSMLTest(precision, m, n);
      break;
    case(36):
      if (argc != 4) {
	PrintMainMenu();
	TearDownGlobalState(); return 0;
      }
      opName = "scal_neg";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      problemInstance.AddDimension(m, "m");
      algPSet = NegateVector(precision, m);
      break;
    default:
      PrintMainMenu();
      TearDownGlobalState();
      return 0;
    }

    problemInstance.SetType(precision);
    problemInstance.SetName(opName);
    ProblemInstanceStats* stats = RunProblemWithRTE(algNum, algPSet, &problemInstance);
    delete stats;
  }

  TearDownGlobalState();
  return 0;
}

#endif //DOLLDLA
