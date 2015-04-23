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

#include <iostream>

#include "benchmark.h"

#if DOLLDLA

#include <assert.h>

#include "benchmarkStats.h"
#include "miscellaneousExamples.h"
#include "problemInstance.h"
#include "problemRunner.h"
#include "singleOperationExamples.h"

using namespace std;

void LGenLevel1ComparisonM(Type type, int m) {
  vector<int> ms;
  for (int i = 1; i <= m; i++) {
    ms.push_back(i);
  }
  LGenLevel1Comparison(type, ms);
}

void LGenLevel1Comparison(Type type, vector<int> ms) {
  cout << "--------------------- LGen Level 1 Comparison -----------------------------\n\n";
  BenchmarkStats benchStats(TypeToStr(type) + "_lgen_level1_comparison");
  for (auto m : ms) {
    RealPSet* test = VMulAddBenchmark(type, m);
    ProblemInstance l1;
    l1.SetName("lgenLevelOneComparison");
    l1.SetType(type);
    l1.AddDimension(m, "m");
    auto pStats = RunProblemWithRTE(1, test, &l1);
    benchStats.AddProblemInstanceStats(pStats);
  }
  benchStats.PrettyPrintStats();
  benchStats.WriteToFiles("benchmarks");
  cout << "\n------------------- LGen Level 1 Comparison ---------------------------\n";
}

void LGenLevel2Comparison(Type type, vector<int> ms, vector<int> ns) {
  assert(ms.size() == ns.size());

  cout << "\n------------------- LGen Level 2 Comparison ---------------------------\n";
  int m, n;
  string benchmarkName = TypeToStr(type) + "_lgenCompL2";
  BenchmarkStats benchStats(benchmarkName);
  for (int i = 0; i < ms.size(); i++) {
    m = ms[i];
    n = ns[i];
    RealPSet* test = LGenCompareL2(type, m, n);
    ProblemInstance lgenCompL2;
    lgenCompL2.SetName("lgenCompL2");
    lgenCompL2.SetType(type);
    lgenCompL2.AddDimension(m, "m");
    lgenCompL2.AddDimension(n, "n");
    auto pStats = RunProblemWithRTE(1, test, &lgenCompL2);
    benchStats.AddProblemInstanceStats(pStats);
  }
  benchStats.PrettyPrintStats();
  benchStats.WriteToFiles("benchmarks");
  cout << "\n------------------- LGen Level 2 Comparison ---------------------------\n";
}

void LGenLevel3Comparison(Type type, vector<int> ms, vector<int> ns, vector<int> ps) {

}

void SVMulAddBenchmark(Type type, int m, int mInc, int numIters) {
  cout << "--------------------- svmul add benchmark -----------------------------\n\n";
  BenchmarkStats benchStats(TypeToStr(type) + "_svmul_add");
  for (int i = 0; i < numIters; i++) {
    RealPSet* test = VMulAddBenchmark(type, m);
    ProblemInstance dotProd;
    dotProd.SetName("svmulAdd");
    dotProd.SetType(type);
    dotProd.AddDimension(m, "m");
    auto pStats = RunProblemWithRTE(1, test, &dotProd);
    benchStats.AddProblemInstanceStats(pStats);
    m += mInc;
  }
  benchStats.PrettyPrintStats();
  benchStats.WriteToFiles("benchmarks");
  cout << "\n------------------- svmul add benchmark ---------------------------\n";
}

void SVMulBenchmark(Type type, VecType vecType, int m, int increment, int numIters) {
  cout << "--------------------- scalar vector multiply benchmark -----------------------------\n\n";
  BenchmarkStats benchStats(TypeToStr(type) + "_svmul");
  for (int i = 0; i < numIters; i++) {
    RealPSet* test = SVMulTest(type, vecType, m);
    ProblemInstance svmul;
    svmul.SetName("svmul");
    svmul.SetType(type);
    svmul.AddDimension(m, "m");
    auto pStats = RunProblemWithRTE(1, test, &svmul);
    benchStats.AddProblemInstanceStats(pStats);
    m += increment;
  }
  benchStats.PrettyPrintStats();
  benchStats.WriteToFiles("benchmarks");
  cout << "\n------------------- end scalar vector multiply benchmark ---------------------------\n";
}

void DotProductBenchmark(Type type, int m, int increment, int numIters) {
  cout << "--------------------- dot product benchmark -----------------------------\n\n";
  BenchmarkStats benchStats(TypeToStr(type) + "_dot_product");
  for (int i = 0; i < numIters; i++) {
    RealPSet* test = DotTest(type, m);
    ProblemInstance dotProd;
    dotProd.SetName("dotProd");
    dotProd.SetType(type);
    dotProd.AddDimension(m, "m");
    auto pStats = RunProblemWithRTE(1, test, &dotProd);
    benchStats.AddProblemInstanceStats(pStats);
    m += increment;
  }
  benchStats.PrettyPrintStats();
  benchStats.WriteToFiles("benchmarks");
  cout << "\n------------------- end dot product benchmark ---------------------------\n";
}

void ColVAddBenchmark(Type type, int m, int increment, int numIters) {
  cout << "--------------------- column vector add benchmark -----------------------------\n\n";
  BenchmarkStats benchStats(TypeToStr(type) + "_col_vector_add");
  for (int i = 0; i < numIters; i++) {
    RealPSet* test = VAddTest(type, ROWVECTOR, m);
    ProblemInstance vadd;
    vadd.SetName("vadd");
    vadd.SetType(type);
    vadd.AddDimension(m, "m");
    auto pStats = RunProblemWithRTE(1, test, &vadd);
    benchStats.AddProblemInstanceStats(pStats);
    m += increment;
  }
  benchStats.PrettyPrintStats();
  benchStats.WriteToFiles("benchmarks");
  cout << "\n------------------- end vector add benchmark ---------------------------\n";
}

void AxpyBenchmark(Type type, VecType vecType) {
  cout << "--------------------- axpy benchmark -----------------------------\n\n";

  int increment = 128;
  int m = 128;
  BenchmarkStats benchStats(TypeToStr(type) + "_" + VecTypeToStr(vecType) + "_axpy");
  for (int i = 0; i < 10; i++) {
    RealPSet* test = Axpy(type, vecType, m);
    ProblemInstance axpy;
    axpy.SetName("axpy");
    axpy.SetType(type);
    axpy.AddDimension(m, "m");
    auto pStats = RunProblemWithRTE(1, test, &axpy);
    benchStats.AddProblemInstanceStats(pStats);
    m += increment;
  }
  benchStats.PrettyPrintStats();
  benchStats.WriteToFiles("benchmarks");

  cout << "\n------------------- axpy benchmark ---------------------------\n";
}

void GemvBenchmark(Type type, bool transpose, int mBase, int mInc, int nBase, int nInc) {
  cout << "--------------------- gemv benchmark -----------------------------\n\n";

  int m = mBase;
  int n = nBase;
  string benchmarkName;
  if (transpose) {
    benchmarkName = TypeToStr(type) + "_transposed_gemv";
  } else {
    benchmarkName = TypeToStr(type) + "_gemv";
  }
  BenchmarkStats benchStats(benchmarkName);
  for (int i = 0; i < 10; i++) {
    RealPSet* test = Gemv(type, transpose, m, n);
    ProblemInstance gemv;
    gemv.SetName("gemv");
    gemv.SetType(type);
    gemv.AddDimension(m, "m");
    gemv.AddDimension(n, "n");
    auto pStats = RunProblemWithRTE(1, test, &gemv);
    benchStats.AddProblemInstanceStats(pStats);
    m += mInc;
    n += nInc;
  }
  benchStats.PrettyPrintStats();
  benchStats.WriteToFiles("benchmarks");
  cout << "\n------------------- gemv benchmark ---------------------------\n";
}

void MAddBenchmark(Type type, int mBase, int mInc, int nBase, int nInc, int numIters) {
  cout << "--------------------- madd benchmark -----------------------------\n\n";

  int m = mBase;
  int n = nBase;
  string benchmarkName = TypeToStr(type) + "_madd";
  BenchmarkStats benchStats(benchmarkName);
  for (int i = 0; i < numIters; i++) {
    RealPSet* test = MAddTest(type, m, n);
    ProblemInstance madd;
    madd.SetName("madd");
    madd.SetType(type);
    madd.AddDimension(m, "m");
    madd.AddDimension(n, "n");
    auto pStats = RunProblemWithRTE(1, test, &madd);
    benchStats.AddProblemInstanceStats(pStats);
    m += mInc;
    n += nInc;
  }
  benchStats.PrettyPrintStats();
  benchStats.WriteToFiles("benchmarks");
  cout << "\n------------------- madd benchmark ---------------------------\n";
}

void LGenCompareL2Benchmark(Type type, int mBase, int mInc, int nBase, int nInc, int numIters) {
  cout << "--------------------- LGen comparison level 2 benchmark -----------------------------\n\n";
  int m = mBase;
  int n = nBase;
  string benchmarkName = TypeToStr(type) + "_lgenCompL2";
  BenchmarkStats benchStats(benchmarkName);
  for (int i = 0; i < numIters; i++) {
    RealPSet* test = LGenCompareL2(type, m, n);
    ProblemInstance lgenCompL2;
    lgenCompL2.SetName("lgenCompL2");
    lgenCompL2.SetType(type);
    lgenCompL2.AddDimension(m, "m");
    lgenCompL2.AddDimension(n, "n");
    auto pStats = RunProblemWithRTE(1, test, &lgenCompL2);
    benchStats.AddProblemInstanceStats(pStats);
    m += mInc;
    n += nInc;
  }
  benchStats.PrettyPrintStats();
  benchStats.WriteToFiles("benchmarks");
  cout << "\n------------------- LGen comparison level 2 benchmark ---------------------------\n";
}

void RunDotProdBenchmarks() {
  DotProductBenchmark(REAL_SINGLE, 128, 128, 10);
  DotProductBenchmark(REAL_DOUBLE, 128, 128, 10);
}

void RunAxpyBenchmarks() {
  AxpyBenchmark(REAL_SINGLE, ROWVECTOR);
  AxpyBenchmark(REAL_SINGLE, COLVECTOR);
  AxpyBenchmark(REAL_DOUBLE, ROWVECTOR);
  AxpyBenchmark(REAL_DOUBLE, COLVECTOR);
}

void RunGemvBenchmarks() {
  GemvBenchmark(REAL_SINGLE, false, 16, 16, 16, 16);
  GemvBenchmark(REAL_DOUBLE, false, 16, 16, 16, 16);
}

void RunAllBenchmarks() {
  cout << "=========================================================================\n";
  cout << "======================== STARTING LLDLA BENCHMARK =======================\n";
  cout << "=========================================================================\n\n";

  RunDotProdBenchmarks();
  RunAxpyBenchmarks();
  RunGemvBenchmarks();

  cout << "\n=========================================================================\n";
  cout << "======================== DONE WITH LLDLA BENCHMARK ======================\n";
  cout << "=========================================================================\n";
}

#endif // DOLLDLA
