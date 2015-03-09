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

#include "benchmarkStats.h"
#include "problemInstance.h"
#include "runBenchmark.h"
#include "singleOperationExamples.h"

using namespace std;

void DotProductBenchmark(Type type, int m, int increment, int numIters) {
  cout << "--------------------- dot product benchmark -----------------------------\n\n";
  BenchmarkStats benchStats(TypeToStr(type) + "_dot_product");
  for (int i = 0; i < numIters; i++) {
    RealPSet* test = DotTest(type, m);
    ProblemInstance dotProd;
    dotProd.SetName("dotProd");
    dotProd.SetType(type);
    dotProd.AddDimension(m, "m");
    auto pStats = RunBenchmark(1, test, &dotProd);
    benchStats.AddProblemInstanceStats(pStats);
    m += increment;
  }
  benchStats.PrettyPrintStats();
  benchStats.WriteToFiles("benchmarks");
  cout << "\n------------------- end dot product benchmark ---------------------------\n";
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
    auto pStats = RunBenchmark(1, test, &axpy);
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
    auto pStats = RunBenchmark(1, test, &gemv);
    benchStats.AddProblemInstanceStats(pStats);
    m += mInc;
    n += nInc;
  }
  benchStats.PrettyPrintStats();
  benchStats.WriteToFiles("benchmarks");
  cout << "\n------------------- gemv benchmark ---------------------------\n";
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
