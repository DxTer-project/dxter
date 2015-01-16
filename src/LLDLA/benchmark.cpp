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


#include "benchmark.h"

#if DOLLDLA

using namespace std;

void DotProductBenchmark() {
  cout << "--------------------- dot product benchmark -----------------------------\n\n";
  int increment = 128;
  int m = 128;
  BenchmarkStats benchStats("singlePrecision dot product");
  for (int i = 0; i < 10; i++) {
    RealPSet* test = DotExample(REAL_SINGLE, m);
    ProblemInstance dotProd;
    dotProd.SetName("dotProdNope");
    dotProd.SetType(REAL_SINGLE);
    dotProd.AddDimension(m, "m");
    ProblemInstanceStats* pStats = RunBenchmark(1, test, &dotProd);
    benchStats.AddProblemInstanceStats(pStats);
    m += increment;
  }
  benchStats.PrettyPrintStats();
  cout << "\n------------------- end dot product benchmark ---------------------------\n";
}

void RunBenchmark() {
  cout << "=========================================================================\n";
  cout << "======================== STARTING LLDLA BENCHMARK =======================\n";
  cout << "=========================================================================\n\n";

  DotProductBenchmark();

  cout << "\n=========================================================================\n";
  cout << "======================== DONE WITH LLDLA BENCHMARK ======================\n";
  cout << "=========================================================================\n";
}

#endif // DOLLDLA
