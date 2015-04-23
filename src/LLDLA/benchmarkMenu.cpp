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

#if DOLLDLA

#include "assert.h"
#include "benchmark.h"
#include "driverUtils.h"

int PromptUserForInt(string valName) {
  cout << valName << ": ";
  int val;
  cin >> val;
  return val;
}

Type PromptUserForType() {
  cout << "Enter datatype (F = single precision float) (D = double precision float): ";
  char typeChar;
  cin >> typeChar;
  auto type = CharToType(typeChar);
  return type;
}

vector<int> FixedIncrementDimension(int start, int end, int inc) {
  vector<int> dim;
  for (int i = start; i <= end; i += inc) {
    dim.push_back(i);
  }
  return dim;
}

vector<int> FixedIncrementCeiling(vector<int> fixedIncs, int divFactor) {
  assert(divFactor != 0);
  vector<int> dim;
  for (auto d : fixedIncs) {
    int res = d / divFactor;
    if (res == 0) {
      dim.push_back(1);
    } else {
      dim.push_back(res);
    }
  }
  return dim;
}

vector<int> FixedValueDimension(int val, int numRepeats) {
  vector<int> fixedDim;
  for (int i = 0; i < numRepeats; i++) {
    fixedDim.push_back(val);
  }
  return fixedDim;
}

vector<int> PromptUserForFixedValueDimension() {
  int val = PromptUserForInt("Constant value");
  int numRepeats = PromptUserForInt("Repeated how many times");
  return FixedValueDimension(val, numRepeats);
}

vector<int> PromptUserForFixedIncrementDimension() {
  int start = PromptUserForInt("Starting value");
  int end = PromptUserForInt("Ending value");
  int inc = PromptUserForInt("Increment");
  return FixedIncrementDimension(start, end, inc);
}

vector<int> PromptUserForFixedIncrementCeiling() {
  auto divFactor = PromptUserForInt("Factor to divide by");
  auto fixedVals = PromptUserForFixedIncrementDimension();
  return FixedIncrementCeiling(fixedVals, divFactor);
}

vector<int> PromptUserForDimension(string dimName) {
  cout << "What kind of values should dimension " << dimName << " have? " << endl;
  cout << "\t0 -> Fixed values i.e. [10, 10, 10, 10]" << endl;
  cout << "\t1 -> Values at fixed increments i.e. [1, 6, 11, 16]" << endl;
  cout << "\t2 -> ceiling of fixed increments divided by a factor" << endl;
  int dimType;
  cin >> dimType;
  switch(dimType) {
  case(0):
    return PromptUserForFixedValueDimension();
  case(1):
    return PromptUserForFixedIncrementDimension();
  case(2):
    return PromptUserForFixedIncrementCeiling();
  default:
    cout << dimType << " is not a valid dimension type" << endl;
    throw;
  }
}

unsigned int GetBenchmarkNumber() {
  cout << "---------------------- BENCHMARKS --------------------------" << endl;
  cout << "     0 -> Run all benchmarks" << endl;
  cout << "     1 -> Dot product benchmark" << endl;
  cout << "     2 -> Axpy benchmarks" << endl;
  cout << "     3 -> Gemv benchmarks" << endl;
  cout << "     4 -> MAdd benchmark" << endl;
  cout << "     5 -> Column VAdd benchmark" << endl;
  cout << "     6 -> SVMul benchmark" << endl;
  cout << "     7 -> y := alpha*x + beta*(y + z)" << endl;
  cout << "     8 -> y := (A + B^T)*x" << endl;
  cout << "     9 -> LGen level 1 comparison (y := alpha*x + beta*(y + z))" << endl;
  cout << "select one of the options listed above: ";
  unsigned int benchmarkOption;
  cin >> benchmarkOption;
  return benchmarkOption;
}

void LGenLevel1ComparisonMenu() {
  Type type;
  vector<int> ms;
  type = PromptUserForType();
  ms = PromptUserForDimension("M");
  LGenLevel1Comparison(type, ms);
}

void LGenLevel2ComparisonMenu() {
  Type type;
  int m, n, mInc, nInc, iters;

  type = PromptUserForType();
  cout << "Starting M value: ";
  cin >> m;
  cout << "M increment: ";
  cin >> mInc;
  cout << "Starting N value: ";
  cin >> n;
  cout << "N increment: ";
  cin >> nInc;
  cout << "Number of iterations: ";
  cin >> iters;
  LGenCompareL2Benchmark(type, m, mInc, n, nInc, iters);
}

void SVMulAddBenchmarkMenu() {
 Type type;
  int m, inc, iters;

  type = PromptUserForType();
  cout << "Starting M value: ";
  cin >> m;
  cout << "M increment: ";
  cin >> inc;
  cout << "Number of iterations: ";
  cin >> iters;
  SVMulAddBenchmark(type, m, inc, iters);
}

void DotProductBenchmarkMenu() {
  Type type;
  int m, inc, iters;

  type = PromptUserForType();
  cout << "Starting M value: ";
  cin >> m;
  cout << "M increment: ";
  cin >> inc;
  cout << "Number of iterations: ";
  cin >> iters;
  DotProductBenchmark(type, m, inc, iters);
}

void MAddBenchmarkMenu() {
  Type type;
  int m, n, mInc, nInc, iters;

  type = PromptUserForType();
  cout << "Starting M value: ";
  cin >> m;
  cout << "M increment: ";
  cin >> mInc;
  cout << "Starting N value: ";
  cin >> n;
  cout << "N increment: ";
  cin >> nInc;
  cout << "Number of iterations: ";
  cin >> iters;
  MAddBenchmark(type, m, mInc, n, nInc, iters);
}

void ColVAddBenchmarkMenu() {
  Type type;
  int m, inc, iters;

  type = PromptUserForType();
  cout << "Starting M value: ";
  cin >> m;
  cout << "M increment: ";
  cin >> inc;
  cout << "Number of iterations: ";
  cin >> iters;
  ColVAddBenchmark(type, m, inc, iters);
}

void SVMulBenchmarkMenu() {
  Type type;
  VecType vecType;
  int m, mInc, iters;

  cout << "Enter vector type (C = column) (R = row): ";
  char vecTypeChar;
  cin >> vecTypeChar;
  vecType = CharToVecType(vecTypeChar);

  type = PromptUserForType();
  cout << "Starting M value: ";
  cin >> m;
  cout << "M increment: ";
  cin >> mInc;
  cout << "Number of iterations: ";
  cin >> iters;
  SVMulBenchmark(type, vecType, m, mInc, iters);
}

void RunBenchmarkNumber(unsigned int num) {
  switch(num) {
  case(0):
    RunAllBenchmarks();
    break;
  case(1):
    DotProductBenchmarkMenu();
    break;
  case(2):
    RunAxpyBenchmarks();
    break;
  case(3):
    RunGemvBenchmarks();
    break;
  case(4):
    MAddBenchmarkMenu();
    break;
  case(5):
    ColVAddBenchmarkMenu();
  case(6):
    SVMulBenchmarkMenu();
    break;
  case(7):
    SVMulAddBenchmarkMenu();
    break;
  case(8):
    LGenLevel2ComparisonMenu();
    break;
  case(9):
    LGenLevel1ComparisonMenu();
    break;
  default:
    cout << "Error: " << num << " is not a valid benchmark number" << endl;
    break;
  }
}

void BenchmarkMenu() {
  unsigned int option = GetBenchmarkNumber();
  RunBenchmarkNumber(option);
}

#endif // DOLLDLA
