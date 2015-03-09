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

#include "benchmark.h"
#include "driverUtils.h"

unsigned int GetBenchmarkNumber() {
  cout << "---------------------- BENCHMARKS --------------------------" << endl;
  cout << "     0 -> Run all benchmarks" << endl;
  cout << "     1 -> Dot product benchmark" << endl;
  cout << "     2 -> Axpy benchmarks" << endl;
  cout << "     3 -> Gemv benchmarks" << endl << endl;
  cout << "select one of the options listed above: ";
  unsigned int benchmarkOption;
  cin >> benchmarkOption;
  return benchmarkOption;
}

void DotProductBenchmarkMenu() {
  Type type;
  int m, inc, iters;

  cout << "Enter datatype (F = single precision float) (D = double precision float): ";
  char typeChar;
  cin >> typeChar;
  type = CharToType(typeChar);
  cout << "Starting M value: ";
  cin >> m;
  cout << "M increment: ";
  cin >> inc;
  cout << "Number of iterations: ";
  cin >> iters;
  DotProductBenchmark(type, m, inc, iters);
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
