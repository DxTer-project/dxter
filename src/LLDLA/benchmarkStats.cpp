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

#include "benchmarkStats.h"

#if DOLLDLA

BenchmarkStats::BenchmarkStats(string name) {
  m_name = new string(name);
  m_problemInstances = new vector<ProblemInstanceStats*>();
}

BenchmarkStats::~BenchmarkStats() {
  delete m_name;
  delete m_problemInstances;
}

void BenchmarkStats::AddProblemInstanceStats(ProblemInstanceStats* stats) {
  m_problemInstances->push_back(stats);
}

void BenchmarkStats::PrettyPrintStats() {
  cout << "----------------------------- " << *m_name << " resuts ----------------------------\n";
  for (auto problemStats : *m_problemInstances) {
    problemStats->PrettyPrintPerformanceStats();
  }
  cout << "----------------------------------------------------------------------------\n";
}

#endif // DOLLDLA
