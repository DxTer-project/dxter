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

#include "problemInstanceStats.h"

#if DOLLDLA

void ProblemInstanceStats::ComputeBestAndWorstImplementations(Type type) {
  double bestAvgFlopsPerCycle = 0.0;
  double bestFlopsPerCycle = 0.0;
  double worstFlopsPerCycle = arch->FlopsPerCycle(m_type);
  for (auto impl : *m_implementationStats) {
    if (impl->GetAvgFlopsPerCycle() > bestAvgFlopsPerCycle) {
      bestAvgFlopsPerCycle = impl->GetAvgFlopsPerCycle();
      m_bestAvgFlopsPerCycleImpl = impl;
    }

    if (impl->GetBestFlopsPerCycle() > bestFlopsPerCycle) {
      bestFlopsPerCycle = impl->GetBestFlopsPerCycle();
      m_bestFlopsPerCycleImpl = impl;
    }

    if (impl->GetWorstFlopsPerCycle() < worstFlopsPerCycle) {
      worstFlopsPerCycle = impl->GetWorstFlopsPerCycle();
      m_worstFlopsPerCycleImpl = impl;
    }
  }
}

vector<ImplementationStats*>* ProblemInstanceStats::ComputeImplementationStats(ImplementationRuntimeMap* impls) {
  vector<ImplementationStats*>* implStats = new vector<ImplementationStats*>();
  for (std::pair<GraphNum, TimeVec> pair : *impls) {
    GraphNum impNum = pair.first;
    TimeVec* times = &(pair.second);
    ImplementationStats* impStats = new ImplementationStats(impNum, m_type, m_cost, times);
    implStats->push_back(impStats);
  }
  return implStats;
}

ProblemInstanceStats::ProblemInstanceStats(ProblemInstance* problemInstance, ImplementationRuntimeMap* impls) {
  m_cost = problemInstance->GetCost();
  m_name = new string(problemInstance->GetName());
  m_type = problemInstance->GetType();
  m_implementationStats = ComputeImplementationStats(impls);
  ComputeBestAndWorstImplementations(problemInstance->GetType());
}

ProblemInstanceStats::~ProblemInstanceStats() {
  for (auto implStats : *m_implementationStats) {
    delete implStats;
  }
  delete m_implementationStats;
}

double ProblemInstanceStats::GetBestAvgFlopsPerCycle() {
  return m_bestAvgFlopsPerCycleImpl->GetAvgFlopsPerCycle();
}

GraphNum ProblemInstanceStats::GetBestAvgFlopsPerCycleImpl() {
  return m_bestAvgFlopsPerCycleImpl->GetNum();
}

void ProblemInstanceStats::PrettyPrintPerformanceStats() {
  cout << "================== PERFORMANCE RESULTS FOR " << *m_name << " =====================\n";
  cout << "\n&&&&&&&&&&&&&&&&&& Problem Summary &&&&&&&&&&&&&&&&&&&" << endl;
  cout << "Datatype             : " << TypeToStr(m_type) << endl;
  cout << "Flop count           : " << m_cost << endl;
  cout << "# of Implementations : " << m_implementationStats->size() << endl;
  cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n" << endl;

  for (auto implStats : *m_implementationStats) {
    implStats->PrettyPrintPerformanceStats();
  }

  cout << "=======================================================================================================\n";
}

string ProblemInstanceStats::CSVLineColumnTitles() {
  return "DUMMY,COLUMN,TITLES";
}

string ProblemInstanceStats::CSVLine() {
  return "123,456,789";
}

#endif // DOLLDLA
