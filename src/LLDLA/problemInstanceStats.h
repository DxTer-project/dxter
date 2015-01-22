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

#ifndef PROBLEM_INSTANCE_STATS_H_
#define PROBLEM_INSTANCE_STATS_H_

#include <memory>

#include "base.h"
#include "implementationStats.h"
#include "problemInstance.h"
#include "runnerUtils.h"

#if DOLLDLA

class ProblemInstanceStats {

 private:
  unique_ptr<string> m_name;
  Type m_type;
  Cost m_cost;

  ImplementationStats* m_bestAvgFlopsPerCycleImpl;
  ImplementationStats* m_bestFlopsPerCycleImpl;
  ImplementationStats* m_worstFlopsPerCycleImpl;

  vector<string*>* m_dimNames;
  vector<int>* m_dimValues;
  vector<unique_ptr<ImplementationStats>> m_implementationStats;

  vector<string*> AllDimNames(ProblemInstance* problemInstance);
  void ComputeBestAndWorstImplementations(Type type);
  void ComputeImplementationStats(ImplementationRuntimeMap* impls);

 public:
  ProblemInstanceStats(ProblemInstance* problemInstance, ImplementationRuntimeMap* impls);
  ~ProblemInstanceStats();

  string CSVLineColumnTitles();
  double GetBestAvgFlopsPerCycle();
  GraphNum GetBestAvgFlopsPerCycleImpl();
  void PrettyPrintPerformanceStats();
  string CSVLine();
  void WriteImplementationCSV(string impCSVPath);
};

#endif // DOLLDLA

#endif // PROBLEM_INSTANCE_STATS_H_
