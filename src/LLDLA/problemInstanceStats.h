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

#include "base.h"
#include "implementationStats.h"
#include "problemInstance.h"

#if DOLLDLA

class ProblemInstanceStats {

 private:
  ImplementationStats* m_bestAvgFlopsPerCycleImpl;
  ImplementationStats* m_bestFlopsPerCycleImpl;
  ImplementationStats* m_worstFlopsPerCycleImpl;

  ProblemInstance* m_problemInstance;
  vector<ImplementationStats*>* m_implementationStats;

  void ComputeBestAndWorstImplementations(Type type);
  vector<ImplementationStats*>* ComputeImplementationStats(ImplementationRuntimeMap* impls);

 public:
  ProblemInstanceStats(ProblemInstance* problemInstance, ImplementationRuntimeMap* impls);
  ~ProblemInstanceStats();

  double GetBestAvgFlopsPerCycle();
  GraphNum GetBestAvgFlopsPerCycleImpl();
  void PrettyPrintPerformanceStats();
};

#endif // DOLLDLA

#endif // PROBLEM_INSTANCE_STATS_H_