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

#include "implementationStats.h"

#if DOLLDLA

void ImplementationStats::ComputeFlopsPerCycle(Type type, Cost flopCost, TimeVec* runtimes) {
  m_meanFlopsPerCycle = 0.0;
  m_bestFlopsPerCycle = 0.0;
  m_worstFlopsPerCycle = arch->FlopsPerCycle(type);

  double numRuns = 0.0;
  double flopsPerCycle;

  for (double runtime : *runtimes) {
    flopsPerCycle = flopCost / runtime;

    m_meanFlopsPerCycle += flopsPerCycle;

    if (m_bestFlopsPerCycle < flopsPerCycle) {
      m_bestFlopsPerCycle = flopsPerCycle;
    }

    if (m_worstFlopsPerCycle > flopsPerCycle) {
      m_worstFlopsPerCycle = flopsPerCycle;
    }

    numRuns += 1.0;
  }

  m_meanFlopsPerCycle = m_meanFlopsPerCycle / numRuns;
}

void ImplementationStats::ComputePercentagesOfPeak(Type type) {
  m_meanPercentOfPeak = (m_meanFlopsPerCycle / arch->FlopsPerCycle(type)) * 100;
  m_worstPercentOfPeak = (m_worstFlopsPerCycle / arch->FlopsPerCycle(type)) * 100;
  m_bestPercentOfPeak = (m_bestFlopsPerCycle / arch->FlopsPerCycle(type)) * 100;
}

void ImplementationStats::ComputePerformanceNumbers(Type type, Cost flopCost, TimeVec* runtimes) {
  ComputeFlopsPerCycle(type, flopCost, runtimes);
  ComputePercentagesOfPeak(type);
}

ImplementationStats::ImplementationStats(GraphNum implNumber, Type type, Cost flopCost, TimeVec* runtimes) {
  m_implNumber = implNumber;
  ComputePerformanceNumbers(type, flopCost, runtimes);
}

GraphNum ImplementationStats::GetNum() {
  return m_implNumber;
}

double ImplementationStats::GetAvgFlopsPerCycle() {
  return m_meanFlopsPerCycle;
}

double ImplementationStats::GetAvgPercentOfPeak() {
  return m_meanPercentOfPeak;
}

double ImplementationStats::GetBestFlopsPerCycle() {
  return m_bestFlopsPerCycle;
}

double ImplementationStats::GetBestPercentOfPeak() {
  return m_bestPercentOfPeak;
}

double ImplementationStats::GetWorstFlopsPerCycle() {
  return m_worstFlopsPerCycle;
}

string ImplementationStats::CSVLine() {
  string line = std::to_string((long long int) GetNum()) + ",";
  line += std::to_string((long double) GetAvgFlopsPerCycle()) + ",";
  line += std::to_string((long double) GetAvgPercentOfPeak()) + ",";
  line += std::to_string((long double) GetBestFlopsPerCycle()) + ",";
  line += std::to_string((long double) GetBestPercentOfPeak()) + ",";
  return line;
}

string ImplementationStats::CSVLineColumnTitles() {
  return "implementation_number,avg_flops_per_cycle,avg_percent_of_peak,best_flops_per_cycle,best_percent_of_peak,";
}

void ImplementationStats::PrettyPrintPerformanceStats()
{
  cout << "IMPLEMENTATION #" << m_implNumber << endl;
  cout << "     Average flops/cycle = " << m_meanFlopsPerCycle << endl;
  cout << "               % of peak = " << m_meanPercentOfPeak << endl << endl;
  cout << "        Best flops/cycle = " << m_bestFlopsPerCycle << endl;
  cout << "               % of peak = " << m_bestPercentOfPeak << endl << endl;
  cout << "       Worst flops/cycle = " << m_worstFlopsPerCycle << endl;
  cout << "               % of peak = " << m_worstPercentOfPeak << endl << endl;

  return;
}

#endif // DOLLDLA
