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

#include "LLDLA.h"

#if DOLLDLA

class ImplementationStats {

 private:
  double
    m_meanFlopsPerCycle,
    m_meanPercentOfPeak,
    m_bestPercentOfPeak,
    m_worstPercentOfPeak,
    m_bestFlopsPerCycle,
    m_worstFlopsPerCycle;

  GraphNum m_implNumber;

  void ComputeFlopsPerCycle(Type type, Cost flopsCost, TimeVec* runtimes);
  void ComputePercentagesOfPeak(Type type);
  void ComputePerformanceNumbers(Type type, Cost flopCost, TimeVec* runtimes);

 public:
  ImplementationStats(GraphNum implNumber, Type type, Cost flopCost, TimeVec* runtimes);

  void PrettyPrintPerformanceStats();

};

#endif // DOLLDLA
