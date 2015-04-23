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

#ifndef ONE_STAGE_TIMING_RESULT_H
#define ONE_STAGE_TIMING_RESULT_H

#include "timingResult.h"

#if DOLLDLA

class OneStageTimingResult : public TimingResult {
 private:
  GraphNum m_implNum;
  vector<double> m_runtimes;

 public:
  virtual ~OneStageTimingResult() {}
  OneStageTimingResult(GraphNum implNum, TimeVec* runtimes);

  virtual GraphNum GetNum() { return m_implNum; }
  virtual double GetAvgFlopsPerCycle() { return 0.0; }
  virtual double GetAvgPercentOfPeak() { return 0.0; }

  virtual string CSVLine() { return ""; }
  virtual string CSVLineColumnTitles() { return ""; }
  virtual void PrettyPrintPerformanceStats() { return; }

  virtual vector<double>* GetTimes() { return &m_runtimes; }
};

#endif // DOLLDLA

#endif // ONE_STAGE_TIMING_RESULT_H
