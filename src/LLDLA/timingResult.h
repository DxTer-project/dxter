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

#ifndef TIMING_RESULT_H
#define TIMING_RESULT_H

#if DOLLDLA

class TimingResult {
 public:
  virtual ~TimingResult() {}
  virtual GraphNum GetNum() {
    LOG_FAIL("TimingResult::GetNum() cannot be called");
    throw;
  }
  virtual double GetAvgFlopsPerCycle() {
    LOG_FAIL("TimingResult::GetAvgFlopsPerCycle() cannot be called");
    throw;
  }
  virtual double GetAvgPercentOfPeak() {
    LOG_FAIL("TimingResult::GetAvgPercentOfPeak() cannot be called");
    throw;
  }

  virtual string CSVLine() {
    LOG_FAIL("TimingResult::CSVLine() cannot be called");
    throw;
  }
  virtual string CSVLineColumnTitles() {
    LOG_FAIL("TimingResult::CSVLineColumnTitles() cannot be called");
    throw;
  }
  virtual void PrettyPrintPerformanceStats() {
    LOG_FAIL("TimingResult::PrettyPrintPerformanceStats() cannot be called");
    throw;
  }
};

#endif // DOLLDLA

#endif // TIMING_RESULT_H
