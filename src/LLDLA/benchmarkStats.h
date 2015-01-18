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

#ifndef BENCHMARK_STATS_H_
#define BENCHMARK_STATS_H_

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "problemInstanceStats.h"
#include "runnerUtils.h"

#if DOLLDLA

class BenchmarkStats {

 private:
  string* m_name;
  vector<ProblemInstanceStats*>* m_problemInstances;

  void CreateAllBenchmarksDirectory(string benchmarkDirName);
  string CreateThisBenchmarksDirectory(string benchmarkDirName);

  string CSVColumnTitles();

  string CreateImplCSVDirectory(string benchmarkDirPath);
  void WriteInstanceDataCSV(string benchmarkDirPath);
  void WriteImplementationDataCSVs(string benchmarkDirPath);

 public:
  BenchmarkStats(string name);
  ~BenchmarkStats();

  void AddProblemInstanceStats(ProblemInstanceStats* stats);
  
  void PrettyPrintStats();
  void WriteToFiles(string dirName);
};

#endif // DOLLDLA

#endif // BENCHMARK_STATS_H_
