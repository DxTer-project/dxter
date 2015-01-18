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
  m_name = new string(NoWhitespace(name));
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

void BenchmarkStats::CreateAllBenchmarksDirectory(string benchmarkDirName) {
  struct stat st = {0};
  if (stat(benchmarkDirName.c_str(), &st) == -1) {
    mkdir(benchmarkDirName.c_str(), 0700);
  }
}

string BenchmarkStats::CreateThisBenchmarksDirectory(string benchmarkDirName) {
  struct stat st = {0};
  string* benchmarkPath = new string(benchmarkDirName + "/" + *m_name + DateAndTimeString());
  if (stat(benchmarkPath->c_str(), &st) == -1) {
    mkdir(benchmarkPath->c_str(), 0700);
  } else {
    cout << "ERROR: " << benchmarkDirName << " already exists!" << endl;
    throw;
  }
  return *benchmarkPath;
}

void BenchmarkStats::WriteToFiles(string dirName) {
  CreateAllBenchmarksDirectory(dirName);
  string benchmarkPath = CreateThisBenchmarksDirectory(dirName);
  WriteInstanceDataCSV(benchmarkPath);
  WriteImplementationDataCSVs(benchmarkPath);
}

string BenchmarkStats::CSVColumnTitles() {
  ProblemInstanceStats* problemInst = m_problemInstances->front();
  return problemInst->CSVLineColumnTitles();
}

void BenchmarkStats::WriteInstanceDataCSV(string benchmarkDirPath) {
  ofstream instanceDataFile(benchmarkDirPath + "/instanceData.csv");
  instanceDataFile << CSVColumnTitles() << endl;
  for (auto problemInstStats : *m_problemInstances) {
    instanceDataFile << problemInstStats->CSVLine() << endl;
  }
  instanceDataFile.close();
}

string BenchmarkStats::CreateImplCSVDirectory(string benchmarkDirPath) {
  struct stat st = {0};
  string* implDirPath = new string(benchmarkDirPath + "/implementationData");
  if (stat(implDirPath->c_str(), &st) == -1) {
    mkdir(implDirPath->c_str(), 0700);
  } else {
    cout << "ERROR: " << implDirPath << " already exists!" << endl;
    throw;
  }
  return *implDirPath;
}

void BenchmarkStats::WriteImplementationDataCSVs(string benchmarkDirPath) {
  string impCSVPath = CreateImplCSVDirectory(benchmarkDirPath);
  for (auto problemInstStats : *m_problemInstances) {
    problemInstStats->WriteImplementationCSV(impCSVPath);
  }
}

#endif // DOLLDLA
