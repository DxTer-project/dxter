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

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "runtimeEvaluation.h"

#if DOLLDLA

#include <assert.h>

#include "oneStageTimingResult.h"

RuntimeEvaluator::RuntimeEvaluator(string evalDirName) {
  m_evalDirName = evalDirName;
}

void RuntimeEvaluator::WriteTestCodeToFile(string executableName, string testCode) {
  string testFileName = executableName + ".c";
  std::ofstream outStream(testFileName);
  if (outStream.is_open()) {
    outStream << testCode << endl;
    outStream.close();
  } else {
    cout << "ERROR: RuntimeEvaluator could not create file " << testFileName << endl;
    LOG_FAIL("replacement for throw call");
    throw;
  }
}

void RuntimeEvaluator::CompileTest(string executableName) {
  string testFileName = executableName + ".c";
  cout << "All implementations written to files\n";
  cout << "Compile string is " << arch->CompileString(executableName, testFileName) << endl;
  int compileRes = system(arch->CompileString(executableName, testFileName).c_str());
  cout << "Compile result = " << std::to_string((long long int) compileRes) << endl;
}

void RuntimeEvaluator::RunTest(string executableName) {
  string runStr = "./" + executableName;
  int runRes = system(runStr.c_str());
  cout << "Run string is " << runStr << endl;
  cout << "Run result = " << std::to_string((long long int) runRes) << endl;
}

void RuntimeEvaluator::CleanUpTest(string executableName) {
  string removeExecutable = "rm -f " + executableName;
  system(removeExecutable.c_str());
}

vector<TimingResult*>* RuntimeEvaluator::ReadTimingData(TimingSetting timingSetting, string dataFileName, int numImpls) {
  cout << "Size of imps = " << std::to_string((long long int) numImpls) << endl;
  cout << "Calling ReadTimeDataFromFile" << endl;

  auto timingResults = ReadTimeDataFromFile(timingSetting, dataFileName, numImpls);
  string removeDataFile = "rm -f " + dataFileName;
  system(removeDataFile.c_str());
  return timingResults;
}

vector<TimingResult*>* RuntimeEvaluator::EvaluateBatch(SanityCheckSetting sanityCheckSetting, TimingSetting timingSetting, RuntimeTest test, vector<pair<GraphNum, ImplInfo>>* impls, string referenceImp) {
  string executableName = m_evalDirName + "/" + test.m_operationName;
  string testCode = test.MakeTestCode(sanityCheckSetting, timingSetting, impls, referenceImp);

  WriteTestCodeToFile(executableName, testCode);
  CompileTest(executableName);
  RunTest(executableName);
  CleanUpTest(executableName);

  auto timeResults = ReadTimingData(timingSetting, test.m_dataFileName, impls->size());
  return timeResults;
}

vector<vector<pair<GraphNum, ImplInfo>>*>* RuntimeEvaluator::OneBatch(map<GraphNum, ImplInfo>* imps) {
  auto currentBatch = new vector<pair<GraphNum, ImplInfo>>();
  for (auto impl : *imps) {
    pair<GraphNum, ImplInfo> newPair(impl.first, impl.second);
    currentBatch->push_back(newPair);
  }
  auto batches = new vector<vector<pair<GraphNum, ImplInfo>>*>();
  batches->push_back(currentBatch);
  return batches;
}

vector<vector<pair<GraphNum, ImplInfo>>*>* RuntimeEvaluator::BreakIntoBatches(map<GraphNum, ImplInfo>* imps, unsigned int batchSize) {
  cout << "Entering BreakIntoBatches" << endl;
  cout << "Size of imps = " << std::to_string((long long int) imps->size()) << endl;
  if (imps->size() <= batchSize) {
    return OneBatch(imps);
  }
  auto batches = new vector<vector<pair<GraphNum, ImplInfo>>*>();
  auto currentBatch = new vector<pair<GraphNum, ImplInfo>>();
  unsigned int i = 0;
  for (auto impl : *imps) {
    if (i == batchSize) {
      batches->push_back(currentBatch);
      currentBatch = new vector<pair<GraphNum, ImplInfo>>();
      i = 0;
    } else {
      i++;
    }
    pair<GraphNum, ImplInfo> newPair(impl.first, impl.second);
    currentBatch->push_back(newPair);
    cout << "Size of current batch = " << std::to_string((long long int) currentBatch->size()) << endl;
  }
  return batches;
}

vector<TimingResult*>* RuntimeEvaluator::EvaluateImplementations(SanityCheckSetting sanityCheckSetting, TimingSetting timingSetting, RuntimeTest test, map<GraphNum, ImplInfo>* imps, string referenceImp) {
  cout << "Entering EvaulateImplementations" << endl;
  auto batchVec = BreakIntoBatches(imps, 100000);
  cout << "size of batchVec = " << std::to_string((long long int) batchVec->size()) << endl;
  auto results = new vector<TimingResult*>();

  for (auto batch : *batchVec) {
    cout << "Evalutating batches" << endl;
    auto batchTimeResults = EvaluateBatch(sanityCheckSetting, timingSetting, test, batch, referenceImp);
    results->insert(results->end(), batchTimeResults->begin(), batchTimeResults->end());
  }

  for (auto batch : *batchVec) {
    delete batch;
  }
  delete batchVec;

  return results;
}

bool RuntimeEvaluator::IsImplementationSeparator(string token) {
  if (token == "#") {
    return true;
  }
  return false;
}

vector<string> RuntimeEvaluator::GetFileTokens(string fileName) {
  std::ifstream dataStream(fileName);
  std::stringstream buffer;
  buffer << dataStream.rdbuf();
  dataStream.close();
  string timeData = buffer.str();
  std::vector<string> runtimeStrings;
  Tokenize(timeData, runtimeStrings, "\n");
  return runtimeStrings;
}

vector<TimingResult*>* RuntimeEvaluator::ReadTimeDataFromFile(TimingSetting timingSetting, string fileName, int numImpls) {
  vector<string> runtimeStrings = GetFileTokens(fileName);
  vector<TimingResult*>* timingResults = new vector<TimingResult*>();
  int i = 1;
  TimeVec* impTimes = new TimeVec();
  bool atEndOfImplTimeData = false;
  for (auto token : runtimeStrings) {
    assert(!(atEndOfImplTimeData && IsImplementationSeparator(token)));

    if (IsImplementationSeparator(token)) {
      atEndOfImplTimeData = true;
    } else if (atEndOfImplTimeData) {
      atEndOfImplTimeData = false;
      int implNum = stoi(token);
      assert(implNum >= 0);
      timingResults->push_back(new OneStageTimingResult(implNum, impTimes));
      i++;
      delete impTimes;
      impTimes = new TimeVec();
    } else {
      impTimes->push_back(std::stod(token));
    }
  }

  if (i != numImpls + 1) {
    cout << "ERROR: In RuntimeEvaluator::ReadDataFromFile i = " << i;
    cout << " numImpls = " << numImpls << endl;
    LOG_FAIL("replacement for throw call");
    throw;
  }
  return timingResults;
}

void RuntimeEvaluator::Tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ") {
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos)
    {
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      lastPos = str.find_first_not_of(delimiters, pos);
      pos = str.find_first_of(delimiters, lastPos);
    }
}

#endif // DOLLDLA
