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

vector<TimingResult*>* RuntimeEvaluator::EvaluateImplementations(SanityCheckSetting sanityCheckSetting, TimingSetting timingSetting, RuntimeTest test, ImplementationMap* imps, string referenceImp) {
  string executableName = m_evalDirName + "/" + test.m_operationName;
  string testCode = test.MakeTestCode(sanityCheckSetting, timingSetting, imps, referenceImp);

  WriteTestCodeToFile(executableName, testCode);
  CompileTest(executableName);
  RunTest(executableName);
  CleanUpTest(executableName);

  auto timeResults = ReadTimingData(timingSetting, test.m_dataFileName, imps->size());
  return timeResults;
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
    if (IsImplementationSeparator(token)) {
      atEndOfImplTimeData = true;
      timingResults->push_back(new OneStageTimingResult(i, impTimes));
      i++;
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
