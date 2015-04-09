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

#ifndef RUNTIME_EVALUATION_H_
#define RUNTIME_EVALUATION_H_

#include <map>
#include <vector>

#include "LLDLA.h"
#include "base.h"
#include "runnerUtils.h"
#include "runtimeTest.h"
#include "timingResult.h"

#if DOLLDLA

class RuntimeEvaluator
{

 protected:
  bool IsImplementationSeparator(string token);
  void CompileTest(string executableName);
  void RunTest(string executableName);
  void CleanUpTest(string executableName);
  vector<TimingResult*>* ReadTimingData(TimingSetting timingSetting, string dataFileName, int numImpls);
  void WriteTestCodeToFile(string testFileName, string testCode);
  void WriteTestCodeToFile(SanityCheckSetting sanityCheckSetting, TimingSetting timingSetting, string executableName, string testFileName);

 public:
  string m_evalDirName;

  RuntimeEvaluator(string evalDirName);
  vector<TimingResult*>* EvaluateImplementations(SanityCheckSetting sanityCheckSetting, TimingSetting timingSetting, RuntimeTest test, ImplementationMap* imps, string referenceImp);

  vector<TimingResult*>* ReadTimeDataFromFile(TimingSetting timingSetting, string fileName, int numImpls);
  void Tokenize(const string& str, vector<string>& tokens, const string& delimiters);
};

#endif // DOLLDLA

#endif // RUNTIME_EVALUATION_H_
