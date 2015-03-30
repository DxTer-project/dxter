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
#include "timingResult.h"

#if DOLLDLA

#include "problemInstanceStats.h"

using namespace std;

enum SanityCheckSetting { CHECKALLBUFFERS, NONE };
enum TimingSetting { ONEPHASETIMING, TWOPHASETIMING };

class RuntimeTest
{
 protected:
  vector<string> m_defines;
  vector<string> m_argNames;
  vector<string> m_headers;
  string m_correctTestFileName;
  vector<string> m_operationArgs;
  vector<string> m_argDeclarations;

  int m_chunkSize;
  Type m_type;

  string ToCStatements(vector<string> lines);
  string CArgList(vector<string> args);
  string MakeImpFuncs(ImplementationMap* imps);
  string ImplementationFunctions(ImplementationMap* imps, string referenceImp);
  string MainFunction();
  string AllBufferSanityChecks(ImplementationMap* imps, string referenceImpName);
  string SanityChecks(SanityCheckSetting sanityCheckSetting, ImplementationMap* imps, string referenceImpName);
  string OnePhaseTimingCode(ImplementationMap* imps, string operationName);
  string TimingCode(TimingSetting timingSetting, ImplementationMap* imps, string operationName);
  string HeadersAndDefines(ImplementationMap* imps);
  string MakeFunc(string funcName, string funcBody);
  string CorrectnessCheck(ImplementationMap* imps, string referenceImpName);
  string AllocateArgBuffers(string postfix);
  string FillBuffersWithRandValues(string postfix);
  string CopyArgBuffersTo(string postfix);
  vector<string> ArgBuffers(string postfix);
  string CheckArgBufferDiffs(string refPostfix, string testPostfix, string testName);
  string TimingLoop(int i);
  string TimingLoops(ImplementationMap* imps);
  void AddIncludes();
  void AddMiscellaneousDefines();

 public:
  int m_minCycles;
  string m_dataFileName;
  string m_operationName;

  RuntimeTest(Type m_type, string operationName, vector<string> argNames, vector<string> argDeclarations, vector<string> defines, int numIterations);
  string MakeTestCode(SanityCheckSetting sanityCheckSetting, TimingSetting timingSetting, ImplementationMap* imps, string referenceImp);
};

class RuntimeEvaluator
{

 private:
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
