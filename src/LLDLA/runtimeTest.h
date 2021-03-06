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

#ifndef RUNTIME_TEST_H
#define RUNTIME_TEST_H

#include "LLDLA.h"

#if DOLLDLA

#include "lldlaUniverse.h"
#include "problemInstanceStats.h"

using namespace std;

enum SanityCheckSetting { CHECKALLBUFFERS, CHECKOUTPUTBUFFERS, NONE };
enum TimingSetting { ONEPHASETIMING, TWOPHASETIMING };

class RuntimeTest {
 protected:
  vector<string> m_defines;
  vector<string> m_argNames;
  vector<string> m_outputNames;
  vector<string> m_headers;
  string m_correctTestFileName;
  vector<string> m_operationArgs;
  vector<string> m_argDeclarations;

  int m_chunkSize;
  Type m_type;

  string ToCStatements(vector<string> lines);
  string CArgList(vector<string> args);
  string MakeImpFuncs(vector<pair<GraphNum, ImplInfo>>* imps);
  string ImplementationFunctions(vector<pair<GraphNum, ImplInfo>>* imps, string referenceImp);
  string SetupFunction();
  string SetupFunctions();
  string MainFunction();
  string SanityCheckBufferAllocation();
  string OutputBufferCorrectnessCheck(vector<pair<GraphNum, ImplInfo>>* imps, string referenceImpName);
  string OutputBufferSanityChecks(vector<pair<GraphNum, ImplInfo>>* imps, string referenceImpName);
  string AllBufferSanityChecks(vector<pair<GraphNum, ImplInfo>>* imps, string referenceImpName);
  string SanityChecks(SanityCheckSetting sanityCheckSetting, vector<pair<GraphNum, ImplInfo>>* imps, string referenceImpName);
  string OnePhaseTimingCode(vector<pair<GraphNum, ImplInfo>>* imps, string operationName);
  string TwoPhaseTimingCode(vector<pair<GraphNum, ImplInfo>>* imps, string operationName);
  string TimingCode(TimingSetting timingSetting, vector<pair<GraphNum, ImplInfo>>* imps, string operationName);
  string HeadersAndDefines(unsigned int numImplementations);
  string MakeFunc(string funcName, string funcBody);
  string CorrectnessCheck(vector<pair<GraphNum, ImplInfo>>* imps, string referenceImpName);
  string AllocateArgBuffers(const vector<string> args, string postfix);
  string FillBuffersWithRandValues(const vector<string> argNames, string postfix);
  string CopyArgBuffersTo(string postfix);
  vector<string> ArgBuffers(string postfix);
  string CheckArgBufferDiffs(const vector<string> args, string refPostfix, string testPostfix, string testName);
  string TimingLoop(unsigned int i);
  string TimingLoops(vector<pair<GraphNum, ImplInfo>>* imps);
  string TwoPhaseTimingLoop(unsigned int i);
  string TwoPhaseTimingLoops(vector<pair<GraphNum, ImplInfo>>* imps);
  void AddIncludes();
  void AddMiscellaneousDefines();

 public:
  int m_minCycles;
  string m_dataFileName;
  string m_operationName;

  RuntimeTest(ProblemInstance* prob, LLDLAUniverse* uni, unsigned int minCycles);
  string MakeTestCode(SanityCheckSetting sanityCheckSetting, TimingSetting timingSetting, vector<pair<GraphNum, ImplInfo>>* imps, string referenceImp);
};

#endif // DOLLDLA

#endif
