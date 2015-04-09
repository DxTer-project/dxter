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

#include "problemInstanceStats.h"

using namespace std;

enum SanityCheckSetting { CHECKALLBUFFERS, NONE };
enum TimingSetting { ONEPHASETIMING, TWOPHASETIMING };

class RuntimeTest {
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
  string SetupFunction();
  string SetupFunctions();
  string MainFunction();
  string AllBufferSanityChecks(unsigned int numImpls, string referenceImpName);
  string SanityChecks(SanityCheckSetting sanityCheckSetting, unsigned int numImpls, string referenceImpName);
  string OnePhaseTimingCode(unsigned int numImpls, string operationName);
  string TwoPhaseTimingCode(unsigned int numImpls, string operationName);
  string TimingCode(TimingSetting timingSetting, unsigned int numImpls, string operationName);
  string HeadersAndDefines(unsigned int numImplementations);
  string MakeFunc(string funcName, string funcBody);
  string CorrectnessCheck(unsigned int numImpls, string referenceImpName);
  string AllocateArgBuffers(string postfix);
  string FillBuffersWithRandValues(string postfix);
  string CopyArgBuffersTo(string postfix);
  vector<string> ArgBuffers(string postfix);
  string CheckArgBufferDiffs(string refPostfix, string testPostfix, string testName);
  string TimingLoop(int i);
  string TimingLoops(unsigned int numImpls);
  string TwoPhaseTimingLoop(int i);
  string TwoPhaseTimingLoops(unsigned int numImpls);
  void AddIncludes();
  void AddMiscellaneousDefines();

 public:
  int m_minCycles;
  string m_dataFileName;
  string m_operationName;

  RuntimeTest(Type m_type, string operationName, vector<string> argNames, vector<string> argDeclarations, vector<string> defines, int numIterations);
  string MakeTestCode(SanityCheckSetting sanityCheckSetting, TimingSetting timingSetting, ImplementationMap* imps, string referenceImp);
};

#endif // DOLLDLA

#endif
