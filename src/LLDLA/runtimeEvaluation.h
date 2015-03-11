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

#if DOLLDLA

using namespace std;

class RuntimeTest
{
 public:
  vector<string> m_defines;
  vector<string> m_argNames;
  vector<string> m_headers;
  string m_operationName;
  string m_dataFileName;
  string m_correctTestFileName;
  vector<string> m_operationArgs;
  vector<string> m_argDeclarations;
  int m_minCycles;
  int m_chunkSize;
  Type m_type;

  RuntimeTest(Type m_type, string operationName, vector<string> argNames, vector<string> argDeclarations, vector<string> defines, int numIterations);
  //  string MakeTestCode(ImplementationMap* imps);
  string MakeTestCodeWithCorrectnessCheck(ImplementationMap* imps, string referenceImp);
  string ToCStatements(vector<string> lines);
  string CArgList(vector<string> args);
  string MakeImpFuncs(ImplementationMap* imps);
  string MakeFunc(string funcName, string funcBody);
  //  string MainFuncCode(ImplementationMap* imps);
  string MainFuncCodeWithCorrectnessCheck(ImplementationMap* imps, string referenceImpName);
  string CorrectnessCheck(ImplementationMap* imps, string referenceImpName);
  string AllocateArgBuffers(string postfix);
  string FillBuffersWithRandValues(string postfix);
  string CopyArgBuffersTo(string postfix);
  vector<string> ArgBuffers(string postfix);
  string CheckArgBufferDiffs(string refPostfix, string testPostfix, string testName);
  string TimingLoop(ImplementationMap* imps);
  void AddIncludes();
  void AddMiscellaneousDefines();
};

class RuntimeEvaluator
{

 private:
  bool IsImplementationSeparator(string token);

 public:
  string m_evalDirName;
  int m_minCycles;

  RuntimeEvaluator(string evalDirName);
  //  ImplementationRuntimeMap EvaluateImplementations(RuntimeTest test, ImplementationMap* imps);
  ImplementationRuntimeMap EvaluateImplementationsWithCorrectnessCheck(RuntimeTest test, ImplementationMap* imps, string referenceImp);
  ImplementationRuntimeMap ReadTimeDataFromFile(string fileName, int numImpls);
  void Tokenize(const string& str, vector<string>& tokens, const string& delimiters);
};

#endif // DOLLDLA

#endif // RUNTIME_EVALUATION_H_
