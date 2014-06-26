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

#include "LLDLA.h"
#include <map>
#include <vector>

using namespace std;

typedef std::map<unsigned int, string> ImplementationMap;
typedef std::map<unsigned int, vector<double>> ImplementationRuntimeMap;
typedef std::pair<unsigned int, string> NumImplementationPair;
typedef std::pair<unsigned int, vector<double>> NumRuntimePair;

class RuntimeTest
{
 public:
  vector<string> m_defines;
  vector<string> m_argNames;
  vector<string> m_headers;
  string m_operationName;
  vector<string> m_operationArgs;
  vector<string> m_argDeclarations;

  RuntimeTest(string operationName, vector<string> argNames, vector<string> argDeclarations, vector<string> defines);
  string MakeTestCode(ImplementationMap imps);
  string ToCStatements(vector<string> lines);
  string CArgList(vector<string> args);
  string MakeImpFuncs(ImplementationMap imps);
  string MainFuncCode(ImplementationMap imps);
  string AllocateArgBuffers();
  string TimingLoop(ImplementationMap imps);
};


class RuntimeEvaluator
{
 public:
  string m_evalDirName;
  string m_dataFileName;
  int m_numIterations;

  RuntimeEvaluator(string evalDirName, int numIterations);
  std::map<unsigned int, vector<double>> EvaluateImplementations(RuntimeTest test, ImplementationMap imps);
  std::map<unsigned int, vector<double>> ReadTimeDataFromFile(int numImpls, int numIters);
  void Tokenize(const string& str, vector<string>& tokens, const string& delimiters);
};
