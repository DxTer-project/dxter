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

class RuntimeEvaluator
{
 public:
  string m_evalDirName;
  string m_driverCode;
  string m_driverFileName;
  string m_operationName;
  string m_functionPrelude;
  int m_numIterations;
  string m_dataFileName;

  RuntimeEvaluator(string evalDirName, string driverFileName, string operationName, string functionPrelude);
  std::map<unsigned int, vector<double>> ImplementationRuntimeMap(ImplementationMap imps);
  std::map<unsigned int, vector<double>> CompileAndRunAllImplementations(ImplementationMap imps);
  vector<double> CompileAndRunImplementation(NumImplementationPair numImp);
  void WriteImplementationHeaderToDriverFile(string impHeaderName);
  void WriteImplementationsToFiles(ImplementationMap imps);
  vector<double> ReadTimeDataFromFile();
  void Tokenize(const string& str, vector<string>& tokens, const string& delimiters);
  void ClearDriverFile();
};

class RuntimeTest
{
 public:
  vector<string> m_defines;
  vector<string> m_headers;
  string m_operationName;
  vector<string> m_operationArgs;
  vector<string> m_argDeclarations;

  RuntimeTest(string operationName, vector<string> argDeclarations, vector<string> defines);
  string MakeTestCode(ImplementationMap imps);
  string ToCStatements(vector<string> lines);
  string CArgList(vector<string> args);
  string MakeImpFuncs(ImplementationMap imps);
  string MainFuncCode(ImplementationMap imps);
  string AllocateArgBuffers();
};
