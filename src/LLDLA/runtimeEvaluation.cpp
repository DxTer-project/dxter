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

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "runtimeEvaluation.h"

RuntimeEvaluator::RuntimeEvaluator(string evalDirName, string driverFileName, string operationName, string functionPrelude)
{
  m_evalDirName = evalDirName;
  m_driverFileName = driverFileName + ".c";
  std::ifstream t(m_evalDirName + "/" + m_driverFileName);
  std:: stringstream buffer;
  buffer << t.rdbuf();
  m_driverCode = buffer.str();
  m_operationName = operationName;
  m_functionPrelude = functionPrelude;
}

void RuntimeEvaluator::CompileAndRunAllImplementations(ImplementationMap imps)
{
  chdir("runtimeEvaluation/");
  ImplementationMap::iterator it;
  for (it = imps.begin(); it != imps.end(); ++it) {
    CompileAndRunImplementation(NumImplementationPair(it->first, it->second));
  }
}

void RuntimeEvaluator::CompileAndRunImplementation(NumImplementationPair numImp)
{
  WriteImplementationHeaderToDriverFile(m_operationName + "_" + std::to_string(numImp.first) + ".h");
  system("pwd");
  const char *cc = ("gcc -mfpmath=sse -msse3 -o " + m_operationName + " " + m_driverFileName + " utils.c").c_str();
  system(cc);
  const char *run_command = ("./" + m_operationName).c_str();
  system(run_command);
  ClearDriverFile();
}

void RuntimeEvaluator::WriteImplementationHeaderToDriverFile(string impHeaderName)
{
  string fileContents = "#include \"" + impHeaderName + "\"\n" + m_driverCode;
  std::ofstream outputFile(m_driverFileName);
  outputFile << fileContents;
  outputFile.close();
}

void RuntimeEvaluator::ClearDriverFile()
{
  std::ofstream outputFile(m_driverFileName);
  outputFile << m_driverCode;
  outputFile.close();
}

void RuntimeEvaluator::WriteImplementationsToFiles(ImplementationMap imps)
{
  ImplementationMap::iterator it;
  for (it = imps.begin(); it != imps.end(); ++it) {
    unsigned int implementationNum = it->first;
    string implementationStr = m_functionPrelude + "\n" + it->second + "\n}\n";
    string fileName = m_operationName + "_" + std::to_string(implementationNum) + ".h";
    std::ofstream outputFile(m_evalDirName + "/" + fileName);
    outputFile << implementationStr;
    outputFile.close();
  }
}
