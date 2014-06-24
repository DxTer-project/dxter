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

  m_numIterations = 10;
  m_dataFileName = "time_data.txt";
}

std::map<unsigned int, vector<double>> RuntimeEvaluator::ImplementationRuntimeMap(ImplementationMap imps)
{
  return CompileAndRunAllImplementations(imps);
}

std::map<unsigned int, vector<double>> RuntimeEvaluator::CompileAndRunAllImplementations(ImplementationMap imps)
{
  chdir("runtimeEvaluation/");
  std::map<unsigned int, vector<double>> runMap;
  ImplementationMap::iterator it;
  for (it = imps.begin(); it != imps.end(); ++it) {
    vector<double> runtimes = CompileAndRunImplementation(NumImplementationPair(it->first, it->second));
    runMap.insert(NumRuntimePair(it->first, runtimes));
  }
  return runMap;
}

vector<double> RuntimeEvaluator::CompileAndRunImplementation(NumImplementationPair numImp)
{
  WriteImplementationHeaderToDriverFile(m_operationName + "_" + std::to_string(numImp.first) + ".h");
  system("pwd");
  const char *cc = ("gcc -mfpmath=sse -msse3 -o " + m_operationName + " " + m_driverFileName + " utils.c").c_str();
  system(cc);
  const char *run_command = ("./" + m_operationName + " > " + m_dataFileName).c_str();
  system(run_command);
  ClearDriverFile();
  return ReadTimeDataFromFile();
}

void RuntimeEvaluator::WriteImplementationHeaderToDriverFile(string impHeaderName)
{
  string fileContents = "#define NUM_ITERATIONS " + std::to_string(m_numIterations) + "\n";
  fileContents = fileContents + "#include \"" + impHeaderName + "\"\n";
  fileContents = fileContents +  m_driverCode;
  std::ofstream outputFile(m_driverFileName);
  outputFile << fileContents;
  outputFile.close();
}

vector<double> RuntimeEvaluator::ReadTimeDataFromFile()
{
  std::ifstream dataStream(m_dataFileName);
  std::stringstream buffer;
  buffer << dataStream.rdbuf();
  string timeData = buffer.str();
  vector<string> runtimeStrings;
  Tokenize(timeData, runtimeStrings, "\n");
  vector<double> runtimes;
  for (std::vector<string>::iterator it = runtimeStrings.begin(); it != runtimeStrings.end(); ++it) {
    runtimes.push_back(stod(*it));
  }
  return runtimes;
}

void RuntimeEvaluator::Tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ")
{
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  string::size_type pos     = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos)
    {
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      lastPos = str.find_first_not_of(delimiters, pos);
      pos = str.find_first_of(delimiters, lastPos);
    }
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
