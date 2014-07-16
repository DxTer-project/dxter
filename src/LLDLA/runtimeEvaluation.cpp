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
nnn    along with DxTer.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "runtimeEvaluation.h"


RuntimeTest::RuntimeTest(string operationName, vector<string> argNames, vector<string> argDeclarations, vector<string> defines, int numIterations, int chunkSize)
{
  m_operationName = operationName;
  m_argNames = argNames;
  m_argDeclarations = argDeclarations;
  m_defines = defines;
  m_numIterations = numIterations;
  m_chunkSize = chunkSize;
  m_headers.push_back("#include \"row_stride_lldla_primitives.h\"");
  m_headers.push_back("#include \"gen_stride_lldla_primitives.h\"");
  m_headers.push_back("#include \"utils.h\"");
  m_headers.push_back("#include <string.h>");
  m_headers.push_back("#include <unistd.h>");
  m_defines.push_back("#define BUF_SIZE 1000000");
  m_defines.push_back("#define NUM_ITERATIONS " + std::to_string(m_numIterations));
  m_defines.push_back("#define CHUNK_SIZE " + std::to_string(m_chunkSize));
  m_defines.push_back("#define MUVALUE 2");
  m_defines.push_back("#define min(a,b) ((a) < (b) ? (a) : (b))");
}

string RuntimeTest::MakeTestCode(ImplementationMap imps)
{
  int numImplementations = imps.size();
  m_defines.push_back("#define NUM_ALGS " + std::to_string(numImplementations));
  string headersAndDefines = ToCStatements(m_headers) + "\n" + ToCStatements(m_defines);
  string implementationFunctions = MakeImpFuncs(imps);
  string driverCode = MainFuncCode(imps);
  string testCode = headersAndDefines + "\n" + implementationFunctions + "\n" + driverCode;
  return testCode;
}

string RuntimeTest::MakeTestCodeWithCorrectnessCheck(ImplementationMap imps, string referenceImp)
{
  int numImplementations = imps.size();
  m_defines.push_back("#define NUM_ALGS " + std::to_string(numImplementations));
  string headersAndDefines = ToCStatements(m_headers) + "\n" + ToCStatements(m_defines);
  string refFuncName = m_operationName + "_test";
  string implementationFunctions = MakeImpFuncs(imps) + MakeFunc(refFuncName, referenceImp);
  string driverCode = MainFuncCodeWithCorrectnessCheck(imps, m_operationName + "_test");
  string testCode = headersAndDefines + "\n" + implementationFunctions + "\n" + driverCode;
  return testCode;
}

string RuntimeTest::MainFuncCode(ImplementationMap imps)
{
  string prototype = "int main() {\n";
  string argBufferAllocation = AllocateArgBuffers("");
  string timingSetup = "\tint i, j, k;\n\tclock_t begin, end;\n\tdouble exec_time;\n";    
  string mainFunc = prototype + argBufferAllocation + "\n" + timingSetup;
  string timingLoop = TimingLoop(imps);
  mainFunc = mainFunc + "\n" + timingLoop + "\n}";
  return mainFunc;
}

string RuntimeTest::MainFuncCodeWithCorrectnessCheck(ImplementationMap imps, string referenceImpName)
{
  string prototype = "int main() {\n";
  string argBufferAllocation = AllocateArgBuffers("") + "\n";
  argBufferAllocation += AllocateArgBuffers("_ref") + "\n";
  argBufferAllocation += AllocateArgBuffers("_test") + "\n";
  argBufferAllocation += FillBuffersWithRandValues("") + "\n";
  string correctnessCheck = argBufferAllocation + CopyArgBuffersTo("_ref") + "\n";
  correctnessCheck += CorrectnessCheck(imps, referenceImpName);
  string timingSetup = "\tint i, j, k;\n\tclock_t begin, end;\n\tdouble exec_time;\n";
  string mainFunc = prototype + correctnessCheck + "\n" + timingSetup;
  string timingLoop = TimingLoop(imps);
  mainFunc = mainFunc + "\n" + timingLoop + "\n}";
  return mainFunc;
}

string RuntimeTest::CorrectnessCheck(ImplementationMap imps, string referenceImpName)
{
  int i;
  string correctnessCheck = "";
  vector<string> argBuffers = ArgBuffers("_ref");
  vector<string> testBuffers = ArgBuffers("_test");
  string callRefImp = "\t" + referenceImpName + "(" + CArgList(argBuffers) + ");\n\n";
  correctnessCheck += callRefImp;
  for (i = 1; i <= imps.size(); i++) {
    correctnessCheck += CopyArgBuffersTo("_test") + "\n\t";
    correctnessCheck += m_operationName + "_" + std::to_string(i) + "(" + CArgList(testBuffers) + ");\n";
    correctnessCheck += CheckArgBufferDiffs("_ref", "_test") + "\n";
  }
  return correctnessCheck;
}

string RuntimeTest::CheckArgBufferDiffs(string refPostfix, string testPostfix)
{
  vector<string>::iterator argIter;
  vector<string> diffChecks;
  for (argIter = m_argNames.begin(); argIter != m_argNames.end(); ++argIter) {
    string refBuf = *argIter + refPostfix;
    string testBuf = *argIter + testPostfix;
    string diffCheck = "\ttest_buffer_diff(BUF_SIZE, " + refBuf + ", " + testBuf + ", \"Sanity Check\");";
    diffChecks.push_back(diffCheck);
  }
  return ToCStatements(diffChecks);
}

vector<string> RuntimeTest::ArgBuffers(string postfix)
{
  vector<string> argBufs;
  vector<string>::iterator argIter;
  for (argIter = m_argNames.begin(); argIter != m_argNames.end(); ++argIter) {
    string name = *argIter;
    argBufs.push_back(name + postfix);
  }
  return argBufs;
}

string RuntimeTest::TimingLoop(ImplementationMap imps)
{
  int i;
  string loopBody = "";
  for (i = 1; i <= imps.size(); i++) {
    string opName = m_operationName + "_" + std::to_string(i);
    //    loopBody += "\tprintf(\"Starting test of " + opName + "\\n\");\n";
    loopBody += "\tfor (j = 0; j < NUM_ITERATIONS; j++) {\n";
    loopBody += "\t\tbegin = clock();\n";
    loopBody += "\t\tfor (k = 0; k < CHUNK_SIZE; k++) {\n";
    loopBody += "\t\t\t" + opName + "(" + CArgList(m_argNames) + ");\n";
    loopBody += "\t\t}\n";
    loopBody += "\t\tend = clock();\n";
    loopBody += "\t\texec_time = (double) (end - begin);\n";
    loopBody += "\t\tchar exec_time_str[100];\n";
    loopBody += "\t\tsprintf(exec_time_str, \"%f\\n\", exec_time);\n";
    loopBody += "\t\tsize_t trash = write(1, exec_time_str, strlen(exec_time_str));\n";
    loopBody += "\t}\n";
  }
  return loopBody;
}


string RuntimeTest::FillBuffersWithRandValues(string postfix)
{
  std::vector<string> bufferFills;
  std::vector<string>::iterator argIter;
  for (argIter = m_argNames.begin(); argIter != m_argNames.end(); ++argIter) {
    string bufName = *argIter;
    bufferFills.push_back("\trand_doubles(BUF_SIZE, " + bufName + ");");
  }
  return ToCStatements(bufferFills);
}

string RuntimeTest::CopyArgBuffersTo(string postfix)
{
  std::vector<string> bufferFills;
  std::vector<string>::iterator argIter;
  for (argIter = m_argNames.begin(); argIter != m_argNames.end(); ++argIter) {
    string srcBuffer = *argIter;
    string destBuffer = *argIter + postfix;
    bufferFills.push_back("\tcopy_buffer(BUF_SIZE, " + srcBuffer + ", " + destBuffer + ");");
  }
  return ToCStatements(bufferFills);
}

string RuntimeTest::AllocateArgBuffers(string postfix)
{
  std::vector<string> argAllocs;
  std::vector<string>::iterator argIter;
  for (argIter = m_argDeclarations.begin(); argIter != m_argDeclarations.end(); ++argIter) {
    string bufName = *argIter + postfix;
    argAllocs.push_back("\t" + bufName + " = alloc_aligned_16(BUF_SIZE * sizeof(double));");
  }
  return ToCStatements(argAllocs);
}
string RuntimeTest::ToCStatements(vector<string> lines)
{
  string cStatements = "";
  std::vector<string>::iterator lineIt;
  for (lineIt = lines.begin(); lineIt != lines.end(); ++lineIt) {
    cStatements = cStatements + *lineIt + "\n";
  }
  return cStatements;
}

string RuntimeTest::CArgList(vector<string> args)
{
  string argList = "";
  std::vector<string>::iterator argsIt;
  int i = 0;
  for (argsIt = args.begin(); argsIt != args.end(); ++argsIt) {
    argList = argList + *argsIt;
    if (i < args.size() - 1) {
      argList = argList + ", ";
    }
    i++;
  }
  return argList;
}

string RuntimeTest::MakeImpFuncs(ImplementationMap imps)
{
  string endOfFuncDec = "(" + CArgList(m_argDeclarations) + ")";
  string allImplementationFuncs = "";
  ImplementationMap::iterator impIt;
  for (impIt = imps.begin(); impIt != imps.end(); ++impIt) {
    string funcName = m_operationName + "_" + std::to_string(impIt->first);
    string funcBody = impIt->second;
    string funcDec = MakeFunc(funcName, funcBody);
    allImplementationFuncs = allImplementationFuncs + funcDec;
  }
  return allImplementationFuncs;
}

string RuntimeTest::MakeFunc(string funcName, string funcBody) {
  string funcDec = "void " + funcName;
  funcDec = funcDec + "(" + CArgList(m_argDeclarations) + ")" + "{\n";
  funcDec = funcDec + funcBody + "return;\n}\n";
  return funcDec;
}

RuntimeEvaluator::RuntimeEvaluator(string evalDirName)
{
  m_evalDirName = evalDirName;
  m_numIterations = 0;
  m_dataFileName = "time_data.txt";
}

std::map<unsigned int, vector<double>> RuntimeEvaluator::EvaluateImplementations(RuntimeTest test, ImplementationMap imps)
{
  m_numIterations = test.m_numIterations;
  string executableName = test.m_operationName;
  string testFileName = executableName + ".c";
  std::ofstream outStream(m_evalDirName + "/" + testFileName);
  outStream << test.MakeTestCode(imps);
  outStream.close(); 
  const char *evalDir = (m_evalDirName + "/").c_str();
  chdir(evalDir);
  string compileStr = "gcc -O3 -mfpmath=sse -msse3 -finline-functions -o " +  executableName;
  compileStr += " " + testFileName + " utils.c";
  system(compileStr.c_str());
  string runStr = "./" + executableName + " > " + m_dataFileName;
  system(runStr.c_str());
  return ReadTimeDataFromFile(imps.size());
}

std::map<unsigned int, vector<double>> RuntimeEvaluator::EvaluateImplementationsWithCorrectnessCheck(RuntimeTest test, ImplementationMap imps, string referenceImp)
{
  m_numIterations = test.m_numIterations;
  string executableName = test.m_operationName;
  string testFileName = executableName + ".c";
  std::ofstream outStream(m_evalDirName + "/" + testFileName);
  outStream << test.MakeTestCodeWithCorrectnessCheck(imps, referenceImp);
  outStream.close(); 
  const char *evalDir = (m_evalDirName + "/").c_str();
  chdir(evalDir);
  string compileStr = "gcc -O3 -mfpmath=sse -msse3 -finline-functions -o " +  executableName;
  compileStr += " " + testFileName + " utils.c";
  system(compileStr.c_str());
  string runStr = "./" + executableName + " > " + m_dataFileName;
  system(runStr.c_str());
  return ReadTimeDataFromFile(imps.size());  
}

std::map<unsigned int, vector<double>> RuntimeEvaluator::ReadTimeDataFromFile(int numImpls)
{
  std::ifstream dataStream(m_dataFileName);
  std::stringstream buffer;
  buffer << dataStream.rdbuf();
  string timeData = buffer.str();
  std::map<unsigned int, vector<double>> runtimeMap;
  std::vector<string> runtimeStrings;
  Tokenize(timeData, runtimeStrings, "\n");
  std::vector<string>::iterator it = runtimeStrings.begin();
  int i, j;
  for (i = 1; i <= numImpls; i++) {
    vector<double> impTimes;
    for (j = 0; j < m_numIterations; j++) {
      impTimes.push_back(std::stod(*it));
      it++;
    }
    runtimeMap.insert(NumRuntimePair(i, impTimes));
  }
  return runtimeMap;
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
