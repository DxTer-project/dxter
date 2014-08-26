/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2014, The University of Texas and Bryan Marker

    DxTer is free software: you can redistribute it and/or modify
n    it under the terms of the GNU General Public License as published by
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

#if DOLLDLA

RuntimeTest::RuntimeTest(Type type, string operationName, vector<string> argNames, vector<string> argDeclarations, vector<string> defines, int numIterations, int chunkSize)
{
  m_type = type;
  m_operationName = operationName;
  m_argNames = argNames;
  m_argDeclarations = argDeclarations;
  m_defines = defines;
  m_numIterations = numIterations;
  m_chunkSize = chunkSize;
  m_dataFileName = m_operationName + "_time_data";
  m_correctTestFileName = m_operationName + "_correctness_test_results";
  AddIncludes();
  AddMiscellaneousDefines();
}

void RuntimeTest::AddIncludes()
{
  m_headers.push_back("#include <immintrin.h>");
  m_headers.push_back("#include \"utils.h\"");
  m_headers.push_back("#include <string.h>");
  m_headers.push_back("#include <unistd.h>");
  return;
}

void RuntimeTest::AddMiscellaneousDefines()
{
  m_defines.push_back("#define BUF_SIZE 1000000");
  m_defines.push_back("#define NUM_ITERATIONS " + std::to_string((long long int) m_numIterations));
  m_defines.push_back("#define CHUNK_SIZE " + std::to_string((long long int) m_chunkSize));
  m_defines.push_back("#define min(a,b) ((a) < (b) ? (a) : (b))");
  m_defines.push_back("#define ALLOC_BUFFER(size) alloc_aligned_16((size))");
  m_defines.push_back("#define MUVALUE " + std::to_string((long long int) arch->VecRegWidth(m_type)) + "\n");
  m_defines.push_back(arch->VecRegTypeDec(m_type));
  if (m_type == REAL_SINGLE) {
    m_defines.push_back("#define NUM_SIZE sizeof(float)");
    m_defines.push_back("#define FILL_WITH_RAND_VALUES(size, buf) rand_floats((size), (buf))");
    m_defines.push_back("#define TEST_BUFFER_DIFF(size, b1, b2, test_name) test_buffer_diff_float((size), (b1), (b2), (test_name))");
    m_defines.push_back("#define COPY_BUFFER(size, b1, b2) copy_buffer_float((size), (b1), (b2))");
  } else {
    m_defines.push_back("#define NUM_SIZE sizeof(double)");
    m_defines.push_back("#define FILL_WITH_RAND_VALUES(size, buf) rand_doubles((size), (buf))");
    m_defines.push_back("#define TEST_BUFFER_DIFF(size, b1, b2, test_name) test_buffer_diff((size), (b1), (b2), (test_name))");
    m_defines.push_back("#define COPY_BUFFER(size, b1, b2) copy_buffer((size), (b1), (b2))");
  }

  return;
}

string RuntimeTest::MakeTestCode(ImplementationMap imps)
{
  int numImplementations = imps.size();
  m_defines.push_back("#define NUM_ALGS " + std::to_string((long long int) numImplementations));
  string headersAndDefines = ToCStatements(m_headers) + "\n" + ToCStatements(m_defines);
  string implementationFunctions = MakeImpFuncs(imps);
  string driverCode = MainFuncCode(imps);
  string testCode = headersAndDefines + "\n" + implementationFunctions + "\n" + driverCode;
  return testCode;
}

string RuntimeTest::MakeTestCodeWithCorrectnessCheck(ImplementationMap imps, string referenceImp)
{
  int numImplementations = imps.size();
  m_defines.push_back("#define NUM_ALGS " + std::to_string((long long int) numImplementations));
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
  string argBufferAllocation = "\tprintf(\"Starting buffer allocation\\n\");\n";
  argBufferAllocation += AllocateArgBuffers("");
  argBufferAllocation += "FILE *" + m_dataFileName + " = fopen(\"" + m_dataFileName + "\", \"w\");\n";
  argBufferAllocation += "\tprintf(\"Done with allocation\\n\");\n";
  string timingSetup = "\tint i, j, k;\n\tstruct timeval begin, end;\n\tdouble exec_time;\n";
  string mainFunc = prototype + argBufferAllocation + "\n" + timingSetup;
  string timingLoop = TimingLoop(imps);
  timingLoop += "\n\tfclose(" + m_dataFileName + ");\n";
  mainFunc = mainFunc + "\n" + timingLoop + "\n}";
  return mainFunc;
}

string RuntimeTest::MainFuncCodeWithCorrectnessCheck(ImplementationMap imps, string referenceImpName)
{
  string prototype = "int main() {\n";
  string argBufferAllocation = "\tprintf(\"Starting buffer allocation\\n\");\n";
  argBufferAllocation += AllocateArgBuffers("") + "\n";
  argBufferAllocation += AllocateArgBuffers("_ref") + "\n";
  argBufferAllocation += AllocateArgBuffers("_test") + "\n";
  argBufferAllocation += FillBuffersWithRandValues("") + "\n";
  argBufferAllocation += "\tFILE *" + m_dataFileName + " = fopen(\"" + m_dataFileName + "\", \"w\");\n";
  argBufferAllocation += "\tprintf(\"Done with allocation\\n\");\n";
  string correctnessCheck = argBufferAllocation + CopyArgBuffersTo("_ref") + "\n";
  correctnessCheck += CorrectnessCheck(imps, referenceImpName);
  string timingSetup = "\tint i, j, k;\n\tstruct timeval begin, end;\n\tdouble exec_time;\n";
  string mainFunc = prototype + correctnessCheck + "\n" + timingSetup;
  string timingLoop = TimingLoop(imps);
  timingLoop += "\n\tfclose(" + m_dataFileName + ");\n";
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
    correctnessCheck += m_operationName + "_" + std::to_string((long long int) i) + "(" + CArgList(testBuffers) + ");\n";
    correctnessCheck += CheckArgBufferDiffs("_ref", "_test",std::to_string((long long int) i)) + "\n";
  }
  return correctnessCheck;
}

string RuntimeTest::CheckArgBufferDiffs(string refPostfix, string testPostfix, string testName)
{
  vector<string>::iterator argIter;
  vector<string> diffChecks;
  for (argIter = m_argNames.begin(); argIter != m_argNames.end(); ++argIter) {
    string refBuf = *argIter + refPostfix;
    string testBuf = *argIter + testPostfix;
    string diffCheck = "\tTEST_BUFFER_DIFF(BUF_SIZE, " + refBuf + ", " + testBuf + ", \"Sanity Check " + testName + "\");";
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
    string opName = m_operationName + "_" + std::to_string((long long int) i);
    loopBody += "\tfor (j = 0; j < NUM_ITERATIONS; j++) {\n";
    loopBody += "\t\tgettimeofday( &begin, NULL );\n";
    loopBody += "\t\tfor (k = 0; k < CHUNK_SIZE; k++) {\n";
    loopBody += "\t\t\t" + opName + "(" + CArgList(m_argNames) + ");\n";
    loopBody += "\t\t}\n";
    loopBody += "\t\tgettimeofday( &end, NULL );\n";
    loopBody += "\t\texec_time = end.tv_sec - begin.tv_sec;\n";
    loopBody += "\t\texec_time += (end.tv_usec - begin.tv_usec) * 1.0e-6;\n";
    loopBody += "\t\tchar exec_time_str[100];\n";
    loopBody += "\t\tsprintf(exec_time_str, \"%f\\n\", exec_time);\n";
    loopBody += "\t\tsize_t trash = fprintf(" + m_dataFileName + ", \"%s\", exec_time_str);\n";
    loopBody += "\t}\n";
    loopBody += "\tprintf(\"Done evaluating " + opName + "\\n\");\n";
  }
  return loopBody;
}


string RuntimeTest::FillBuffersWithRandValues(string postfix)
{
  std::vector<string> bufferFills;
  std::vector<string>::iterator argIter;
  for (argIter = m_argNames.begin(); argIter != m_argNames.end(); ++argIter) {
    string bufName = *argIter;
    bufferFills.push_back("\tFILL_WITH_RAND_VALUES(BUF_SIZE, " + bufName + ");");
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
    bufferFills.push_back("\tCOPY_BUFFER(BUF_SIZE, " + srcBuffer + ", " + destBuffer + ");");
  }
  return ToCStatements(bufferFills);
}

string RuntimeTest::AllocateArgBuffers(string postfix)
{
  std::vector<string> argAllocs;
  std::vector<string>::iterator argIter;
  for (argIter = m_argDeclarations.begin(); argIter != m_argDeclarations.end(); ++argIter) {
    string bufName = *argIter + postfix;
    argAllocs.push_back("\t" + bufName + " = ALLOC_BUFFER(BUF_SIZE * sizeof(NUM_SIZE));");
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
    string funcName = m_operationName + "_" + std::to_string((long long int) impIt->first);
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
}

ImplementationRuntimeMap RuntimeEvaluator::EvaluateImplementations(RuntimeTest test, ImplementationMap imps)
{
  m_numIterations = test.m_numIterations;
  string executableName = test.m_operationName;
  string testFileName = executableName + ".c";
  std::ofstream outStream(m_evalDirName + "/" + testFileName);
  outStream << test.MakeTestCode(imps);
  outStream.close();
  cout << "All implementations written to files\n";
  const char *evalDir = (m_evalDirName + "/").c_str();
  chdir(evalDir);
  system(arch->CompileString(executableName, testFileName).c_str());
  cout << "Compiled\n";
  string runStr = "./" + executableName;
  system(runStr.c_str());
  string removeExecutable = "rm -f " + executableName;
  system(removeExecutable.c_str());
  return ReadTimeDataFromFile(test.m_dataFileName, imps.size());
}

ImplementationRuntimeMap RuntimeEvaluator::EvaluateImplementationsWithCorrectnessCheck(RuntimeTest test, ImplementationMap imps, string referenceImp)
{
  m_numIterations = test.m_numIterations;
  string executableName = test.m_operationName;
  string testFileName = executableName + ".c";
  std::ofstream outStream(m_evalDirName + "/" + testFileName);
  outStream << test.MakeTestCodeWithCorrectnessCheck(imps, referenceImp);
  outStream.close(); 
  cout << "All implementations written to files\n";
  const char *evalDir = (m_evalDirName + "/").c_str();
  chdir(evalDir);
  system(arch->CompileString(executableName, testFileName).c_str());
  string runStr = "./" + executableName;
  system(runStr.c_str());
  string removeExecutable = "rm -f " + executableName;
  system(removeExecutable.c_str());
  return ReadTimeDataFromFile(test.m_dataFileName, imps.size());  
}

ImplementationRuntimeMap RuntimeEvaluator::ReadTimeDataFromFile(string fileName, int numImpls)
{
  std::ifstream dataStream(fileName);
  std::stringstream buffer;
  buffer << dataStream.rdbuf();
  string timeData = buffer.str();
  ImplementationRuntimeMap runtimeMap;
  std::vector<string> runtimeStrings;
  Tokenize(timeData, runtimeStrings, "\n");
  std::vector<string>::iterator it = runtimeStrings.begin();
  int i, j;
  for (i = 1; i <= numImpls; i++) {
    TimeVec impTimes;
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

#endif // DOLLDLA
