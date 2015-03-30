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

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "runtimeEvaluation.h"

#if DOLLDLA

#include "oneStageTimingResult.h"

RuntimeTest::RuntimeTest(Type type, string operationName, vector<string> argNames, vector<string> argDeclarations, vector<string> defines, int minCycles)
{
  m_type = type;
  m_operationName = operationName;
  m_argNames = argNames;
  m_argDeclarations = argDeclarations;
  m_defines = defines;
  m_minCycles = minCycles;
  m_dataFileName = m_operationName + "_time_data";
  m_correctTestFileName = m_operationName + "_correctness_test_results";
  AddIncludes();
  AddMiscellaneousDefines();
}

void RuntimeTest::AddIncludes()
{
  m_headers.push_back("#include <immintrin.h>");
  m_headers.push_back("#include \"rdtsc.h\"");
  m_headers.push_back("#include <string.h>");
  m_headers.push_back("#include \"utils.h\"");
  m_headers.push_back("#include <unistd.h>");
  return;
}

void RuntimeTest::AddMiscellaneousDefines()
{
  cout << "\n\n\n\n\n\n\nAdding misc defines\n\n\n\n\n";
  m_defines.push_back("#define BUF_SIZE 1000000");
  m_defines.push_back("#define MIN_CYCLES " + std::to_string((long long int) m_minCycles));
  m_defines.push_back("#define min(a,b) ((a) < (b) ? (a) : (b))");
  m_defines.push_back("#define ALLOC_BUFFER(size) alloc_aligned_16((size))");
  m_defines.push_back("#define MUVALUE " + std::to_string((long long int) arch->VecRegWidth(m_type)) + "\n");
  m_defines.push_back(arch->VecRegTypeDec(m_type));
  if (m_type == REAL_SINGLE) {
    cout << "Datatype is real single\n";
    m_defines.push_back("#define NUM_SIZE sizeof(float)");
    m_defines.push_back("#define FILL_WITH_RAND_VALUES(size, buf) rand_floats((size), (buf))");
    m_defines.push_back("#define TEST_BUFFER_DIFF(size, b1, b2, test_name) test_buffer_diff_float((size), (b1), (b2), (test_name))");
    m_defines.push_back("#define COPY_BUFFER(size, b1, b2) copy_buffer_float((size), (b1), (b2))");
  } else if (m_type == REAL_DOUBLE) {
    cout << "Datatype is real double\n";
    m_defines.push_back("#define NUM_SIZE sizeof(double)");
    m_defines.push_back("#define FILL_WITH_RAND_VALUES(size, buf) rand_doubles((size), (buf))");
    m_defines.push_back("#define TEST_BUFFER_DIFF(size, b1, b2, test_name) test_buffer_diff((size), (b1), (b2), (test_name))");
    m_defines.push_back("#define COPY_BUFFER(size, b1, b2) copy_buffer((size), (b1), (b2))");
  } else {
    cout << "Error: Unsupported type for runtime test data\n";
    LOG_FAIL("replacement for throw call");
    throw;
  }
  return;
}

string RuntimeTest::HeadersAndDefines(ImplementationMap* imps) {
  int numImplementations = imps->size();
  m_defines.push_back("#define NUM_ALGS " + std::to_string((long long int) numImplementations));
  string headersAndDefines = ToCStatements(m_headers) + "\n" + ToCStatements(m_defines);
  return headersAndDefines;
}

string RuntimeTest::ImplementationFunctions(ImplementationMap* imps, string referenceImp) {
  string refFuncName = m_operationName + "_test";
  string implementationFunctions = MakeImpFuncs(imps) + MakeFunc(refFuncName, referenceImp);
  return implementationFunctions;
}

string RuntimeTest::MainFunction() {
  string mainFunc = "";
  mainFunc += "int main() {\n";
  mainFunc += "\tsanity_check_implementations();\n";
  mainFunc += "\ttime_implementations();\n";
  mainFunc += "}";
  return mainFunc;
}

string RuntimeTest::SanityChecks(SanityCheckSetting sanityCheckSetting, ImplementationMap* imps, string referenceImpName) {
  if (sanityCheckSetting == CHECKALLBUFFERS) {
    return AllBufferSanityChecks(imps, referenceImpName);
  } else {
    cout << "SanityCheckSetting NONE is not yet supported" << endl;
    LOG_FAIL("replacement for throw call");
    throw;
  }
}

string RuntimeTest::AllBufferSanityChecks(ImplementationMap* imps, string referenceImpName) {
  string prototype = "void sanity_check_implementations() {\n";
  string argBufferAllocation = "\tprintf(\"Starting buffer allocation\\n\");\n";
  argBufferAllocation += AllocateArgBuffers("") + "\n";
  argBufferAllocation += AllocateArgBuffers("_ref") + "\n";
  argBufferAllocation += AllocateArgBuffers("_test") + "\n";
  argBufferAllocation += FillBuffersWithRandValues("") + "\n";
  argBufferAllocation += "\tprintf(\"Done with allocation\\n\");\n";
  string correctnessCheck = prototype + argBufferAllocation + CopyArgBuffersTo("_ref") + "\n";
  correctnessCheck += CorrectnessCheck(imps, referenceImpName);
  correctnessCheck += "}\n";
  return correctnessCheck;
}

string RuntimeTest::TimingCode(TimingSetting timingSetting, ImplementationMap* imps, string operationName) {
  if (timingSetting == ONEPHASETIMING) {
    return OnePhaseTimingCode(imps, operationName);
  } else {
    cout << "TWOPHASETIMING not supported yet" << endl;
    LOG_FAIL("replacement for throw call");
    throw;
  }
}

string RuntimeTest::OnePhaseTimingCode(ImplementationMap* imps, string operationName) {
  string prototype = "void time_implementations() {\n";
  string argBufferAllocation = "\tprintf(\"Starting buffer allocation\\n\");\n";
  argBufferAllocation += AllocateArgBuffers("") + "\n";
  argBufferAllocation += FillBuffersWithRandValues("") + "\n";
  argBufferAllocation += "\tFILE *" + m_dataFileName + " = fopen(\"" + m_dataFileName + "\", \"w\");\n";
  argBufferAllocation += "\tprintf(\"Done with allocation\\n\");\n";

  string timingSetup = "\tint i, j, k;\n\tlong long start_time, end_time, exec_time, total_cycles;\n";
  string timingFunc = prototype + argBufferAllocation + "\n" + timingSetup;
  string timingLoop = TimingLoops(imps);
  timingLoop += "\n\tfclose(" + m_dataFileName + ");\n";
  timingFunc += "\n" + timingLoop + "\n}\n";
  return timingFunc;
}

string RuntimeTest::MakeTestCode(SanityCheckSetting sanityCheckSetting, TimingSetting timingSetting, ImplementationMap* imps, string referenceImp) {
  string hds = HeadersAndDefines(imps);
  string impFuncs = ImplementationFunctions(imps, referenceImp);
  string sanityCheckFunc = SanityChecks(sanityCheckSetting, imps, m_operationName + "_test");
  string timingFunc = TimingCode(timingSetting, imps, m_operationName + "_test");
  string mainFunc = MainFunction();
  string testCode = hds + "\n" + impFuncs + "\n" + sanityCheckFunc + "\n" + timingFunc + "\n" + mainFunc;
  return testCode;
}

string RuntimeTest::CorrectnessCheck(ImplementationMap* imps, string referenceImpName)
{
  string correctnessCheck = "";
  vector<string> argBuffers = ArgBuffers("_ref");
  vector<string> testBuffers = ArgBuffers("_test");
  string callRefImp = "\t" + referenceImpName + "(" + CArgList(argBuffers) + ");\n\n";
  correctnessCheck += callRefImp;
  unsigned int i;
  for (i = 1; i <= imps->size(); i++) {
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

vector<string> RuntimeTest::ArgBuffers(string postfix) {
  vector<string> argBufs;
  vector<string>::iterator argIter;
  for (argIter = m_argNames.begin(); argIter != m_argNames.end(); ++argIter) {
    string name = *argIter;
    argBufs.push_back(name + postfix);
  }
  return argBufs;
}

string RuntimeTest::TimingLoop(int i) {
  string loopBody = "";
  string opName = m_operationName + "_" + std::to_string((long long int) i);
  loopBody += "\ttotal_cycles = 0;\n";
  loopBody += "\twhile (total_cycles < MIN_CYCLES) {\n";
  loopBody += "\t\tstart_time = rdtsc();\n";
  loopBody += "\t\t" + opName + "(" + CArgList(m_argNames) + ");\n";
  loopBody += "\t\tend_time = rdtsc();\n";
  loopBody += "\t\texec_time = end_time - start_time;\n";
  loopBody += "\t\ttotal_cycles += exec_time;\n";
  loopBody += "\t\tchar exec_time_str[100];\n";
  loopBody += "\t\tsprintf(exec_time_str, \"%lld\\n\", exec_time);\n";
  loopBody += "\t\tsize_t trash = fprintf(" + m_dataFileName + ", \"%s\", exec_time_str);\n";
  loopBody += "\t}\n";
  loopBody += "\tfprintf(" + m_dataFileName + ", \"#\\n\");\n";
  loopBody += "\tprintf(\"Done evaluating " + opName + "\\n\");\n";
  return loopBody;
}

string RuntimeTest::TimingLoops(ImplementationMap* imps) {
  unsigned int i;
  string loopBody = "";
  for (i = 1; i <= imps->size(); i++) {
    loopBody += TimingLoop(i);
  }
  return loopBody;
}

string RuntimeTest::FillBuffersWithRandValues(string postfix) {
  std::vector<string> bufferFills;
  std::vector<string>::iterator argIter;
  for (argIter = m_argNames.begin(); argIter != m_argNames.end(); ++argIter) {
    string bufName = *argIter;
    bufferFills.push_back("\tFILL_WITH_RAND_VALUES(BUF_SIZE, " + bufName + ");");
  }
  return ToCStatements(bufferFills);
}

string RuntimeTest::CopyArgBuffersTo(string postfix) {
  std::vector<string> bufferFills;
  for (auto srcBuffer : m_argNames) {
    string destBuffer = srcBuffer + postfix;
    bufferFills.push_back("\tCOPY_BUFFER(BUF_SIZE, " + srcBuffer + ", " + destBuffer + ");");
  }
  return ToCStatements(bufferFills);
}

string RuntimeTest::AllocateArgBuffers(string postfix) {
  std::vector<string> argAllocs;
  std::vector<string>::iterator argIter;
  for (argIter = m_argDeclarations.begin(); argIter != m_argDeclarations.end(); ++argIter) {
    string bufName = *argIter + postfix;
    argAllocs.push_back("\t" + bufName + " = ALLOC_BUFFER(BUF_SIZE * sizeof(NUM_SIZE));");
  }
  return ToCStatements(argAllocs);
}

string RuntimeTest::ToCStatements(vector<string> lines) {
  string cStatements = "";
  for (auto line : lines) {
    cStatements = cStatements + line + "\n";
  }
  return cStatements;
}

string RuntimeTest::CArgList(vector<string> args) {
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

string RuntimeTest::MakeImpFuncs(ImplementationMap* imps) {
  string endOfFuncDec = "(" + CArgList(m_argDeclarations) + ")";
  string allImplementationFuncs = "";
  ImplementationMap::iterator impIt;
  for (impIt = imps->begin(); impIt != imps->end(); ++impIt) {
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

RuntimeEvaluator::RuntimeEvaluator(string evalDirName) {
  m_evalDirName = evalDirName;
}

void RuntimeEvaluator::WriteTestCodeToFile(string executableName, string testCode) {
  string testFileName = executableName + ".c";
  std::ofstream outStream(testFileName);
  if (outStream.is_open()) {
    outStream << testCode << endl;
    outStream.close();
  } else {
    cout << "ERROR: RuntimeEvaluator could not create file " << testFileName << endl;
    LOG_FAIL("replacement for throw call");
    throw;
  }
}

void RuntimeEvaluator::CompileTest(string executableName) {
  string testFileName = executableName + ".c";
  cout << "All implementations written to files\n";
  cout << "Compile string is " << arch->CompileString(executableName, testFileName) << endl;
  int compileRes = system(arch->CompileString(executableName, testFileName).c_str());
  cout << "Compile result = " << std::to_string((long long int) compileRes) << endl;
}

void RuntimeEvaluator::RunTest(string executableName) {
  string runStr = "./" + executableName;
  int runRes = system(runStr.c_str());
  cout << "Run string is " << runStr << endl;
  cout << "Run result = " << std::to_string((long long int) runRes) << endl;
}

void RuntimeEvaluator::CleanUpTest(string executableName) {
  string removeExecutable = "rm -f " + executableName;
  system(removeExecutable.c_str());
}

vector<TimingResult*>* RuntimeEvaluator::ReadTimingData(TimingSetting timingSetting, string dataFileName, int numImpls) {
  cout << "Size of imps = " << std::to_string((long long int) numImpls) << endl;
  cout << "Calling ReadTimeDataFromFile" << endl;

  auto timingResults = ReadTimeDataFromFile(timingSetting, dataFileName, numImpls);
  string removeDataFile = "rm -f " + dataFileName;
  system(removeDataFile.c_str());
  return timingResults;
}

vector<TimingResult*>* RuntimeEvaluator::EvaluateImplementations(SanityCheckSetting sanityCheckSetting, TimingSetting timingSetting, RuntimeTest test, ImplementationMap* imps, string referenceImp) {
  string executableName = m_evalDirName + "/" + test.m_operationName;
  string testCode = test.MakeTestCode(sanityCheckSetting, timingSetting, imps, referenceImp);

  WriteTestCodeToFile(executableName, testCode);
  CompileTest(executableName);
  RunTest(executableName);
  CleanUpTest(executableName);

  auto timeResults = ReadTimingData(timingSetting, test.m_dataFileName, imps->size());
  return timeResults;
}

bool RuntimeEvaluator::IsImplementationSeparator(string token) {
  if (token == "#") {
    return true;
  }
  return false;
}

vector<TimingResult*>* RuntimeEvaluator::ReadTimeDataFromFile(TimingSetting timingSetting, string fileName, int numImpls) {
  std::ifstream dataStream(fileName);
  std::stringstream buffer;
  buffer << dataStream.rdbuf();
  dataStream.close();
  string timeData = buffer.str();
  std::vector<string> runtimeStrings;
  Tokenize(timeData, runtimeStrings, "\n");
  cout << "Num impls " << std::to_string((long long int) numImpls) << endl;
  cout << "Number of tokens " << std::to_string((long long int) runtimeStrings.size()) << endl;
  vector<TimingResult*>* timingResults = new vector<TimingResult*>();
  int i = 1;
  TimeVec* impTimes = new TimeVec();
  for (auto token : runtimeStrings) {
    if (IsImplementationSeparator(token)) {
      timingResults->push_back(new OneStageTimingResult(i, impTimes));
      i++;
      impTimes = new TimeVec();
    } else {
      impTimes->push_back(std::stod(token));
    }
  }

  if (i != numImpls + 1) {
    cout << "ERROR: In RuntimeEvaluator::ReadDataFromFile i = " << i;
    cout << " numImpls = " << numImpls << endl;
    LOG_FAIL("replacement for throw call");
    throw;
  }
  return timingResults;
}

void RuntimeEvaluator::Tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ") {
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos)
    {
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      lastPos = str.find_first_not_of(delimiters, pos);
      pos = str.find_first_of(delimiters, lastPos);
    }
}

#endif // DOLLDLA
