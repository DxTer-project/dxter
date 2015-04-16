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

#include "runtimeTest.h"

#if DOLLDLA

#include "avx.h"

RuntimeTest::RuntimeTest(ProblemInstance* prob, LLDLAUniverse* uni, unsigned int minCycles) {
  m_type = prob->GetType();
  m_operationName = prob->GetName();
  m_argNames = uni->m_argNames;
  m_outputNames = uni->m_outputNames;
  m_argDeclarations = uni->m_declarationVectors;
  m_defines = uni->m_constantDefines;
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

string RuntimeTest::HeadersAndDefines(unsigned int numImplementations) {
  m_defines.push_back("#define NUM_ALGS " + std::to_string((long long int) numImplementations));
  string headersAndDefines = ToCStatements(m_headers) + "\n" + ToCStatements(m_defines);
  for (auto ext : *arch->SupportedExtensions()) {
    cout << "EXTENSION NAME = " + ext->SetupFuncName() << endl;
    headersAndDefines += ext->GlobalDeclarations();
  }
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
  mainFunc += "\tset_up_test();\n";
  mainFunc += "\tsanity_check_implementations();\n";
  mainFunc += "\ttime_implementations();\n";
  mainFunc += "}";
  return mainFunc;
}

string RuntimeTest::SanityChecks(SanityCheckSetting sanityCheckSetting, unsigned int numImpls, string referenceImpName) {
  if (sanityCheckSetting == CHECKALLBUFFERS) {
    return AllBufferSanityChecks(numImpls, referenceImpName);
  } else if (sanityCheckSetting == CHECKOUTPUTBUFFERS) {
    return OutputBufferSanityChecks(numImpls, referenceImpName);
  } else {
    cout << "SanityCheckSetting NONE is not yet supported" << endl;
    LOG_FAIL("replacement for throw call");
    throw;
  }
}

string RuntimeTest::SanityCheckBufferAllocation() {
  string argBufferAllocation = "void sanity_check_implementations() {\n";
  argBufferAllocation += "\tprintf(\"Starting buffer allocation\\n\");\n";
  argBufferAllocation += AllocateArgBuffers(m_argDeclarations, "") + "\n";
  argBufferAllocation += AllocateArgBuffers(m_argDeclarations, "_ref") + "\n";
  argBufferAllocation += AllocateArgBuffers(m_argDeclarations, "_test") + "\n";
  argBufferAllocation += FillBuffersWithRandValues(m_argNames, "") + "\n";
  argBufferAllocation += "\tprintf(\"Done with allocation\\n\");\n";
  return argBufferAllocation;
}

string RuntimeTest::OutputBufferSanityChecks(unsigned int numImpls, string referenceImpName) {
  auto argBufferAllocation = SanityCheckBufferAllocation();
  string correctnessCheck = argBufferAllocation + CopyArgBuffersTo("_ref") + "\n";
  correctnessCheck += OutputBufferCorrectnessCheck(numImpls, referenceImpName);
  correctnessCheck += "}\n";
  return correctnessCheck;
}

string RuntimeTest::AllBufferSanityChecks(unsigned int numImpls, string referenceImpName) {
  auto argBufferAllocation = SanityCheckBufferAllocation();
  string correctnessCheck = argBufferAllocation + CopyArgBuffersTo("_ref") + "\n";
  correctnessCheck += CorrectnessCheck(numImpls, referenceImpName);
  correctnessCheck += "}\n";
  return correctnessCheck;
}

string RuntimeTest::TimingCode(TimingSetting timingSetting, unsigned int numImpls, string operationName) {
  if (timingSetting == ONEPHASETIMING) {
    return OnePhaseTimingCode(numImpls, operationName);
  } else {
    return TwoPhaseTimingCode(numImpls, operationName);
  }
}

string RuntimeTest::OnePhaseTimingCode(unsigned int numImpls, string operationName) {
  string prototype = "void time_implementations() {\n";
  string argBufferAllocation = "\tprintf(\"Starting buffer allocation\\n\");\n";
  argBufferAllocation += AllocateArgBuffers(m_argDeclarations, "") + "\n";
  argBufferAllocation += FillBuffersWithRandValues(m_argNames, "") + "\n";
  argBufferAllocation += "\tFILE *" + m_dataFileName + " = fopen(\"" + m_dataFileName + "\", \"w\");\n";
  argBufferAllocation += "\tprintf(\"Done with allocation\\n\");\n";

  string timingSetup = "\tint i, j, k;\n\tlong long start_time, end_time, exec_time, total_cycles;\n";
  string timingFunc = prototype + argBufferAllocation + "\n" + timingSetup;
  string timingLoop = TimingLoops(numImpls);
  timingLoop += "\n\tfclose(" + m_dataFileName + ");\n";
  timingFunc += "\n" + timingLoop + "\n}\n";
  return timingFunc;
}

string RuntimeTest::TwoPhaseTimingCode(unsigned int numImpls, string operationName) {
  string prototype = "void time_implementations() {\n";
  string argBufferAllocation = "\tprintf(\"Starting buffer allocation\\n\");\n";
  argBufferAllocation += AllocateArgBuffers(m_argDeclarations, "") + "\n";
  argBufferAllocation += FillBuffersWithRandValues(m_argNames, "") + "\n";
  argBufferAllocation += "\tFILE *" + m_dataFileName + " = fopen(\"" + m_dataFileName + "\", \"w\");\n";
  argBufferAllocation += "\tprintf(\"Done with allocation\\n\");\n";

  string timingSetup = "\tint i, j, k;\n\tlong long start_time, end_time, exec_time, num_runs, total_cycles, avg_exec_time;\n";
  timingSetup += "\tchar exec_time_str[100];\n";
  timingSetup += "\tsize_t trash;\n";
  string timingFunc = prototype + argBufferAllocation + "\n" + timingSetup;
  string timingLoop = TwoPhaseTimingLoops(numImpls);
  timingLoop += "\n\tfclose(" + m_dataFileName + ");\n";
  timingFunc += "\n" + timingLoop + "\n}\n";
  return timingFunc;
}

string RuntimeTest::SetupFunction() {
  string decl = "void set_up_test() {\n";
  for (auto extension : *arch->SupportedExtensions()) {
    decl += "\t" + extension->SetupFuncName() + "();\n";
  }
  decl += "}\n";
  return decl;
}

string RuntimeTest::SetupFunctions() {
  string setupFuncs = "";
  for (auto ext : *arch->SupportedExtensions()) {
    setupFuncs += ext->SetupFunc() + "\n";
  }
  return setupFuncs;
}
string RuntimeTest::MakeTestCode(SanityCheckSetting sanityCheckSetting, TimingSetting timingSetting, ImplementationMap* imps, string referenceImp) {
  unsigned int numImpls = imps->size();
  string hds = HeadersAndDefines(numImpls);
  string extSetupFuncs = SetupFunctions();
  hds += "\n" + SetupFunctions();
  string impFuncs = ImplementationFunctions(imps, referenceImp);
  string setupFunc = SetupFunction();
  string sanityCheckFunc = SanityChecks(sanityCheckSetting, numImpls, m_operationName + "_test");
  string timingFunc = TimingCode(timingSetting, imps->size(), m_operationName + "_test");
  string mainFunc = MainFunction();
  string testCode = hds + "\n" + setupFunc + "\n" + impFuncs + "\n" + sanityCheckFunc + "\n" + timingFunc + "\n" + mainFunc;
  return testCode;
}

string RuntimeTest::OutputBufferCorrectnessCheck(unsigned int numImpls, string referenceImpName)
{
  string correctnessCheck = "";
  vector<string> argBuffers = ArgBuffers("_ref");
  vector<string> testBuffers = ArgBuffers("_test");
  string callRefImp = "\t" + referenceImpName + "(" + CArgList(argBuffers) + ");\n\n";
  correctnessCheck += callRefImp;
  unsigned int i;
  for (i = 1; i <= numImpls; i++) {
    correctnessCheck += CopyArgBuffersTo("_test") + "\n\t";
    correctnessCheck += m_operationName + "_" + std::to_string((long long int) i) + "(" + CArgList(testBuffers) + ");\n";
    correctnessCheck += CheckArgBufferDiffs(m_outputNames, "_ref", "_test",std::to_string((long long int) i)) + "\n";
  }
  return correctnessCheck;
}

string RuntimeTest::CorrectnessCheck(unsigned int numImpls, string referenceImpName)
{
  string correctnessCheck = "";
  vector<string> argBuffers = ArgBuffers("_ref");
  vector<string> testBuffers = ArgBuffers("_test");
  string callRefImp = "\t" + referenceImpName + "(" + CArgList(argBuffers) + ");\n\n";
  correctnessCheck += callRefImp;
  unsigned int i;
  for (i = 1; i <= numImpls; i++) {
    correctnessCheck += CopyArgBuffersTo("_test") + "\n\t";
    correctnessCheck += m_operationName + "_" + std::to_string((long long int) i) + "(" + CArgList(testBuffers) + ");\n";
    correctnessCheck += CheckArgBufferDiffs(m_argNames, "_ref", "_test",std::to_string((long long int) i)) + "\n";
  }
  return correctnessCheck;
}

string RuntimeTest::CheckArgBufferDiffs(const vector<string> args, string refPostfix, string testPostfix, string testName)
{
  vector<string>::iterator argIter;
  vector<string> diffChecks;
  for (auto argName : args) {
    string refBuf = argName + refPostfix;
    string testBuf = argName + testPostfix;
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

string RuntimeTest::TimingLoops(unsigned int numImpls) {
  unsigned int i;
  string loopBody = "";
  for (i = 1; i <= numImpls; i++) {
    loopBody += TimingLoop(i);
  }
  return loopBody;
}

string RuntimeTest::TwoPhaseTimingLoops(unsigned int numImpls) {
  unsigned int i;
  string loopBody = "";
  for (i = 1; i <= numImpls; i++) {
    loopBody += TwoPhaseTimingLoop(i);
  }
  return loopBody;
}

string RuntimeTest::TwoPhaseTimingLoop(int i) {
  string loopBody = "";
  string opName = m_operationName + "_" + std::to_string((long long int) i);
  loopBody += "\ttotal_cycles = 0;\n";
  loopBody += "\tnum_runs = 0;\n";
  loopBody += "\twhile (total_cycles < MIN_CYCLES) {\n";
  loopBody += "\t\tnum_runs++;\n";
  loopBody += "\t\tstart_time = rdtsc();\n";
  loopBody += "\t\t" + opName + "(" + CArgList(m_argNames) + ");\n";
  loopBody += "\t\tend_time = rdtsc();\n";
  loopBody += "\t\texec_time = end_time - start_time;\n";
  loopBody += "\t\ttotal_cycles += exec_time;\n";
  loopBody += "\t}\n";
  loopBody += "\tfor (j = 0; j < 1; j++) {\n";
  loopBody += "\t\tstart_time = rdtsc();\n";
  loopBody += "\t\tfor (i = 0; i < num_runs; i++) {\n";
  loopBody += "\t\t\t" + opName + "(" + CArgList(m_argNames) + ");\n";
  loopBody += "\t\t}\n";
  loopBody += "\t\tend_time = rdtsc();\n";
  loopBody += "\t\texec_time = end_time - start_time;\n";
  loopBody += "\t\tavg_exec_time = exec_time / num_runs;\n";
  loopBody += "\t\tmemset(exec_time_str, 0, 100);\n";
  loopBody += "\t\tsprintf(exec_time_str, \"%lld\\n\", avg_exec_time);\n";
  loopBody += "\t\ttrash = fprintf(" + m_dataFileName + ", \"%s\", exec_time_str);\n";
  loopBody += "\t}\n";
  loopBody += "\tfprintf(" + m_dataFileName + ", \"#\\n\");\n";
  loopBody += "\tprintf(\"Done evaluating " + opName + "\\n\");\n";
  return loopBody;
}

string RuntimeTest::FillBuffersWithRandValues(const vector<string> argNames, string postfix) {
  std::vector<string> bufferFills;
  std::vector<string>::iterator argIter;
  for (auto bufName : argNames) {
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

string RuntimeTest::AllocateArgBuffers(const vector<string> args, string postfix) {
  std::vector<string> argAllocs;
  for (auto argName : args) {
    string bufName = argName + postfix;
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

#endif // DOLLDLA
