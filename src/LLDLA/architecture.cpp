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

#include "architecture.h"
#include "base.h"

#if DOLLDLA

#include "logging.h"

int Architecture::VecRegWidth(Type type)
{
  //  cout << "Getting reg width\n";
  if (type == REAL_SINGLE) {
    //    cout << "Type is single\n";
    return SVecRegWidth();
  } else if (type == REAL_DOUBLE) {
    //    cout << "Type is double\n";
    return DVecRegWidth();
  } else {
    cout << "Error: VecRegWidth bad type\n";
    LOG_FAIL("Replacement for call to throw;");
    throw;
  }
}

string Architecture::VecRegTypeDec(Type type)
{
  if (type == REAL_SINGLE) {
    return SVecRegTypeDec();
  } else if (type == REAL_DOUBLE) {
    return DVecRegTypeDec();
  } else {
    cout << "Error: VecRegWidth bad type\n";
    LOG_FAIL("Replacement for call to throw;");
    throw;
  }
}

string Architecture::TypeName(Type type)
{
  if (type == REAL_SINGLE) {
    return STypeName();
  } else if (type == REAL_DOUBLE) {
    return DTypeName();
  } else {
    cout << "Error: VecRegWidth bad type\n";
    LOG_FAIL("Replacement for call to throw;");
    throw;
  }
}

string Architecture::AddCode(Type type, string operand1, string operand2, string result)
{
  if (type == REAL_SINGLE) {
    return SAddCode(operand1, operand2, result);
  } else if (type == REAL_DOUBLE) {
    return DAddCode(operand1, operand2, result);
  } else {
    cout << "Error: VecRegWidth bad type\n";
    LOG_FAIL("Replacement for call to throw;");
    throw;
  }
}

string Architecture::MulCode(Type type, string operand1, string operand2, string result)
{
  if (type == REAL_SINGLE) {
    return SMulCode(operand1, operand2, result);
  } else if (type == REAL_DOUBLE) {
    return DMulCode(operand1, operand2, result);
  } else {
    cout << "Error: VecRegWidth bad type\n";
    LOG_FAIL("Replacement for call to throw;");
    throw;
  }
}

string Architecture::FMACode(Type type, string operand1, string operand2, string operand3, string result)
{
  if (type == REAL_SINGLE) {
    return SFMACode(operand1, operand2, operand3, result);
  } else if (type == REAL_DOUBLE) {
    return DFMACode(operand1, operand2, operand3, result);
  } else {
    cout << "Error: VecRegWidth bad type\n";
    LOG_FAIL("Replacement for call to throw;");
    throw;
  }
}

string Architecture::StridedLoad(Type type, string memPtr, string receivingLoc, string stride)
{
  if (type == REAL_SINGLE) {
    return SStridedLoad(memPtr, receivingLoc, stride);
  } else if (type == REAL_DOUBLE) {
    return DStridedLoad(memPtr, receivingLoc, stride);
  } else {
    cout << "Error: VecRegWidth bad type\n";
    LOG_FAIL("Replacement for call to throw;");
    throw;
  }
}

string Architecture::PackedLoad(Type type, string memPtr, string receivingLoc, string stride, int residual)
{
  if (type == REAL_SINGLE) {
    return SPackedLoad(memPtr, receivingLoc, stride, residual);
  } else if (type == REAL_DOUBLE) {
    return DPackedLoad(memPtr, receivingLoc, stride, residual);
  } else {
    cout << "Error: VecRegWidth bad type\n";
    LOG_FAIL("Replacement for call to throw;");
    throw;
  }
}

string Architecture::ContiguousLoad(Type type, string memPtr, string receivingLoc)
{
  if (type == REAL_SINGLE) {
    return SContiguousLoad(memPtr, receivingLoc);
  } else if (type == REAL_DOUBLE) {
    return DContiguousLoad(memPtr, receivingLoc);
  } else {
    cout << "Error: VecRegWidth bad type\n";
    LOG_FAIL("Replacement for call to throw;");
    throw;
  }
}

string Architecture::DuplicateLoad(Type type, string memPtr, string receivingLoc)
{
  if (type == REAL_SINGLE) {
    return SDuplicateLoad(memPtr, receivingLoc);
  } else if (type == REAL_DOUBLE) {
    return DDuplicateLoad(memPtr, receivingLoc);
  } else {
    cout << "Error: VecRegWidth bad type\n";
    LOG_FAIL("Replacement for call to throw;");
    throw;
  }
}

string Architecture::AccumCode(Type type, string memPtr, string startingLoc)
{
  if (type == REAL_SINGLE) {
    return SAccumCode(memPtr, startingLoc);
  } else if (type == REAL_DOUBLE) {
    return DAccumCode(memPtr, startingLoc);
  } else {
    cout << "Error: VecRegWidth bad type\n";
    LOG_FAIL("Replacement for call to throw;");
    throw;
  }
}

string Architecture::ContiguousStore(Type type, string memPtr, string startingLoc)
{
  if (type == REAL_SINGLE) {
    return SContiguousStore(memPtr, startingLoc);
  } else if (type == REAL_DOUBLE) {
    return DContiguousStore(memPtr, startingLoc);
  } else {
    cout << "Error: VecRegWidth bad type\n";
    LOG_FAIL("Replacement for call to throw;");
    throw;
  }
}

string Architecture::StridedStore(Type type, string memPtr, string startingLoc, string stride)
{
  if (type == REAL_SINGLE) {
    return SStridedStore(memPtr, startingLoc, stride);
  } else if (type == REAL_DOUBLE) {
    return DStridedStore(memPtr, startingLoc, stride);
  } else {
    cout << "Error: VecRegWidth bad type\n";
    LOG_FAIL("Replacement for call to throw;");
    throw;
  }
}

string Architecture::UnpackStore(Type type, string memPtr, string startingLoc, string stride, int residual) {
  if (type == REAL_SINGLE) {
    return SUnpackStore(memPtr, startingLoc, stride, residual);
  } else if (type == REAL_DOUBLE) {
    return DUnpackStore(memPtr, startingLoc, stride, residual);
  } else {
    cout << "Error: VecRegWidth bad type\n";
    LOG_FAIL("Replacement for call to throw;");
    throw;
  }
}

string Architecture::ZeroVar(Type type, string varName)
{
  if (type == REAL_SINGLE) {
    return SZeroVar(varName);
  } else if (type == REAL_DOUBLE) {
    return DZeroVar(varName);
  } else {
    cout << "Error: VecRegWidth bad type\n";
    LOG_FAIL("Replacement for call to throw;");
    throw;
  }
}

double Architecture::FlopsPerCycle(Type type)
{
  if (type == REAL_SINGLE) {
    return SFlopsPerCycle();
  } else if (type == REAL_DOUBLE) {
    return DFlopsPerCycle();
  } else {
    cout << "Error: VecRegWidth bad type\n";
    LOG_FAIL("Replacement for call to throw;");
    throw;
  }
}

string Architecture::SPackedLoad(string memPtr, string receivingLoc, string stride, int residual) {
  string loadCode = SZeroVar(receivingLoc);
  for (int i = 0; i < residual; i++) {
    string indStr = std::to_string((long long int) i);
    loadCode += receivingLoc + ".f[" + indStr + "] = *(" + memPtr + " + " + indStr + " * " + stride + "); ";
  }
  return loadCode;
}

string Architecture::DPackedLoad(string memPtr, string receivingLoc, string stride, int residual) {
  string loadCode = DZeroVar(receivingLoc);
  for (int i = 0; i < residual; i++) {
    string indStr = std::to_string((long long int) i);
    loadCode += receivingLoc + ".d[" + indStr + "] = *(" + memPtr + " + " + indStr + " * " + stride + "); ";
  }
  return loadCode;
}

string Architecture::SUnpackStore(string memPtr, string startingLoc, string stride, int residual) {
  string storeCode = "";
  for (int i = 0; i < residual; i++) {
    auto resInd = std::to_string((long long int) i);
    storeCode += "*(" + memPtr + " + " + resInd + " * " + stride + ") = " + startingLoc + ".f[" + resInd + "]; ";
  }
  return storeCode;
}

string Architecture::DUnpackStore(string memPtr, string startingLoc, string stride, int residual) {
  string storeCode = "";
  for (int i = 0; i < residual; i++) {
    auto resInd = std::to_string((long long int) i);
    storeCode += "*(" + memPtr + " + " + resInd + " * " + stride + ") = " + startingLoc + ".d[" + resInd + "]; ";
  }
  return storeCode;
}

string AMDEngSample::CompileString(string executableName, string testFileName)
{
  string compileStr = "gcc -O3 -mavx -march=native -mfma -finline-functions -funroll-loops -o ";
  compileStr += executableName + " " + testFileName + " runtimeEvaluation/utils.c";
  return compileStr;
}

double AMDEngSample::SFlopsPerCycle()
{
  return 16.0;
}

double AMDEngSample::DFlopsPerCycle()
{
  return 8.0;
}

double AMDEngSample::CyclesPerSecond()
{
  return 3.7e9;
}

int AMDEngSample::SVecRegWidth()
{
  return 4;
}

string AMDEngSample::SVecRegTypeDec()
{
  return "typedef union {\n\t__m128 v;\n\tfloat f[4];\n} svec_reg;\n";
}

string AMDEngSample::STypeName()
{
  return "svec_reg";
}

string AMDEngSample::SAddCode(string operand1, string operand2, string result)
{
  return result + ".v = _mm_add_ps( " + operand1 + ".v , " + operand2 + ".v );\n";
}

string AMDEngSample::SMulCode(string operand1, string operand2, string result)
{
  return result + ".v = _mm_mul_ps( " + operand1 + ".v , " + operand2 + ".v );\n";
}

string AMDEngSample::SFMACode(string operand1, string operand2, string operand3, string result)
{
  return result + ".v = _mm_fmadd_ps( " + operand1 + ".v , " + operand2 + ".v , " + operand3 + ".v );\n";
}

string AMDEngSample::SAccumCode(string memPtr, string startingLoc)
{
  return "*" + memPtr + " += "
    + startingLoc + ".f[0] + "
    + startingLoc + ".f[1] + "
    + startingLoc + ".f[2] + "
    + startingLoc + ".f[3];\n";
}

string AMDEngSample::SContiguousLoad(string memPtr, string receivingLoc)
{
  return receivingLoc + ".v = _mm_load_ps( " + memPtr + " );\n";
}

string AMDEngSample::SDuplicateLoad(string memPtr, string receivingLoc)
{
  return receivingLoc + ".v = _mm_load_ps1( " + memPtr + " );\n";
}

string AMDEngSample::SStridedLoad(string memPtr, string receivingLoc, string stride)
{
  return receivingLoc + ".f[0] = *(" + memPtr + "); "
    + receivingLoc + ".f[1] = *(" + memPtr + " + " + stride + " ); "
    + receivingLoc + ".f[2] = *(" + memPtr + " + 2 * " + stride + " ); "
    + receivingLoc + ".f[3] = *(" + memPtr + " + 3 * " + stride + " );\n";
}

string AMDEngSample::SContiguousStore(string memPtr, string startingLoc)
{
  return "_mm_store_ps( " + memPtr + ", " + startingLoc + ".v );\n";
}

string AMDEngSample::SStridedStore(string memPtr, string startingLoc, string stride)
{
  return "*" + memPtr + " = " + startingLoc + ".f[0]; "
    + "*(" + memPtr + " + " + stride + ") = " + startingLoc + ".f[1]; "
    + "*(" + memPtr + " + 2 * " + stride + ") = " + startingLoc + ".f[2]; "
    + "*(" + memPtr + " + 3 * " + stride + ") = " + startingLoc + ".f[3];\n";
}

string AMDEngSample::SZeroVar(string varName)
{
  return varName + ".v = _mm_setzero_ps();\n";
}

int AMDEngSample::DVecRegWidth()
{
  return 2;
}

string AMDEngSample::DVecRegTypeDec()
{
  return "typedef union {\n\t__m128d v;\n\tdouble d[2];\n} dvec_reg;\n";
}

string AMDEngSample::DTypeName()
{
  return "dvec_reg";
}

string AMDEngSample::DAddCode(string operand1, string operand2, string result)
{
  return result + ".v = _mm_add_pd( " + operand1 + ".v , " + operand2 + ".v );\n";
}

string AMDEngSample::DMulCode(string operand1, string operand2, string result)
{
  return result + ".v = _mm_mul_pd( " + operand1 + ".v , " + operand2 + ".v );\n";
}

string AMDEngSample::DFMACode(string operand1, string operand2, string operand3, string result)
{
  return result + ".v = _mm_fmadd_pd( " + operand1 + ".v , " + operand2 + ".v , " + operand3 + ".v );\n";
}

string AMDEngSample::DAccumCode(string memPtr, string startingLoc)
{
  return "*" + memPtr + " += "
    + startingLoc + ".d[0] + "
    + startingLoc + ".d[1];\n";
}

string AMDEngSample::DContiguousLoad(string memPtr, string receivingLoc)
{
  return receivingLoc + ".v = _mm_load_pd( " + memPtr + " );\n";
}

string AMDEngSample::DStridedLoad(string memPtr, string receivingLoc, string stride)
{
  return receivingLoc + ".d[0] = *(" + memPtr + "); "
    + receivingLoc + ".d[1] = *(" + memPtr + " + " + stride + " );\n";
}

string AMDEngSample::DDuplicateLoad(string memPtr, string receivingLoc)
{
  return receivingLoc + ".v = _mm_loaddup_pd( " + memPtr + " );\n";
}

string AMDEngSample::DContiguousStore(string memPtr, string startingLoc)
{
  return "_mm_store_pd( " + memPtr + ", " + startingLoc + ".v );\n";
}

string AMDEngSample::DStridedStore(string memPtr, string startingLoc, string stride)
{
  return "*" + memPtr + " = " + startingLoc + ".d[0]; "
    + "*(" + memPtr + " + " + stride + ") = " + startingLoc + ".d[1];\n";
}

string AMDEngSample::DZeroVar(string varName)
{
  return varName + ".v = _mm_setzero_pd();\n";
}

string Stampede::CompileString(string executableName, string testFileName)
{
  string compileStr = "icc -O3 -xhost -ip -ipo -fargument-noalias-global -o ";
  compileStr += executableName + " " + testFileName + " runtimeEvaluation/utils.c";
  return compileStr;
}

double Stampede::CyclesPerSecond()
{
  return 2.7e9;
}

double Stampede::DFlopsPerCycle()
{
  return 8.0;
}

double Stampede::SFlopsPerCycle()
{
  return 16.0;
}

int Stampede::SVecRegWidth()
{
  return 8;
}

string Stampede::SVecRegTypeDec()
{
  return "typedef union {\n\t__m256 v;\n\tfloat f[8];\n} svec_reg;\n";
}

string Stampede::STypeName()
{
  return "svec_reg";
}

string Stampede::SAddCode(string operand1, string operand2, string result)
{
  return result + ".v = _mm256_add_ps( " + operand1 + ".v , " + operand2 + ".v );\n";
}

string Stampede::SMulCode(string operand1, string operand2, string result)
{
  return result + ".v = _mm256_mul_ps( " + operand1 + ".v , " + operand2 + ".v );\n";
}

string Stampede::SFMACode(string operand1, string operand2, string operand3, string result)
{
  return result + ".v = _mm256_add_ps( _mm256_mul_ps( " + operand1 + ".v, " + operand2 + ".v ), " + operand3 + ".v );\n";
}

string Stampede::SAccumCode(string memPtr, string startingLoc)
{
  return "*" + memPtr + " += "
    + startingLoc + ".f[0] + "
    + startingLoc + ".f[1] + "
    + startingLoc + ".f[2] + "
    + startingLoc + ".f[3] + "
    + startingLoc + ".f[4] + "
    + startingLoc + ".f[5] + "
    + startingLoc + ".f[6] + "
    + startingLoc + ".f[7];\n";

}

string Stampede::SContiguousLoad(string memPtr, string receivingLoc)
{
  return receivingLoc + ".v = _mm256_loadu_ps( " + memPtr + " );\n";
}

string Stampede::SDuplicateLoad(string memPtr, string receivingLoc)
{
  return receivingLoc + ".v = _mm256_broadcast_ss( " + memPtr + " );\n";
}

string Stampede::SStridedLoad(string memPtr, string receivingLoc, string stride)
{
  return receivingLoc + ".f[0] = *(" + memPtr + "); "
    + receivingLoc + ".f[1] = *(" + memPtr + " + " + stride + " ); "
    + receivingLoc + ".f[2] = *(" + memPtr + " + 2 * " + stride + " ); "
    + receivingLoc + ".f[3] = *(" + memPtr + " + 3 * " + stride + " ); "
    + receivingLoc + ".f[4] = *(" + memPtr + " + 4 * " + stride + " ); "
    + receivingLoc + ".f[5] = *(" + memPtr + " + 5 * " + stride + " ); "
    + receivingLoc + ".f[6] = *(" + memPtr + " + 6 * " + stride + " ); "
    + receivingLoc + ".f[7] = *(" + memPtr + " + 7 * " + stride + " );\n";
}

string Stampede::SContiguousStore(string memPtr, string startingLoc)
{
  return "_mm256_storeu_ps( " + memPtr + ", " + startingLoc + ".v );\n";
}

string Stampede::SStridedStore(string memPtr, string startingLoc, string stride)
{
  return "*" + memPtr + " = " + startingLoc + ".f[0]; "
    + "*(" + memPtr + " + " + stride + ") = " + startingLoc + ".f[1]; "
    + "*(" + memPtr + " + 2 * " + stride + ") = " + startingLoc + ".f[2]; "
    + "*(" + memPtr + " + 3 * " + stride + ") = " + startingLoc + ".f[3]; "
    + "*(" + memPtr + " + 4 * " + stride + ") = " + startingLoc + ".f[4]; "
    + "*(" + memPtr + " + 5 * " + stride + ") = " + startingLoc + ".f[5]; "
    + "*(" + memPtr + " + 6 * " + stride + ") = " + startingLoc + ".f[6]; "
    + "*(" + memPtr + " + 7 * " + stride + ") = " + startingLoc + ".f[7];\n";
}

string Stampede::SZeroVar(string varName)
{
  return varName + ".v = _mm256_setzero_ps();\n";
}

int Stampede::DVecRegWidth()
{
  return 4;
}

string Stampede::DVecRegTypeDec()
{
  return "typedef union {\n\t__m256d v;\n\tdouble d[4];\n} dvec_reg;\n";
}

string Stampede::DTypeName()
{
  return "dvec_reg";
}

string Stampede::DAddCode(string operand1, string operand2, string result)
{
  return result + ".v = _mm256_add_pd( " + operand1 + ".v , " + operand2 + ".v );\n";
}

string Stampede::DMulCode(string operand1, string operand2, string result)
{
  return result + ".v = _mm256_mul_pd( " + operand1 + ".v , " + operand2 + ".v );\n";
}

string Stampede::DFMACode(string operand1, string operand2, string operand3, string result)
{
  return result + ".v = _mm256_add_pd( _mm256_mul_pd( " + operand1 + ".v, " + operand2 + ".v ), " + operand3 + ".v );\n";
}

string Stampede::DAccumCode(string memPtr, string startingLoc)
{
  return "*" + memPtr + " += "
    + startingLoc + ".d[0] + "
    + startingLoc + ".d[1] + "
    + startingLoc + ".d[2] + "
    + startingLoc + ".d[3];\n";

}

string Stampede::DContiguousLoad(string memPtr, string receivingLoc)
{
  return receivingLoc + ".v = _mm256_loadu_pd( " + memPtr + " );\n";
}

string Stampede::DStridedLoad(string memPtr, string receivingLoc, string stride)
{
  return receivingLoc + ".d[0] = *(" + memPtr + "); "
    + receivingLoc + ".d[1] = *(" + memPtr + " + " + stride + " ); "
    + receivingLoc + ".d[2] = *(" + memPtr + " + 2 * " + stride + " ); "
    + receivingLoc + ".d[3] = *(" + memPtr + " + 3 * " + stride + " );\n";
}

string Stampede::DDuplicateLoad(string memPtr, string receivingLoc)
{
  return receivingLoc + ".v = _mm256_broadcast_sd( " + memPtr + " );\n";
}

string Stampede::DContiguousStore(string memPtr, string startingLoc)
{
  return "_mm256_storeu_pd( " + memPtr + ", " + startingLoc + ".v );\n";
}

string Stampede::DStridedStore(string memPtr, string startingLoc, string stride)
{
  return "*" + memPtr + " = " + startingLoc + ".d[0]; "
    + "*(" + memPtr + " + " + stride + ") = " + startingLoc + ".d[1]; "
    + "*(" + memPtr + " + 2 * " + stride + ") = " + startingLoc + ".d[2]; "
    + "*(" + memPtr + " + 3 * " + stride + ") = " + startingLoc + ".d[3];\n";
}

string Stampede::DZeroVar(string varName)
{
  return varName + ".v = _mm256_setzero_pd();\n";
}

string HaswellMacbook::CompileString(string executableName, string testFileName)
{
  string compileStr = "clang -O3 -mavx -march=native -mfma -funroll-loops -o ";
  compileStr += executableName + " " + testFileName + " runtimeEvaluation/utils.c";
  return compileStr;
}

double HaswellMacbook::CyclesPerSecond()
{
  return 1.4e9;
}

double HaswellMacbook::DFlopsPerCycle()
{
  return 16.0;
}

double HaswellMacbook::SFlopsPerCycle()
{
  return 32.0;
}

string HaswellMacbook::SFMACode(string operand1, string operand2, string operand3, string result)
{
  return result + ".v = _mm256_fmadd_ps( " + operand1 + ".v, " + operand2 + ".v, " + operand3 + ".v );\n";
}
string HaswellMacbook::DFMACode(string operand1, string operand2, string operand3, string result)
{
  return result + ".v = _mm256_fmadd_pd( " + operand1 + ".v, " + operand2 + ".v, " + operand3 + ".v );\n";
}

#endif // DOLLDLA
