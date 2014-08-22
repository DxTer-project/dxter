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

#include "architecture.h"
#include "base.h"

#if DOLLDLA

int Architecture::VecRegWidth(Type type)
{
  if (type == REAL_SINGLE) {
    return SVecRegWidth();
  } else if (type == REAL_DOUBLE) {
    return DVecRegWidth();
  } else {
    cout << "Error: VecRegWidth bad type\n";
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
    throw;
  }
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
  return receivingLoc + ".f[0] = *(" + memPtr + ");\n"
    + receivingLoc + ".f[1] = *(" + memPtr + " + " + stride + " );\n"
    + receivingLoc + ".f[2] = *(" + memPtr + " + 2 * " + stride + " );\n"
    + receivingLoc + ".f[3] = *(" + memPtr + " + 3 * " + stride + " );\n";
}

string AMDEngSample::SContiguousStore(string memPtr, string startingLoc)
{
  return "_mm_store_ps( " + memPtr + ", " + startingLoc + ".v );\n";
}

string AMDEngSample::SStridedStore(string memPtr, string startingLoc, string stride)
{
  return "*" + memPtr + " = " + startingLoc + ".f[0];\n"
    + "*(" + memPtr + " + " + stride + ") = " + startingLoc + ".f[1];\n"
    + "*(" + memPtr + " + 2 * " + stride + ") = " + startingLoc + ".f[2];\n"
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
  return receivingLoc + ".d[0] = *(" + memPtr + ");\n"
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
  return "*" + memPtr + " = " + startingLoc + ".d[0];\n"
    + "*(" + memPtr + " + " + stride + ") = " + startingLoc + ".d[1];\n";
}

string AMDEngSample::DZeroVar(string varName)
{
  return varName + ".v = _mm_setzero_pd();\n";
}

#endif // DOLLDLA
