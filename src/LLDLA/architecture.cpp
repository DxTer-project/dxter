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

int AMDEngSample::SVecRegWidth()
{
  return 4;
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
  return "*" + memPtr + " = "
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
  return receivingLoc + ".f[0] = *(" + receivingLoc + ");\n"
    + receivingLoc + ".f[1] = *(" + receivingLoc + " + stride);\n"
    + receivingLoc + ".f[2] = *(" + receivingLoc + " + 2 * stride);\n"
    + receivingLoc + ".f[3] = *(" + receivingLoc + " + 3 * stride);\n";
}

string AMDEngSample::SContiguousStore(string memPtr, string startingLoc)
{
  return "_mm_store_ps( " + memPtr + ", " + startingLoc + ".v );\n";
}

string AMDEngSample::SStridedStore(string memPtr, string startingLoc, string stride)
{
  return "*" + memPtr + " = " + startingLoc + ".f[0];\n"
    + "*(" + memPtr + " + stride) = " + startingLoc + ".f[1];\n"
    + "*(" + memPtr + " + 2 * stride) = " + startingLoc + ".f[2];\n"
    + "*(" + memPtr + " + 3 * stride) = " + startingLoc + ".f[3];\n";
}

int AMDEngSample::DVecRegWidth()
{
  return 2;
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
  return "*" + memPtr + " = "
    + startingLoc + ".d[0] + "
    + startingLoc + ".d[1];\n";
}

string AMDEngSample::DContiguousLoad(string memPtr, string receivingLoc)
{
  return receivingLoc + ".v = _mm_load_pd( " + memPtr + " );\n";
}

string AMDEngSample::DStridedLoad(string memPtr, string receivingLoc, string stride)
{
  return receivingLoc + ".d[0] = *(" + receivingLoc + ")"
    + receivingLoc + ".d[1] = *(" + receivingLoc + " + stride)";
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
    + "*(" + memPtr + " + stride) = " + startingLoc + ".d[1];\n";
}


int AMDEngSample::VecRegWidth(Type type)
{
  if (type == REAL_SINGLE) {
    return SVecRegWidth();
  } else if (type == REAL_DOUBLE) {
    return DVecRegWidth();
  } else {
    cout << "Error: AMDEngSample VecRegWidth bad type\n";
    throw;
  }
}

#endif // DOLLDLA
