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

#include "contigCopy.h"

#if DOLLDLA

void ContiguousCopy::PrintCode(IndStream& out) {
  out.Indent();
  string inputName = GetInputName(0).m_name;
  string inputSize = "(" + InputDataType(0).m_numRowsVar + " * " + InputDataType(0).m_numColsVar + ")";

  string receivingName = GetInputName(1).m_name;

  if (Input(0)->GetDataType() == REAL_SINGLE) {
    inputSize += "* sizeof(float)";
  } else {
    inputSize += "* sizeof(double)";
  }
  *out << "memcpy( (void*) " << receivingName << ", (void*) " << inputName << ", ";
  *out << inputSize << " );\n";
}

void ContiguousCopy::Prop() {
  if (!IsValidCost(m_cost)) {
    Copy::Prop();
    if (!InputIsContiguous(0) || !InputIsContiguous(1)) {
      cout << "ERROR: Contiguous copy on non-contiguous operands" << endl;
      LOG_FAIL("Replacement for call to throw;");
    }
    m_cost = ZERO;
  }
}

#endif // DOLLDLA
