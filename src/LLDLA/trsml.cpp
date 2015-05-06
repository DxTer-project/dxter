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

#include "trsml.h"

#if DOLLDLA

void TRSML::PrintCode(IndStream& out) {
  out.Indent();
  if (GetDataType() == REAL_DOUBLE) {
    *out << "simple_trsml(";
  } else {
    throw;
  }
  *out << InputDataType(0).m_numRowsVar << ", " <<
    InputDataType(0).m_numColsVar << ", " <<
    GetInputName(0).str() << ", " <<
    InputDataType(0).m_rowStrideVar << ", " <<
    InputDataType(0).m_colStrideVar << ", " <<
    GetInputName(1).str() << ", " <<
    InputDataType(1).m_rowStrideVar << ", " <<
    InputDataType(1).m_colStrideVar << ");\n";
  return;
}

void TRSML::SanityCheckInputDimensions() {
  if (GetInputNumRows(0) != GetInputNumCols(0)) {
    LOG_FAIL("Bad dimensions in TRSML");
    throw;
  }
  if (GetInputNumRows(0) != GetInputNumRows(1)) {
    LOG_FAIL("Bad dimensions in TRSML");
    throw;
  }
}

void TRSML::Prop() {
  if (!IsValidCost(m_cost)) {
    SanityCheckInputDimensions();
    m_cost = ZERO;
  }
}

Phase TRSML::MaxPhase() const {
  switch (m_layer)
    { 
    case(ABSLAYER):
      return LLDLALOOPPHASE;
    case(LLDLAMIDLAYER):
      return LLDLAPRIMPHASE;
    case (LLDLAPRIMITIVELAYER):
      return NUMPHASES; 
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }
}

Node* TRSML::BlankInst() {
  throw;
}

void TRSML::Duplicate(const Node* orig, bool shallow, bool possMerging) {
  throw;
}

NodeType TRSML::GetType() const {
  return "TRSML";
}

#endif // DOLLDLA
