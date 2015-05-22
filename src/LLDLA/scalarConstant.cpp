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

#include "scalarConstant.h"

#if DOLLDLA

ScalarConstant::ScalarConstant(string name, Type dataType, double value)
  : LocalInput::LocalInput(name, 1, 1, 1, 1, dataType) {
  m_value = value;
}

void ScalarConstant::Prop() {
  if (!IsValidCost(m_cost)) {
    m_cost = ZERO;
  }
  return;
}

void ScalarConstant::PrintCode(IndStream& out) {
  Type dataType = m_dataTypeInfo.m_type;
  string typeName;
  string varName = m_varName.m_name;
  string byteArray = varName + "_name";
  if (dataType == REAL_SINGLE) {
    typeName = "float ";
  } else {
    typeName = "double ";
  }

  out.Indent();
  *out << typeName << " " << varName << "[1];\n";
  out.Indent();
  *out << varName << "[0] = " << std::to_string((long double) m_value) << ";\n";

  return;
}

void ScalarConstant::AddVariables(VarSet& set) const {
  return;
}

NodeType ScalarConstant::GetType() const {
  return "Scalar constant " + LayerNumToStr(GetLayer());
}


#endif // DOLLDLA
