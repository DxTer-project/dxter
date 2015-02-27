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

#include "localInput.h"

#if DOLLDLA

void LocalInput::Prop() {
  if (!IsValidCost(m_cost)) {
    InputNode::Prop();

    if (!IsContiguous()) {
      cout << "ERROR: Local inputs must be contiguous in memory" << endl;
      throw;
    }
    m_cost = 0;
  }
}

void LocalInput::PrintCode(IndStream& out) {
  Type dataType = m_dataTypeInfo.m_type;
  //  string size = "";
  if (dataType == REAL_SINGLE) {
    out.Indent();
    *out << "float ";
    //    size += "sizeof(float)";
  } else {
    out.Indent();
    *out << "double ";
    //    size += "sizeof(double)";
  }
  string size = m_dataTypeInfo.m_numRowsVar + " * " + m_dataTypeInfo.m_numColsVar + " + 1";
  //  *out << m_varName.m_name << " = alloc_aligned_16( " << size << " );\n" << endl;
  *out << m_varName.m_name << "_array[" << size << "]" << " = {0};" << endl;
  if (dataType == REAL_SINGLE) {
    out.Indent();
    *out << "float* ";
  } else {
    out.Indent();
    *out << "double* ";
  }
  *out << m_varName.m_name << " = (((unsigned long long)" << m_varName.m_name << "_array) % 16) != 0 ? &(" << m_varName.m_name << "_array[1]) : &(" << m_varName.m_name << "_array[0]);" << endl;
}

NodeType LocalInput::GetType() const {
  return "LocalInput " + LayerNumToStr(GetLayer());
}

void LocalInput::AddVariables(VarSet& set) const {
  string typeName;
  if (GetDataType() == REAL_SINGLE) {
    typeName = "float* ";
  } else if (GetDataType() == REAL_DOUBLE) {
    typeName = "double* ";
  } else {
    cout << "Unsupported datatype: " << GetDataType() << endl;
  }
  string inputVarDecl = typeName + m_varName.m_name + ";";

  string uint = "unsigned int ";
  string nRowsVarName = uint + m_dataTypeInfo.m_numRowsVar + ";";
  string nColsVarName = uint + m_dataTypeInfo.m_numRowsVar + ";";
  string rowStrideVarName = uint + m_dataTypeInfo.m_rowStrideVar + ";";
  string colStrideVarName = uint + m_dataTypeInfo.m_colStrideVar + ";";

  Var inputNameVar(DirectVarDeclType, inputVarDecl, GetDataType());
  Var nRowsVar(DirectVarDeclType, nRowsVarName, GetDataType());
  Var nColsVar(DirectVarDeclType, nColsVarName, GetDataType());
  Var rowStrideVar(DirectVarDeclType, rowStrideVarName, GetDataType());
  Var colStrideVar(DirectVarDeclType, colStrideVarName, GetDataType());

  set.insert(inputNameVar);
  set.insert(nRowsVar);
  set.insert(nColsVar);
  set.insert(rowStrideVar);
  set.insert(colStrideVar);
}
#endif // DOLLDLA
