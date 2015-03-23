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
  string typeName;
  string varName = m_varName.m_name;
  string byteArray = varName + "_name";
  if (dataType == REAL_SINGLE) {
    typeName = "float";
  } else {
    typeName = "double";
  }
  string size = m_dataTypeInfo.m_numRowsVar + " * " + m_dataTypeInfo.m_numColsVar + " + 32";
  out.Indent();
  *out << typeName << " " << varName << "[" << size << "*sizeof(" << typeName << ")] __attribute__((aligned(32))) = {0};" << endl;

  /*  string arrName = m_varName.m_name + "_char_array";
  string arrPtrName = arrName + "_array_ptr";
  out.Indent();
  *out << typeName << " " << arrName << "[" << size << "*sizeof(" << typeName << ") + 32" << "] = {0};" << endl;
  out.Indent();
  *out << "unsigned char* " << arrPtrName << " = (unsigned char*)" << arrName << ";" << endl;
  string shiftStr = "32 - (((unsigned int) " + arrPtrName + ") % 32)";
  out.Indent();
  *out << arrPtrName << " = " << arrPtrName << " + " << shiftStr << ";" << endl;
  out.Indent();
  *out << m_varName.m_name << " = (" << typeName << "*)" << arrPtrName << ";" << endl;*/
  //  *out << m_varName.m_name << " = alloc_aligned_32(sizeof(" + typeName + ")*" + size + ");";
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

  string uint = "const unsigned int ";
  string nRowsVarName = uint + m_dataTypeInfo.m_numRowsVar + " = ";
  nRowsVarName += std::to_string((long long int) m_msize.OnlyEntry()) + ";";
  string nColsVarName = uint + m_dataTypeInfo.m_numColsVar + " = ";
  nColsVarName += std::to_string((long long int) m_nsize.OnlyEntry()) + ";";
  string rowStrideVarName = uint + m_dataTypeInfo.m_rowStrideVar + " = ";
  rowStrideVarName += std::to_string((long long int) m_dataTypeInfo.m_rowStrideVal) + ";";
  string colStrideVarName = uint + m_dataTypeInfo.m_colStrideVar + " = ";
  colStrideVarName += std::to_string((long long int) m_dataTypeInfo.m_colStrideVal) + ";";

  Var inputNameVar(DirectVarDeclType, inputVarDecl, GetDataType());
  Var nRowsVar(DirectVarDeclType, nRowsVarName, GetDataType());
  Var nColsVar(DirectVarDeclType, nColsVarName, GetDataType());
  Var rowStrideVar(DirectVarDeclType, rowStrideVarName, GetDataType());
  Var colStrideVar(DirectVarDeclType, colStrideVarName, GetDataType());

  //  set.insert(inputNameVar);
  set.insert(nRowsVar);
  set.insert(nColsVar);
  set.insert(rowStrideVar);
  set.insert(colStrideVar);
}

#endif // DOLLDLA
