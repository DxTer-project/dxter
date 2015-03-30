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

#include "pack.h"

#if DOLLDLA

Pack::Pack(Layer layer) {
  m_layer = layer;
}

const DataTypeInfo& Pack::DataType(ConnNum num) const {
  if (num != 0) {
    LOG_FAIL("replacement for throw call");
  }
  return InputDataType(1);
}

bool Pack::Overwrites(const Node* input, ConnNum num) const {
  const NodeConn* conn = m_inputs[1];
  return conn->m_n == input && conn->m_num == num;
}

void Pack::PrintCode(IndStream& out) {
  string toPack = GetInputName(0).m_name;
  string packInto = GetInputName(1).m_name;

  string toPackM = InputDataType(0).m_numRowsVar;
  string toPackN = InputDataType(0).m_numColsVar;

  string toPackRowStride = InputDataType(0).m_rowStrideVar;
  string toPackColStride = InputDataType(0).m_colStrideVar;

  string packIntoM = InputDataType(1).m_numRowsVar;
  string packIntoN = InputDataType(1).m_numColsVar;

  string packIntoRowStride = InputDataType(1).m_rowStrideVar;
  string packIntoColStride = InputDataType(1).m_colStrideVar;

  out.Indent();
  if (Input(0)->GetDataType() == REAL_SINGLE) {
    *out << "pack_float( ";
  } else {
    *out << "pack_double( ";
  }

  *out << toPack << ", " << packInto << ", " <<
    toPackM << ", " << toPackN << ", " <<
    toPackRowStride << ", " << toPackColStride;
  *out << ", " <<
    packIntoM << ", " << packIntoN << ", " <<
    packIntoRowStride << ", " << packIntoColStride << " );" << endl;
}

void Pack::SanityCheckReceivingDimension() {
  if (IsInputScalar(1)) {
    cout << "Error: Trying to pack into a scalar" << endl;
    LOG_FAIL("replacement for throw call");
  }
}

void Pack::Prop() {
  if (!IsValidCost(m_cost)) {
    DLAOp<2, 1>::Prop();

    SanityCheckReceivingDimension();

    m_cost = GetInputM(1)->SumProds11(*GetInputN(1));
  }
}

int Pack::PackM() {
  return (int) GetInputNumRows(0);
}

int Pack::PackN() {
  return (int) GetInputNumCols(0);
}

#endif // DOLLDLA
