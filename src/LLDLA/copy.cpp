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

#include <assert.h>
#include "copy.h"

#if DOLLDLA

Copy::Copy(Layer layer) {
  m_layer = layer;
}

const DataTypeInfo& Copy::DataType(ConnNum num) const {
  assert(num == 0);
  return InputDataType(1);
}

bool Copy::Overwrites(const Node* input, ConnNum num) const {
  const NodeConn* conn = m_inputs[1];
  return conn->m_n == input && conn->m_num == num;
}

void Copy::Prop() {
  if (!IsValidCost(m_cost)) {
    DLAOp<2, 1>::Prop();
    if (!InputsAreSameSize(0, 1)) {
      cout << "Input(0)" << endl;
      cout << Input(0)->GetNodeClass() << endl;
      cout << "Input(1)" << endl;
      cout << Input(1)->GetNodeClass() << endl;
      cout << "ERROR: Copy operands do not have the same size" << endl;
      cout << "Input 0 # rows = " << GetInputNumRows(0) << endl;
      cout << "Input 0 # cols = " << GetInputNumCols(0) << endl;
      cout << "Input 1 # rows = " << GetInputNumRows(1) << endl;
      cout << "Input 1 # cols = " << GetInputNumCols(1) << endl;
      this->m_poss->PrintTransVecUp();

      LOG_FAIL("Copy::Prop failed, inputs are not the same size");
      throw;
    }
    m_cost = GetInputM(0)->SumProds11(*GetInputN(0));
  }
}

void Copy::PrintCode(IndStream& out) {
  const DataTypeInfo& src = InputDataType(0);
  const DataTypeInfo& dst = InputDataType(1);


  out.Indent();
  if (Input(0)->GetDataType() == REAL_SINGLE) {
    *out << "copy_float(";
  } else {
    *out << "copy_double(";
  }

  *out << src.m_numRowsVar << ", " << src.m_numColsVar << ", ";
  *out << GetInputName(0).m_name << ", ";
  *out << src.m_rowStrideVar << ", " << src.m_colStrideVar << ", ";
  *out << GetInputName(1).m_name << ", ";
  *out << dst.m_rowStrideVar << ", " << dst.m_colStrideVar << ");" << endl;
}

#endif // DOLLDLA
