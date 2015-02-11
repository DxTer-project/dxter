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
  return DataType(num);
}

bool Copy::Overwrites(const Node* input, ConnNum num) const {
  const NodeConn* conn = m_inputs[1];
  return conn->m_n == input && conn->m_num == num;
}

void Copy::Prop() {
  if (!IsValidCost(m_cost)) {
    DLAOp<2, 1>::Prop();
    if (!InputDataType(0).IsSameSizeAs(InputDataType(1))) {
      cout << "ERROR: Copy operands do not have the same size" << endl;
      throw;
    }
    m_cost = ZERO;
  }
}

void Copy::PrintCode(IndStream& out) {
  cout << "ERROR: Copy is must be refined into a specific implementation. ";
  cout << "It cannot be implemented" << endl;
  throw;
}

#endif // DOLLDLA
