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

#include "unpack.h"

#if DOLLDLA

Unpack::Unpack(Layer layer) {
  m_layer = layer;
}

const DataTypeInfo& Unpack::DataType(ConnNum num) const {
  if (num != 0) {
    throw;
  }
  return InputDataType(1);
}

bool Unpack::Overwrites(const Node* input, ConnNum num) const {
  const NodeConn* conn = m_inputs[1];
  return conn->m_n == input && conn->m_num == num;
}

void Unpack::PrintCode(IndStream& out) {
  throw;
}

void Unpack::Prop() {
  if (!IsValidCost(m_cost)) {
    DLAOp<2, 1>::Prop();
    m_cost = ZERO;
  }
}

#endif // DOLLDLA
