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

#include "scalarMulVerticalPackToMultipleOfMu.h"

#if DOLLDLA

#include "packingUtils.h"

bool ScalarMulVerticalPackToMultipleOfMu::CanApply(const Node* node) const {
  const DLANode* dlaNode = static_cast<const DLANode*>(node);
  if (dlaNode->GetLayer() != m_fromLayer) {
    return false;
  }
  if (dlaNode->GetNodeClass() != m_nodeTypeName) {
    return false;
  }
  return !(dlaNode->InputMIsMultipleOfVecRegWidth(1))
    && dlaNode->GetInputNumRows(1) > 1
    && dlaNode->GetInputNumCols(1) > 0;
}

void ScalarMulVerticalPackToMultipleOfMu::Apply(Node* node) const {
  //  Pack* packA = PackToMultipleOf(m_toLayer, node->Input(1), node->InputConnNum(1), node, 1, DIMM, node->GetVecRegWidth());

  
  return;
}

#endif // DOLLDLA
