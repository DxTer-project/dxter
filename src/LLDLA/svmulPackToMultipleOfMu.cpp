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

#include "svmulPackToMultipleOfMu.h"

#if DOLLDLA

#include "horizontalUnpack.h"
#include "packingUtils.h"
#include "verticalUnpack.h"
#include "SVMul.h"

bool SVMulVerticalPackToMultipleOfMu::CanApply(const Node* node) const {
  if (node->GetNodeClass() == SVMul::GetClass()) {
    const SVMul* svmul = static_cast<const SVMul*>(node);
    if (svmul->GetLayer() != m_fromLayer) {
      return false;
    }
    return !(svmul->InputMIsMultipleOfVecRegWidth(1))
      && svmul->GetInputNumRows(1) > 1
      && svmul->GetInputNumCols(1) == 1;
  }
  throw;
}

void SVMulVerticalPackToMultipleOfMu::Apply(Node* node) const {
  Pack* packA = PackToMultipleOf(m_toLayer, node->Input(1), node->InputConnNum(1), node, 1, DIMM, node->GetVecRegWidth());

  SVMul* newSVMul = new SVMul(ABSLAYER);
  newSVMul->AddInputs(4,
		      node->Input(0), node->InputConnNum(0),
		      packA, 0);

  Unpack* unpack = new VerticalUnpack(ABSLAYER);
  unpack->AddInputs(4,
		    newSVMul, 0,
		    node->Input(1), node->InputConnNum(1));

  node->m_poss->AddUp(node->m_poss->m_possNodes, unpack, false, true);
  node->RedirectChildren(unpack, 0);
  node->m_poss->DeleteChildAndCleanUp(node);

  return;
}

bool SVMulHorizontalPackToMultipleOfMu::CanApply(const Node* node) const {
  if (node->GetNodeClass() == SVMul::GetClass()) {
    const SVMul* svmul = static_cast<const SVMul*>(node);
    if (svmul->GetLayer() != m_fromLayer) {
      return false;
    }
    return !(svmul->InputNIsMultipleOfVecRegWidth(1))
      && svmul->GetInputNumCols(1) > 1
      && svmul->GetInputNumRows(1) == 1;
  }
  throw;
}

void SVMulHorizontalPackToMultipleOfMu::Apply(Node* node) const {
  Pack* packA = PackToMultipleOf(m_toLayer, node->Input(1), node->InputConnNum(1), node, 1, DIMN, node->GetVecRegWidth());

  SVMul* newSVMul = new SVMul(ABSLAYER);
  newSVMul->AddInputs(4,
		      node->Input(0), node->InputConnNum(0),
		      packA, 0);

  Unpack* unpack = new HorizontalUnpack(ABSLAYER);
  unpack->AddInputs(4,
		    newSVMul, 0,
		    node->Input(1), node->InputConnNum(1));

  node->m_poss->AddUp(node->m_poss->m_possNodes, unpack, false, true);
  node->RedirectChildren(unpack, 0);
  node->m_poss->DeleteChildAndCleanUp(node);

  return;
}

#endif // DOLLDLA
