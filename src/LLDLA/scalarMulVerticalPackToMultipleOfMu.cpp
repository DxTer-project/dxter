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
#include "smmul.h"
#include "verticalUnpack.h"

bool ScalarMulVerticalPackToMultipleOfMu::CanApply(const Node* node) const {
  if (node->GetNodeClass() == SMMul::GetClass()) {
    const SMMul* smmul = static_cast<const SMMul*>(node);
    if (smmul->GetLayer() != m_fromLayer) {
      return false;
    }
    return !(smmul->InputMIsMultipleOfVecRegWidth(1))
      && smmul->GetInputNumRows(1) > 1
      && smmul->GetInputNumCols(1) > 0;
  }
  throw;
}

void ScalarMulVerticalPackToMultipleOfMu::Apply(Node* node) const {
  Pack* packA = PackToMultipleOf(m_toLayer, node->Input(1), node->InputConnNum(1), node, 1, DIMM, node->GetVecRegWidth());

  SMMul* newSMMul = new SMMul(ABSLAYER);
  newSMMul->AddInputs(4,
		      node->Input(0), node->InputConnNum(0),
		      packA, 0);

  Unpack* unpack = new VerticalUnpack(ABSLAYER);
  unpack->AddInputs(4,
		    newSMMul, 0,
		    node->Input(1), node->InputConnNum(1));

  node->m_poss->AddUp(node->m_poss->m_possNodes, unpack, false, true);
  node->RedirectChildren(unpack, 0);
  node->m_poss->DeleteChildAndCleanUp(node);

  
  return;
}

#endif // DOLLDLA
