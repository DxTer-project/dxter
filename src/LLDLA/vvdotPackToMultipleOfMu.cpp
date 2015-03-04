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

#include "vvdotPackToMultipleOfMu.h"

#if DOLLDLA

#include "packingUtils.h"
#include "vvdot.h"

bool VVDotPackToMultipleOfMu::CanApply(const Node* node) const {
  const VVDot* dot = static_cast<const VVDot*>(node);
  if (dot->GetLayer() != m_fromLayer) {
    return false;
  }

  return !(dot->InputNIsMultipleOfVecRegWidth(0))
    && dot->GetInputNumRows(0) > 0;
}

void VVDotPackToMultipleOfMu::Apply(Node* node) const {
  cout << "Applying vvdot pack" << endl;
  //  throw;
  Pack* packX = PackToMultipleOf(m_toLayer, node->Input(0), node->InputConnNum(0), node, 0, DIMN, node->GetVecRegWidth());

  Pack* packY = PackToMultipleOf(m_toLayer, node->Input(1), node->InputConnNum(1), node, 1, DIMM, node->GetVecRegWidth());

  VVDot* newDot = new VVDot(m_toLayer);
  newDot->AddInputs(6,
		    packX, 0,
		    packY, 0,
		    node->Input(2), node->InputConnNum(2));

  //  Unpack* unpack = PackBinarySymmetricOperation(m_toLayer, node, DIMM, node->GetVecRegWidth());
  node->m_poss->AddUp(node->m_poss->m_possNodes, newDot, false, true);
  node->RedirectChildren(newDot, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

#endif // DOLLDLA
