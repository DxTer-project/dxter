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

#include "vvdotSplitToMainAndResidual.h"

#if DOLLDLA

#include "base.h"
#include "packingUtils.h"
#include "vvdot.h"

VVDotSplitToMainAndResidual::VVDotSplitToMainAndResidual(Layer fromLayer, Layer toLayer) {
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
}

bool VVDotSplitToMainAndResidual::CanApply(const Node* node) const {
  if (node->GetNodeClass() == VVDot::GetClass()) {
    const VVDot* vvdot = static_cast<const VVDot*>(node);
    if (vvdot->GetLayer() != m_fromLayer) {
      return false;
    }
    return !(vvdot->InputNIsMultipleOfVecRegWidth(0))
      && (vvdot->GetInputNumCols(0) > vvdot->GetVecRegWidth());
  }

  LOG_FAIL("replacement for throw call");
  throw;
}

void VVDotSplitToMainAndResidual::Apply(Node* node) const {
  auto partX = PartitionIntoMainAndResidual(m_toLayer, node->Input(0), node->InputConnNum(0), node, 0, DIMN, node->GetVecRegWidth());
  auto partY = PartitionIntoMainAndResidual(m_toLayer, node->Input(1), node->InputConnNum(1), node, 1, DIMM, node->GetVecRegWidth());

  auto mainDot = new VVDot(m_toLayer);
  mainDot->AddInputs(6,
		     partX, 0,
		     partY, 0,
		     node->Input(2), node->InputConnNum(2));

  auto residualDot = new VVDot(m_toLayer);
  residualDot->AddInputs(6,
			 partX, 1,
			 partY, 1,
			 mainDot, 0);

  node->m_poss->AddNode(partX);
  node->m_poss->AddNode(partY);
  node->m_poss->AddNode(mainDot);
  node->m_poss->AddNode(residualDot);
  
  node->RedirectChildren(residualDot, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;
}

#endif // DOLLDLA
