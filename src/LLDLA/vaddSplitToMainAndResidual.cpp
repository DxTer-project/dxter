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

#include "vaddSplitToMainAndResidual.h"

#if DOLLDLA

#include "base.h"
#include "packingUtils.h"
#include "vadd.h"

VAddSplitToMainAndResidual::VAddSplitToMainAndResidual(Layer fromLayer, Layer toLayer, VecType vecType) {
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_vecType = vecType;
}

bool VAddSplitToMainAndResidual::CanApply(const Node* node) const {
  if (node->GetNodeClass() == VAdd::GetClass()) {
    const VAdd* vadd = static_cast<const VAdd*>(node);
    if (vadd->GetLayer() != m_fromLayer) {
      return false;
    }
    if (vadd->GetVecType() != m_vecType) {
      return false;
    }
    if (m_vecType == ROWVECTOR) {
      return !(vadd->InputNIsMultipleOfVecRegWidth(0))
	&& (vadd->GetInputNumCols(0) > vadd->GetVecRegWidth());
    } else {
      return !(vadd->InputMIsMultipleOfVecRegWidth(0))
	&& (vadd->GetInputNumRows(0) > vadd->GetVecRegWidth());
    }
  }
  LOG_FAIL("replacement for throw call");
  throw;
}

void VAddSplitToMainAndResidual::Apply(Node* node) const {
  DimName splitDim = m_vecType == ROWVECTOR ? DIMN : DIMM;
  auto recombine = SplitBinarySymmetricOperationIntoMainAndResidual(m_toLayer, node, splitDim, node->GetVecRegWidth());
  node->m_poss->AddUp(node->m_poss->m_possNodes, recombine, false, true);
  node->RedirectChildren(recombine, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;
}

#endif // DOLLDLA
