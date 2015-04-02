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

#include "vaddPackResidualToVRW.h"

#if DOLLDLA

#include "packingUtils.h"
#include "unpack.h"
#include "vadd.h"

bool VAddPackResidualToVRW::CanApply(const Node* node) const {
  if (node->GetNodeClass() == VAdd::GetClass()) {
    const VAdd* vadd = static_cast<const VAdd*>(node);

    if (vadd->GetLayer() != m_fromLayer) {
      return false;
    }
    if (vadd->GetVecType() != m_vecType) {
      return false;
    }

    if (m_vecType == ROWVECTOR) {
      return vadd->GetInputNumCols(0) < vadd->GetVecRegWidth();
    } else {
      return vadd->GetInputNumRows(0) < vadd->GetVecRegWidth();
    }
  }
  LOG_FAIL("Bad node type in VAddPackResidualToVRW");
  throw;
}

void VAddPackResidualToVRW::Apply(Node* node) const {
  DimName packDim = m_vecType == ROWVECTOR ? DIMN : DIMM;
  Unpack* unpack = PackBinarySymmetricOperation(m_toLayer, node, packDim, node->GetVecRegWidth());
  node->m_poss->AddUp(node->m_poss->m_possNodes, unpack, false, true);
  node->RedirectChildren(unpack, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;
}

#endif // DOLLDLA
