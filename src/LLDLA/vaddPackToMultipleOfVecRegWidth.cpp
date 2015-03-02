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

#include "vaddPackToMultipleOfVecRegWidth.h"

#if DOLLDLA

#include "packingUtils.h"
#include "vadd.h"


VAddPackToMultipleOfVecRegWidth::VAddPackToMultipleOfVecRegWidth(Layer fromLayer, Layer toLayer, DimName dim) {
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_dim = dim;
  m_vecType = dim == DIMM ? COLVECTOR : ROWVECTOR;
}

bool VAddPackToMultipleOfVecRegWidth::CanApply(const Node* node) const {
  if (node->GetNodeClass() == VAdd::GetClass()) {
    VAdd* vadd = static_cast<VAdd*>(node);
    
    if (vadd->GetVecType() == ROWVECTOR) {
      return !(vadd->InputNIsMultipleOfVecRegWidth(0))
	&& vadd->GetInputNumRows(0) == 1
	&& (((int) vadd->GetInputNumCols(0)) % vadd->GetVecRegWidth() != 0);
    } else {
      return !(vadd->InputMIsMultipleOfVecRegWidth(0))
	&& vadd->GetInputNumCols(0) == 1
	&& (((int) vadd->GetInputNumRows(0)) % vadd->GetVecRegWidth() != 0);
    }
    return false;
  }
  throw;
}

void VAddPackToMultipleOfVecRegWidth::Apply(Node* node) const {
  cout << "Applying VAddPackToMultipleOfVecRegWidth" << endl;
  DLANode* dlaNode = static_cast<DLANode*>(node);
  cout << "Node input 0 # rows = " << dlaNode->GetInputNumRows(0) << endl;
  cout << "Node input 0 # cols = " << dlaNode->GetInputNumCols(0) << endl;
  cout << "Node input 1 # rows = " << dlaNode->GetInputNumRows(1) << endl;
  cout << "Node input 1 # cols = " << dlaNode->GetInputNumCols(1) << endl;
  cout << "Calling PackBinarySymmetricOperation" << endl;
  Unpack* unpack = PackBinarySymmetricOperation(node, m_dim, node->GetVecRegWidth());
  cout << "Done with PackBinarySymmetricOperation" << endl;
  node->m_poss->AddUp(node->m_poss->m_possNodes, unpack, false, true);
  node->RedirectChildren(unpack, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;
}

#endif // DOLLDLA
