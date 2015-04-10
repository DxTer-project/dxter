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

#include "vrwVAddToRegArith.h"

#if DOLLDLA

#include "regArith.h"
#include "regLoadStore.h"

VRWVAddToRegArith::VRWVAddToRegArith(Layer fromLayer, Layer toLayer) {
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
}

string VRWVAddToRegArith::GetType() const {
  return "VRWVAddToRegArith";
}

bool VRWVAddToRegArith::IsExactSizeVAdd(const VAdd* vadd) const {
  return vadd->InputIsVRW(0) && vadd->InputIsVRW(1);
}

bool VRWVAddToRegArith::CanApply(const Node* node) const {
  if (node->GetNodeClass() == VAdd::GetClass()) {
    const VAdd* vadd = static_cast<const VAdd*>(node);
    return IsExactSizeVAdd(vadd);
  }
  LOG_FAIL("Bad class for node in VRWVAddToRegArith::CanApply");
  throw;
}

void VRWVAddToRegArith::Apply(Node* node) const {
  LoadToRegs* loadX = new LoadToRegs();
  loadX->AddInput(node->Input(0), node->InputConnNum(0));

  LoadToRegs* loadY = new LoadToRegs();
  loadY->AddInput(node->Input(1), node->InputConnNum(1));

  Add* add = new Add();
  add->AddInput(loadX, 0);
  add->AddInput(loadY, 0);

  StoreFromRegs* storeToY = new StoreFromRegs();
  storeToY->AddInput(add, 0);
  storeToY->AddInput(node->Input(1), node->InputConnNum(1));

  node->m_poss->AddNode(loadX);
  node->m_poss->AddNode(loadY);
  node->m_poss->AddNode(add);
  node->m_poss->AddNode(storeToY);

  node->RedirectChildren(storeToY, 0);
  node->m_poss->DeleteChildAndCleanUp(node);

  return;
}

#endif // DOLLDLA
