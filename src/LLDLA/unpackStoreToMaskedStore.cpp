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

#include "unpackStoreToMaskedStore.h"

#if DOLLDLA

#include "maskedStore.h"
#include "regLoadStore.h"

bool UnpackStoreToMaskedStore::CanApply(const Node* node) const {
  if (node->GetNodeClass() == UnpackStoreFromRegs::GetClass()) {
    auto store = static_cast<const UnpackStoreFromRegs*>(node);
    auto rs = store->InputDataType(1).m_rowStride;
    auto cs = store->InputDataType(1).m_colStride;
    return (store->IsInputColVector(1) && IsUnitStride(rs)) ||
      (store->IsInputRowVector(1) && IsUnitStride(cs));
  }
  throw;
}

void UnpackStoreToMaskedStore::Apply(Node* node) const {
  auto maskedStore = new MaskedStore();
  maskedStore->AddInputs(4,
			 node->Input(0), node->InputConnNum(0),
			 node->Input(1), node->InputConnNum(1));

  node->RedirectChildren(maskedStore, 0);

  node->m_poss->AddNode(maskedStore);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;
}

#endif // DOLLDLA
