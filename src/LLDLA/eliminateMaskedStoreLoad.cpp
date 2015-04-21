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

#include "eliminateMaskedStoreLoad.h"

#if DOLLDLA

#include "maskedLoad.h"
#include "maskedStore.h"
#include "regLoadStore.h"

bool EliminateMaskedStoreLoad::CanApply(const Node* node) const {
  if (node->GetNodeClass() == MaskedStore::GetClass()) {
    if (node->NumChildrenOfOutput(0) == 2) {
      if ((node->Child(0)->GetNodeClass() == PackedLoadToRegs::GetClass() ||
	   node->Child(0)->GetNodeClass() == MaskedLoad::GetClass())
	  && (node->Child(1)->GetNodeClass() == UnpackStoreFromRegs::GetClass()
	      || node->Child(1)->GetNodeClass() == MaskedStore::GetClass())) {
	return true;
      } else {
	return false;
      }
    }
    return false;
  }
  throw;
}

void EliminateMaskedStoreLoad::Apply(Node* node) const {
  auto superfluousLoad = node->Child(0);
  auto finalMaskedStore = node->Child(1);

  auto newFinalMaskedStore = new UnpackStoreFromRegs();
  newFinalMaskedStore->AddInputs(4,
			   finalMaskedStore->Input(0), finalMaskedStore->InputConnNum(0),
			   node->Input(1), node->InputConnNum(1));

  finalMaskedStore->RedirectChildren(newFinalMaskedStore, 0);

  node->m_poss->AddNode(newFinalMaskedStore);
  node->m_poss->DeleteChildAndCleanUp(finalMaskedStore);

  superfluousLoad->RedirectChildren(node->Input(0), node->InputConnNum(0));

  node->m_poss->DeleteChildAndCleanUp(superfluousLoad);
  return;
}

#endif // DOLLDLA
