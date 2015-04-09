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

#include "eliminatePackedStoreLoad.h"

#if DOLLDLA

#include "maskedLoad.h"
#include "regLoadStore.h"

bool EliminatePackedStoreLoad::CanApply(const Node* node) const {
  if (node->GetNodeClass() == UnpackStoreFromRegs::GetClass()) {
    if (node->NumChildrenOfOutput(0) == 2) {
      if ((node->Child(0)->GetNodeClass() == PackedLoadToRegs::GetClass() ||
	   node->Child(0)->GetNodeClass() == MaskedLoad::GetClass())
	  && node->Child(1)->GetNodeClass() == UnpackStoreFromRegs::GetClass()) {
	return true;
      } else {
	/*	cout << node->Child(0)->GetNodeClass() << endl;
	cout << node->Child(0)->Child(0)->GetNodeClass() << endl;
	cout << node->Child(0)->Child(0)->Child(0)->GetNodeClass() << endl;
	cout << node->Child(0)->Child(0)->Child(0)->Child(0)->GetNodeClass() << endl;
	cout << node->Child(0)->Child(0)->Child(0)->Child(0)->Child(0)->GetNodeClass() << endl;*/
	return false;
      }
    }
    return false;
  }
  throw;
}

void EliminatePackedStoreLoad::Apply(Node* node) const {
  auto superfluousLoad = node->Child(0);
  auto finalPackedStore = node->Child(1);

  auto newFinalPackedStore = new UnpackStoreFromRegs();
  newFinalPackedStore->AddInputs(4,
			   finalPackedStore->Input(0), finalPackedStore->InputConnNum(0),
			   node->Input(1), node->InputConnNum(1));

  finalPackedStore->RedirectChildren(newFinalPackedStore, 0);

  node->m_poss->AddNode(newFinalPackedStore);
  node->m_poss->DeleteChildAndCleanUp(finalPackedStore);

  superfluousLoad->RedirectChildren(node->Input(0), node->InputConnNum(0));

  node->m_poss->DeleteChildAndCleanUp(superfluousLoad);
  return;
}

#endif // DOLLDLA
