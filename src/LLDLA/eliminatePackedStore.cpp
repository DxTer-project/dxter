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

#include "eliminatePackedStore.h"

#if DOLLDLA

#include "maskedLoad.h"
#include "regLoadStore.h"

bool EliminatePackedStore::CanApply(const Node* node) const {
  if (node->GetNodeClass() == UnpackStoreFromRegs::GetClass()) {
    if (node->NumChildrenOfOutput(0) == 1) {
      /*	cout << node->Child(0)->GetNodeClass() << endl;
	cout << node->Child(0)->Child(0)->GetNodeClass() << endl;
	cout << node->Child(0)->Child(0)->Child(0)->GetNodeClass() << endl;
	cout << node->Child(0)->Child(0)->Child(0)->Child(0)->GetNodeClass() << endl;
	cout << node->Child(0)->Child(0)->Child(0)->Child(0)->Child(0)->GetNodeClass() << endl << endl;*/
      if (node->Child(0)->GetNodeClass() == PackedLoadToRegs::GetClass() ||
	  node->Child(0)->GetNodeClass() == MaskedLoad::GetClass()) {
	return true;
      } else {
	return false;
      }
    }
    return false;
  }
  throw;
}

void EliminatePackedStore::Apply(Node* node) const {
  auto superfluousLoad = node->Child(0);

  superfluousLoad->RedirectChildren(node->Input(0), node->InputConnNum(0));

  node->m_poss->DeleteChildAndCleanUp(superfluousLoad);
  return;
}

#endif // DOLLDLA
