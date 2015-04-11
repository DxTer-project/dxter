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

#include "eliminateStoreLoadOut.h"

#if DOLLDLA

#include "regLoadStore.h"

bool EliminateStoreLoadOut::HasNoOutputs(const CombineSingleIter* comb) const {
  if (comb->NumChildrenOfOutput(0) > 1) {
    cout << "More than one child of combine" << endl;
    return false;
  } else {
    cout << "One child of combine" << endl;
    cout << "Child type is" << endl;
    cout << comb->Child(0)->GetNodeClass() << endl;
    cout << "Child of combine has children? ";
    cout << (comb->Child(0)->m_children.empty() ? "NO" : "YES") << endl;
    return (comb->Child(0)->GetNodeClass() == CombineSingleIter::GetClass()) &&
      comb->Child(0)->m_children.empty();
  }
}

bool EliminateStoreLoadOut::CanApply(const Node* node) const {
  if (node->GetNodeClass() == StoreFromRegs::GetClass()) {
    //    cout << "store from regs found" << endl;
    if (node->NumChildrenOfOutput(0) == 2) {
      //      cout << "with 2 children" << endl;
      if (node->Child(1)->GetNodeClass() == LoadToRegs::GetClass()
	  && node->Child(0)->GetNodeClass() == CombineSingleIter::GetClass()) {
	auto comb = static_cast<const CombineSingleIter*>(node->Child(0));
	return HasNoOutputs(comb);
      } else {
	return false;
      }
    }
    return false;
  }
  throw;
}

void EliminateStoreLoadOut::Apply(Node* node) const {
  auto combine = node->Child(0);
  auto load = node->Child(1);

  load->RedirectChildren(node->Input(0), node->InputConnNum(0));
  combine->ChangeInput2Way(combine->Input(1), combine->InputConnNum(1), node->Input(1), node->InputConnNum(1));

  node->m_poss->DeleteChildAndCleanUp(load);
  return;
}

#endif // DOLLDLA
