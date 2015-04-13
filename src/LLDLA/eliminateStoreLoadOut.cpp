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
    return false;
  } else {
    Node* n = comb->Child(0);
    while (!n->m_children.empty()) {
      n = n->Child(0);
      if (n->GetNodeClass() != CombineSingleIter::GetClass()) {
	return false;
      }
      if (n->NumChildrenOfOutput(0) > 1) {
	return false;
      }
    }
    return true;
  }
}

bool EliminateStoreLoadOut::CanApply(const Node* node) const {
  if (node->GetNodeClass() == StoreFromRegs::GetClass()) {
    if (node->NumChildrenOfOutput(0) == 2) {
      if (node->OutputHasChildOfClass(0, LoadToRegs::GetClass()) &&
	  node->OutputHasChildOfClass(0, CombineSingleIter::GetClass())) {
	auto comb = static_cast<const CombineSingleIter*>(node->ChildOfOutputWithClass(0, CombineSingleIter::GetClass()));
	return HasNoOutputs(comb);
      }
      return false;
    }
    return false;
  }
  throw;
}

void EliminateStoreLoadOut::Apply(Node* node) const {
  auto combine = node->ChildOfOutputWithClass(0, CombineSingleIter::GetClass());
  auto load = node->ChildOfOutputWithClass(0, LoadToRegs::GetClass());

  load->RedirectChildren(node->Input(0), node->InputConnNum(0));
  combine->ChangeInput2Way(combine->Input(1), combine->InputConnNum(1), node->Input(1), node->InputConnNum(1));

  node->m_poss->DeleteChildAndCleanUp(load);
  return;
}

#endif // DOLLDLA
