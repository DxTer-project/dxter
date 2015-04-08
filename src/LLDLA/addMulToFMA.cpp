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

#include "addMulToFMA.h"

#if DOLLDLA

#include "regArith.h"

bool AddMulToFMA::CanApply(const Node* node) const {
  if (node->GetNodeClass() != Mul::GetClass()) {
    throw;
  }

  return node->NumChildrenOfOutput(0) == 1 &&
    node->Child(0)->GetNodeClass() == Add::GetClass() &&
    node->Child(0)->Input(0) == node &&
    node->Child(0)->Input(1) != node &&
    node->Child(0)->Input(1) != node->Input(0);
}

void AddMulToFMA::Apply(Node* node) const {
  auto fma = new FMAdd();
  fma->AddInputs(6,
		 node->Input(0), node->InputConnNum(0),
		 node->Input(1), node->InputConnNum(1),
		 node->Child(0)->Input(1), node->Child(0)->InputConnNum(1));

  node->Child(0)->RedirectChildren(fma, 0);
  node->m_poss->AddNode(fma);
  node->m_poss->DeleteChildAndCleanUp(node->Child(0));
  return;
}

#endif // DOLLDLA
