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

#include "eliminateRecombine.h"

#if DOLLDLA

#include "partition.h"
#include "recombine.h"

bool EliminateRecombine::PartitionsAreIdentical(const Node* node) const {
  auto firstRecombine = static_cast<const Recombine*>(node);
  auto part = static_cast<const Partition*>(firstRecombine->Child(0));
  if (*part->GetM(0) == firstRecombine->GetInputNumRows(0)
      && *part->GetN(0) == firstRecombine->GetInputNumCols(0)
      && *part->GetM(1) == firstRecombine->GetInputNumRows(1)
      && *part->GetN(1) == firstRecombine->GetInputNumCols(1)) {
    return true;
  } else {
    return false;
  }
}

bool EliminateRecombine::OutputIsSuperfluous(const Node* node) const {
  if (node->NumChildrenOfOutput(0) == 1) {
    if (node->Child(0)->GetNodeClass() == Partition::GetClass()) {
      return PartitionsAreIdentical(node);
    } else {
      return false;
    }
  } else {
    return false;
  }
}

bool EliminateRecombine::CanApply(const Node* node) const {
  if (node->GetNodeClass() == Recombine::GetClass()) {
    return OutputIsSuperfluous(node);
  }
  throw;
}

void EliminateRecombine::Apply(Node* node) const {
  auto superfluousPart = node->Child(0);

  superfluousPart->RedirectChildren(0, node->Input(0), node->InputConnNum(0));
  superfluousPart->RedirectChildren(1, node->Input(1), node->InputConnNum(1));

  node->m_poss->DeleteChildAndCleanUp(superfluousPart);
  return;
}

#endif // DOLLDLA
