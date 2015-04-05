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

#include "residualVAddToRegArith.h"

#if DOLLDLA

#include "regArith.h"
#include "regLoadStore.h"

ResidualVAddToRegArith::ResidualVAddToRegArith(Layer fromLayer, Layer toLayer) {
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
}

string ResidualVAddToRegArith::GetType() const {
  return "ResidualVAddToRegArith";
}

bool ResidualVAddToRegArith::IsResidualVAdd(const VAdd* vadd) const {
  if (vadd->GetVecType() == ROWVECTOR) {
    return *vadd->GetInputM(0) == ONE &&
      vadd->GetInputNumCols(0) < vadd->GetVecRegWidth();
  } else {
    return *vadd->GetInputN(0) == ONE &&
      vadd->GetInputNumRows(0) < vadd->GetVecRegWidth();
  }
}

bool ResidualVAddToRegArith::CanApply(const Node* node) const {
  if (node->GetNodeClass() == VAdd::GetClass()) {
    const VAdd* vadd = static_cast<const VAdd*>(node);
    return IsResidualVAdd(vadd);
  }
  LOG_FAIL("Bad class for node in ResidualVAddToRegArith::CanApply");
  throw;
}

void ResidualVAddToRegArith::Apply(Node* node) const {
  auto loadX = new PackedLoadToRegs();
  loadX->AddInput(node->Input(0), node->InputConnNum(0));

  auto loadY = new PackedLoadToRegs();
  loadY->AddInput(node->Input(1), node->InputConnNum(1));

  auto add = new Add();
  add->AddInput(loadX, 0);
  add->AddInput(loadY, 0);

  auto storeToY = new UnpackStoreFromRegs();
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
