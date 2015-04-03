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

#include "vrwVVDotToRegArith.h"

#if DOLLDLA

#include "regArith.h"
#include "regLoadStore.h"

VRWVVDotToRegArith::VRWVVDotToRegArith(Layer fromLayer, Layer toLayer) {
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
}

string VRWVVDotToRegArith::GetType() const {
  return "VRWVVDotToRegArith";
}

bool VRWVVDotToRegArith::IsExactSizeVVDot(const VVDot* vvdot) const {
    return *vvdot->GetInputM(0) == ONE &&
      vvdot->GetInputNumCols(0) == vvdot->GetVecRegWidth() &&
      *vvdot->GetInputN(1) == ONE &&
      vvdot->GetInputNumRows(1) == vvdot->GetVecRegWidth();
}

bool VRWVVDotToRegArith::CanApply(const Node* node) const {
  if (node->GetNodeClass() == VVDot::GetClass()) {
    const VVDot* vvdot = static_cast<const VVDot*>(node);
    return IsExactSizeVVDot(vvdot);
  }
  LOG_FAIL("Bad class for node in VRWVVDotToRegArith::CanApply");
  throw;
}

void VRWVVDotToRegArith::Apply(Node* node) const {
  auto loadX = new LoadToRegs();
  loadX->AddInput(node->Input(0), node->InputConnNum(0));

  auto loadY = new LoadToRegs();
  loadY->AddInput(node->Input(1), node->InputConnNum(1));

  auto accum = new TempVecReg();
  accum->AddInput(node->Input(2), node->InputConnNum(2));

  auto fma = new FMAdd();
  fma->AddInputs(6,
		 loadX, 0,
		 loadY, 0,
		 accum, 0);

  auto accumInC = new AccumReg();
  accumInC->AddInput(fma, 0);
  accumInC->AddInput(node->Input(2), node->InputConnNum(2));

  node->m_poss->AddNode(loadX);
  node->m_poss->AddNode(loadY);
  node->m_poss->AddNode(accum);
  node->m_poss->AddNode(fma);
  node->m_poss->AddNode(accumInC);

  node->RedirectChildren(accumInC, 0);
  node->m_poss->DeleteChildAndCleanUp(node);

  return;
}

#endif // DOLLDLA
