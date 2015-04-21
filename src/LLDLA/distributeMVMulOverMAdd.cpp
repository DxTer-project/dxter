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

#include "distributeMVMulOverMAdd.h"

#if DOLLDLA

#include "madd.h"
#include "mvmul.h"

string DistributeMVMulOverMAdd::GetType() const {
  return "DistributeMVMulOverMAdd";
}

bool DistributeMVMulOverMAdd::CanApply(const Node* node) const {
  if (node->GetNodeClass() == MVMul::GetClass()) {
    auto mvmul = static_cast<const MVMul*>(node);
    if (mvmul->GetLayer() != m_fromLayer) {
      return false;
    }
    auto in0 = mvmul->Input(0);
    if (in0->GetNodeClass() == MAdd::GetClass() && in0->NumChildrenOfOutput(0) == 1) {      
      return true;
    } else {
      return false;
    }
  }
  throw;
}

void DistributeMVMulOverMAdd::Apply(Node* node) const {
  auto madd = static_cast<MAdd*>(node->Input(0));

  auto newMVMul1 = new MVMul(m_toLayer);
  newMVMul1->AddInputs(6,
		       madd->Input(0), madd->InputConnNum(0),
		       node->Input(1), node->InputConnNum(1),
		       node->Input(2), node->InputConnNum(2));

  auto newMVMul2 = new MVMul(m_toLayer);
  newMVMul2->AddInputs(6,
		       madd->Input(1), madd->InputConnNum(1),
		       node->Input(1), node->InputConnNum(1),
		       newMVMul1, 0);

  node->m_poss->AddNode(newMVMul1);
  node->m_poss->AddNode(newMVMul2);

  node->RedirectChildren(newMVMul2, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;
}

#endif // DOLLDLA
