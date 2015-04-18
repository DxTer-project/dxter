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

#include "residualSVMulAddToRegArith.h"

#if DOLLDLA

#include "regArith.h"
#include "regLoadStore.h"
#include "svmulAdd.h"

string ResidualSVMulAddToRegArith::GetType() const {
  switch(m_vecType)
    {
    case(ROWVECTOR):
      return "ResidualSVMulAddToRegArith - row vector";
    case(COLVECTOR):
      return "ResidualSVMulAddToRegArith - column vector";
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }
}

bool ResidualSVMulAddToRegArith::CanApply(const Node* node) const {
  if (node->GetNodeClass() == SVMulAdd::GetClass()) {
    const SVMulAdd* svmulAdd = static_cast<const SVMulAdd*>(node);
    if (svmulAdd->GetLayer() != m_fromLayer) {
      return false;
    }
    if (svmulAdd->GetVecType() != m_vecType) {
      return false;
    }
    if (m_vecType == ROWVECTOR) {
      return svmulAdd->GetInputNumCols(1) < svmulAdd->GetVecRegWidth();
    } else {
      return svmulAdd->GetInputNumRows(1) < svmulAdd->GetVecRegWidth();
    }
  }
  LOG_FAIL("bad type in SVMulAddPackToRegArith::CanApply");
  throw;
}

void ResidualSVMulAddToRegArith::Apply(Node* node) const {
  auto svmul = static_cast<SVMulAdd*>(node);

  auto dup = new DuplicateRegLoad();
  dup->AddInput(svmul->Input(0), svmul->InputConnNum(0));

  auto loadA = new PackedLoadToRegs();
  loadA->AddInput(svmul->Input(1), svmul->InputConnNum(1));

  auto mul = new Mul();
  mul->AddInput(dup, 0);
  mul->AddInput(loadA, 0);

  auto loadB = new PackedLoadToRegs();
  loadB->AddInput(svmul->Input(2), svmul->InputConnNum(2));

  auto add = new Add();
  add->AddInputs(4,
		 mul, 0,
		 loadB, 0);

  auto storeVec = new UnpackStoreFromRegs();
  storeVec->AddInput(add, 0);
  storeVec->AddInput(svmul->Input(2), svmul->InputConnNum(2));

  svmul->m_poss->AddNode(dup);
  svmul->m_poss->AddNode(loadA);
  svmul->m_poss->AddNode(loadB);
  svmul->m_poss->AddNode(mul);
  svmul->m_poss->AddNode(add);
  svmul->m_poss->AddNode(storeVec);

  svmul->RedirectChildren(storeVec, 0);
  svmul->m_poss->DeleteChildAndCleanUp(svmul);

  return;
}

#endif // DOLLDLA
