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

#include "vrwSVMulToRegArith.h"

#if DOLLDLA

#include "regArith.h"
#include "regLoadStore.h"

VRWSVMulToRegArith::VRWSVMulToRegArith(Layer fromLayer, Layer toLayer) {
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
}

string VRWSVMulToRegArith::GetType() const {
  return "VRWSVMulToRegArith";
}

bool VRWSVMulToRegArith::IsExactSizeSVMul(const SVMul* vadd) const {
  if (vadd->GetVecType() == ROWVECTOR) {
    return *vadd->GetInputM(1) == ONE &&
      vadd->GetInputNumCols(1) == vadd->GetVecRegWidth();
  } else {
    return *vadd->GetInputN(1) == ONE &&
      vadd->GetInputNumRows(1) == vadd->GetVecRegWidth();
  }
}

bool VRWSVMulToRegArith::CanApply(const Node* node) const {
  if (node->GetNodeClass() == SVMul::GetClass()) {
    const SVMul* svmul = static_cast<const SVMul*>(node);
    return IsExactSizeSVMul(svmul);
  }
  LOG_FAIL("Bad class for node in VRWSVMulToRegArith::CanApply");
  throw;
}

void VRWSVMulToRegArith::Apply(Node* node) const {
  SVMul* svmul = static_cast<SVMul*>(node);

  DuplicateRegLoad* dup = new DuplicateRegLoad();
  dup->AddInput(svmul->Input(0), svmul->InputConnNum(0));

  LoadToRegs* loadA = new LoadToRegs();
  loadA->AddInput(svmul->Input(1), svmul->InputConnNum(1));

  Mul* mul = new Mul();
  mul->AddInput(dup, 0);
  mul->AddInput(loadA, 0);

  StoreFromRegs* storeVec = new StoreFromRegs();
  storeVec->AddInput(mul, 0);
  storeVec->AddInput(node->Input(1), node->InputConnNum(1));

  node->m_poss->AddNode(dup);
  node->m_poss->AddNode(loadA);
  node->m_poss->AddNode(mul);
  node->m_poss->AddNode(storeVec);

  node->RedirectChildren(storeVec, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

#endif // DOLLDLA
