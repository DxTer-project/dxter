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

#include "vrwSVMulAddToRegArith.h"

#if DOLLDLA

#include "regArith.h"
#include "regLoadStore.h"

VRWSVMulAddToRegArith::VRWSVMulAddToRegArith(Layer fromLayer, Layer toLayer) {
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
}

string VRWSVMulAddToRegArith::GetType() const {
  return "VRWSVMulAddToRegArith";
}

bool VRWSVMulAddToRegArith::IsExactSizeSVMulAdd(const SVMulAdd* svmul) const {
  if (svmul->InputIsVRW(1)) {
    return true;
  } else {
    cout << "Input 1 to SVMulAdd is not VRW" << endl;
    cout << "Num rows: " << svmul->GetInputNumRows(1) << endl;
    cout << "Num cols: " << svmul->GetInputNumCols(1) << endl;
    svmul->GetInputM(1)->Print();
    svmul->GetInputN(1)->Print();
    return false;
  }
}

bool VRWSVMulAddToRegArith::CanApply(const Node* node) const {
  if (node->GetNodeClass() == SVMulAdd::GetClass()) {
    cout << "Checking SVMulAdd VRW conversion" << endl;
    const SVMulAdd* svmul = static_cast<const SVMulAdd*>(node);
    return IsExactSizeSVMulAdd(svmul);
  }
  LOG_FAIL("Bad class for node in VRWSVMulAddToRegArith::CanApply");
  throw;
}

void VRWSVMulAddToRegArith::Apply(Node* node) const {
  cout << "Applying VRWSVMulAddToRegArith" << endl;
  auto svmul = static_cast<SVMulAdd*>(node);

  auto dup = new DuplicateRegLoad();
  dup->AddInput(svmul->Input(0), svmul->InputConnNum(0));

  auto loadA = new LoadToRegs();
  loadA->AddInput(svmul->Input(1), svmul->InputConnNum(1));

  auto mul = new Mul();
  mul->AddInput(dup, 0);
  mul->AddInput(loadA, 0);

  auto loadB = new LoadToRegs();
  loadB->AddInput(svmul->Input(2), svmul->InputConnNum(2));

  auto add = new Add();
  add->AddInputs(4,
		 mul, 0,
		 loadB, 0);

  auto storeVec = new StoreFromRegs();
  storeVec->AddInput(add, 0);
  storeVec->AddInput(node->Input(2), node->InputConnNum(2));

  node->m_poss->AddNode(dup);
  node->m_poss->AddNode(loadA);
  node->m_poss->AddNode(loadB);
  node->m_poss->AddNode(mul);
  node->m_poss->AddNode(add);
  node->m_poss->AddNode(storeVec);

  node->RedirectChildren(storeVec, 0);
  node->m_poss->DeleteChildAndCleanUp(node);

  return;
}

#endif // DOLLDLA
