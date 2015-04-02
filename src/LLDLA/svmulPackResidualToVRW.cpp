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

#include "svmulPackResidualToVRW.h"

#if DOLLDLA

#include "horizontalUnpack.h"
#include "packingUtils.h"
#include "vadd.h"
#include "verticalUnpack.h"
#include "svmul.h"

string SVMulPackResidualToVRW::GetType() const {
  switch(m_vecType)
    {
    case(ROWVECTOR):
      return "SVMulPackResidualToVRW - row vector";
    case(COLVECTOR):
      return "SVMulPackResidualToVRW - column vector";
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }
}

bool SVMulPackResidualToVRW::CanApply(const Node* node) const {
  if (node->GetNodeClass() == SVMul::GetClass()) {
    const SVMul* svmul = static_cast<const SVMul*>(node);
    if (svmul->GetLayer() != m_fromLayer) {
      return false;
    }
    if (svmul->GetVecType() != m_vecType) {
      return false;
    }
    if (m_vecType == ROWVECTOR) {
      return svmul->GetInputNumCols(1) < svmul->GetVecRegWidth();
    } else {
      return svmul->GetInputNumRows(1) < svmul->GetVecRegWidth();
    }
  }
  LOG_FAIL("bad type in SVMulPackToVRW::CanApply");
  throw;
}

void SVMulPackResidualToVRW::Apply(Node* node) const {
  cout << "Applying SVMulPackResidualToVRW" << endl;
  SVMul* svmul = static_cast<SVMul*>(node);
  svmul->GetInputM(1)->Print();
  svmul->GetInputN(1)->Print();
  cout << VecTypeToString(svmul->GetVecType()) << endl;
  DimName dim = m_vecType == ROWVECTOR ? DIMN : DIMM;
  Pack* packA = PackToMultipleOf(m_toLayer, node->Input(1), node->InputConnNum(1), node, 1, dim, node->GetVecRegWidth());

  SVMul* newSVMul = new SVMul(m_toLayer);
  newSVMul->AddInputs(4,
		      node->Input(0), node->InputConnNum(0),
		      packA, 0);

  Unpack* unpack;
  if (m_vecType == COLVECTOR) {
    unpack = new VerticalUnpack(m_toLayer);
  } else {
    unpack = new HorizontalUnpack(m_toLayer);
  }
  unpack->AddInputs(4,
		    newSVMul, 0,
		    node->Input(1), node->InputConnNum(1));

  node->m_poss->AddUp(node->m_poss->m_possNodes, unpack, false, true);
  node->RedirectChildren(unpack, 0);
  node->m_poss->DeleteChildAndCleanUp(node);

  return;
}

#endif // DOLLDLA
