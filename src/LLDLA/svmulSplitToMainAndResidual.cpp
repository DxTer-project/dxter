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

#include "svmulSplitToMainAndResidual.h"

#if DOLLDLA

#include "packingUtils.h"
#include "svmul.h"

SVMulSplitToMainAndResidual::SVMulSplitToMainAndResidual(Layer fromLayer, Layer toLayer, VecType vecType) {
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_vecType = vecType;
}

bool SVMulSplitToMainAndResidual::CanApply(const Node* node) const {
  if (node->GetNodeClass() == SVMul::GetClass()) {
    const SVMul* svmul = static_cast<const SVMul*>(node);
    if (svmul->GetLayer() != m_fromLayer) {
      return false;
    }
    if (svmul->GetVecType() != m_vecType) {
      return false;
    }
    if (m_vecType == ROWVECTOR) {
      return !(svmul->InputNIsMultipleOfVecRegWidth(1))
	&& (svmul->GetInputNumCols(1) > svmul->GetVecRegWidth());
    } else {
      return !(svmul->InputMIsMultipleOfVecRegWidth(1))
	&& (svmul->GetInputNumRows(1) > svmul->GetVecRegWidth());
    }
  }
  throw;
}

void SVMulSplitToMainAndResidual::Apply(Node* node) const {
  DimName splitDim = m_vecType == ROWVECTOR ? DIMN : DIMM;
  Dir partDir = splitDim == DIMN ? HORIZONTAL : VERTICAL;
  auto partVec = PartitionIntoMainAndResidual(m_toLayer, node->Input(1), node->InputConnNum(1), node, 1, splitDim, node->GetVecRegWidth());

  auto mainSVMul = new SVMul(m_toLayer);
  mainSVMul->AddInputs(4,
		       node->Input(0), node->InputConnNum(0),
		       partVec, 0);

  auto residualSVMul = new SVMul(m_toLayer);
  residualSVMul->AddInputs(4,
		       node->Input(0), node->InputConnNum(0),
		       partVec, 1);

  auto rec = new Recombine(m_toLayer, partDir);
  rec->AddInputs(6,
		 mainSVMul, 0,
		 residualSVMul, 0,
		 node->Input(1), node->InputConnNum(1));

  node->m_poss->AddNode(partVec);
  node->m_poss->AddNode(mainSVMul);
  node->m_poss->AddNode(residualSVMul);
  node->m_poss->AddNode(rec);

  node->RedirectChildren(rec, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

#endif // DOLLDLA
