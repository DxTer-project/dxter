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

#include "svmulAddSplitToMainAndResidual.h"

#if DOLLDLA

#include "packingUtils.h"
#include "svmulAdd.h"

SVMulAddSplitToMainAndResidual::SVMulAddSplitToMainAndResidual(Layer fromLayer, Layer toLayer, VecType vecType) {
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_vecType = vecType;
}

bool SVMulAddSplitToMainAndResidual::CanApply(const Node* node) const {
  if (node->GetNodeClass() == SVMulAdd::GetClass()) {
    const SVMulAdd* svmul = static_cast<const SVMulAdd*>(node);
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

void SVMulAddSplitToMainAndResidual::Apply(Node* node) const {
  auto svmulAdd = static_cast<SVMulAdd*>(node);
  DimName splitDim = m_vecType == ROWVECTOR ? DIMN : DIMM;
  Dir partDir = splitDim == DIMN ? HORIZONTAL : VERTICAL;

  auto partVecA = PartitionIntoMainAndResidual(m_toLayer, svmulAdd->Input(1), svmulAdd->InputConnNum(1), svmulAdd, 1, splitDim, svmulAdd->GetVecRegWidth());
  auto partVecB = PartitionIntoMainAndResidual(m_toLayer, svmulAdd->Input(2), svmulAdd->InputConnNum(2), svmulAdd, 2, splitDim, svmulAdd->GetVecRegWidth());

  auto mainSVMulAdd = new SVMulAdd(m_toLayer, svmulAdd->GetVecType());
  mainSVMulAdd->AddInputs(6,
			  svmulAdd->Input(0), svmulAdd->InputConnNum(0),
			  partVecA, 0,
			  partVecB, 0);

  auto residualSVMulAdd = new SVMulAdd(m_toLayer, svmulAdd->GetVecType());
  residualSVMulAdd->AddInputs(6,
			      svmulAdd->Input(0), svmulAdd->InputConnNum(0),
			      partVecA, 1,
			      partVecB, 1);

  auto rec = new Recombine(m_toLayer, partDir);
  rec->AddInputs(6,
		 mainSVMulAdd, 0,
		 residualSVMulAdd, 0,
		 svmulAdd->Input(2), svmulAdd->InputConnNum(2));

  svmulAdd->m_poss->AddNode(partVecA);
  svmulAdd->m_poss->AddNode(partVecB);
  svmulAdd->m_poss->AddNode(mainSVMulAdd);
  svmulAdd->m_poss->AddNode(residualSVMulAdd);
  svmulAdd->m_poss->AddNode(rec);

  svmulAdd->RedirectChildren(rec, 0);
  svmulAdd->m_poss->DeleteChildAndCleanUp(svmulAdd);
}

#endif // DOLLDLA
