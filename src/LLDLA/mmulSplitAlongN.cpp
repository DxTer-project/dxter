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

#include "mmul.h"
#include "mmulSplitAlongN.h"
#include "packingUtils.h"

#if DOLLDLA

bool MMulSplitAlongN::CanApply(const Node* node) const {
  if (node->GetNodeClass() == Gemm::GetClass()) {
    auto gemm = static_cast<const Gemm*>(node);
    if (gemm->GetLayer() != m_fromLayer) {
      return false;
    }
    return !(gemm->GetInputN(1)->EvenlyDivisibleBy(gemm->GetVecRegWidth())) &&
      *(gemm->GetInputN(1)) > gemm->GetVecRegWidth();
  }
  return false;
}

void MMulSplitAlongN::Apply(Node* node) const {
  auto oldGemm = static_cast<Gemm*>(node);

  auto partB = PartitionIntoMainAndResidual(m_toLayer, oldGemm->Input(1), oldGemm->InputConnNum(1), oldGemm, 1, DIMN, oldGemm->GetVecRegWidth());
  auto partC = PartitionIntoMainAndResidual(m_toLayer, oldGemm->Input(2), oldGemm->InputConnNum(2), oldGemm, 2, DIMN, oldGemm->GetVecRegWidth());

  auto mainGemm = new Gemm(ABSLAYER, NORMAL, NORMAL, COEFONE, COEFONE, oldGemm->GetDataType());
  mainGemm->AddInputs(6,
		      oldGemm->Input(0), oldGemm->InputConnNum(0),
		      partB, 0,
		      partC, 0);

  auto residualGemm = new Gemm(ABSLAYER, NORMAL, NORMAL, COEFONE, COEFONE, oldGemm->GetDataType());
  residualGemm->AddInputs(6,
			  oldGemm->Input(0), oldGemm->InputConnNum(0),
			  partB, 1,
			  partC, 1);			  

  auto rec = new Recombine(m_toLayer, VERTICAL);
  rec->AddInputs(6,
		 mainGemm, 0,
		 residualGemm, 0,
		 oldGemm->Input(2), node->InputConnNum(2));

  node->m_poss->AddNode(partB);
  node->m_poss->AddNode(partC);
  node->m_poss->AddNode(mainGemm);
  node->m_poss->AddNode(residualGemm);
  node->m_poss->AddNode(rec);

  node->RedirectChildren(rec, 0);
  node->m_poss->DeleteChildAndCleanUp(oldGemm);

  return;
}

#endif // DOLLDLA
