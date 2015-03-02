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

#include "copyToContigCopy.h"

#if DOLLDLA

#include "contigCopy.h"
#include "copy.h"

CopyToContigCopy::CopyToContigCopy(Layer fromLayer, Layer toLayer) {
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
}

bool CopyToContigCopy::CanApply(const Node* node) const {
  if (node->GetNodeClass() == Copy::GetClass()) {
    const Copy* copy = static_cast<const Copy*>(node);
    return copy->InputIsContiguous(0) &&
      copy->InputIsContiguous(1);
  }
  throw;
}

void CopyToContigCopy::Apply(Node* node) const {
  Copy* copy = static_cast<Copy*>(node);
  auto contigCopy = new ContiguousCopy(m_toLayer);
  contigCopy->AddInputs(4,
			copy->Input(0), copy->InputConnNum(0),
			copy->Input(1), copy->InputConnNum(1));

  copy->m_poss->AddNode(contigCopy);
  copy->RedirectChildren(contigCopy, 0);
  copy->m_poss->DeleteChildAndCleanUp(copy);
}

#endif // DOLLDLA
