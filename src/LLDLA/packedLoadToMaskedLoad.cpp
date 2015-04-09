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

#include "packedLoadToMaskedLoad.h"

#if DOLLDLA

#include "maskedLoad.h"
#include "regLoadStore.h"

bool PackedLoadToMaskedLoad::CanApply(const Node* node) const {
  if (node->GetNodeClass() == PackedLoadToRegs::GetClass()) {
    auto load = static_cast<const LoadToRegs*>(node);
    auto rs = load->InputDataType(0).m_rowStride;
    auto cs = load->InputDataType(0).m_colStride;
    return (load->IsInputColVector(0) && IsUnitStride(rs)) ||
      (load->IsInputRowVector(0) && IsUnitStride(cs));
  }
  throw;
}

void PackedLoadToMaskedLoad::Apply(Node* node) const {
  auto maskedLoad = new MaskedLoad();
  maskedLoad->AddInput(node->Input(0), node->InputConnNum(0));

  node->RedirectChildren(maskedLoad, 0);

  node->m_poss->AddNode(maskedLoad);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;
}

#endif // DOLLDLA
