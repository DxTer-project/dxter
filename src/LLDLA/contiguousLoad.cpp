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

#include "contiguousLoad.h"

#if DOLLDLA

#include "costModel.h"

void ContiguousLoad::SanityCheckInputDimensions() {
  if (!(IsInputRowVector(0) || IsInputColVector(0))) {
    throw;
  }
}

void ContiguousLoad::Prop() {
  if (!IsValidCost(m_cost)) {
    DLANode::Prop();
    SanityCheckInputDimensions();
    m_cost = costModel->ContigVecLoadCost();
  }
}

void ContiguousLoad::PrintCode(IndStream& out) {
  string toLoadName = GetInputNameStr(0);
  string loadStr = GetNameStr(0);

  out.Indent();
  *out << arch->ContiguousLoad(GetDataType(), toLoadName, loadStr);
}

#endif // DOLLDLA
