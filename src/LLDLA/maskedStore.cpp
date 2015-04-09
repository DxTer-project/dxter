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

#include "maskedStore.h"

#if DOLLDLA

#include "costModel.h"

void MaskedStore::SanityCheckInputDimensions() {
  if (!(IsInputRowVector(1) || IsInputColVector(1)) || !InputIsContiguous(1)) {
    throw;
  }
  if (GetInputNumRows(1) >= Input(1)->GetVecRegWidth() ||
      GetInputNumCols(1) >= Input(1)->GetVecRegWidth()) {
    throw;
  }
}

void MaskedStore::Prop() {
  if (!IsValidCost(m_cost)) {
    SanityCheckInputDimensions();
    m_cost = costModel->ContigVecLoadCost();
  }
}

void MaskedStore::PrintCode(IndStream& out) {
  string buffer = GetInputNameStr(1);
  string regName = GetInputNameStr(0);
  string assignOpName, maskName;

  unsigned int residualSize;
  if (GetInputNumRows(0) == 1) {
    residualSize = GetInputNumCols(0);
  } else {
    residualSize = GetInputNumRows(0);
  }

  string residualInd = std::to_string((long long int) (residualSize - 1));

  if (GetDataType() == REAL_SINGLE) {
    assignOpName = "_mm256_maskstore_ps";
    maskName = "lsm0" + residualInd;
  } else if (GetDataType() == REAL_DOUBLE) {
    assignOpName = "_mm256_maskstore_pd";
    maskName = "llsm0" + residualInd;
  } else {
    throw;
  }
  out.Indent();
  *out << assignOpName + "( " + buffer + " , " + maskName + ".v , " + regName + ".v );\n" << endl;
}

#endif // DOLLDLA
