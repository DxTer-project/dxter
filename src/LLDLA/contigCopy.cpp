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

#include "contigCopy.h"

#if DOLLDLA

void ContiguousCopy::PrintCode(IndStream& out) {
  
}

void ContiguousCopy::Prop() {
  if (!IsValidCost(m_cost)) {
    Copy::Prop();
    if (!InputDataType(0).IsContiguous() || !InputDataType(1).IsContiguous()) {
      cout << "ERROR: Contiguous copy on non-contiguous operands" << endl;
      throw;
    }
    m_cost = ZERO;
  }
}

#endif // DOLLDLA
