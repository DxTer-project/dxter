/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2014, The University of Texas and Bryan Marker

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

#include "LLDLA.h"
#include "partition.h"

#if DOLLDLA

Partition::Partition(Layer layer, Dir partType, Size partStart, Size totalSize)
{
  if (totalSize <= partStart) {
    cout << "Error in Partition: totalSize <= partStart" << endl;
    throw;
  }
  m_layer = layer;
  m_partType = partType;
  m_startSizes = new Sizes(partStart);
  m_endSizes = new Sizes(totalSize - partStart);

  startName = GetInputName(0);
  endName = GetInputName(0);
  if (partType == HORIZONTAL) {
    endName.m_name += "_HORIZONTAL_PARTITION";
  } else {
    endName.m_name += "_VERTICAL_PARTITION";
  }
  return;
}

void Partition::Prop()
{
  if (!IsValidCost(m_cost)) {
  DLANode::Prop();
    m_cost = 0;
  }
  return;
}

Name Partition::GetName(ConnNum num) const
{
  if (num > 1) {
    cout << "num > 1 in Partition::GetName" << endl;
    throw;
  }

  if (num == 0) {
    return startName;
  } else {
    return endName;
  }
}

void Partition::AddVariables(VarSet &set) const
{
  string varDecl;
  if (m_partType == HORIZONTAL) {
    varDecl = arch->TypeName(GetDataType()) + " " +  GetInputNameStr(0) + "_HORIZONTAL_PARTITION;";
  } else {
    varDecl = arch->TypeName(GetDataType()) + " " +  GetInputNameStr(0) + "_VERTICAL_PARTITION;";
  }
  Var var(DirectVarDeclType, varDecl, GetDataType());
  set.insert(var);
}

const Sizes* Partition::GetM(ConnNum num) const {
  if (num > 1) {
    throw;
  }

  if (m_partType == HORIZONTAL) {
    return GetInputM(num);
  } else {
    if (num == 0) {
      return m_startSizes;
    } else {
      return m_endSizes;
    }
  }
  throw;
}

const Sizes* Partition::GetN(ConnNum num) const {
  if (num > 1) {
    throw;
  }

  if (m_partType == VERTICAL) {
    return GetInputN(num);
  } else {
    if (num == 0) {
      return m_startSizes;
    } else {
      return m_endSizes;
    }
  }
}

#endif // DOLLDLA
