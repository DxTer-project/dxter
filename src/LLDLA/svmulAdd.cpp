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

#include "svmulAdd.h"

#if DOLLDLA

string SVecTypeToString(VecType vType) {
  if (vType == COLVECTOR) {
    return "COLVECTOR";
  } else {
    return "ROWVECTOR";
  }
}

void SVMulAdd::Prop() {
  if (!IsValidCost(m_cost)) {
    
    m_cost = 0;
  }
}

SVMulAdd::SVMulAdd(Layer layer, VecType vecType)
{
  m_layer = layer;
  m_vecType = vecType;
}

Node* SVMulAdd::BlankInst()
{
  return new SVMulAdd(LLDLAPRIMITIVELAYER, COLVECTOR);
}

Node* SVMulAdd::GetNewInst()
{
  return new SVMulAdd(m_layer, GetVecType());
}

void SVMulAdd::Duplicate(const Node *orig, bool shallow, bool possMerging) {
  DLAOp<3, 1>::Duplicate(orig, shallow, possMerging);
  const SVMulAdd* vadd = static_cast<const SVMulAdd*>(orig);
  m_layer = vadd->m_layer;
  m_vecType = vadd->GetVecType();
}

NodeType SVMulAdd::GetType() const
{
  return "SVMulAdd" + LayerNumToStr(GetLayer()) + SVecTypeToString(GetVecType());
}

VecType SVMulAdd::GetVecType() const {
  return m_vecType;
}

Phase SVMulAdd::MaxPhase() const
{
  switch (m_layer)
    { 
    case(ABSLAYER):
      return LLDLALOOPPHASE;
    case(LLDLAMIDLAYER):
      return LLDLAPRIMPHASE;
    case (LLDLAPRIMITIVELAYER):
      return NUMPHASES; 
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }
}


#endif // DOLLDLA
