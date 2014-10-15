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

#include "LLDLATranspose.h"
#include "DLAOp.h"
#include "loopSupport.h"
#include "LLDLA.h"

#if DOLLDLA

LLDLATranspose::LLDLATranspose(Layer layer)
{
  m_layer = layer;
}

void LLDLATranspose::PrintCode(IndStream &out)
{
  return;
}

Phase LLDLATranspose::MaxPhase() const
{
  return NUMPHASES;
  /*  switch (m_layer)
    { 
    case(ABSLAYER):
      return LLDLALOOPPHASE;
    case(LLDLAMIDLAYER):
      return LLDLAPRIMPHASE;
    case (LLDLAPRIMITIVELAYER):
      return NUMPHASES; 
    default:
      throw;
      }*/
}

NodeType LLDLATranspose::GetType() const
{
  return "LLDLATranspose";
}

Node* LLDLATranspose::BlankInst()
{
  return new LLDLATranspose(ABSLAYER);
}

void LLDLATranspose::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<1,1>::Duplicate(orig, shallow, possMerging);
  const LLDLATranspose *rhs = (LLDLATranspose*)orig;
  m_layer = rhs->m_layer;
  return;
}

void LLDLATranspose::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<1, 1>::Prop();

    m_cost = ZERO;
  }
  return;
}

LLDLATransposeLowerLayer::LLDLATransposeLowerLayer(Layer toLayer, Layer fromLayer)
{
  m_toLayer = toLayer;
  m_fromLayer = fromLayer;
}

bool LLDLATransposeLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == LLDLATranspose::GetClass()) {
    const LLDLATranspose *trans = (LLDLATranspose*) node;
    if (trans->GetLayer() != m_fromLayer) {
      return false;
    }
    return true;
  }
  cout << "Error: Applying LLDLATransposeLowerLayer to non LLDLATranspose node\n";
  throw;
}

void LLDLATransposeLowerLayer::Apply(Node *node) const
{
  LLDLATranspose* trans = (LLDLATranspose*) node;
  trans->SetLayer(m_toLayer);
}

string LLDLATransposeLowerLayer::GetType() const
{
  return "LLDLATranspose lower layer " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer);
}

const DataTypeInfo& LLDLATranspose::DataType(ConnNum num) const
{
  return m_info;
}

void LLDLATranspose::ClearDataTypeCache()
{
  m_info.m_rowStride = BADSTRIDE;
}

void LLDLATranspose::BuildDataTypeCache()
{
  const DataTypeInfo &in = InputDataType(0);

  m_info.m_rowStride = in.m_colStride;
  m_info.m_colStride = in.m_rowStride;
  m_info.m_numRowsVar = in.m_numColsVar;
  m_info.m_numColsVar = in.m_numRowsVar;
  m_info.m_rowStrideVar = in.m_colStrideVar;
  m_info.m_colStrideVar = in.m_rowStrideVar;
}

const Sizes* LLDLATranspose::GetN(ConnNum num) const
{
  if (num > 0) {
    throw;
  }
  return GetInputM(0);
}

const Sizes* LLDLATranspose::GetM(ConnNum num) const
{
  if (num > 0) {
    throw;
  }
  return GetInputN(0);
}

#endif // DOLLDLA
