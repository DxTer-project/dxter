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

#include "recombine.h"

#if DOLLDLA

Recombine::Recombine(Layer layer, Dir partType)
{
  m_layer = layer;
  m_partType = partType;
}

void Recombine::SanityCheckDimensions() {
  if ((m_partType == VERTICAL) &&
      (GetInputNumRows(0) + GetInputNumRows(1) != GetInputNumRows(2))) {
    LOG_FAIL("Bad dimension in VERTICAL Recombine::SanityCheckDimensions");
    throw;
  }
  if ((m_partType == HORIZONTAL) &&
      (GetInputNumCols(0) + GetInputNumCols(1) != GetInputNumCols(2))) {
    LOG_FAIL("Bad dimension in HORIZONTAL Recombine::SanityCheckDimensions");
    throw;
  }

}

void Recombine::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3, 1>::Prop();
    SanityCheckDimensions();
    m_cost = 0;
  }

  return;
}

void Recombine::PrintCode(IndStream &out)
{
  return;
}

Name Recombine::GetName(ConnNum num) const
{
  if (num != 0) {
    cout << "Error: Invalid connection number for Recombine::GetName" << endl;
    LOG_FAIL("replacement for throw call");
  }
  return GetInputName(2);
}

const Sizes* Recombine::GetM(ConnNum num) const {
  if (num != 0) {
    cout << "Error: Invalid connection number for Recombine::GetM" << endl;
    LOG_FAIL("replacement for throw call");
  }
  return GetInputM(2);
}

const Sizes* Recombine::GetN(ConnNum num) const {
  if (num != 0) {
    cout << "Error: Invalid connection number for Recombine::GetN" << endl;
    LOG_FAIL("replacement for throw call");
  }
  return GetInputN(2);
}

const Type Recombine::GetDataType() const {
  return Input(2)->GetDataType();
}

void Recombine::Duplicate(const Node* orig, bool shallow, bool possMerging)
{
  DLAOp<3, 1>::Duplicate(orig, shallow, possMerging);
  const Recombine* rec = static_cast<const Recombine*>(orig);
  m_layer = rec->m_layer;
  m_partType = rec->m_partType;
}

const DataTypeInfo& Recombine::DataType(ConnNum num) const
{
  if (num != 0) {
    cout << "Error: num > 0 in Recombine::DataType\n";
    LOG_FAIL("replacement for throw call");
  }
  return InputDataType(2);
}

RecombineLowerLayer::RecombineLowerLayer(Layer fromLayer, Layer toLayer)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
}

bool RecombineLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Recombine::GetClass()) {
    const Recombine *recombine = static_cast<const Recombine*>(node);
    if (recombine->GetLayer() != m_fromLayer) {
      return false;
    }
    return true;
  }
  cout << "Error: Applying RecombineLowerLayer to non-Recombine node\n";
  cout << "Node class is: " << node->GetClass() << endl;
  cout << "GetNodeClass is: " << node->GetNodeClass() << endl;
  LOG_FAIL("replacement for throw call");
}

void RecombineLowerLayer::Apply(Node *node) const
{
  Recombine *recombine = static_cast<Recombine*>(node);
  recombine->SetLayer(m_toLayer);
}

string RecombineLowerLayer::GetType() const
{
  return "Recombine lower layer " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer);
}

#endif // DOLLDLA
