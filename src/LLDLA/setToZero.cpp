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

#include "LLDLA.h"
#include "setToZero.h"

#if DOLLDLA

SetToZero::SetToZero(Layer layer) {
  m_layer = layer;
}

NodeType SetToZero::GetType() const {
  return "SetToZero " + LayerNumToStr(GetLayer());
}

Node* SetToZero::BlankInst() {
  return new SetToZero(ABSLAYER);
}

Node* SetToZero::GetNewInst() {
  return BlankInst();
}

Phase SetToZero::MaxPhase() const {
  switch (m_layer)
    { 
    case(ABSLAYER):
      return LLDLALOOPPHASE;
    case(LLDLAMIDLAYER):
      return LLDLAPRIMPHASE;
    case (LLDLAPRIMITIVELAYER):
      return NUMPHASES; 
    default:
      throw;
    }
}

void SetToZero::Prop() {
  if (!IsValidCost(m_cost)) {
    DLAOp<1, 1>::Prop();
    
    if (!InputDataType(0).IsContiguous()) {
      cout << "ERROR: Currently SetToZero does not support general stride operands" << endl;
      cout << "Input 0 row stride = " << InputDataType(0).m_rowStrideVal << endl;
      cout << "Input 0 col stride = " << InputDataType(0).m_colStrideVal << endl;
      throw;
    }

    m_cost = GetInputM(0)->SumProds11(*GetInputN(0));
  }
}

void SetToZero::NaiveZeroPrintout(IndStream& out) {
  string inputName = GetInputName(0).m_name;
  string inputM = InputDataType(0).m_numRowsVar;
  string inputN = InputDataType(0).m_numColsVar;
  string rowStride = InputDataType(0).m_rowStrideVar;
  string colStride = InputDataType(0).m_colStrideVar;

  out.Indent();
  if (Input(0)->GetDataType() == REAL_SINGLE) {
    *out << "set_to_zero_float( ";
  } else {
    *out << "set_to_zero_double( ";
  }

  *out << inputName << ", " <<
    inputM << ", " <<
    inputN << ", " <<
    rowStride << ", " <<
    colStride << " );\n" << endl;
}

void SetToZero::MemsetPrintout(IndStream& out) {
  out.Indent();
  string inputName = GetInputName(0).m_name;
  string inputSize = "(" + InputDataType(0).m_numRowsVar + " * " + InputDataType(0).m_numColsVar + ")";
  if (Input(0)->GetDataType() == REAL_SINGLE) {
    inputSize += "* sizeof(float)";
  } else {
    inputSize += "* sizeof(double)";
  }
  *out << "memset( " << inputName << ", 0, " << inputSize << " );\n" << endl;
}

void SetToZero::PrintCode(IndStream& out) {
  if (m_layer == ABSLAYER) {
    NaiveZeroPrintout(out);
  } else {
    MemsetPrintout(out);
  }
}

ClassType SetToZero::GetNodeClass() const {
  return GetClass();
}

ClassType SetToZero::GetClass() {
  return "SetToZero";
}

Name SetToZero::GetName(ConnNum num) const {
  if (num != 0) {
    cout << "ERROR: Bad ConnNum in SetToZero::GetName" << endl;
    throw;
  }
  return GetInputName(0);
}

bool SetToZero::IsReadOnly() const {
  return false;
}

bool SetToZero::Overwrites(const Node* input, ConnNum num) const {
  return true;
}

bool SetToZero::IsDataDependencyOfInput() const {
  return true;
}

SetToZeroLowerLayer::SetToZeroLowerLayer(Layer fromLayer, Layer toLayer) {
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
}

string SetToZeroLowerLayer::GetType() const {
  return "SetToZeroLowerLayer from " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer);
}

bool SetToZeroLowerLayer::CanApply(const Node* node) const {
  if (node->GetNodeClass() == SetToZero::GetClass()) {
    const SetToZero* stz = (SetToZero*) node;
    if (stz->GetLayer() != m_fromLayer) {
      return false;
    }
    return true;
  }
  throw;
}

void SetToZeroLowerLayer::Apply(Node* node) const {
  SetToZero* stz = (SetToZero*) node;
  stz->SetLayer(m_toLayer);
}

bool SetToZeroLowerLayer::IsRef() const {
  return true;
}

#endif // DOLLDLA
