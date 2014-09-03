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


#include "layers.h"
#if DOTENSORS
#include "base.h"
#include "xpay.h"
#include "string.h"
#include "helperNodes.h"

using namespace std;

Xpay::Xpay(Layer layer, Coef alpha)
  : m_alpha(alpha)
{ 
  SetLayer(layer);
}

NodeType Xpay::GetType() const
{
  return "Xpay " + LayerNumToStr(GetLayer());
}

void Xpay::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<2,1>::Duplicate(orig, shallow, possMerging);
  const Xpay *axpy = (Xpay*)orig;
  m_alpha = axpy->m_alpha;
}

void Xpay::FlattenCore(ofstream &out) const
{
  DLAOp<2,1>::FlattenCore(out);
  WRITE(m_alpha);
}

void Xpay::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<2,1>::UnflattenCore(in, info);
  READ(m_alpha);
}

Phase Xpay::MaxPhase() const 
{
  if (m_layer == DMLAYER)
    throw;
  else if (m_layer == SMLAYER)
    return NUMPHASES;
  else
    throw;
}

void Xpay::PrintCode(IndStream &out)
{
  out.Indent();
  *out << "Xpay( ";
  *out << GetInputName(0).str() << ", ";
  out << m_alpha;
  *out << ", " << GetInputName(1).str() << " );\n";
}

void Xpay::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();

    if (m_layer == ABSLAYER)
      throw;
    if (m_layer == DMLAYER)
      throw;
    else if (m_layer == SMLAYER) {
      m_cost = 3 * TotalNumberOfLocalElements(0);
      if (InputDataType(0).m_dist != InputDataType(1).m_dist)
	throw;
    }
    else {
      cout << LayerNumToStr(m_layer) << endl;
      throw;
    }
  }
}

#endif // DOTENSORS
