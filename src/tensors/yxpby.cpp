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
#include "yxpby.h"
#include "string.h"
#include "helperNodes.h"

using namespace std;

YxpBy::YxpBy(Layer layer, Coef beta)
  : m_beta(beta)
{ 
  SetLayer(layer);
}

NodeType YxpBy::GetType() const
{
  return "YxpBy " + LayerNumToStr(GetLayer());
}

void YxpBy::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<2,1>::Duplicate(orig, shallow, possMerging);
  const YxpBy *axpy = (YxpBy*)orig;
  m_beta = axpy->m_beta;
}

void YxpBy::FlattenCore(ofstream &out) const
{
  DLAOp<2,1>::FlattenCore(out);
  WRITE(m_beta);
}

void YxpBy::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<2,1>::UnflattenCore(in, info);
  READ(m_beta);
}

Phase YxpBy::MaxPhase() const 
{
  if (m_layer == DMLAYER)
    throw;
  else if (m_layer == SMLAYER)
    return NUMPHASES;
  else
    throw;
}

void YxpBy::PrintCode(IndStream &out)
{
  out.Indent();
  *out << "YxpBy( ";
  *out << GetInputName(0).str() << ", ";
  out << m_beta;
  *out << ", " << GetInputName(1).str() << " );\n";
}

void YxpBy::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();

    if (m_layer == ABSLAYER)
      m_cost = 2 * TotalNumberOfElements(0);
    if (m_layer == DMLAYER)
      m_cost = 2 * TotalNumberOfElements(0);
    else if (m_layer == SMLAYER) {
      m_cost = 2 * TotalNumberOfLocalElements(0);
      if (InputDataType(0).m_dist != InputDataType(1).m_dist)
	throw;

      if (m_beta == COEFZERO)
	throw;
      if (m_beta == COEFNEGONE) {
	if (GetInputNameStr(0) == GetInputNameStr(1))
	  throw;
      }
    }
    else {
      cout << LayerNumToStr(m_layer) << endl;
      throw;
    }
  }
}

#endif // DOTENSORS
