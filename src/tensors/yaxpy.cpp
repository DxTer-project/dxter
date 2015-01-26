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


#include "layers.h"
#if DOTENSORS
#include "base.h"
#include "yaxpy.h"
#include "string.h"
#include "helperNodes.h"

using namespace std;

Yaxpy::Yaxpy(Layer layer, Coef alpha)
  : m_alpha(alpha)
{ 
  SetLayer(layer);
}

void Yaxpy::AlignInfo(string &align,
                          DimVec &alignModes,
                          DimVec &alignModesSrc)
{
  align= GetInputNameStr(0);
  Dim numDims = InputNumDims(0);
  DimVec tmp;
  IdentDimVec(numDims, tmp);
  alignModes = tmp;
  alignModesSrc = tmp;
}

NodeType Yaxpy::GetType() const
{
  return "Yaxpy " + LayerNumToStr(GetLayer());
}

void Yaxpy::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<2,1>::Duplicate(orig, shallow, possMerging);
  const Yaxpy *axpy = (Yaxpy*)orig;
  m_alpha = axpy->m_alpha;
}

void Yaxpy::FlattenCore(ofstream &out) const
{
  DLAOp<2,1>::FlattenCore(out);
  WRITE(m_alpha);
}

void Yaxpy::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<2,1>::UnflattenCore(in, info);
  READ(m_alpha);
}

Phase Yaxpy::MaxPhase() const 
{
  if (m_layer == DM1LAYER || m_layer == DM2LAYER)
    throw;
  else if (m_layer == SMLAYER)
    return NUMPHASES;
  else
    throw;
}

void Yaxpy::PrintCode(IndStream &out)
{
  out.Indent();
  *out << "YAxpy( ";
  out << m_alpha;
  *out << ", " << GetInputName(0).str() << ", "
       << GetInputName(1).str() << " );\n";
}

void Yaxpy::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();

    if (m_layer == ABSLAYER)
      m_cost = 2 * TotalNumberOfElements(0);
    else if (m_layer == DM1LAYER || m_layer == DM2LAYER)
      m_cost = 2 * TotalNumberOfElements(0);
    else if (m_layer == SMLAYER) {
      m_cost = 2 * TotalNumberOfLocalElements(0);
      if (InputDataType(0).GetEffectiveDist() != InputDataType(1).GetEffectiveDist())
	throw;

      if (m_alpha == COEFZERO)
	throw;
      if (m_alpha == COEFNEGONE) {
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
