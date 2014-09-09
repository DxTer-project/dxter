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
#include "zaxpby.h"
#include "string.h"
#include "helperNodes.h"

using namespace std;


ZAxpBy::ZAxpBy(Layer layer, Coef alpha, Coef beta)
  : m_alpha(alpha), m_beta(beta)
{ 
  if (layer!=SMLAYER) {
    cout << "no refinements for ZAxpBy, so use SM\n";
    throw;
  }
  SetLayer(layer);
}


NodeType ZAxpBy::GetType() const
{
  return "ZAxpBy " + LayerNumToStr(GetLayer());
}

void ZAxpBy::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<3,1>::Duplicate(orig, shallow, possMerging);
  const ZAxpBy *axpy = (ZAxpBy*)orig;
  m_alpha = axpy->m_alpha;
  m_beta = axpy->m_beta;
}

void ZAxpBy::FlattenCore(ofstream &out) const
{
  DLAOp<3,1>::FlattenCore(out);
  WRITE(m_alpha);
  WRITE(m_beta);
}

void ZAxpBy::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<3,1>::UnflattenCore(in, info);
  READ(m_alpha);
  READ(m_beta);
}

Phase ZAxpBy::MaxPhase() const 
{
  if (m_layer == DMLAYER)
    return DPTENSORPHASE;
  else if (m_layer == SMLAYER)
    return NUMPHASES;
  else
    throw;
}
/*
bool ZAxpBy::ShouldCullDP() const 
{
#if DODPTENSORPHASE
  return m_layer == DMLAYER;
#else
  throw;
#endif
}
*/
bool ZAxpBy::DoNotCullDP() const 
{
#if DODPTENSORPHASE
  return m_layer == DMLAYER;
#else
  throw;
#endif
}

void ZAxpBy::PrintCode(IndStream &out)
{
  out.Indent();
  *out << "ZAxpBy( ";
  out << m_alpha;
  *out << ", " << GetInputName(0).str() << ", ";
  out << m_beta;
  *out << ", " << GetInputName(1).str() << ", "
       << GetInputName(2).str() << " );\n";
}

void ZAxpBy::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3,1>::Prop();

    if (m_layer == ABSLAYER || m_layer == DMLAYER) {
      m_cost = 3 * TotalNumberOfElements(0);
      Dim numDims = InputNumDims(0);
      if (InputNumDims(1) != numDims || InputNumDims(2) != numDims)
	throw;
      for (Dim dim = 0; dim < numDims; ++dim) {
	if (*InputLen(0,dim) != *InputLen(1,dim))
	  throw;
	if (*InputLen(0,dim) != *InputLen(2,dim))
	  throw;
      }
    }
    else if (m_layer == SMLAYER) {
      m_cost = 3 * TotalNumberOfLocalElements(0);
      Dim numDims = InputNumDims(0);
      if (InputNumDims(1) != numDims || InputNumDims(2) != numDims)
	throw;
      for (Dim dim = 0; dim < numDims; ++dim) {
	const DataTypeInfo &in0Type = InputDataType(0);
	const DataTypeInfo &in1Type = InputDataType(1);
	const DataTypeInfo &in2Type = InputDataType(2);
	if (in0Type.m_dist.m_dists[dim] != in1Type.m_dist.m_dists[dim])
	  throw;
	if (in0Type.m_dist.m_dists[dim] != in2Type.m_dist.m_dists[dim])
	  throw;
	
	if (*InputLen(0,dim) != *InputLen(1,dim))
	  throw;
	if (*InputLen(0,dim) != *InputLen(2,dim))
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
