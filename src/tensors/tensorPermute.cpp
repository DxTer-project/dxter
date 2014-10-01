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

#include "tensorPermute.h"

#if DOTENSORS

Permute::Permute(string start, string end, Layer layer)
 : m_permutation(start,end)
{
  SetLayer(layer);
}

Permute::Permute(const Permutation &permutation, Layer layer)
{
  SetLayer(layer);
  m_permutation = permutation;
}

void Permute::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig,shallow, possMerging);
  const Permute *origNode = (Permute*)orig;
  m_info = origNode->m_info;
  m_permutation = origNode->m_permutation;
}

NodeType Permute::GetType() const 
{
  return LayerNumToStr(GetLayer()) + "perm" + m_permutation.Str();
}

void Permute::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLANode::Prop();

    if (m_inputs.size() != 1) {
      cout << "m_inputs.size() != 1\n";
      throw;
    }
    
    if (m_children.empty())
      throw;
 
    Input(0)->Prop();
    
    switch (GetLayer()) 
      {
      case (ABSLAYER):
	{
	  m_cost = ZERO;
	  break;
	}
      case (SMLAYER):
	{
	  m_cost = 1;
	  if (InputLocalLen(0,0)->NumSizes() != 1)
	    throw;
	  const unsigned int numDims = InputNumDims(0);
	  if (m_permutation.Size() != numDims)
	    throw;
	  for (Dim dim = 0; dim < numDims; ++dim)
	    m_cost *= (*InputLocalLen(0,dim))[0];
	  m_cost *= PSIW + PSIR;
	  break;
	}
      default:
	throw;
    }
  }
}

Phase Permute::MaxPhase() const
{
  switch (GetLayer()) 
    {
    case (ABSLAYER):
      return DPTENSORPHASE;
      break;
    case (SMLAYER):
      return NUMPHASES;
      break;
    default:
      throw;
    }
}

const Dim Permute::NumDims(ConnNum num) const
{
  if (num > 0)
    throw;
  return InputNumDims(0);
}

const Sizes* Permute::Len(ConnNum num, Dim dim) const
{
  if (num > 0)
    throw;
  if (dim >= m_permutation.Size())
    throw;
  return InputLen(0,m_permutation.Map(dim));
}

const Sizes* Permute::LocalLen(ConnNum num, Dim dim) const
{
  if (num > 0)
    throw;
  if (dim >= m_permutation.Size())
    throw;
  return InputLocalLen(0,m_permutation.Map(dim));
}



void Permute::BuildDataTypeCache()
{
  const DataTypeInfo &inputType = InputDataType(0);
  if (m_permutation.Size() != inputType.m_dist.m_numDims)
    throw;
  m_info.m_dist = inputType.m_dist;
  for(Dim dim = 0; dim < inputType.m_dist.m_numDims; ++dim) {
    m_info.m_dist.m_dists[dim] = inputType.m_dist.m_dists[m_permutation.Map(dim)];
  }
  if (DistTypeEqual(m_info.m_dist, inputType.m_dist))
    throw;
}

void Permute::FlattenCore(ofstream &out) const
{
  throw;
  //  out << m_permutation << endl;
}
   
void Permute::UnflattenCore(ifstream &in, SaveInfo &info)
{
  throw;
  //  getline(in, m_permutation);
}

Name Permute::GetName(ConnNum num) const
{
  if (num > 0)
    throw;
  Name name = GetInputName(0);
  name.m_type = m_info.m_dist;
  name.m_name += "_p";
  return name;
}

void Permute::PrintCode(IndStream &out)
{  
  //Reflect in AddVars
  out.Indent();
  throw;

}

void Permute::AddVariables(VarSet &set) const
{
  DLANode::AddVariables(set);
  Var var(PermutationVarType, m_permutation.m_permutation);
  set.insert(var);  
}

bool LowerPermute::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Permute::GetClass())
    throw;
  const Permute *perm = (Permute*)node;
  return (perm->GetLayer() == ABSLAYER);
}


void LowerPermute::Apply(Node *node) const
{
  Permute *perm = (Permute*)node;
  perm->SetLayer(SMLAYER);
}

#endif //DOTENSORS
