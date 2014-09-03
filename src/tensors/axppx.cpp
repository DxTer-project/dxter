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
#include "axppx.h"
#include "string.h"
#include "helperNodes.h"

using namespace std;

Axppx::Axppx(Layer layer, Coef alpha, Coef beta, string start, string end)
  : m_alpha(alpha), m_beta(beta)
{ 
  SetLayer(layer);
  if (start.length() != end.length())
    throw;
  string::iterator iter = end.begin();
  for(; iter != end.end(); ++iter) {
    m_permutation.push_back(start.find(*iter));
  }
}

Axppx::Axppx(Layer layer, Coef alpha, Coef beta)
  : m_alpha(alpha), m_beta(beta)
{ 
  SetLayer(layer);
}

Axppx::Axppx(Layer layer, Coef alpha, Coef beta, const DimVec &perm)
  : m_alpha(alpha), m_beta(beta)
{ 
  SetLayer(layer);
  m_permutation = perm;
}

NodeType Axppx::GetType() const
{
  return "Axppx " + LayerNumToStr(GetLayer());
}

void Axppx::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<3,1>::Duplicate(orig, shallow, possMerging);
  const Axppx *axpy = (Axppx*)orig;
  m_alpha = axpy->m_alpha;
  m_beta = axpy->m_beta;
  m_permutation = axpy->m_permutation;
}

void Axppx::FlattenCore(ofstream &out) const
{
  DLAOp<3,1>::FlattenCore(out);
  WRITE(m_alpha);
  WRITE(m_alpha);
  throw;
  //m_permutation
}

void Axppx::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<3,1>::UnflattenCore(in, info);
  throw;
  //m_permutation
}

Phase Axppx::MaxPhase() const 
{
  if (m_layer == DMLAYER)
    return DPTENSORPHASE;
  else if (m_layer == SMLAYER)
    return NUMPHASES;
  else
    throw;
}
/*
bool Axppx::ShouldCullDP() const 
{
#if DODPTENSORPHASE
  return m_layer == DMLAYER;
#else
  throw;
#endif
}
*/
bool Axppx::DoNotCullDP() const 
{
#if DODPTENSORPHASE
  return m_layer == DMLAYER;
#else
  throw;
#endif
}

void Axppx::PrintCode(IndStream &out)
{
  out.Indent();
  *out << "Axppx( ";
  out << m_alpha;
  *out << ", " << GetInputName(0).str() << ", ";
  out << m_beta;
  *out << ", " << PermutationVarName(m_permutation)
       << ", " << GetInputName(1).str() << ", "
       << GetInputName(2).str() << " );\n";
}

void Axppx::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3,1>::Prop();

    if (m_layer == ABSLAYER)
      m_cost = 0;
    if (m_layer == DMLAYER)
      m_cost = 0;
    else if (m_layer == SMLAYER) {
      m_cost = 3 * TotalNumberOfLocalElements(0);
      /*
      DistType tmp;
      tmp.SetToDefault(InputDataType(2).m_dist.m_numDims);
      if (tmp != DataType(0).m_dist)
	throw;
      */
    }
    else {
      cout << LayerNumToStr(m_layer) << endl;
      throw;
    }
  }
}



bool DistAxppxToDefaultLocalAxppx::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Axppx::GetClass())
    return false;
  if (((Axppx*)node)->GetLayer() != DMLAYER)
    return false;
  return true;
}

void DistAxppxToDefaultLocalAxppx::Apply(Node *node) const
{
  Axppx *orig = (Axppx*)node;
  Axppx *newAxppx = new Axppx(SMLAYER, orig->m_alpha, orig->m_beta, orig->m_permutation);

  newAxppx->AddInput(node->Input(0),node->InputConnNum(0));

  if (!orig->m_permutation.empty()) {
    const DataTypeInfo &inputType = orig->InputDataType(0);
    DistType newType = inputType.m_dist;
    for(Dim dim = 0; dim < inputType.m_dist.m_numDims; ++dim) {
      newType.m_dists[dim] = inputType.m_dist.m_dists[orig->m_permutation[dim]];
    }
    RedistNode *redist = new RedistNode(newType);
    redist->AddInput(node->Input(1),node->InputConnNum(1));
    newAxppx->AddInput(redist, 0);
    node->m_poss->AddNode(redist);
  }
  else
    newAxppx->AddInput(node->Input(1),node->InputConnNum(1));
  
  newAxppx->AddInput(node->Input(2),node->InputConnNum(2));

  node->m_poss->AddNode(newAxppx);
  
  node->RedirectChildren(newAxppx,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

#endif // DOTENSORS
