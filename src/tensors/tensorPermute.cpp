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

#include "tensorPermute.h"

#if DOTENSORS

#include "loopTunnel.h"
#include "tempVarNode.h"
#include "helperNodes.h"
#include "tensorRedist.h"

Permute::Permute(string start, string end, Layer layer)
  : m_permutation(start,end),
    m_zero(false)
{
  if (start.empty())
    throw;
  SetLayer(layer);
}

Permute::Permute(const Permutation &permutation, Layer layer)
  : m_zero(false)
{
  if (!permutation.Size())
    throw;
  SetLayer(layer);
  m_permutation = permutation;
}

void Permute::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig,shallow, possMerging);
  const Permute *origNode = (Permute*)orig;
  m_info = origNode->m_info;
  m_permutation = origNode->m_permutation;
  m_zero = origNode->m_zero;
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
      cout << m_inputs.size() << endl;
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
	  if (m_permutation.Size() != InputNumDims(0))
	    throw;
	  m_cost = (PSIW + PSIR) * TotalNumberOfLocalElements(0);
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
  return InputLen(0,m_permutation.MapFinishToStart(dim));
}

const Sizes* Permute::LocalLen(ConnNum num, Dim dim) const
{
  if (num > 0)
    throw;
  if (dim >= m_permutation.Size())
    throw;
  return InputLocalLen(0,m_permutation.MapFinishToStart(dim));
}

void Permute::BuildDataTypeCache()
{
  Node *in = Input(0);
  m_info = in->DataType(InputConnNum(0));
  m_info.SetPerm(m_info.GetPerm().ComposeWith(m_permutation));
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
  name.m_type = m_info.GetDist();
  name.m_permutation = m_info.GetPerm();
  //  name.m_name += "_perm" + m_info.GetPerm().Str();
  return name;
}

void Permute::PrintCode(IndStream &out)
{  
  //Reflect in AddVars
  if (m_zero) {
    out.Indent();
    *out << "tempShape = " << GetInputNameStr(0) << ".Shape();\n";
    out.Indent();
    *out << GetNameStr(0) << ".ResizeTo( tempShape );\n";
    out.Indent();
    *out << "Scal( ";
    out << COEFZERO;
    *out << ", "
	 << GetNameStr(0) << " );\n";
  }
  else {
    out.Indent();
    *out << "Permute( " << GetInputNameStr(0) << ", " << GetNameStr(0) << " );\n";
  }
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

bool MovePermuteIntoTempVarNode::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Permute::GetClass())
    throw;
  const Permute *perm = (Permute*)node;
  if (perm->m_zero)
    return false;
  if (perm->Input(0)->GetNodeClass() == TempVarNode::GetClass())
    return perm->Input(0)->m_children.size() == 1;
  return false;
}


void MovePermuteIntoTempVarNode::Apply(Node *node) const
{
  Permute *perm = (Permute*)node;
  TempVarNode *temp = (TempVarNode*)(perm->Input(0));
  temp->m_info.SetPerm(temp->m_info.GetPerm().ComposeWith(perm->m_permutation));
  perm->RedirectAllChildren(temp);
  perm->m_poss->DeleteChildAndCleanUp(perm);
}

bool MovePermuteIntoRedist::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Permute::GetClass())
    throw;
  const Permute *perm = (Permute*)node;
  if (perm->m_zero)
    return false;
  if (perm->Input(0)->GetNodeClass() == RedistNode::GetClass())
    return perm->Input(0)->m_children.size() == 1;
  return false;
}


void MovePermuteIntoRedist::Apply(Node *node) const
{
  Permute *perm = (Permute*)node;
  RedistNode *temp = (RedistNode*)(perm->Input(0));

  RedistNode *redist = new RedistNode(temp->m_info.GetDist(), temp->m_info.GetPerm().ComposeWith(perm->m_permutation),
				      temp->m_align, temp->m_alignModes, temp->m_alignModesSrc);
  
  redist->AddInput(temp->Input(0), temp->InputConnNum(0));

  perm->m_poss->AddNode(redist);
  perm->RedirectAllChildren(redist);
  perm->m_poss->DeleteChildAndCleanUp(perm);
}


bool PermuteLoopHoist::CanApply(const Node *node) const
{
  if (CurrPhase != FINALOPTPHASE)
    return false;
  if (node->GetNodeClass() != LoopTunnel::GetClass())
    throw;
  if (!node->IsTunnel(SETTUNIN))
    return false;
  bool foundOne = false;
  Permutation perm;
  const LoopTunnel *setTun = (LoopTunnel*)node;
  //should have inlined
  if (setTun->m_pset->IsShadow())
    throw;
  if (setTun->m_children.size() == 0)
    throw;
  NodeConnVecConstIter setIter = setTun->m_children.begin();
  for(; setIter != setTun->m_children.end(); ++setIter) {
    const LoopTunnel *possTun = (LoopTunnel*)((*setIter)->m_n);
    bool foundHere = false;
    bool foundSomethingElse = false;
    NodeConnVecConstIter possIter = possTun->m_children.begin();
    for(; possIter != possTun->m_children.end(); ++possIter) {
      const Node *child = (*possIter)->m_n;
      if (child->GetNodeClass() == Permute::GetClass()) {
	const Permute *thisPerm = (Permute*)child;
	if (foundHere || foundOne) {
	  if (perm != thisPerm->m_permutation)
	    throw;
	}
	else 
	  perm = thisPerm->m_permutation;
	foundHere = true;
      }
      else if (!child->IsTunnel(POSSTUNOUT))
	foundSomethingElse = true;
    }
    if (foundHere && foundSomethingElse) {
      throw;
    }
    foundOne |= foundHere;
  }
  return foundOne;
}

void PermuteLoopHoist::Apply(Node *node) const
{
  if (node->GetNodeClass() != LoopTunnel::GetClass())
    throw;
  if (!node->IsTunnel(SETTUNIN))
    throw;
  LoopTunnel *setTun = (LoopTunnel*)node;
  if (setTun->m_pset->GetPosses().size() > 1)
    throw;
  LoopTunnel *possTun = (LoopTunnel*)(setTun->Child(0));
  Permute *perm = NULL;
  NodeConnVecConstIter possIter = possTun->m_children.begin();
  for(; possIter != possTun->m_children.end(); ++possIter) {
    const Node *child = (*possIter)->m_n;
    if (child->GetNodeClass() == Permute::GetClass()) {
      perm = (Permute*)child;
      break;
    }
    else if (!child->IsTunnel(POSSTUNOUT))
      throw;
  }
  if (!perm)
    throw;
  Permute *newPermute = new Permute(perm->m_permutation, perm->GetLayer());
  newPermute->AddInput(setTun->Input(0), setTun->InputConnNum(0));
  setTun->ChangeInput2Way(setTun->Input(0), setTun->InputConnNum(0), newPermute, 0);
  setTun->m_poss->AddNode(newPermute);
  newPermute->BuildDataTypeCache();

  setTun->ClearDataTypeCache();

  perm->RedirectChildren(possTun, 0);
  perm->m_poss->DeleteChildAndCleanUp(perm);


  LoopTunnel *possTunOut = possTun->GetMatchingOutTun();
  if (possTunOut->Input(0)->GetNodeClass() == Permute::GetClass()) {
    Permute *outPermute = (Permute*)(possTunOut->Input(0));
    Permute *newOutPermute = new Permute(outPermute->m_permutation, outPermute->GetLayer());
    LoopTunnel *setTunOut = (LoopTunnel*)(possTunOut->Child(0));
    if (setTunOut->m_children.size() == 0)
      throw;
    setTunOut->RedirectChildren(newOutPermute, 0);
    newOutPermute->AddInput(setTunOut, 0);
    //sanity check
    if (setTunOut->m_poss != setTun->m_poss)
      throw;
    setTunOut->m_poss->AddNode(newOutPermute);
    outPermute->RedirectChildren(outPermute->Input(0), outPermute->InputConnNum(0));
    outPermute->m_poss->DeleteChildAndCleanUp(outPermute);
    newOutPermute->BuildDataTypeCache();
  }
}

bool CombinePermutations::CanApply(const Node *node) const
{
  const Permute *perm = (Permute*)node;
  NodeConnVecConstIter iter = perm->m_children.begin();
  for(; iter != perm->m_children.end(); ++iter) {
    const Node *child = (*iter)->m_n;
    if (child->GetNodeClass() == Permute::GetClass())
      return true;
  }
  return false;
}

void CombinePermutations::Apply(Node *node) const
{
  Permute *perm = (Permute*)node;
  for(int i = 0; i < (int)perm->m_children.size(); ++i) {
    Node *tmp = perm->Child(i);
    if (tmp->GetNodeClass() == Permute::GetClass()) {
      Permute *child = (Permute*)tmp;
      Permutation composition = perm->m_permutation.ComposeWith(child->m_permutation);
      if (composition.HasPerm()) {
	Permute *newPerm = new Permute(composition, perm->GetLayer());
	newPerm->AddInput(perm->Input(0), perm->InputConnNum(0));
	child->RedirectChildren(newPerm, 0);
	child->m_poss->AddNode(newPerm);
	newPerm->BuildDataTypeCache();
      }
      else {
	child->RedirectChildren(perm->Input(0), perm->InputConnNum(0));
      }
      bool needToBreak = false;
      if (perm->m_children.size() == 1)
	needToBreak = true;
      child->m_poss->DeleteChildAndCleanUp(child);
      if (needToBreak)
	break;
      else
	--i;
    }
  }
  perm->m_poss->BuildDataTypeCache();
}

bool CombineScaleAndPermutation::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != ScaleNode::GetClass())
    throw;
  const ScaleNode *scale = (ScaleNode*)node;
  if (scale->m_val != COEFZERO)
    return false;
  if (scale->m_children.size() != 1)
    return false;
  return scale->Child(0)->GetNodeClass() == Permute::GetClass();
}

void CombineScaleAndPermutation::Apply(Node *node) const
{
  ScaleNode *scale = (ScaleNode*)node;
  Permute *perm = (Permute*)(node->Child(0));
  Permute *newPerm = new Permute(perm->m_permutation, perm->GetLayer());
  newPerm->m_info = perm->m_info;
  newPerm->m_zero = true;
  newPerm->AddInput(scale->Input(0), scale->InputConnNum(0));
  perm->RedirectChildren(newPerm);
  perm->m_poss->AddNode(newPerm);
  perm->m_poss->DeleteChildAndCleanUp(perm);
}

#endif //DOTENSORS
