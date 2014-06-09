/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2013, The University of Texas and Bryan Marker

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



#include "DLANode.h"
#include "poss.h"
#include "elemRedist.h"
#include <sstream>

DLANode::DLANode()
{
  m_cost = -1;
  m_layer = BADLAYER;
}

DLANode::DLANode(Layer layer)
{
  m_cost = -1;
  m_layer = layer;
}

#if DOELEM
const DistType& DLANode::GetDistType(unsigned int num) const
{ 
  return MC_MR;
}
#elif DOTENSORS
const DistType& DLANode::GetDistType(unsigned int num) const
{ 
  throw;
  //  return DEFAULTDISTTYPE;
}
#endif

void DLANode::ClearBeforeProp()
{
  m_cost = -1;
}


void DLANode::SanityCheck()
{
  Node::SanityCheck();
  NodeConnVecIter iter1 = m_inputs.begin();
  for( ; iter1 != m_inputs.end(); ++iter1 ) {
    if (!(*iter1)->m_n->IsDLA()) {
      cout << "input not DLA" << endl;
      throw;
    }
  }
  
  NodeConnVecConstIter iter2 = m_children.begin();
  for( ; iter2 != m_children.end(); ++iter2 ) {
    if (!(*iter2)->m_n->IsDLA()) {
      cout << "output not DLA" << endl;
      throw;
    }
  }
}

#if DODM
const DistType& DLANode::InputDistType(unsigned int num) const
{
  DLANode *in = (DLANode*)Input(num);
  unsigned int inNum = InputConnNum(num);
  return in->GetDistType(inNum);
}
#endif

#if DOLLDLA
Stride DLANode::InputRowStride(unsigned int num) const
{
  const NodeConn *conn = m_inputs[num];
  DLANode *in = (DLANode*)(conn->m_n);
  unsigned int inNum = conn->m_num;
  return in->RowStride(inNum);
}

Stride DLANode::InputColStride(unsigned int num) const
{
  const NodeConn *conn = m_inputs[num];
  DLANode *in = (DLANode*)(conn->m_n);
  unsigned int inNum = conn->m_num;
  return in->ColStride(inNum);
}
#endif //DOLLDLA

string DLANode::GetCostStr()
{
  std::stringstream str;
  str << m_cost;
  return str.str();
}

bool DLANode::HasProped() const
{
  return IsValidCost(m_cost);
}

void DLANode::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode *origDLA = (DLANode*)orig;
  m_cost = origDLA->m_cost;
  m_layer = origDLA->m_layer;
  
  Node::Duplicate(orig, shallow, possMerging);
}

void DLANode::FlattenCore(ofstream &out) const
{
  WRITE(m_layer);
}

void DLANode::UnflattenCore(ifstream &in, SaveInfo &info)
{
  READ(m_layer);
}

#if TWOD
const Sizes* DLANode::GetInputM(unsigned int num) const
{
  if (num >= m_inputs.size()) {
    cout << "bad size 2\n";
    cout << num << " >= " << m_inputs.size() << endl;
    throw;
  }
  const DLANode *in = (DLANode*)Input(num);
  unsigned int inNum = InputConnNum(num);
  return in->GetM(inNum);
}

const Sizes* DLANode::GetInputN(unsigned int num) const 
{
  if (num >= m_inputs.size()) {
    cout << "bad size 4\n";
    throw;
  }
  const NodeConn *conn = m_inputs[num];
  const DLANode *in = (DLANode*)(conn->m_n);
  unsigned int inNum = conn->m_num;
  return in->GetN(inNum);
}

const Sizes* DLANode::InputLocalM(unsigned int num) const
{
  if (num >= m_inputs.size()) {
    cout << "bad size 2\n";
    cout << num << " >= " << m_inputs.size() << endl;
    throw;
  }
  DLANode *in = (DLANode*)Input(num);
  unsigned int inNum = InputConnNum(num);
  return in->LocalM(inNum);
}

const Sizes* DLANode::InputLocalN(unsigned int num) const 
{
  if (num >= m_inputs.size()) {
    cout << "bad size 4\n";
    throw;
  }
  DLANode *in = (DLANode*)Input(num);
  unsigned int inNum = InputConnNum(num);
  return in->LocalN(inNum);
}

#else
const Sizes* DLANode::InputLen(unsigned int num, Dim dim) const
{
  if (num >= m_inputs.size()) {
    cout << "bad size 2\n";
    cout << num << " >= " << m_inputs.size() << endl;
    throw;
  }
  const DLANode *in = (DLANode*)Input(num);
  unsigned int inNum = InputConnNum(num);
  return in->Len(inNum, dim);
}


const Sizes* DLANode::InputLocalLen(unsigned int num, Dim dim) const
{
  if (num >= m_inputs.size()) {
    cout << "bad size 2\n";
    cout << num << " >= " << m_inputs.size() << endl;
    throw;
  }
  DLANode *in = (DLANode*)Input(num);
  unsigned int inNum = InputConnNum(num);
  return in->LocalLen(inNum, dim);
}


const Dim DLANode::InputNumDims(unsigned int num) const
{
  if (num >= m_inputs.size()) {
    cout << "bad size 2\n";
    cout << num << " >= " << m_inputs.size() << endl;
    throw;
  }
  DLANode *in = (DLANode*)Input(num);
  unsigned int inNum = InputConnNum(num);
  return in->NumDims(inNum);
}


#endif

#if DOELEM
DLANode* DLANode::FindNonRedistParent(unsigned int num)
{
  unsigned int trash;
  return FindNonRedistParent(num, trash);
}

DLANode* DLANode::FindNonRedistParent(unsigned int num, unsigned int &parentNum)
{
  Node *parent = Input(num);
  parentNum = InputConnNum(num);
  while (parent) {
    DLANode *ddlaNode = (DLANode*)parent;
    if (ddlaNode->GetNodeClass() == RedistNode::GetClass()) {
      parent = ddlaNode->Input(0);
      parentNum = InputConnNum(0);
    }
    else
      return ddlaNode;
  }
  return NULL;
}
#endif

#if TWOD
bool DLANode::IsScalar(unsigned int num) const
{
  return GetM(num)->AllOnes() && GetN(num)->AllOnes();
}
#endif

DLANode* DLANode::FindSideEffectingUser(unsigned int num)
{
  NodeConnVecIter iter = m_children.begin();
  for(; iter != m_children.end(); ++iter) {
    if ((*iter)->m_num == num) {
      DLANode *child = (DLANode*)(*iter)->m_n;
#if DOELEM
      if (child->GetNodeClass() == RedistNode::GetClass()) {
        return child->FindSideEffectingUser(0);
      }
      else 
#elif DOTENSORS
	cout << "CheckThisSpot\n";
      throw;
      //	CheckThisSpot();
#endif
	{
        if (child->Overwrites(this, num))
          return child;
      }
    }
  }
  return NULL;
}

#if DOBLIS
void DLANode::UpdateInnerPackingMultiple(PackSize size)
{
  cout << GetNodeClass() << endl;
  cout << endl;
  for(unsigned int i = 0; i < m_inputs.size(); ++i)
    cout << i << " " << Input(i)->GetNodeClass() << endl;
  cout << endl;
  for(unsigned int i = 0; i < m_children.size(); ++i)
    cout << i << " " << Child(i)->GetNodeClass() << endl;

  throw;
}
#endif

#if DOELEM
void DLACullDP(Poss *poss, bool &cullIfPossible, bool &doNotCull)
{
  doNotCull = false;
  cullIfPossible = false;
  NodeVecIter iter = poss->m_possNodes.begin();
  for(; iter != poss->m_possNodes.end(); ++iter) {
    Node *node = *iter;
    if (node->IsDLA()) {
      DLANode *ddla = (DLANode*)node;
      if (ddla->DoNotCullDP()) {
        doNotCull = true;
        return;
      }
      if (!cullIfPossible) {
        cullIfPossible = ddla->ShouldCullDP();
      }
    }
  }
}

void DLACullRO(Poss *poss, bool &cullIfPossible, bool &doNotCull)
{
  cullIfPossible = false;
  doNotCull = false;
}
#endif

#if TWOD
void DLACullLA(Poss *poss, bool &cullIfPossible, bool &doNotCull)
{
  cullIfPossible = false;
  doNotCull = false;
}

#if DOSQM || DOSM
void DLACullSR(Poss *poss, bool &cullIfPossible, bool &doNotCull)
{
  doNotCull = false;
  cullIfPossible = false;
  NodeVecIter iter = poss->m_possNodes.begin();
  for(; iter != poss->m_possNodes.end(); ++iter) {
    Node *node = *iter;
    if (node->IsDLA()) {
      DLANode *ddla = (DLANode*)node;
      if (CurrPhase == SR1PHASE) {
	if (ddla->GetLayer() == S1LAYER) {
	  if (ddla->m_hasRefined) {
	    cullIfPossible = true;
	  }
	  else {
	    doNotCull = true;
	    return;
	  }
	}
      }
      if (CurrPhase == SR2PHASE) {
	if (ddla->GetLayer() == S2LAYER) {
	  if (ddla->m_hasRefined) {
	    cullIfPossible = true;
	  }
	  else {
	    doNotCull = true;
	    return;
	  }
	}
      }
      if (CurrPhase == SR3PHASE) {
	if (ddla->GetLayer() == S3LAYER) {
	  if (ddla->m_hasRefined) {
	    cullIfPossible = true;
	  }
	  else {
	    doNotCull = true;
	    return;
	  }
	}
      }
    }
  }
}
#endif
#endif


#if DOTENSORS
void TenCullDP(Poss *poss, bool &cullIfPossible, bool &doNotCull)
{
  doNotCull = false;
  cullIfPossible = false;
}

void TenCullRO(Poss *poss, bool &cullIfPossible, bool &doNotCull)
{
  cullIfPossible = false;
  doNotCull = false;
}
#endif

#if DOLLDLA
void LLDLACull(Poss *poss, bool &cullIfPossible, bool &doNotCull)
{
  doNotCull = false;
  cullIfPossible = false;
}
#endif
