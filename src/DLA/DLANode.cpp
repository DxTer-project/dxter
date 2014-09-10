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

void DLANode::Prop()
{
  //  Node::Prop();
  CheckConnections();
}

void DLANode::ClearBeforeProp()
{
  m_cost = -1;
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
const Sizes* DLANode::GetInputM(ConnNum num) const
{
  if (num >= m_inputs.size()) {
    cout << "bad size 2\n";
    cout << num << " >= " << m_inputs.size() << endl;
    throw;
  }
  const DLANode *in = (DLANode*)Input(num);
  ConnNum inNum = InputConnNum(num);
  return in->GetM(inNum);
}

const Sizes* DLANode::GetInputN(ConnNum num) const 
{
  if (num >= m_inputs.size()) {
    cout << "bad size 4\n";
    throw;
  }
  const NodeConn *conn = m_inputs[num];
  const DLANode *in = (DLANode*)(conn->m_n);
  ConnNum inNum = conn->m_num;
  return in->GetN(inNum);
}

#if DODM
const Sizes* DLANode::InputLocalM(ConnNum num) const
{
  if (num >= m_inputs.size()) {
    cout << "bad size 2\n";
    cout << num << " >= " << m_inputs.size() << endl;
    throw;
  }
  DLANode *in = (DLANode*)Input(num);
  ConnNum inNum = InputConnNum(num);
  return in->LocalM(inNum);
}

const Sizes* DLANode::InputLocalN(ConnNum num) const 
{
  if (num >= m_inputs.size()) {
    cout << "bad size 4\n";
    throw;
  }
  DLANode *in = (DLANode*)Input(num);
  ConnNum inNum = InputConnNum(num);
  return in->LocalN(inNum);
}
#endif //DODM

#else
const Sizes* DLANode::InputLen(ConnNum num, Dim dim) const
{
  if (num >= m_inputs.size()) {
    cout << "bad size 2\n";
    cout << num << " >= " << m_inputs.size() << endl;
    throw;
  }
  const DLANode *in = (DLANode*)Input(num);
  ConnNum inNum = InputConnNum(num);
  return in->Len(inNum, dim);
}


const Sizes* DLANode::InputLocalLen(ConnNum num, Dim dim) const
{
  if (num >= m_inputs.size()) {
    cout << "bad size 2\n";
    cout << num << " >= " << m_inputs.size() << endl;
    throw;
  }
  DLANode *in = (DLANode*)Input(num);
  ConnNum inNum = InputConnNum(num);
  return in->LocalLen(inNum, dim);
}


const Dim DLANode::InputNumDims(ConnNum num) const
{
  if (num >= m_inputs.size()) {
    cout << "bad size 2\n";
    cout << num << " >= " << m_inputs.size() << endl;
    throw;
  }
  DLANode *in = (DLANode*)Input(num);
  ConnNum inNum = InputConnNum(num);
  return in->NumDims(inNum);
}

Size DLANode::TotalNumberOfLocalElements(ConnNum num) const
{
  Dim numDims = NumDims(num);
  Size totSize = 0;
  const unsigned int totNumIters = LocalLen(num,0)->NumSizes();
  for (unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
    Size temp = 1;
    for (Dim dim = 0; dim < numDims; ++dim) {
      temp *= (*LocalLen(num,dim))[iteration];
    }
    totSize += temp;
  }
  return totSize;
}

Size DLANode::TotalNumberOfElements(ConnNum num) const
{
  Dim numDims = NumDims(num);
  Size totSize = 0;
  const unsigned int totNumIters = LocalLen(num,0)->NumSizes();
  for (unsigned int iteration = 0; iteration < totNumIters; ++iteration) {
    Size temp = 1;
    for (Dim dim = 0; dim < numDims; ++dim) {
      temp *= (*Len(num,dim))[iteration];
    }
    totSize += temp;
  }
  return totSize;
}


#endif

#if DOELEM
DLANode* DLANode::FindNonRedistParent(ConnNum num)
{
  ConnNum trash;
  return FindNonRedistParent(num, trash);
}

DLANode* DLANode::FindNonRedistParent(ConnNum num, ConnNum &parentNum)
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
bool DLANode::IsScalar(ConnNum num) const
{
  return GetM(num)->AllOnes() && GetN(num)->AllOnes();
}
#endif

DLANode* DLANode::FindSideEffectingUser(ConnNum num)
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
  for(ConnNum i = 0; i < m_inputs.size(); ++i)
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
	  if (ddla->HasRefined()) {
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
	  if (ddla->HasRefined()) {
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
	  if (ddla->HasRefined()) {
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


#if DODM

DataTypeInfo::DataTypeInfo()
{

}


DataTypeInfo::DataTypeInfo(DistType dist)
  :m_dist(dist)
{
}


DataTypeInfo::DataTypeInfo(const DataTypeInfo &rhs)
{
  m_dist = rhs.m_dist;
}

DataTypeInfo& DataTypeInfo::operator=(const DataTypeInfo &rhs)
{
  m_dist = rhs.m_dist;
  return *this;
}

#endif

#if TWOD
bool DLANode::IsInputRowVector(ConnNum num) const
{
  return *GetInputM(num) == 1;
}

bool DLANode::IsInputColVector(ConnNum num) const
{
  return *GetInputN(num) == 1;
}

bool DLANode::IsInputScalar(ConnNum num) const
{
  return *GetInputM(num) == 1 &&
    *GetInputN(num) == 1;
}

#if DOLLDLA
/*bool DLANode::InputIsMuByMu(ConnNum num) {
  return *GetInputM(num) == LLDLA_MU &&
    *GetInputN(num) == LLDLA_MU;
    }*/
#endif //DOLLDLA
#endif //TWOD
