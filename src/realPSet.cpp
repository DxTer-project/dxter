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



#include "base.h"
#include "transform.h"
#include "realPSet.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "loopSupport.h"

extern unsigned int M_phase;

#define FORMSETSEARLY 0


RealPSet::RealPSet()
  : m_functionality(), m_mergeLeft(NULL), m_mergeRight(NULL)
{
}

RealPSet::RealPSet(Poss *poss)
: m_mergeLeft(NULL), m_mergeRight(NULL)
{
  Init(poss);
}

void RealPSet::Init(Poss *poss)
{
  m_functionality = poss->GetFunctionalityString();

  if (m_functionality.empty()) {
    cout << "starting PSet without functionality\n";
    throw;
  }
  if (IsLoop()) {
    //    Loop *loop = (Loop*)this;
    m_functionality += (char)((dynamic_cast<const LoopInterface*>(this))->GetBSSize().GetSize());
  }
  /*
  //Make single tunnels with multiple inputs/outputs into individual tunnels
  //Poss mergin with multiple intput/output tunnels is very buggy
  poss->ExpandTunnels();
  */
  
  //Go through the input tunnels of the poss, create a set tunnel,
  // change the inputs from connecting to poss tunnels to set tunnels
  for(unsigned int i = 0; i < poss->m_inTuns.size(); ++i) {
    Tunnel *possTun = (Tunnel*)(poss->InTun(i));
    if (!possTun->IsTunnel(POSSTUNIN)) {
      cout << "bad poss tunnel\n";
      throw;
    }
    Tunnel *setTun = possTun->GetSetTunnel();
    for(ConnNum j = 0; j < possTun->m_inputs.size(); ++j) {
      NodeConn *conn = possTun->InputConn(j);
      setTun->AddInput(conn->m_n, conn->m_num);
    }
    while(possTun->m_inputs.size()) {
      NodeConn *conn = possTun->InputConn(0);
      conn->m_n->RemoveChild(possTun, conn->m_num);
      delete conn;
      possTun->m_inputs.erase(possTun->m_inputs.begin());
    }
    m_inTuns.push_back(setTun);
    setTun->m_pset = this;
  }
  
  for(unsigned int i = 0; i < poss->m_outTuns.size(); ++i) {
    Tunnel *possTun = (Tunnel*)(poss->OutTun(i));
    Tunnel *setTun = possTun->GetSetTunnel();
    for(unsigned int j = 0; j < possTun->m_children.size(); ++j) {
      NodeConn *conn = possTun->m_children[j];
      conn->m_n->ChangeInput1Way(possTun,conn->m_num,setTun,conn->m_num);
      delete conn;
    }
    possTun->m_children.clear();
    m_outTuns.push_back(setTun);
    setTun->m_pset = this;
  }
  
  AddPoss(poss);
}

void RealPSet::UpdateRealPSetPointers(RealPSet *oldPtr, RealPSet *newPtr)
{
  //  cout << "updating " << this << endl;
  if (m_mergeLeft == oldPtr) {
    m_mergeLeft = newPtr;
    if (newPtr == NULL) {
      m_mergeRight = NULL;
      m_leftInMap.clear();
      m_rightInMap.clear();
      m_leftOutMap.clear();
      m_rightOutMap.clear();
    }
    return;
  }
  if (m_mergeRight == oldPtr) {
    m_mergeRight = newPtr;
    if (newPtr == NULL) {
      m_mergeLeft = NULL;
      m_leftInMap.clear();
      m_rightInMap.clear();
      m_leftOutMap.clear();
      m_rightOutMap.clear();
    }
    return;
  }

  if (m_mergeMap.empty()) {
    if (!newPtr)
      return;
    else 
      throw;
  }
  PSetMapIter setIter = m_mergeMap.begin();
  for(; setIter != m_mergeMap.end(); ++setIter) {
    if (setIter->first == oldPtr) {
      if (newPtr == NULL) {
	m_mergeMap.erase(setIter);	
	return;
      }
      RealPSet *second = setIter->second;
      m_mergeMap.erase(setIter);
      m_mergeMap.insert(PSetMapPair(newPtr, second));
      return;
    }
    if (setIter->second == oldPtr) {
      if (newPtr == NULL) {
	m_mergeMap.erase(setIter);	
	return;
      }
      RealPSet *first = setIter->first;
      m_mergeMap.erase(setIter);
      m_mergeMap.insert(PSetMapPair(first, newPtr));
      return;
    }
  }
  throw;
}

void RealPSet::Migrate()
{
  ShadowPSet *shadowToReplace = NULL;
  for(unsigned int shadowNum = 0; shadowNum < m_shadows.size(); ++shadowNum) {
    BasePSet *tmp = m_shadows[shadowNum];
    if (!(tmp->m_flags & SETISDELETINGFLAG)) {
      shadowToReplace = (ShadowPSet*)tmp;
      break;
    }
  }

  //  cout << "migrating " << this << endl;
  if (!shadowToReplace) {
    //    cout << "deleting, no shadows " << this << endl;
#if PRINTTRACKING
    cout << "deleting, no shadows " << this << endl;
#endif
    if (m_mergeLeft) {
      //      cout << "updating left\n";
      m_mergeLeft->UpdateRealPSetPointers(this, NULL); 
      m_mergeLeft = NULL;
    }
    if (m_mergeRight) {
      //      cout << "updating right\n";
      m_mergeRight->UpdateRealPSetPointers(this, NULL);
      m_mergeLeft = NULL;
    }
    //    cout << "map's size : " << m_mergeMap.size() << endl;
    PSetMapIter setIter = m_mergeMap.begin();
    for(; setIter != m_mergeMap.end(); ++setIter) {
      //      cout << "updating first\n";
      setIter->first->UpdateRealPSetPointers(this, NULL);
      //      cout << "updating second\n";
      setIter->second->UpdateRealPSetPointers(this, NULL);
    }
    m_mergeMap.clear();
    return;
  }

  //  m_shadows.erase(m_shadows.begin());


  BasePSet *newBase = GetNewInst();
  if (!newBase->IsReal())
    throw;

  RealPSet *newSet = (RealPSet*)newBase;

  //  cout << "migrating " << this << " to shadow " << shadowToReplace << ", replaced with " << newSet << endl;

  if (shadowToReplace->IsLoop() != IsLoop())
    throw;

  if (shadowToReplace->m_realPSet != this) {
    cout << "shadow doesn't point to me\n";
    throw;
  }


  newSet->m_functionality = m_functionality;
  newSet->m_shadows.swap(m_shadows);



  
#if PRINTTRACKING
  cout << "migrating " << this << " to shadow " << shadowToReplace << ", replaced with " << newSet << endl;
#endif


  if (m_mergeLeft)
    m_mergeLeft->UpdateRealPSetPointers(this, newSet);
  if (m_mergeRight)
    m_mergeRight->UpdateRealPSetPointers(this, newSet);
  PSetMapIter setIter = m_mergeMap.begin();
  for(; setIter != m_mergeMap.end(); ++setIter) {
    setIter->first->UpdateRealPSetPointers(this, newSet);
    setIter->second->UpdateRealPSetPointers(this, newSet);
  }


  newSet->m_leftInMap.swap(m_leftInMap);
  newSet->m_rightInMap.swap(m_rightInMap);
  newSet->m_leftOutMap.swap(m_leftOutMap);
  newSet->m_rightOutMap.swap(m_rightOutMap);
  newSet->m_mergeMap.swap(m_mergeMap);
  newSet->m_mergeLeft = m_mergeLeft;
  newSet->m_mergeRight = m_mergeRight;

  PSetVecIter iter = newSet->m_shadows.begin();
  for(; iter != newSet->m_shadows.end(); ++iter) {
    ((ShadowPSet*)(*iter))->m_realPSet = newSet;
  }

  newSet->m_inTuns.swap(shadowToReplace->m_inTuns);
  if (newSet->m_inTuns.size() != m_inTuns.size())
    throw;
	
  NodeVecIter nodeIter = newSet->m_inTuns.begin();
  NodeVecIter oldSetNodeIter = m_inTuns.begin();
  for( ; nodeIter != newSet->m_inTuns.end(); ++nodeIter, ++oldSetNodeIter) {
    Tunnel *tun = (Tunnel*)(*nodeIter);
    Tunnel *oldTun = (Tunnel*)(*oldSetNodeIter);
    oldTun->RemoveAllChildren2Way();
    tun->MigrateFromOldTun(oldTun);
    tun->m_pset = newSet;
    NodeConnVecIter connIter = tun->m_children.begin();
    for(; connIter != tun->m_children.end(); ++connIter) {
      Node *child = (*connIter)->m_n;
      ((Tunnel*)child)->m_pset = newSet;
    }
  }
  m_inTuns.clear();

  newSet->m_outTuns.swap(shadowToReplace->m_outTuns);
  if (newSet->m_outTuns.size() != m_outTuns.size())
    throw;

  nodeIter = newSet->m_outTuns.begin();
  oldSetNodeIter = m_outTuns.begin();
  for( ; nodeIter != newSet->m_outTuns.end(); ++nodeIter, ++oldSetNodeIter) {
    Node *tun = *nodeIter;
    ((Tunnel*)tun)->m_pset = newSet;
    (*oldSetNodeIter)->RemoveAllInputs2Way();
    NodeConnVecIter connIter = tun->m_inputs.begin();
    for(; connIter != tun->m_inputs.end(); ++connIter) {
      Node *in = (*connIter)->m_n;
      ((Tunnel*)in)->m_pset = newSet;
    }
  }
  m_outTuns.clear();
  
  newSet->m_posses.swap(m_posses);
  PossMMapIter possIter = newSet->m_posses.begin();
  for( ; possIter != newSet->m_posses.end(); ++possIter) {
    Poss *poss = possIter->second;
    if (newSet->m_inTuns.size() != poss->m_inTuns.size())
      throw;
    for(unsigned int i = 0; i < newSet->m_inTuns.size(); ++i) {
      Tunnel *setTun = (Tunnel*)(newSet->m_inTuns[i]);
      Tunnel *possTun = (Tunnel*)(poss->m_inTuns[i]);
      if (!possTun->m_inputs.empty()) {
	cout << possTun << endl;
	cout << possTun->m_inputs.size() << endl;
	cout << possTun->m_inputs[0] << endl;
	cout << possTun << endl;
	cout << possTun->GetNodeClass() << endl;
	cout << possTun->GetType() << endl;
	cout.flush();
	cout << "is " << possTun->Input(0)->GetNodeClass() << endl;
	throw;
      }
      possTun->AddInput(setTun,0);
    }
    for(unsigned int i = 0; i < newSet->m_outTuns.size(); ++i) {
      Tunnel *setTun = (Tunnel*)(newSet->m_outTuns[i]);
      Tunnel *possTun = (Tunnel*)(poss->m_outTuns[i]);
      if (!possTun->m_children.empty())
	throw;
      setTun->AddInput(possTun,0);
    }
    poss->m_pset = newSet;
  }
  shadowToReplace->m_ownerPoss->AddPSet(newSet, true);
  PSetVecIter psetIter = shadowToReplace->m_ownerPoss->m_sets.begin();
  for(; psetIter != shadowToReplace->m_ownerPoss->m_sets.end(); ++psetIter) {
    if (*psetIter == shadowToReplace) {
      shadowToReplace->m_ownerPoss->m_sets.erase(psetIter);
      break;
    }
  }
  newSet->ClearDeletingRecursively();
  delete shadowToReplace;
  //  cout << "building on " << newSet->m_ownerPoss << endl;
  //  newSet->m_ownerPoss->BuildDataTypeCache();
}

RealPSet::~RealPSet()
{
  if (!(m_flags & SETTOPLEVELFLAG) && (m_inTuns.empty() || m_outTuns.empty()))
    throw;
  SetDeletingRecursively();
  Migrate();
  if (!m_inTuns.empty()) {
    NodeVecIter iter = m_inTuns.begin();
    for(; iter != m_inTuns.end(); ++iter) {
      if (!(*iter)->m_inputs.empty())
	throw;
      (*iter)->RemoveAllChildren2Way();
    }
    m_inTuns.clear();
  }
  if (!m_outTuns.empty()) {
    NodeVecIter iter = m_outTuns.begin();
    for(; iter != m_outTuns.end(); ++iter) {
      if (!(*iter)->m_children.empty())
	throw;
      (*iter)->RemoveAllInputs2Way();
    }
    m_outTuns.clear();
  }
  if (!m_shadows.empty()) {
    PSetVecIter iter = m_shadows.begin();
    for( ; iter != m_shadows.end(); ++iter) {
      ((ShadowPSet*)*iter)->m_realPSet = NULL;
    }
    m_shadows.clear();
  }
  PossMMapIter iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter) {
    delete (*iter).second;
  }
}

void RealPSet::AddPossesOrDispose(PossMMap &mmap, PossMMap *added)
{
  if (m_functionality.empty())
    throw;
  PossMMapIter newIter = mmap.begin();
  for( ; newIter != mmap.end(); ++newIter) {
    Poss *poss = (*newIter).second;
    poss->RemoveConnectionToSet();

    if (poss->m_inTuns.size() != m_inTuns.size()) {
      cout << "wrong number\n";
      throw;
    }
    
    for(unsigned int i = 0; i < poss->m_inTuns.size(); ++i) {
      NodeConn *conn = new NodeConn(InTun(i),0);
      poss->InTun(i)->m_inputs.push_back(conn);
    }
    
    bool existing = false;
    PossMMapRangePair pair = m_posses.equal_range(poss->GetHash());
    for( ; !existing && pair.first != pair.second; ++pair.first) {
      if (*((*(pair.first)).second) == *poss) {
        delete poss;
        existing = true;
      }
    }
    
    if (!existing) {
      poss->RemoveConnectionToSet();
      AddPoss(poss);
      if (added)
	added->insert(PossMMapPair(poss->GetHash(),poss));
    }
  }
}

void RealPSet::AddPoss(Poss *poss)
{
  if (m_functionality.empty()) {
    if (m_posses.size())
      throw;
    else {
      m_functionality = poss->GetFunctionalityString();
      if (IsLoop()) {
	//	Loop *loop = (Loop*)this;
	m_functionality += (char)((dynamic_cast<const LoopInterface*>(this))->GetBSSize().GetSize());
      }
      if (m_functionality.empty())
	throw;
    }
  }
  if (m_inTuns.size() != poss->m_inTuns.size()) {
    cout << "New poss doesn't have same number of inputs\n";
    throw;
  }
  if (m_outTuns.size() != poss->m_outTuns.size()) {
    cout << "New poss doesn't have same number of outputs\n";
    cout << m_outTuns.size() << " != " << poss->m_outTuns.size() << endl;
    throw;
  }
  
  for(unsigned int i = 0; i < m_inTuns.size(); ++i) {
    Node *setTun = InTun(i);
    Node *possTun = poss->InTun(i);
    if (!possTun->m_inputs.empty()) {
      cout << "(!possTun->m_inTuns.empty()\n";
      throw;
    }
    possTun->AddInput(setTun,0);
  }
  
  for(unsigned int i = 0; i < m_outTuns.size(); ++i) {
    Node *setTun = OutTun(i);
    Node *possTun = poss->OutTun(i);
    if (!possTun->m_children.empty()) {
      cout << "!possTun->m_outTuns.empty()\n";
      cout << possTun->GetType() << " has child "
	   << possTun->Child(0)->GetType() << endl;
      cout << possTun << endl;
      throw;
    }
    setTun->AddInput(possTun,0);
  }
  
  m_posses.insert(PossMMapPair(poss->GetHash(),poss));
  poss->m_pset = this;
}

bool RealPSet::operator==(const BasePSet &rhs) const
{
  if (rhs.IsShadow()) {
    return *this == *(rhs.GetReal());
    //    return (*this)==(*(((ShadowPSet&)rhs).m_realPSet));
  }
  if (!rhs.IsReal())
    throw;
  const RealPSet &realRhs = (RealPSet&)rhs;
  if (m_inTuns.size() != realRhs.m_inTuns.size()
      || m_outTuns.size() != realRhs.m_outTuns.size())
    return false;
  if (GetFunctionalityString() != realRhs.GetFunctionalityString()) {
    return false;
  }
  else {
    if (IsLoop()) {
      if (!realRhs.IsLoop())
        return false;
      else {
	const LoopInterface *loop1 = dynamic_cast<const LoopInterface*>(this);
	const LoopInterface *loop2 = dynamic_cast<const LoopInterface*>(&rhs);
	if (((RealLoop*)this)->IsUnrolled() != dynamic_cast<const RealLoop*>(&rhs)->IsUnrolled())
	  return false;
#if TWOD
	if (loop1->GetDimName() != loop2->GetDimName())
	  return false;
#endif
	if (loop1->GetBSSize() != loop2->GetBSSize())
	  return false;
	if (loop1->GetBS() == 0) {
	  cout << "BS bs\n";
	  cout << loop1->GetBSSize().GetSize() << endl;
	  throw;
	}
        for (unsigned int i = 0; i < m_inTuns.size(); ++i) {
          const LoopTunnel *tun1 = (LoopTunnel*)(m_inTuns[i]);
          const LoopTunnel *tun2 = (LoopTunnel*)(realRhs.m_inTuns[i]);
#if DOBLIS
          if (((Loop*)(tun1->m_pset))->m_comm != ((Loop*)(tun2->m_pset))->m_comm)
            return false;
#endif
          if (tun1->m_statTL != tun2->m_statTL
              || tun1->m_statBL != tun2->m_statBL
              || tun1->m_statTR != tun2->m_statTR
              || tun1->m_statBR != tun2->m_statBR)
            return false;
          if (tun1->IsSplit()) {
            if (tun2->GetNodeClass() == tun1->GetNodeClass()) {
#if TWOD
              if (((SplitBase*)tun1)->m_dir != ((SplitBase*)tun2)->m_dir)
                return false;
#else
	      if (((SplitBase*)tun1)->m_partDim != ((SplitBase*)tun2)->m_partDim)
                return false;
#endif
	      if (tun1->GetNodeClass() == SplitUnrolled::GetClass()) {
		if (((SplitUnrolled*)tun1)->m_unrollFactor != ((SplitUnrolled*)tun2)->m_unrollFactor)
		  return false;
	      }
            }
            else
              return false;
          }
          else if (tun2->IsSplit()) {
            return false;
          }
        }
      }
    }
#if DOBLIS
    else if (IsCritSect()) {
      if (!realRhs.IsCritSect())
	return false;
    }
#endif
    return true;
  }
}

Cost RealPSet::Prop()
{
  if(m_flags & SETHASPROPEDFLAG)
    return m_cost;

  if (m_posses.empty()) {
    cout << "I'm empty\n";
    throw;
  }

  if (m_functionality.empty()) {
    cout << m_posses.size() << endl;
    (*(m_posses.begin())).second->PrintSetConnections();
    throw;
  }

  if (m_mergeLeft && !m_mergeRight)
    throw;
  if (m_mergeRight && !m_mergeLeft)
    throw;
  if (m_mergeLeft) {
    if (m_leftInMap.empty()
	|| m_rightInMap.empty()
	|| m_leftOutMap.empty()
	|| m_rightOutMap.empty())
      {
	throw;
      }

    const int inSize = m_inTuns.size();
    const int outSize = m_outTuns.size();
    vector<int>::iterator mapIter;
    for(mapIter = m_leftInMap.begin(); mapIter != m_leftInMap.end(); ++mapIter) {
      if (*mapIter >= inSize)
	throw;
    }
    for(mapIter = m_rightInMap.begin(); mapIter != m_rightInMap.end(); ++mapIter) {
      if (*mapIter >= inSize) {
	cout << "mapping to " << *mapIter << " out of " << inSize << endl;
	throw;
      }
    }
    for(mapIter = m_leftOutMap.begin(); mapIter != m_leftOutMap.end(); ++mapIter) {
      if (*mapIter >= outSize)
	throw;
    }
    for(mapIter = m_rightOutMap.begin(); mapIter != m_rightOutMap.end(); ++mapIter) {
      if (*mapIter >= outSize) {
	cout << *mapIter << " >= " << outSize << endl;
	throw;
      }
    }
  }
  else {
    if (!m_leftInMap.empty()
	|| !m_rightInMap.empty()
	|| !m_leftOutMap.empty()
	|| !m_rightOutMap.empty())
      {
	throw;
      }
  }


  //BAM Par + check for > 1
  for (unsigned int i = 0; i < m_inTuns.size(); ++i) {
    Node *in = InTun(i);
    for (unsigned int j = 0; j < in->m_children.size(); ++j) {
      Node *child = in->m_children[j]->m_n;
      if (child->m_inputs.size() != 1) {
        cout << "child->m_inputs.size() != 1\n";
        throw;
      }
      if (child->Input(0) != in) {
        cout << "child->m_inputs[0]->m_n != in\n";
        throw;
      }
    }
  }
  
  for (unsigned int i = 0; i < m_outTuns.size(); ++i) {
    Node *out = m_outTuns[i];
    for (unsigned int j = 0; j < out->m_inputs.size(); ++j) {
      Node *parent = out->Input(j);
      if (parent->m_children.size() != 1) {
        cout << "parent->m_children.size() != 1\n";
        throw;
      }
      if (parent->m_children[0]->m_n != out) {
        cout << "parent->m_children[0]->m_n != out\n";
        throw;
      }
    }
  }
  
  if (!IsTopLevel() && !m_ownerPoss) {
    cout << "no owner\n";
    throw;
  }


  for(unsigned int i = 0; i < m_inTuns.size(); ++i) {
    InTun(i)->Prop();
  }
  PossMMapIter iter;
  int j = 0;
#if !USESHADOWS
#pragma omp parallel private(j,iter)
#endif
  {
    iter = m_posses.begin();
    j = 0;
    int size = m_posses.size();
#if !USESHADOWS
#pragma omp for schedule(static) 
#endif
    for (int i = 0; i < size; ++i) {
      if (j > i) {
	cout << "uhoh\n";
	throw;
      }
      while (j < i) {
	++iter;
	++j;
      }
      
      (*iter).second->Prop();
    }
  }


  iter = m_posses.begin();
  while (iter != m_posses.end()) {
    Poss *poss = (*iter).second;
    if (!poss->IsSane()) {
      RemoveAndDeletePoss(poss, false);
      m_posses.erase(iter);
      iter = m_posses.begin();
      if (!m_posses.size()) {
	cout << "Ran out of posses in set " << this << endl;
	throw;
      }
    }
    else {
      ++iter;
    }
  }

  //  m_cost = -1;
  iter = m_posses.begin();
  do {
    Poss *poss = (*iter).second;
    if (m_cost <= 0)
      m_cost = poss->m_cost;
    else if (m_cost > poss->m_cost)
      m_cost = poss->m_cost;
    if (poss->GetHash() != (*iter).first) {
      m_posses.erase(iter);
	  
      bool existing = false;
      PossMMapRangePair pair = m_posses.equal_range(poss->GetHash());
      for( ; !existing && pair.first != pair.second; ++pair.first) {
	if (*((*(pair.first)).second) == *poss) {
	  RemoveAndDeletePoss(poss, false);
	  existing = true;
	}
      }
    
      if (!existing)
	m_posses.insert(PossMMapPair(poss->GetHash(),poss));
      
      iter = m_posses.begin();
    }
    else {
      ++iter;
    }
  } while (iter != m_posses.end());


  PSetVecIter shadowIter = m_shadows.begin();
  for(; shadowIter != m_shadows.end(); ++shadowIter) {
    if (((ShadowPSet*)(*shadowIter))->m_realPSet != this)
      throw;
    if ((*shadowIter)->IsLoop() != IsLoop()) {
      cout << "shadow and real loop satus don't agree\n";
      throw;
    }
  }
  
  m_flags |= SETHASPROPEDFLAG;
  return m_cost;
}

void RealPSet::Cull(Phase phase)
{
  PossMMapIter iter2 = m_posses.begin();
  for(; iter2 != m_posses.end(); ++iter2) 
    (*iter2).second->m_flags &= ~POSSISSANEFLAG;

  int j;
  PossMMapIter iter;
  
#pragma omp parallel private(j,iter)
  {
    iter = m_posses.begin();
    j = 0;
    int size = m_posses.size();
#pragma omp for schedule(static) 
    for (int i = 0; i < size; ++i) {
      if (i < j) {
	cout << "uhoh";
	cout.flush();
	throw;
      }
      while (j < i) {
	++iter;
	++j;
      }
      Poss *poss = (*iter).second;
      if (poss->IsSane())
	throw;
      poss->m_flags |= POSSISSANEFLAG;
      poss->Cull(phase);
    }
  }

  iter = m_posses.begin();
  while (iter != m_posses.end()) {
    Poss *poss = (*iter).second;
    if (!poss->IsSane()) {
      RemoveAndDeletePoss(poss, false);
      m_posses.erase(iter);
      iter = m_posses.begin();
      if (!m_posses.size()) {
        cout << "Ran out of posses in set " << this << endl;
        throw;
      }
    }
    else {
      ++iter;
    }
  }
}

void RealPSet::CullWorstPerformers(double percentToCull, int ignoreThreshold)
{
  if (m_posses.size() > ignoreThreshold) {
    SortedPossQueue queue;
    PossMMapIter iter = m_posses.begin();
    for(; iter != m_posses.end(); ++iter) {
      queue.push(iter->second);
    }
    if (queue.size() != m_posses.size())
      throw;
    int numToKeep = ceil((1.0-percentToCull) * m_posses.size());
    for(; numToKeep > 0; --numToKeep)
      queue.pop();
    while(!queue.empty()) {
      RemoveAndDeletePoss(queue.top(), true);
      queue.pop();
    }
  }
  PossMMapIter iter;
  int j = 0;
#pragma omp parallel private(j,iter)
  {
    iter = m_posses.begin();
    j = 0;
    int size = m_posses.size();
#pragma omp for schedule(static) 
    for (int i = 0; i < size; ++i) {
      if (j > i) {
	cout << "uhoh\n";
	throw;
      }
      while (j < i) {
	++iter;
	++j;
      }
      iter->second->CullWorstPerformers(percentToCull, ignoreThreshold);
    }
  }
}

void RealPSet::RemoveAndDeletePoss(Poss *poss, bool removeFromMyList)
{
  if (removeFromMyList && m_posses.size() <= 1) {
    if (m_posses.size()) {
      poss->ForcePrint();
    }
    throw;
  }
  if (removeFromMyList) {
    PossMMapIter possIter = m_posses.begin();
    bool found = false;
    for(; !found && possIter != m_posses.end(); ++possIter) {
      if ((*possIter).second == poss) {
        m_posses.erase(possIter);
        found = true;
	break;
      }
    }
    if (!found)
      throw;
  }
  
  if (poss->m_inTuns.size() != m_inTuns.size()) {
    cout << "(poss->m_inTuns.size() != m_inTuns.size())\n";
    throw;
  }
  for (unsigned int i = 0; i < m_inTuns.size(); ++i)
    InTun(i)->RemoveChild(poss->InTun(i),0);
  
  if (poss->m_outTuns.size() != m_outTuns.size()) {
    cout << "(poss->m_outTuns.size() != m_outTuns.size())\n";
    throw;
  }
  for (unsigned int i = 0; i < m_outTuns.size(); ++i)
    OutTun(i)->RemoveInput(poss->OutTun(i),0);
  poss->m_pset = (RealPSet*)(0xDEADBEEF);
  delete poss;
}

void RealPSet::ClearBeforeProp()
{
  BasePSet::ClearBeforeProp();

  m_cost = -1;
  
  PossMMapIter iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter)
    (*iter).second->ClearBeforeProp();
}


bool RealPSet::TakeIter(const TransMap &transMap,
                    const TransMap &simplifiers)
{
  bool newOne = false;
  PossMMap actuallyAdded;
  
#ifdef _OPENMP
  static omp_lock_t lock;
  static bool inited = false;
  if (!inited) {
    omp_init_lock(&lock);
    inited = true;
  }
#endif
  
  PossMMap mmap;
  PossMMapIter iter;
  int j = 0;

#pragma omp parallel private(j,iter)
  {
    iter = m_posses.begin();
    j = 0;
    int size = m_posses.size();
#pragma omp for schedule(static) 
    for (int i = 0; i < size; ++i) {
      while (j < i) {
	++iter;
	++j;
      }
      Poss *poss = (*iter).second;
      if (!poss->m_fullyExpanded) {
	bool didSomething;
	PossMMap newPosses;
	didSomething = poss->TakeIter(transMap, simplifiers, newPosses);
	if (didSomething && (!newOne || newPosses.size() > 0)) {
#ifdef _OPENMP
	  omp_set_lock(&lock);
#endif
	  newOne = true;
	  PossMMapIter newPossesIter = newPosses.begin();
	  for(; newPossesIter != newPosses.end(); ++newPossesIter) {
	    if (!AddPossToMMap(mmap, (*newPossesIter).second, (*newPossesIter).second->GetHash()))
	      delete (*newPossesIter).second;
	  }
#ifdef _OPENMP
	  omp_unset_lock(&lock);
#endif
	}
      }
    }
  }
  //have to add these at the end or we'd be adding posses while iterating
  // over the posses
  // BAM: Or do I?
  if (mmap.size()) {
    if (mmap.size() > 100) {
      cout << "\t\tAdding " << mmap.size() << " posses\n";
      //      cout << "\t\t" << m_functionality << endl;
    }
    AddPossesOrDispose(mmap, &actuallyAdded);
    if (mmap.size() > 100)
      cout << "\t\tDone adding ( " << actuallyAdded.size() << " actually added )\n";
    PossMMapIter added = actuallyAdded.begin();
    for(; added != actuallyAdded.end(); ++added) {
      (*added).second->BuildDataTypeCache();
#if FORMSETSEARLY
      if (CurrPhase == DPTENSORPHASE)
       	(*added).second->FormSets(SUMSCATTERTENSORPHASE);
#endif
    }
  }
  return newOne;
}

void RealPSet::Duplicate(const BasePSet *orig, NodeMap &map, bool possMerging, bool useShadows)
{
  BasePSet::Duplicate(orig, map, possMerging,useShadows);
  if (!orig->IsReal())
    throw;
  const RealPSet *real = (RealPSet*)orig;
  m_functionality = real->m_functionality;
  if (m_functionality.empty())
    throw;
  PossMMapConstIter iter2 = real->m_posses.begin();
  for( ; iter2 != real->m_posses.end(); ++iter2) {
    const Poss *oldPoss = (*iter2).second;
    Poss *newPoss = new Poss;
    newPoss->Duplicate(oldPoss, map, possMerging, useShadows);
    m_posses.insert(PossMMapPair(newPoss->GetHash(),newPoss));
    newPoss->m_pset = this;
  }
  
}

void RealPSet::PatchAfterDuplicate(NodeMap &map)
{
  PossMMapIter iter2 = m_posses.begin();
  for( ; iter2 != m_posses.end(); ++iter2) {
    (*iter2).second->PatchAfterDuplicate(map);
    if ((*iter2).first != (*iter2).second->GetHash()) {
      cout << "different hash in PatchAfterDuplicate\n";
      throw;
    }

  }
}

void RealPSet::CombineAndRemoveTunnels()
{
  if (!m_shadows.empty())
    throw;
  for(unsigned int inIdx1 = 0; inIdx1 < m_inTuns.size(); ++inIdx1) {
    for(unsigned int inIdx2 = inIdx1+1; inIdx2 < m_inTuns.size(); ++inIdx2) {
      Node *setInput1 = InTun(inIdx1);
      Node *setInput2 = InTun(inIdx2);
      if (setInput1->m_inputs.size() > 1 || setInput2->m_inputs.size() > 1) {
        cout << "setInput1->m_inputs.size() > 1 || setInput2->m_inputs.size() > 1\n";
        cout << setInput1->m_inputs.size() << endl;
        cout << setInput2->m_inputs.size() << endl;
        throw;
      }
      if (setInput1->m_inputs.size()
          && setInput2->m_inputs.size()
          && setInput1->Input(0) == setInput2->Input(0)
          && setInput1->InputConnNum(0) == setInput2->InputConnNum(0))
	{
	  if (setInput1->IsLoopTunnel() && setInput2->IsLoopTunnel()) {
	    const LoopTunnel *tun1 = (LoopTunnel*)(setInput1);
	    const LoopTunnel *tun2 = (LoopTunnel*)(setInput2);
	    if ((tun1->IsSplit() && !tun2->IsSplit()) || (!tun1->IsSplit() && tun2->IsSplit()))
	      {
		continue;
	      }
	    if (tun1->IsSplit()) {
	      if (tun1->GetNodeClass() != tun2->GetNodeClass())
		continue;
	      if (tun1->GetNodeClass() == SplitUnrolled::GetClass()) {
		if (((SplitUnrolled*)tun1)->m_unrollFactor != ((SplitUnrolled*)tun2)->m_unrollFactor)
		  continue;
	      }
#if TWOD
	      if (((SplitBase*)setInput1)->m_dir != ((SplitBase*)setInput2)->m_dir)
		continue;
#else
	      if (((SplitBase*)setInput1)->m_partDim != ((SplitBase*)setInput2)->m_partDim)
		continue;
#endif
	    }
	  }
	  if (setInput1->m_children.size() != setInput2->m_children.size()) {
	    cout << setInput1->m_children.size() << " vs. " << setInput2->m_children.size() << endl;
	    cout << setInput1 << endl;
	    cout << setInput2 << endl;
	    setInput2->PrintChildren();
	    throw;
	  }
	  NodeConnVecIter connIter1 = setInput1->m_children.begin();
	  NodeConnVecIter connIter2 = setInput2->m_children.begin();
	  for(; connIter1 != setInput1->m_children.end(); ++connIter1, ++connIter2) {
	    Node *possIn1 = (*connIter1)->m_n;
	    Node *possIn2 = (*connIter2)->m_n;
	    possIn2->RedirectAllChildren(possIn1);
	  }
	  if (connIter2 != setInput2->m_children.end()) {
	    cout << "connIter2 != setInput2->m_children.end()\n";
	    throw;
	  }
	  setInput2->Input(0)->RemoveChild(setInput2,setInput2->InputConnNum(0));
	  delete setInput2->InputConn(0);
	  setInput2->m_inputs.clear();
	}
    }
  }
  
  for(unsigned int i = 0; i < m_inTuns.size(); ++i) {
    Node *tun = InTun(i);
    if (tun->m_inputs.empty()) {
      NodeConnVecIter inChildIter = tun->m_children.begin();
      for(; inChildIter != tun->m_children.end(); ++inChildIter) {
        Node *child = (*inChildIter)->m_n;
        if (child->m_children.size()) {
          cout << "setInTun " << i << " " << tun << " " << tun->GetType() << endl;//<< " " << tun->GetName(0).str() << endl;
          cout << "child/possInTun " << child << " " << child->GetType() << endl;
          cout << "poss tunnel connected to input set tunnel without inputs has children\n";
          for (unsigned int j = 0; j < m_inTuns.size(); ++j) {
            cout << "other " << j << " = " << InTun(j)->GetType() << " " << InTun(j)->GetName(0).str() << endl;
          }
          //throw;
        }
        if (!child->IsTunnel()) {
          cout << "!child->IsTunnel()\n";
          throw;
        }
        child->m_poss->DeleteNode(child);
        delete *inChildIter;
      }
      tun->m_children.clear();
      delete tun;
      m_inTuns.erase(m_inTuns.begin()+i);
      for (unsigned int j = 0; j < m_leftInMap.size(); ++j) {
	int val = m_leftInMap[j];
	if (val == (int)i)
	  m_leftInMap[j] = -1;
	if (val > (int)i)
	  m_leftInMap[j] = val-1;
      }
      for (unsigned int j = 0; j < m_rightInMap.size(); ++j) {
	int val = m_rightInMap[j];
	if (val == (int)i)
	  m_rightInMap[j] = -1;
	else if (val > (int)i)
	  m_rightInMap[j] = val-1;
      }
      --i;
    }
  }
  
  for(unsigned int i = 0; i < m_outTuns.size(); ++i) {
    Node *tun = OutTun(i);
    if (tun->IsLoopTunnel()) {
      if (!((LoopTunnel*)tun)->IsConst()) {
        continue;
      }
    }
    if (tun->m_children.empty()) {
      bool skip = false;
      Node *inTun = tun->Input(0);
      if (inTun->IsLoopTunnel()) {
        Node *source = inTun->Input(inTun->m_inputs.size()-1);
        bool found = false;
        for (unsigned int j = 0; !found && j < source->m_children.size(); ++j) {
          NodeConn *conn = source->m_children[j];
          if (conn->m_num == source->NumOutputs()-1) {
            if (conn->m_n != inTun)
              found = true;
          }
        }
        if (!found)
          skip = true;
      }
      if (!skip) {
        NodeConnVecIter iter = tun->m_inputs.begin();
        for( ; iter != tun->m_inputs.end(); ++iter) {
          Node *possOut = (*iter)->m_n;
          delete possOut->m_children[0];
          possOut->m_children.clear();
          possOut->m_poss->DeleteChildAndCleanUp(possOut, true, true);
          delete *iter;
        }
        tun->m_inputs.clear();
        delete tun;
        m_outTuns.erase(m_outTuns.begin()+i);
	
	for(unsigned int j = 0; j < m_leftOutMap.size(); ++j) {
	  int val = m_leftOutMap[j];
	  if (val == (int)i)
	    m_leftOutMap[j] = -1;
	  else if (val > (int)i)
	    m_leftOutMap[j] = val-1;
	}

	for(unsigned int j = 0; j < m_rightOutMap.size(); ++j) {
	  int val = m_rightOutMap[j];
	  if (val == (int)i)
	    m_rightOutMap[j] = -1;
	  else if (val > (int)i)
	    m_rightOutMap[j] = val-1;
	}

        --i;
        continue;
      }
    }
    /*
    if (tun->IsLoopTunnel() && (((LoopTunnel*)tun)->IsCombine())) {
      for(unsigned int j = i+1; j < m_outTuns.size(); ++j) {
        Node *tun2 = OutTun(j);
        if (tun2->GetNodeClass() == tun->GetNodeClass()) {
          if (!((LoopTunnel*)tun)->IsConst()) {
            continue;
          }
          Node *in1 = tun->Input(0);
          Node *in2 = tun2->Input(0);
          if (in1->Input(in1->m_inputs.size()-1) !=
              in2->Input(in2->m_inputs.size()-1)) {
            continue;
          }
          SplitBase *split1 = (SplitBase*)tun;
          SplitBase *split2 = (SplitBase*)tun2;
#if TWOD
          if (split1->m_dir != split2->m_dir) {
            continue;
          }
#else
          if (split1->m_partDim != split2->m_partDim) {
            continue;
          }
#endif
	  if (split1->GetNodeClass() != split2->GetNodeClass())
	    continue;
	  if (split1->GetNodeClass() == SplitUnrolled::GetClass()) {
	    if (((SplitUnrolled*)split1)->m_unrollFactor != ((SplitUnrolled*)split1)->m_unrollFactor)
	      continue;
	  }
          split2->RedirectAllChildren(split1);
          NodeConnVecIter iter = split2->m_inputs.begin();
          for( ; iter != split2->m_inputs.end(); ++iter) {
            Node *possOut = (*iter)->m_n;
            delete possOut->m_children[0];
            possOut->m_children.clear();
            possOut->m_poss->DeleteChildAndCleanUp(possOut);
            delete *iter;
          }
          split2->m_inputs.clear();
          delete split2;
          m_outTuns.erase(m_outTuns.begin()+j);
          --j;
          --i;
        }
      }
    }
    */
  }
}

void RealPSet::Simplify(const TransMap &simplifiers, bool recursive)
{
  //BAM par
  PossMMapIter iter;
  int j = 0;
#pragma omp parallel private(j,iter)
  {
    iter = m_posses.begin();
    j = 0;
    int size = m_posses.size();
#pragma omp for schedule(static) 
    for (int i = 0; i < size; ++i) {
      while (j < i) {
	++iter;
	++j;
      }
      Poss *poss = (*iter).second;
      poss->Simplify(simplifiers, recursive);
    }
  }
  iter = m_posses.begin();
  do {
    Poss *poss = (*iter).second;
    if (poss->GetHash() != (*iter).first) {
      m_posses.erase(iter);

      bool existing = false;
      PossMMapRangePair pair = m_posses.equal_range(poss->GetHash());
      for( ; !existing && pair.first != pair.second; ++pair.first) {
	if (*((*(pair.first)).second) == *poss) {
	  RemoveAndDeletePoss(poss, false);
	  existing = true;
	}
      }

      if (!existing)
	m_posses.insert(PossMMapPair(poss->GetHash(),poss));
      else {
	if (m_posses.size() == 0)
	  throw;
      }

      iter = m_posses.begin();
    }
    else
      ++iter;
  } while (iter != m_posses.end());
}

void RealPSet::ClearFullyExpanded()
{
  PossMMapIter iter = m_posses.begin();
  for( ;iter != m_posses.end(); ++iter)
    (*iter).second->ClearFullyExpanded();
}


bool RealPSet::MergePosses(const TransMap &simplifiers, CullFunction cullFunc)
{
  /*
    Here's the idea:
    If one or more of the posses in this PSet has PSets on them
    that can be merged or a single PSet that can be removed (i.e.
    the posses in that set can be inlined), then do it now.  We want
    to do this from the bottom-up, though.
  */

#ifdef _OPENMP
  static omp_lock_t lock;
  static bool inited = false;
  if (!inited) {
    omp_init_lock(&lock);
    inited = true;
  }
#endif
  bool didMerge = false;
  int numPosses = (int)(m_posses.size());
  
  if (numPosses > 1) {
    PossMMap mmap;
    PossMMapIter iter;
    int j = 0;
#if !USESHADOWS
#pragma omp parallel private(j,iter)
#endif
    {
      iter = m_posses.begin();
      j = 0;
	int size = m_posses.size();
#if !USESHADOWS
#pragma omp for schedule(static) 
#endif
      for (int i = 0; i < size; ++i) {
	while (j < i) {
	  ++iter;
	  ++j;
	}
	Poss *poss = (*iter).second;
        PossMMap newPosses;
	if (poss->MergePosses(newPosses, simplifiers, cullFunc)) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
	    didMerge = true;

	    if (!newPosses.empty()) {
#ifdef _OPENMP
	      omp_set_lock(&lock);
#endif
	      PossMMapIter newPossesIter = newPosses.begin();
	      for(; newPossesIter != newPosses.end(); ++newPossesIter) {
		if (!AddPossToMMap(mmap, (*newPossesIter).second, (*newPossesIter).first))
		  delete (*newPossesIter).second;
	      }
#ifdef _OPENMP
	      omp_unset_lock(&lock);
#endif
	    }
	}
      }
    }
    PossMMapIter mmapIter = mmap.begin();
    for(; mmapIter != mmap.end(); ++mmapIter)
      (*mmapIter).second->BuildDataTypeCache();
    AddPossesOrDispose(mmap);
  }
  else {
    PossMMapIter iter = m_posses.begin();
    for(; iter != m_posses.end() && !didMerge; ++iter) {
      PossMMap newPosses;
      Poss *poss = (*iter).second;
      if (poss->MergePosses(newPosses, simplifiers, cullFunc)) {
	didMerge = true;
        PossMMapIter mmapIter = newPosses.begin();
        for(; mmapIter != newPosses.end(); ++mmapIter)
          (*mmapIter).second->BuildDataTypeCache();
        AddPossesOrDispose(newPosses);
      }
    }
  }
  if (!didMerge) {
    PossMMapIter iter = m_posses.begin();
    while (iter != m_posses.end()) {
      PossMMap newPosses;
      Poss *poss = (*iter).second;
      //This poss has only a single PSet on it
      // Make new posses for me, each containing the posses in
      // poss->m_posses[0];
      if (poss->m_sets.size() == 1 && 
	  poss->m_sets[0]->IsTransparent()
	  && m_ownerPoss
	  && m_ownerPoss->m_sets.size() > 1)
        {
	  //	  didMerge = true;
	  InlinePoss(poss, newPosses);
	  m_posses.erase(iter);
          PossMMapIter mapIter = newPosses.begin();
          for(; mapIter != newPosses.end(); ++mapIter)
            (*mapIter).second->BuildDataTypeCache();
          AddPossesOrDispose(newPosses);
	  iter = m_posses.begin();
        }
      else {
	++iter;
      }
    }
  }
  return didMerge;
}


void RealPSet::FormSets(unsigned int phase)
{
  PossMMap temp = m_posses;
  m_posses.clear();
  PossMMapIter iter = temp.begin();
  while (iter != temp.end()) {
    Poss *poss = (*iter).second;
    poss->FormSets(phase);
    m_posses.insert(PossMMapPair(poss->GetHash(),poss));
    ++iter;
  }
}


GraphNum RealPSet::TotalCount() const
{
  GraphNum tot = 0;
  PossMMapConstIter iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter)
    tot += (*iter).second->TotalCount();
  return tot;
}

void RealPSet::InlinePoss(Poss *inliningPoss, PossMMap &newPosses)
{
#if PRINTTRACKING
  cout << "inlining " << this << endl;
#endif
  if (inliningPoss->m_sets[0]->IsShadow()) {
    inliningPoss->ReplaceShadowSetWithReal(0);
  }

  RealPSet *pset = (RealPSet*)(inliningPoss->m_sets[0]);
  
  NodeIntMap tunnelNumMap;
  
  NodeMap setTunnels;
  int i = 0;
  NodeVecIter iter2 = pset->m_inTuns.begin();
  for(; iter2 != pset->m_inTuns.end(); ++iter2) {
    setTunnels[*iter2] = *iter2;
    tunnelNumMap[*iter2] = i;
    ++i;
  }
  i = 0;
  iter2 = pset->m_outTuns.begin();
  for(; iter2 != pset->m_outTuns.end(); ++iter2) {
    setTunnels[*iter2] = *iter2;
    tunnelNumMap[*iter2] = i;
    ++i;
  }
  
  iter2 = m_inTuns.begin();
  for(; iter2 != m_inTuns.end(); ++iter2) {
    setTunnels[*iter2] = *iter2;
  }
  iter2 = m_outTuns.begin();
  for(; iter2 != m_outTuns.end(); ++iter2) {
    setTunnels[*iter2] = *iter2;
  }
  
  PossMMapIter setIter = pset->m_posses.begin();
  for (; setIter != pset->m_posses.end(); ++setIter) {
    Poss *currPoss = (*setIter).second;
    NodeMap map = setTunnels;
    
    //create a new poss for this that has the same poss nodes as poss and currPoss
    //redirect the poss inputs on currposs to take as input the newPoss's poss inputs
    //  or the poss nodes coppied from poss
    //redirect the poss outputs similarly
    
    Poss *newPoss = new Poss();
    //    cout << "newposs " << newPoss << endl;
    
    newPoss->m_parent = inliningPoss->m_num;
    
    NodeVecIter iter = inliningPoss->m_possNodes.begin();
    for( ; iter != inliningPoss->m_possNodes.end(); ++iter) {
      if (!(*iter)->IsTunnel(SETTUNIN) && !(*iter)->IsTunnel(SETTUNOUT)) {
        Node *newNode = (*iter)->GetNewInst();
        //        cout << "On outer poss, creating newNode " << newNode << " for " << *iter << endl;
        newNode->Duplicate(*iter, false,true);
        newPoss->m_possNodes.push_back(newNode);
        newNode->m_poss = newPoss;
        map[*iter] = newNode;
      }
      else {
        map[*iter] = *iter;
        //        cout << *iter << " is set tunnel on outer poss\n";
      }
    }
    
    iter = currPoss->m_possNodes.begin();
    for( ; iter != currPoss->m_possNodes.end(); ++iter) {
      if (!(*iter)->IsTunnel(POSSTUNIN) && !(*iter)->IsTunnel(POSSTUNOUT)) {
        Node *newNode = (*iter)->GetNewInst();
        //        cout << "On inlining poss, creating newNode " << newNode << " for " << *iter << endl;
        newNode->Duplicate(*iter, false,true);
        newPoss->m_possNodes.push_back(newNode);
        newNode->m_poss = newPoss;
        map[*iter] = newNode;
      }
      else {
        map[*iter] = *iter;
        //        cout << *iter << " is poss tunnel on inlining poss\n";
      }
    }

    {
      PSetVecIter setIter = currPoss->m_sets.begin();
      for(; setIter != currPoss->m_sets.end(); ++setIter) {
#if USESHADOWS
        BasePSet *newSet = (*setIter)->GetNewShadow();
        newSet->Duplicate(*setIter, map, true, true);
#else
        BasePSet *newSet = (*setIter)->GetNewInst();
        newSet->Duplicate(*setIter, map, true, false);
#endif
        newPoss->m_sets.push_back(newSet);
        newSet->m_ownerPoss = newPoss;
      }
      setIter = newPoss->m_sets.begin();
      for(; setIter != newPoss->m_sets.end(); ++setIter) {
        (*setIter)->PatchAfterDuplicate(map);
      }
    }
    
    newPoss->PatchAfterDuplicate(map);


    NodeSet cleanupNodes;
    
    iter = newPoss->m_possNodes.begin();
    for(; iter != newPoss->m_possNodes.end(); ++iter) {
      
      Node *node = *iter;

      NodeConnVecIter connIter = node->m_inputs.begin();
      for(; connIter != node->m_inputs.end(); ++connIter) {
        NodeConn *conn = *connIter;
        NodeIntMapIter mapIter = tunnelNumMap.find(conn->m_n);
        if (mapIter != tunnelNumMap.end()) {
	  if (!conn->m_n->IsTunnel(SETTUNOUT))
	    throw;
          Node *newParent = map[currPoss->OutTun(mapIter->second)->Input(conn->m_num)];
          unsigned int newNum = currPoss->OutTun(mapIter->second)->InputConnNum(conn->m_num);
	  if (newParent->IsTunnel(POSSTUNIN)) {
	    newParent = newParent->Input(0); // SetTunIn;
	    newNum = newParent->InputConnNum(0);
	    newParent = map[newParent->Input(0)];
	    if (!newParent)
	      throw;
	  }
	  else {
	    for(unsigned int i = 0; i < newParent->m_children.size(); ++i) {
	      if (newParent->Child(i)->IsTunnel(POSSTUNOUT)) {
		delete newParent->m_children[i];
		newParent->m_children.erase(newParent->m_children.begin()+i);
		--i;
	      }
	    }
	  }
          newParent->AddChild(node, newNum);
          conn->m_num = newNum;
          conn->m_n = newParent;
        }
      }

      unsigned int size = node->m_children.size();
      for (unsigned int i = 0; i < size; ++i) {
        NodeConn *conn = node->m_children[i];
        NodeIntMapIter mapIter = tunnelNumMap.find(conn->m_n);
        if (mapIter != tunnelNumMap.end()) {
          Node *setInput = conn->m_n;
          node->m_children.erase(node->m_children.begin()+i);
          --i;
          --size;
          unsigned int setInputNum;
          for(setInputNum = 0; setInputNum < setInput->m_inputs.size(); ++setInputNum) {
            if (map[setInput->Input(setInputNum)] == node)
              break;
          }
          if (setInputNum >= setInput->m_inputs.size()) {
            cout << "didn't find setInput that is the child\n";
            throw;
          }
          
          Node *inTun = currPoss->InTun(mapIter->second);
          if (inTun->Input(0) != setInput) {
            cout <<"inTun->Input(0) != setInput\n";
            throw;
          }

          connIter = inTun->m_children.begin();
          for(; connIter != inTun->m_children.end(); ++connIter) {
            if ((*connIter)->m_num == setInputNum) {
	      if (!(*connIter)->m_n->IsTunnel(POSSTUNOUT)) {
		Node *child = map[(*connIter)->m_n];
		child->ChangeInput1Way(inTun, setInputNum, node, conn->m_num);
	      }
            }
          }
          
          delete conn;
        }
	else if (conn->m_n->IsTunnel(POSSTUNOUT) &&
		 conn->m_n->m_poss == currPoss) {
	  if (conn->m_n->Child(0)->m_children.empty()) {
	    node->m_children.erase(node->m_children.begin()+i);
	    --i;
	    --size;
	    delete conn;
	    if (node->m_children.empty())
	      cleanupNodes.insert(node);
	  }
	  //else handled in above code when updating inputs that are set tun out
	}
	else {
	  if (node->IsTunnel(POSSTUNIN) &&
	      conn->m_n->IsTunnel(SETTUNIN)) {
	    throw;
	  }
	}
      }
    }
    


    NodeSetIter cleanupIter = cleanupNodes.begin();
    for(; cleanupIter != cleanupNodes.end(); ++cleanupIter) {
      Node *node = *cleanupIter;
      if (!node->m_children.empty())
	throw;
      newPoss->DeleteChildAndCleanUp(node,false,false,true);
    }

    iter = inliningPoss->m_inTuns.begin();
    for(; iter != inliningPoss->m_inTuns.end(); ++iter)  {
      Node *node = map[*iter];
      if (!node) {
        cout << "!node in dup\n";
        throw;
      }
      newPoss->m_inTuns.push_back(node);
      if (!(*iter)->m_poss) {
        cout << "!(*iter)->m_poss for " << *iter << endl;
        throw;
      }
    }
    
    iter = inliningPoss->m_outTuns.begin();
    for(; iter != inliningPoss->m_outTuns.end(); ++iter) {
      Node *node = map[*iter];
      if (!(*iter)->m_poss) {
        cout << "!(*iter)->m_poss\n";
        throw;
      }
      if (!node) {
        cout << "!node in dup\n";
        throw;
      }
      newPoss->m_outTuns.push_back(node);
    }

    iter = newPoss->m_outTuns.begin();
    for(; iter != newPoss->m_outTuns.end(); ++iter) {
      Node *node = *iter;
      for(unsigned int i = 0; i < node->m_inputs.size(); ++i) {
	Node *in = node->Input(i);
	if (tunnelNumMap.find(in) != tunnelNumMap.end())
	  throw;
	else if (in->m_poss != node->m_poss)
	  throw;	  
      }
    }

    iter = newPoss->m_possNodes.begin();
    for(; iter != newPoss->m_possNodes.end(); ++iter) {
      Node *node = *iter;
      for(unsigned int i = 0; i < node->m_children.size(); ++i) {
	Node *child = node->Child(i);
	if (node->m_poss != child->m_poss) {
	  if (node->IsTunnel() && child->IsTunnel()) {
	    Tunnel *nodeTun = (Tunnel*)node;
	    Tunnel *childTun = (Tunnel*)child;
	    if ((nodeTun->m_tunType == POSSTUNOUT && childTun->m_tunType == SETTUNOUT))
	      continue;
	    //the following can happen if currPoss has a set
	    if ((nodeTun->m_tunType == SETTUNIN && childTun->m_tunType == POSSTUNIN))
	      continue;
	  }
	  cout << node << " " << node->GetType() << endl;
	  cout << "child " << i << " is " << child << " " << child->GetType() << endl;
	  throw;
	}
      }
    }
    
    newPoss->m_transVec = inliningPoss->m_transVec;
    newPoss->m_transVec.insert(newPoss->m_transVec.end(),currPoss->m_transVec.begin(),currPoss->m_transVec.end());
    AddPossToMMap(newPosses, newPoss, newPoss->GetHash());
  }
  
  for (unsigned int i = 0; i < m_inTuns.size(); ++i) {
    InTun(i)->RemoveChild(inliningPoss->InTun(i),0);
  }
  for (unsigned int i = 0; i < m_outTuns.size(); ++i) {
    OutTun(i)->RemoveInput(inliningPoss->OutTun(i),0);
  }
  delete inliningPoss;
}

void RealPSet::Cull(CullFunction cullFunc)
{
  PossMMapIter iter = m_posses.begin();
  GraphNum i = 0;
  if (iter == m_posses.end()){
    cout << "starting with nothing\n";
    throw;
  }
  while (iter != m_posses.end()) {
    Poss *poss = (*iter).second;
    bool cullIfPossible, doNotCull;
    cullFunc(poss, cullIfPossible, doNotCull);
    if (cullIfPossible && !doNotCull) {
      RemoveAndDeletePoss(poss, false);
      m_posses.erase(iter);
      iter = m_posses.begin();
      for (GraphNum j = 0; j < i; ++j)
	++iter;
    }
    else {
      ++i;
      ++iter;
    }
  }
  if (m_posses.size() == 0) {
    cout << "Ran out of posses\n";
    throw;
  }
}


void RealPSet::FlattenCore(ofstream &out) const
{
  unsigned int size = m_posses.size();
  WRITE(size);
  PossMMapConstIter iter2 = m_posses.begin();
  for(; iter2 != m_posses.end(); ++iter2) {
    WRITE((*iter2).first);
    WRITE((*iter2).second);
  }
  iter2 = m_posses.begin();
  for(; iter2 != m_posses.end(); ++iter2)
    (*iter2).second->Flatten(out);
  WRITE(END);
}

void RealPSet::UnflattenCore(ifstream &in, SaveInfo &info)
{
  unsigned int size;
  READ(size);
  for(GraphNum i = 0; i < size; ++i) {
    Poss *newPoss = new Poss;
    Poss *oldPoss;
    size_t hash;
    READ(hash);
    READ(oldPoss);
    (*(info.possMap))[oldPoss] = newPoss;
    m_posses.insert(PossMMapPair(hash,newPoss));
  }
  PossMMapIter iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter) {
    (*iter).second->Unflatten(in, info);
  }
  char tmp;
  READ(tmp);
  if (tmp != END)
    throw;
  iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter) {
    if ((*iter).second->GetHash() != (*iter).first) {
      cout << "not same hash while reading\n";
      throw;
    }
  }
}


void RealPSet::BuildDataTypeCache()
{
  PossMMapIter iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter)
    (*iter).second->BuildDataTypeCache();
}

void RealPSet::ClearDataTypeCache()
{
  PossMMapIter iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter)
    (*iter).second->ClearDataTypeCache();
}


#if DOBLIS
bool RealPSet::RemoveParallelization(Comm comm)
{
#if DOBLIS
  if (IsLoop()) {
    Loop *loop = (Loop*)this;
    if (loop->IsParallel()) {
      if (comm == CORECOMM) {
        return true;
      }
      else if (!CommAllowedWithin(comm, loop->m_comm)) {
        return true;
      }
    }
  }
#endif
  
  GraphNum i = 0;
  PossMMapIter iter = m_posses.begin();
  while (iter != m_posses.end()) {
    bool found = false;
    Poss *poss = (*iter).second;
    PSetVecIter setIter = poss->m_sets.begin();
    for(; !found && setIter != poss->m_sets.end(); ++setIter) {
      PSet *set = *setIter;
      if (set->RemoveParallelization(comm)) {
        found = true;
      }
    }
    if (!found) {
      NodeVecIter nodeIter = poss->m_possNodes.begin();
      for(; !found && nodeIter != poss->m_possNodes.end(); ++nodeIter) {
        if ((*nodeIter)->IsParallel()) {
          Comm parComm = (*nodeIter)->ParallelComm();
          if ((comm == CORECOMM) || !CommAllowedWithin(comm, parComm)) {
            if ((*nodeIter)->RemoveParallelization())
              found = true;
          }
        }
      }
    }
    if (found) {
      if (m_posses.size() <= 1) {
        return true;
      }
      RemoveAndDeletePoss(poss, true);
      iter = m_posses.begin();
      for(GraphNum j = 0; j < i; ++j)
	++iter;
    }
    else {
      ++iter;
      ++i;
    }
  }
  return false;
}

Comm PSet::ParallelismWithinCurrentPosses() const
{
#if DOBLIS
  Comm comm = CORECOMM;
  if (IsLoop()) {
    const Loop *loop = (Loop*)this;
    if (loop->m_comm != CORECOMM) {
      return loop->m_comm;
    }
  }
  const Poss *currPoss = GetCurrPoss();
  PSetVecConstIter iter = currPoss->m_sets.begin();
  for(; iter != currPoss->m_sets.end(); ++iter) {
    const PSet *set = *iter;
    comm = MaxComm(comm,set->ParallelismWithinCurrentPosses());
  }
  NodeVecConstIter iter2 = currPoss->m_possNodes.begin();
  for(; iter2 != currPoss->m_possNodes.end(); ++iter2) {
    const Node *node = *iter2;
    if (node->IsParallel())
      comm = MaxComm(comm, node->ParallelComm());
  }
  return comm;
#else
  return CORECOMM;
#endif
}
#endif //DOBLIS

const string& RealPSet::GetFunctionalityString() const
{
  if (m_functionality.empty()) {
    cout << m_posses.size() << endl;
    throw;
  }
  else
    return m_functionality;
}


ShadowPSet* RealPSet::GetNewShadow()
{
  ShadowPSet *shadow = new ShadowPSet;
  shadow->m_realPSet = this;
  if (shadow->IsLoop() != IsLoop())
    throw;
  m_shadows.push_back(shadow);
  return shadow;
}

void RealPSet::RemoveShadow(ShadowPSet *shadow)
{
  PSetVecIter iter = m_shadows.begin();
  for (; iter != m_shadows.end(); ++iter) {
    if (*iter == shadow) {
      m_shadows.erase(iter);
      return;
    }
  }
  throw;
}

void RealPSet::SetInTunsAsPrinted()
{
  NodeVecIter iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter)
    (*iter)->SetPrinted();  
}


ShadowPSet* RealPSet::GetNewShadowDup(Poss *poss)
{
  ShadowPSet *shadow = GetNewShadow();
  if (IsLoop() != shadow->IsLoop())
    throw;

  poss->AddPSet(shadow, true);
  NodeVecIter iter = m_inTuns.begin();
  for( ; iter != m_inTuns.end(); ++iter) {
    Tunnel *tun = (Tunnel*)((*iter)->GetNewInst());
    tun->Duplicate(*iter,true,false);
    tun->m_pset = shadow;
    shadow->m_inTuns.push_back(tun);
    poss->AddNode(tun);
  }
  iter = m_outTuns.begin();
  for( ; iter != m_outTuns.end(); ++iter) {
    Tunnel *tun = (Tunnel*)((*iter)->GetNewInst());
    tun->Duplicate(*iter,true,false);
    tun->m_pset = shadow;
    shadow->m_outTuns.push_back(tun);
    poss->AddNode(tun);
  }
  return shadow;
}


RealPSet* RealPSet::HasMergedWith(RealPSet *set, bool checkOtherOrder)
{
  PSetMapIter iter = m_mergeMap.find(set);
  if (iter != m_mergeMap.end())
    return iter->second;

  if (checkOtherOrder)
    return set->HasMergedWith(this,false);
  else
    return NULL;
}

bool RealPSet::RemoveLoops(bool *doneSomething)
{
  PossMMapIter iter2 = m_posses.begin();
  while (iter2 != m_posses.end()) {
    Poss *poss = iter2->second;
    bool doneSomething2 = false;
    if (poss->RemoveLoops(&doneSomething2)) {
      *doneSomething = true;
      if (m_posses.size() == 1)
	return true;
      RemoveAndDeletePoss(poss, true);
      iter2 = m_posses.begin();
    }
    else if (doneSomething2) {
      *doneSomething = true;
      if (poss->GetHash() != iter2->first) {
	m_posses.erase(iter2);
	m_posses.insert(PossMMapPair(poss->GetHash(),poss));
	iter2 = m_posses.begin();
      }
      else
	++iter2;
    }
    else {
      ++iter2;
    }
  }
  return m_posses.empty();
}

void RealPSet::SetDeletingRecursively()
{
  m_flags |= SETISDELETINGFLAG;
  PossMMapIter iter = m_posses.begin();
  for( ; iter != m_posses.end(); ++iter) {
    iter->second->SetDeletingRecursively();
  }

}

void RealPSet::ClearDeletingRecursively()
{
  m_flags &= ~SETISDELETINGFLAG;
  PossMMapIter iter = m_posses.begin();
  for( ; iter != m_posses.end(); ++iter) {
    iter->second->ClearDeletingRecursively();
  }
}
