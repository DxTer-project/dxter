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
#include "pset.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "loopSupport.h"

extern unsigned int M_phase;

PSet::PSet()
  : m_isTopLevel(false), m_ownerPoss(NULL)
{
}

PSet::PSet(Poss *poss)
  : m_isTopLevel(false), m_ownerPoss(NULL)
{
  //Make single tunnels with multiple inputs/outputs into individual tunnels
  //Poss mergin with multiple intput/output tunnels is very buggy
  poss->ExpandTunnels();
  
  //Go through the input tunnels of the poss, create a set tunnel,
  // change the inputs from connecting to poss tunnels to set tunnels
  for(unsigned int i = 0; i < poss->m_inTuns.size(); ++i) {
    PossTunnel *possTun = (PossTunnel*)(poss->InTun(i));
    if (!possTun->IsPossTunnel(POSSTUNIN)) {
      cout << "bad poss tunnel\n";
      throw;
    }
    PossTunnel *setTun = possTun->GetSetTunnel();
    for(unsigned int j = 0; j < possTun->m_inputs.size(); ++j) {
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
    PossTunnel *possTun = (PossTunnel*)(poss->OutTun(i));
    PossTunnel *setTun = possTun->GetSetTunnel();
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

PSet::~PSet()
{
  PossMMapIter iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter) {
    delete (*iter).second;
  }
}

void PSet::AddPossesOrDispose(PossMMap &mmap, PossMMap *added)
{
  PossMMapIter newIter = mmap.begin();
  for( ; newIter != mmap.end(); ++newIter) {
    Poss *poss = (*newIter).second;
    poss->RemoveConnectionToSet();
    
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

void PSet::AddPoss(Poss *poss)
{
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
    Node *possTun = (PossTunnel*)(poss->InTun(i));
    if (!possTun->m_inputs.empty()) {
      cout << "(!possTun->m_inTuns.empty()\n";
      throw;
    }
    possTun->m_inputs.clear();
    possTun->AddInput(setTun,0);
  }
  
  for(unsigned int i = 0; i < m_outTuns.size(); ++i) {
    Node *setTun = OutTun(i);
    Node *possTun = (PossTunnel*)(poss->OutTun(i));
    if (!possTun->m_children.empty()) {
      cout << "!possTun->m_outTuns.empty()\n";
      cout << possTun->GetType() << " has child "
	   << possTun->Child(0)->GetType() << endl;
      cout << possTun << endl;
      throw;
    }
    possTun->m_children.clear();
    setTun->AddInput(possTun,0);
  }
  
  m_posses.insert(PossMMapPair(poss->GetHash(),poss));
  poss->m_pset = this;
}


bool PSet::operator==(const Poss &rhs) const
{
  cout << "PSet::operator==(const Poss &rhs) const Not defined\n";
  return false;
}

bool PSet::operator==(const PSet &rhs) const
{
  if (m_inTuns.size() != rhs.m_inTuns.size()
      || m_outTuns.size() != rhs.m_outTuns.size()
      || m_posses.size() != rhs.m_posses.size())
    return false;
  else {
    if (IsLoop()) {
      if (!rhs.IsLoop())
        return false;
      else {
        for (unsigned int i = 0; i < m_inTuns.size(); ++i) {
          const LoopTunnel *tun1 = (LoopTunnel*)(m_inTuns[i]);
          const LoopTunnel *tun2 = (LoopTunnel*)(rhs.m_inTuns[i]);
#if DOBLIS
          if (((Loop*)(tun1->m_pset))->m_comm != ((Loop*)(tun2->m_pset))->m_comm)
            return false;
#endif
          if (tun1->m_statTL != tun2->m_statTL
              || tun1->m_statBL != tun2->m_statBL
              || tun1->m_statTR != tun2->m_statTR
              || tun1->m_statBR != tun2->m_statBR)
            return false;
          if (tun1->GetNodeClass() == Split::GetClass()) {
            if (tun2->GetNodeClass() == Split::GetClass()) {
#if TWOD
              if (((Split*)tun1)->m_dir != ((Split*)tun2)->m_dir)
                return false;
#else
	      if (((Split*)tun1)->m_partDim != ((Split*)tun2)->m_partDim)
                return false;
#endif
            }
            else
              return false;
          }
          else if (tun2->GetNodeClass() == Split::GetClass()) {
            return false;
          }
        }
      }
    }
#if DOBLIS
    else if (IsCritSect()) {
      if (!rhs.IsCritSect())
	return false;
    }
#endif
    if ((*(m_posses.begin())).second->GetHash() != (*(rhs.m_posses.begin())).second->GetHash())
      return false;
    //BAM Really, instead of doing all of this comparison,
    //    it would be better to keep track of what PSets implement
    //    and compare that.  If they implement the same thing, they'll
    //    eventually (with enough iterations) be the same
    /*    for (unsigned int i = 0; i < m_posses.size(); ++i)
	  if (!(*(m_posses[i]) == *(rhs.m_posses[i])))
	  return false;*/
    return true;
  }
}

void PSet::Prop()
{
  if(m_hasProped)
    return;

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
  
  if (!m_isTopLevel && !m_ownerPoss) {
    cout << "no owner\n";
    throw;
  }




  for(unsigned int i = 0; i < m_inTuns.size(); ++i) {
    InTun(i)->Prop();
  }
  PossMMapIter iter;
  int j = 0;
#pragma omp parallel private(j,iter)
  {
    iter = m_posses.begin();
    j = 0;
    int size = m_posses.size();
    //BAM par
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
      (*iter).second->Prop();
    }
  }
  iter = m_posses.begin();
  while (iter != m_posses.end()) {
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
      
      iter = m_posses.begin();
    }
    else {
      ++iter;
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
  m_hasProped = true;
}

void PSet::Cull(Phase phase)
{
  PossMMapIter iter2 = m_posses.begin();
  for(; iter2 != m_posses.end(); ++iter2) 
    (*iter2).second->m_isSane = false;

  int j;
  PossMMapIter iter;
  
#pragma omp parallel private(j,iter)
  {
    iter = m_posses.begin();
    j = 0;
    int size = m_posses.size();
    //BAM par
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
      if (poss->m_isSane)
	throw;
      poss->m_isSane = true;
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

Node* PSet::InTun(unsigned int num) const
{
  if (num >= m_inTuns.size())
    throw;
  return m_inTuns[num];
}

Node* PSet::OutTun(unsigned int num) const
{
  if (num >= m_outTuns.size())
    throw;
  return m_outTuns[num];
}

void PSet::RemoveAndDeletePoss(Poss *poss, bool removeFromMyList)
{
  if (m_posses.size() <= 1) {
    if (m_posses.size()) {
      poss->ForcePrint();
    }
    throw;
  }
  if (removeFromMyList) {
    cout << "m_posses has " << m_posses.size() << endl;
    PossMMapIter possIter = m_posses.begin();
    bool found = false;
    for(; !found && possIter != m_posses.end(); ++possIter) {
      cout << "another iteration\n";
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
  poss->m_pset = (PSet*)(0xDEADBEEF);
  delete poss;
}

void PSet::ClearBeforeProp()
{
  m_hasProped = false;
  for (unsigned int i = 0; i < m_inTuns.size(); ++i)
    InTun(i)->ClearBeforeProp();
  for (unsigned int i = 0; i < m_outTuns.size(); ++i)
    OutTun(i)->ClearBeforeProp();
  
  PossMMapIter iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter)
    (*iter).second->ClearBeforeProp();
}


bool PSet::TakeIter(const TransMap &transMap,
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
  
  int size = m_posses.size();
  
  if (size > 1) {
    PossMMap mmap;
    PossMMapIter iter;
    int j = 0;

#pragma omp parallel private(j,iter)
    {
      iter = m_posses.begin();
      j = 0;
      int size = m_posses.size();
      //BAM par
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
      if (mmap.size() > 50)
        cout << "\t\tAdding " << mmap.size() << " posses\n";
      AddPossesOrDispose(mmap, &actuallyAdded);
      if (mmap.size() > 50)
        cout << "\t\tDone adding ( " << actuallyAdded.size() << " actually added )\n";
      PossMMapIter added = actuallyAdded.begin();
      for(; added != actuallyAdded.end(); ++added) {
        (*added).second->BuildDataTypeCache();
      }
    }
  }
  else {
    Poss *poss = (*(m_posses.begin())).second;
    PossMMap newPosses;
    newOne = poss->TakeIter(transMap, simplifiers, newPosses);
    if (newPosses.size()) {
      if (newPosses.size() > 50)
        cout << "\t\tAdding " << newPosses.size() << " posses\n";
      AddPossesOrDispose(newPosses, &actuallyAdded);
      if (newPosses.size() > 50)
        cout << "\t\tDone adding\n";
      PossMMapIter added = actuallyAdded.begin();
      for(; added != actuallyAdded.end(); ++added) {
        (*added).second->BuildDataTypeCache();
      }
    }
  }
  return newOne;
}

bool PSet::GlobalSimplification(const TransMap &globalSimplifiers, const TransMap &simplifiers)
{
  bool didSomething = false;
  PossMMapIter iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter) {
    if ((*iter).second->GlobalSimplification(globalSimplifiers, simplifiers)) {
      didSomething = true;
      if ((*iter).first != (*iter).second->GetHash()) {
	cout << "GlobalSimplification different hash\n";
	cout << (*iter).first  << " vs " << (*iter).second->GetHash() << endl;
	throw;
      }
    }
  }
  return didSomething;
}


void PSet::Duplicate(const PSet *orig, NodeMap &map, bool possMerging)
{
  m_isTopLevel = orig->m_isTopLevel;
  NodeVecConstIter iter  = orig->m_inTuns.begin();
  for (; iter != orig->m_inTuns.end(); ++iter) {
    PossTunnel *tun = (PossTunnel*)(map[*iter]);
    //expect set tunnel for this set to be duplicated
    // as part of the owning poss's duplication
    if (!tun)
      throw;
    m_inTuns.push_back(tun);
    tun->m_pset = this;
  }
  iter  = orig->m_outTuns.begin();
  for (; iter != orig->m_outTuns.end(); ++iter) {
    PossTunnel *tun = (PossTunnel*)(map[*iter]);
    if (!tun)
      throw;
    m_outTuns.push_back(tun);
    tun->m_pset = this;
  }
  
  PossMMapConstIter iter2 = orig->m_posses.begin();
  for( ; iter2 != orig->m_posses.end(); ++iter2) {
    const Poss *oldPoss = (*iter2).second;
    Poss *newPoss = new Poss;
    newPoss->Duplicate(oldPoss, map, possMerging);
    m_posses.insert(PossMMapPair(newPoss->GetHash(),newPoss));
    newPoss->m_pset = this;
  }
  
}

void PSet::PatchAfterDuplicate(NodeMap &map)
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

void PSet::CombineAndRemoveTunnels()
{
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
	    ClassType type1 = setInput1->GetNodeClass();
	    ClassType type2 = setInput2->GetNodeClass();
	    if ((type1 == Split::GetClass() && type2 != Split::GetClass())
		|| (type1 != Split::GetClass() && type2 == Split::GetClass()))
	      {
		continue;
	      }
	    if (type1 == Split::GetClass()) {
#if TWOD
	      if (((Split*)setInput1)->m_dir != ((Split*)setInput2)->m_dir)
		continue;
#else
	      if (((Split*)setInput1)->m_partDim != ((Split*)setInput2)->m_partDim)
		continue;
#endif
	    }
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
        if (!child->IsPossTunnel()) {
          cout << "!child->IsPossTunnel()\n";
          throw;
        }
        child->m_poss->DeleteNode(child);
        delete *inChildIter;
      }
      tun->m_children.clear();
      delete tun;
      m_inTuns.erase(m_inTuns.begin()+i);
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
        --i;
        continue;
      }
    }
    if (tun->GetNodeClass() == Combine::GetClass()) {
      for(unsigned int j = i+1; j < m_outTuns.size(); ++j) {
        Node *tun2 = OutTun(j);
        if (tun2->GetNodeClass() == Combine::GetClass()) {
          if (!((LoopTunnel*)tun)->IsConst()) {
            continue;
          }
          Node *in1 = tun->Input(0);
          Node *in2 = tun2->Input(0);
          if (in1->Input(in1->m_inputs.size()-1) !=
              in2->Input(in2->m_inputs.size()-1)) {
            continue;
          }
          Split *split1 = (Split*)tun;
          Split *split2 = (Split*)tun2;
#if TWOD
          if (split1->m_dir != split2->m_dir) {
            continue;
          }
#else
          if (split1->m_partDim != split2->m_partDim) {
            continue;
          }
#endif
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
  }
  
}

void PSet::Simplify(const TransMap &simplifiers)
{
  //BAM par
  PossMMapIter iter;
  int j = 0;
#pragma omp parallel private(j,iter)
  {
    iter = m_posses.begin();
    j = 0;
    int size = m_posses.size();
    //BAM par
#pragma omp for schedule(static) 
    for (int i = 0; i < size; ++i) {
      while (j < i) {
	++iter;
	++j;
      }
      Poss *poss = (*iter).second;
      poss->Simplify(simplifiers);
    }
  }
}

// void PSet::RemoveDups()
// {
// #ifdef _OPENMP
//   static omp_lock_t lock;
//   static bool inited = false;
//   if (!inited) {
//     omp_init_lock(&lock);
//     inited = true;
//   }
// #endif
  
//   for (unsigned int i = 0; i < m_posses.size(); ++i) {
//     Poss *poss = m_posses[i];
//     std::set<unsigned int> dups;
//     int size = m_posses.size();
// #pragma omp parallel
//     {
// #pragma omp for
//       for (int j = i+1; j < size; ++j) {
//         if (*poss == *m_posses[j]) {
// #ifdef _OPENMP
//           omp_set_lock(&lock);
// #endif
//           dups.insert(j);
// #ifdef _OPENMP
//           omp_unset_lock(&lock);
// #endif
//         }
//       }
//     }
//     unsigned int count = 0;
//     std::set<unsigned int>::iterator iter = dups.begin();
//     for(; iter != dups.end(); ++iter, ++count) {
//       RemoveAndDeletePoss(m_posses[*iter - count], false);
//       m_posses.erase(m_posses.begin() + *iter - count);
//     }
//   }
// }

void PSet::ClearFullyExpanded()
{
  PossMMapIter iter = m_posses.begin();
  for( ;iter != m_posses.end(); ++iter)
    (*iter).second->ClearFullyExpanded();
}

bool FoundPossUp(Node *node, const PSet *set, NodeVec &queue)
{
  NodeVecIter checkIter = queue.begin();
  for(; checkIter != queue.end(); ++checkIter) {
    if (*checkIter == node) {
      cout << "recursion on node " << node << " " << node->GetNodeClass() << endl;
      throw;
    }
  }
  queue.push_back(node);
  if (node->IsPossTunnel(POSSTUNOUT) || node->IsPossTunnel(POSSTUNIN)) {
    queue.pop_back();
    return false;
  }
  else if (node->IsPossTunnel(SETTUNOUT)) {
    const PossTunnel *tunOut = (PossTunnel*)node;
    if (tunOut->m_pset == set) {
      queue.pop_back();
      return true;
    }
    else {
      const PSet *foundSet = tunOut->m_pset;
      NodeVecConstIter iter = foundSet->m_inTuns.begin();
      for(; iter != foundSet->m_inTuns.end(); ++iter) {
        if (FoundPossUp(*iter, set, queue)) {
          queue.pop_back();
          return true;
        }
      }
    }
  }
  else {
    NodeConnVecConstIter iter = node->m_inputs.begin();
    for (; iter != node->m_inputs.end(); ++iter) {
      if (FoundPossUp((*iter)->m_n, set, queue)) {
        queue.pop_back();
        return true;
      }
    }
  }
  queue.pop_back();
  return false;
}

bool NothingBetween(const PSet *left, const PSet *right)
{
  NodeVecConstIter iter = right->m_inTuns.begin();
  for(; iter != right->m_inTuns.end(); ++iter) {
    Node *input = *iter;
    NodeConnVecConstIter iter2 = input->m_inputs.begin();
    for (; iter2 != input->m_inputs.end(); ++iter2) {
      if ((*iter2)->m_n->IsPossTunnel(SETTUNOUT)) {
        PossTunnel *tunOut = (PossTunnel*)((*iter2)->m_n);
        if (tunOut->m_pset != left) {
          NodeVec queue;
          if (FoundPossUp(tunOut,left,queue))
            return false;
          if (!queue.empty()) {
            cout << "queue not empty!\n";
            throw;
          }
        }
      }
      else {
        NodeVec queue;
        if (FoundPossUp((*iter2)->m_n,left, queue))
          return false;
        if (!queue.empty()) {
          cout << "queue not empty!\n";
          throw;
        }
      }
    }
  }
  return true;
}

bool ShouldMerge(const PSet *set1, const PSet *set2)
{
  unsigned int i, j, k;
  for(i = 0; i < set1->m_inTuns.size(); ++i) {
    const Node *in = set1->m_inTuns[i];
    for(j = 0; j < in->m_inputs.size(); ++j) {
      const Node *inInput = in->Input(j);
      if (inInput->IsPossTunnel()) {
        if (((PossTunnel*)inInput)->m_pset == set2)
          return true;
      }
      for(k = 0; k < inInput->m_children.size(); ++k) {
        const Node *child = inInput->Child(k);
        if (child->IsPossTunnel()) {
          if (((PossTunnel*)child)->m_pset == set2)
            return true;
        }
      }
    }
  }
  for(i = 0; i < set1->m_outTuns.size(); ++i) {
    const Node *out = set1->m_outTuns[i];
    for(j = 0; j < out->m_children.size(); ++j) {
      const Node *child = out->Child(j);
      if (child->IsPossTunnel()) {
        if (((PossTunnel*)child)->m_pset == set2)
          return true;
      }
    }
  }
  return false;
}

bool PSet::CanMerge(PSet *pset) const
{
  bool nothingBetween = NothingBetween(this, pset) && NothingBetween(pset, this);
  if (!nothingBetween)
    return false;
  return ShouldMerge(this, pset);
}

bool PSet::MergePosses(const TransMap &simplifiers, CullFunction cullFunc)
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
    //BAM par
    PossMMapIter iter;
    int j = 0;
#pragma omp parallel private(j,iter)
    {
      iter = m_posses.begin();
	j = 0;
	int size = m_posses.size();
	//BAM par
#pragma omp for schedule(static) 
      for (int i = 0; i < size; ++i) {
	while (j < i) {
	  ++iter;
	  ++j;
	}
	Poss *poss = (*iter).second;
        PossMMap newPosses;
        if (poss->MergePosses(newPosses, simplifiers, cullFunc)) {
#ifdef _OPENMP
          omp_set_lock(&lock);
#endif
          didMerge = true;
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
          InlinePoss(poss, newPosses);
	  m_posses.erase(iter);
	  iter = m_posses.begin();
          PossMMapIter mapIter = newPosses.begin();
          for(; mapIter != newPosses.end(); ++mapIter)
            (*mapIter).second->BuildDataTypeCache();
          AddPossesOrDispose(newPosses);
        }
      else {
	++iter;
      }
    }
  }
  return didMerge;
}

void PSet::FormSets(unsigned int phase)
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

void PSet::ClearPrinted()
{
  (*m_currPoss).second->ClearPrinted();
}

unsigned int PSet::TotalCount() const
{
  unsigned int tot = 0;
  PossMMapConstIter iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter)
    tot += (*iter).second->TotalCount();
  return tot;
}

void PSet::InlinePoss(Poss *inliningPoss, PossMMap &newPosses)
{
  PSet *pset = inliningPoss->m_sets[0];
  
  //  cout << "inlining " << inliningPoss << endl;
  
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
    //    cout << "currPoss " << currPoss << endl;
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
      if (!(*iter)->IsPossTunnel(SETTUNIN) && !(*iter)->IsPossTunnel(SETTUNOUT)) {
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
      if (!(*iter)->IsPossTunnel(POSSTUNIN) && !(*iter)->IsPossTunnel(POSSTUNOUT)) {
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
        PSet *newSet = (*setIter)->GetNewInst();
        newSet->Duplicate(*setIter, map, true);
        newPoss->m_sets.push_back(newSet);
        newSet->m_ownerPoss = newPoss;
      }
      setIter = newPoss->m_sets.begin();
      for(; setIter != newPoss->m_sets.end(); ++setIter) {
        (*setIter)->PatchAfterDuplicate(map);
      }
    }
    
    newPoss->PatchAfterDuplicate(map);
    
    iter = newPoss->m_possNodes.begin();
    for(; iter != newPoss->m_possNodes.end(); ++iter) {
      
      Node *node = *iter;
      NodeConnVecIter connIter = node->m_inputs.begin();
      for(; connIter != node->m_inputs.end(); ++connIter) {
        NodeConn *conn = *connIter;
        NodeIntMapIter mapIter = tunnelNumMap.find(conn->m_n);
        if (mapIter != tunnelNumMap.end()) {
          Node *newParent = map[currPoss->OutTun(mapIter->second)->Input(conn->m_num)];
          unsigned int newNum = currPoss->OutTun(mapIter->second)->InputConnNum(conn->m_num);
          for(unsigned int i = 0; i < newParent->m_children.size(); ++i) {
            if (newParent->Child(i)->IsPossTunnel(POSSTUNOUT)) {
              delete newParent->m_children[i];
              newParent->m_children.erase(newParent->m_children.begin()+i);
              --i;
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
              Node *child = map[(*connIter)->m_n];
              child->ChangeInput1Way(inTun, setInputNum, node, conn->m_num);
            }
          }
          
          delete conn;
        }
      }
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

void PSet::RemoveInTun(Node *tun)
{
  NodeVecIter iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    if (*iter == tun) {
      m_inTuns.erase(iter);
      return;
    }
  }
  throw;
}

void PSet::RemoveOutTun(Node *tun)
{
  NodeVecIter iter = m_outTuns.begin();
  for(; iter != m_outTuns.end(); ++iter) {
    if (*iter == tun) {
      m_outTuns.erase(iter);
      return;
    }
  }
  throw;
}

void PSet::ClearCurrPoss()
{
  m_currPoss = m_posses.begin();
  PossMMapIter iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter) {
    (*iter).second->ClearCurrPoss();
  }
}

bool PSet::IncrementCurrPoss()
{
  if (m_currPoss == m_posses.end())
    throw;
  if ((*m_currPoss).second->IncrementCurrPoss()) {
    m_currPoss++;
    if (m_currPoss == m_posses.end()) {
      m_currPoss = m_posses.begin();
      return true;
    }
  }
  return false;
}

Cost PSet::EvalCurrPoss(TransConstVec &transList)
{
  //Any changes should be reflected in Loop::EvalCurrPoss()
  if (m_currPoss == m_posses.end()) {
    throw;
  }
  return (*m_currPoss).second->EvalCurr(transList);
}

Cost PSet::EvalAndSetBest()
{
  Cost optCost = -1;
  PossMMapIter best;
  PossMMapIter iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter) {
    Poss *poss = (*iter).second;
    Cost tmp = poss->EvalAndSetBest();
    if (optCost < 0 || tmp < optCost) {
      optCost = tmp;
      best = iter;
    }
  }
  m_currPoss = best;
  return optCost;
}

void PSet::PrintCurrPoss(IndStream &out, unsigned int &graphNum)
{
  out.Indent();
  *out << "//**** (out of " << m_posses.size() << ")\n";
  
  if (m_currPoss == m_posses.end()) {
    throw;
  }
  ++out;
  (*m_currPoss).second->Print(out, graphNum);
  --out;
  
  out.Indent();
  *out << "//****\n";
}

Poss* PSet::GetCurrPoss() const
{
  if (m_currPoss == m_posses.end()) {
    throw;
  }
  return (*m_currPoss).second;
}


void PSet::Cull(CullFunction cullFunc)
{
  PossMMapIter iter = m_posses.begin();
  unsigned int i = 0;
  while (iter != m_posses.end()) {
    Poss *poss = (*iter).second;
    bool cullIfPossible, doNotCull;
    cullFunc(poss, cullIfPossible, doNotCull);
    if (cullIfPossible && !doNotCull) {
      RemoveAndDeletePoss(poss, false);
      m_posses.erase(iter);
      iter = m_posses.begin();
      for (unsigned int j = 0; j < i; ++j)
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

void PSet::GetCurrTransVec(TransVec &transVec) const
{
  if (m_currPoss == m_posses.end()) {
    throw;
  }
  (*m_currPoss).second->GetCurrTransVec(transVec);
}


void PSet::FormSetAround()
{
  Poss *owner = m_ownerPoss;
  PSet *newSet = new PSet;
  Poss *newPoss = new Poss;
  
  newPoss->m_pset = newSet;
  newSet->m_posses.insert(PossMMapPair(newPoss->GetHash(),newPoss));
  
  newSet->m_ownerPoss = owner;
  owner->m_sets.push_back(newSet);
  
  
  
  NodeVecIter iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    PossTunnel *tun = (PossTunnel*)(*iter);
    PossTunnel *newSetTun = new PossTunnel(SETTUNIN);
    PossTunnel *newPossTun = new PossTunnel (POSSTUNIN);
    newPoss->AddNode(newPossTun);
    newPoss->m_inTuns.push_back(newPossTun);
    newPossTun->AddInput(newSetTun,0);
    if (tun->m_inputs.size() != 1) {
      throw;
    }
    NodeConn *in = tun->m_inputs[0];
    newSetTun->AddInput(in->m_n, in->m_num);
    tun->ChangeInput2Way(in->m_n, in->m_num, newPossTun, 0);
    newSet->m_inTuns.push_back(newSetTun);
    newSetTun->m_pset = newSet;
    owner->AddNode(newSetTun);
    owner->RemoveFromGraphNodes(tun);
    tun->m_poss=NULL;
    newPoss->AddNode(tun);
  }
  
  
  iter = m_outTuns.begin();
  for(; iter != m_outTuns.end(); ++iter) {
    PossTunnel *tun = (PossTunnel*)(*iter);
    PossTunnel *newSetTun = new PossTunnel(SETTUNOUT);
    PossTunnel *newPossTun = new PossTunnel (POSSTUNOUT);
    newPoss->m_outTuns.push_back(newPossTun);
    newPoss->AddNode(newPossTun);
    newSetTun->AddInput(newPossTun,0);
    tun->RedirectAllChildren(newSetTun);
    newPossTun->AddInput(tun,0);
    
    newSet->m_outTuns.push_back(newSetTun);
    newSetTun->m_pset = newSet;
    owner->AddNode(newSetTun);
    owner->RemoveFromGraphNodes(tun);
    tun->m_poss=NULL;
    newPoss->AddNode(tun);
  }
  
  bool found = false;
  PSetVecIter iter2 = owner->m_sets.begin();
  for(; !found && iter2 != owner->m_sets.end(); ++iter2) {
    if (*iter2 == this) {
      owner->m_sets.erase(iter2);
      found = true;
    }
  }
  if (!found)
    throw;
  newPoss->m_sets.push_back(this);
  m_ownerPoss = newPoss;
}


void PSet::Flatten(ofstream &out) const
{
  WRITE(START);
  WRITE(m_isTopLevel);
  FlattenCore(out);
  unsigned int size;
  if (m_isTopLevel) {
    FullyFlatten(m_inTuns, out);
    FullyFlatten(m_outTuns, out);
  }
  else {
    size = m_inTuns.size();
    WRITE(size);
    NodeVecConstIter iter = m_inTuns.begin();
    for(; iter != m_inTuns.end(); ++iter)
      WRITE(*iter);
    WRITE(END);
    size = m_outTuns.size();
    WRITE(size);
    iter = m_outTuns.begin();
    for(; iter != m_outTuns.end(); ++iter) {
      WRITE(*iter);
    }
  }
  WRITE(END);
  WRITE(END);
  WRITE(m_ownerPoss);
  size = m_posses.size();
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

void PSet::Unflatten(ifstream &in, SaveInfo &info)
{
  char tmp;
  READ(tmp);
  if (tmp != START)
    throw;
  READ(m_isTopLevel);
  UnflattenCore(in,info);
  unsigned int size;
  if (m_isTopLevel) {
    FullyUnflatten(m_inTuns, in, info);
    FullyUnflatten(m_outTuns, in, info);
  }
  else {
    READ(size);
    for(unsigned int i = 0; i < size; ++i) {
      Node *tun;
      READ(tun);
      Swap(&tun,info.nodeMap);
      m_inTuns.push_back(tun);
    }
    READ(tmp);
    if (tmp != END)
      throw;
    READ(size);
    for(unsigned int i = 0; i < size; ++i) {
      Node *tun;
      READ(tun);
      Swap(&tun,info.nodeMap);
      m_outTuns.push_back(tun);
    }
  }
  READ(tmp);
  if (tmp != END)
    throw;
  READ(tmp);
  if (tmp != END)
    throw;
  READ(m_ownerPoss);
  if (!m_isTopLevel)
    Swap(&m_ownerPoss, info.possMap);
  READ(size);
  for(unsigned int i = 0; i < size; ++i) {
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
  READ(tmp);
  if (tmp != END)
    throw;
  if (m_isTopLevel) {
    NodeVecIter iter2 = m_inTuns.begin();
    for(; iter2 != m_inTuns.end(); ++iter2)
      (*iter2)->PatchAfterDuplicate(*(info.nodeMap));
    iter2 = m_outTuns.begin();
    for(; iter2 != m_outTuns.end(); ++iter2)
      (*iter2)->PatchAfterDuplicate(*(info.nodeMap));
  }
  iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter) {
    if ((*iter).second->GetHash() != (*iter).first) {
      cout << "not same hash while reading\n";
      throw;
    }
  }
}

void PSet::AddCurrPossVars(VarSet &set) const
{
  GetCurrPoss()->AddCurrPossVars(set);
}

void PSet::BuildDataTypeCache()
{
  PossMMapIter iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter)
    (*iter).second->BuildDataTypeCache();
}

void PSet::ClearDataTypeCache()
{
  PossMMapIter iter = m_posses.begin();
  for(; iter != m_posses.end(); ++iter)
    (*iter).second->ClearDataTypeCache();
}

bool PSet::CanPrint() const
{
  NodeVecConstIter iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    const Node *in = *iter;
    if (!in->CanPrintCode())
      return false;
  }
  return true;
}

#if DOBLIS
bool PSet::RemoveParallelization(Comm comm)
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
  
  unsigned int i = 0;
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
      for(unsigned int j = 0; j < i; ++j)
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
