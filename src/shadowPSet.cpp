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
  : m_isTopLevel(false), m_ownerPoss(NULL), m_functionality()
{
}

PSet::PSet(Poss *poss)
  : m_isTopLevel(false), m_ownerPoss(NULL)
{
  m_functionality = poss->GetFunctionalityString();

  if (m_functionality.empty()) {
    cout << "starting PSet without functionality\n";
    throw;
  }
  if (IsLoop()) {
    Loop *loop = (Loop*)this;
    m_functionality += (char)(loop->m_bsSize);
  }
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

bool PSet::operator==(const PSet &rhs) const
{
  if (m_inTuns.size() != rhs.m_inTuns.size()
      || m_outTuns.size() != rhs.m_outTuns.size())
    return false;
  if (GetFunctionalityString() != rhs.GetFunctionalityString()) {
    return false;
  }
  else {
    if (IsLoop()) {
      if (!rhs.IsLoop())
        return false;
      else {
#if TWOD
	if (((Loop*)this)->GetDimName() != ((Loop*)(&rhs))->GetDimName())
	  return false;
#endif
	if (((Loop*)this)->m_bsSize != ((Loop*)(&rhs))->m_bsSize)
	  return false;
	if (BSSizeToSize(((Loop*)this)->m_bsSize) == 0) {
	  cout << "BS bs\n";
	  cout << ((Loop*)this)->m_bsSize << endl;
	  throw;
	}
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
      if (!rhs.IsCritSect())
	return false;
    }
#endif
    return true;
  }
}

void ShadowPSet::Prop()
{
  if(m_hasProped)
    return;

  if (!m_realPSet)
    throw;

  //BAM Par + check for > 1
  for (unsigned int i = 0; i < m_inTuns.size(); ++i) {
    Node *in = InTun(i);
    if (!in->IsPossTunnel(SHADOWTUNIN))
      throw;
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
    if (!out->IsPossTunnel(SHADOWTUNOUT))
      throw;
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

  m_realPSet->Prop();

  m_hasProped = true;
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


GraphNum ShadowPSet::TotalCount() const
{
  return m_realPSet->TotalCount();
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



void PSet::GetCurrTransVec(TransVec &transVec) const
{
  if (m_currPoss == m_posses.end()) {
    throw;
  }
  (*m_currPoss).second->GetCurrTransVec(transVec);
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


#if DOBLIS
bool PSet::RemoveParallelization(Comm comm)
{
  throw;
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

const string& PSet::GetFunctionalityString() const
{
  return m_realPSet->GetFunctionalityString();
}
