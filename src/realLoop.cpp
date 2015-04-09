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



#include "loopSupport.h"
#include "elemRedist.h"
#include <cmath>
#include <climits>
#include "pack.h"
#include "critSect.h"
#include "blis.h"

int RealLoop::M_currLabel = 0;

RealLoop::RealLoop()
: IntLoop<RealPSet>(), 
#if TWOD
 m_dim(BADDIM),
#endif
  m_type(UNKNOWNLOOP),
#if DOBLIS
 m_comm(CORECOMM),
#endif
  m_currIter(0)
{
  AssignNewLabel();
  m_bsSize = BadBS;
}

RealLoop::RealLoop(LoopType type)
: 
#if TWOD
  m_dim(BADDIM),
#endif
  m_type(type)
#if DOBLIS
, m_comm(CORECOMM)
#endif
  , m_currIter(0)

{
#if DOELEM
  if (m_type == ELEMLOOP)
    m_bsSize = ElemBS;
  else
#endif
    m_bsSize = BadBS;
  AssignNewLabel();
}

RealLoop::RealLoop(LoopType type, Poss *poss, BSSize bsSize)
  :
  IntLoop<RealPSet>(poss),
#if TWOD
 m_dim(BADDIM),
#endif
 m_type(type), m_bsSize(bsSize)
#if DOBLIS
, m_comm(CORECOMM)
#endif
  , m_currIter(0)
{
  unsigned int i;
  for(i = 0; i < poss->m_inTuns.size(); ++i) {
    LoopTunnel *tun = (LoopTunnel*)(poss->m_inTuns[i]);
    if (tun->IsSplit()) {
      SplitBase *split = (SplitBase*)tun;
      if (split->m_isControlTun)
        split->m_isControlTun = true;
    }
  }
  AssignNewLabel();
  
  bool foundControl = false;
  //  NodeVecIter iter = m_inTuns.begin();
  //  for(; iter != m_inTuns.end(); ++iter) {
  //    Node *in = *iter;
  for(auto in : m_inTuns) {
    if (!in->IsLoopTunnel()) {
      cout << "non loop tunnel on loop!\n";
      LOG_FAIL("replacement for throw call");
    }
    if (((LoopTunnel*)in)->IsSplit())
    {
      SplitBase *split = (SplitBase*)in;
      
      if (split->m_isControlTun) {
        if (foundControl) {
          LOG_FAIL("replacement for throw call");
	} else {
          foundControl = true;
	}
      }
    }
  }
  if (!foundControl) {
    cout << "ERROR: Could not find loop control\n";
    LOG_FAIL("replacement for throw call");
  }
}

void RealLoop::Duplicate(const BasePSet *orig, NodeMap &map, bool possMerging, bool useShadows)
{
  IntLoop<RealPSet>::Duplicate(orig,map,possMerging,useShadows);
  const RealLoop *loop = (RealLoop*)orig;
  m_label = loop->m_label;
  m_bsSize = loop->m_bsSize;
  m_type = loop->m_type;
#if DOBLIS
  m_comm = loop->m_comm;
#endif
#if TWOD
  m_dim = loop->m_dim;
#endif
  if (GetBS() == 0) {
    cout << "duplicating a loop with zero blocksize\n";
    LOG_FAIL("replacement for throw call");
  }
}

void RealLoop::AssignNewLabel()
{
  m_label.insert(M_currLabel++);
}

void RealLoop::SetBS(BSSize size)
{
  m_bsSize = size;
}

int RealLoop::GetBS() const
{
  return (int)(m_bsSize.GetSize());
}

void RealLoop::FlattenCore(ofstream &out) const
{
  WRITE(m_type);
  WRITE(m_bsSize);
#if DOBLIS
  WRITE(m_comm);
#endif
#if TWOD
  WRITE(m_dim);
#endif
  unsigned int size = m_label.size();
  WRITE(size);
  IntSetConstIter iter = m_label.begin();
  for(; iter != m_label.end(); ++iter)
    WRITE(*iter);
}


void RealLoop::UnflattenCore(ifstream &in, SaveInfo &info)
{
  READ(m_type);
  READ(m_bsSize);
#if DOBLIS
  READ(m_comm);
#endif
#if TWOD
  READ(m_dim);
#endif
  unsigned int size;
  READ(size);
  for(unsigned int i = 0; i < size; ++i) {
    int tmp;
    READ(tmp);
    m_label.insert(tmp);
  }
}

void RealLoop::FlattenStatic(ofstream &out)
{
  WRITE(M_currLabel);
}

void RealLoop::UnflattenStatic(ifstream &in)
{
  READ(M_currLabel);
}

void RealLoop::StartFillingTunnels() {
  for (auto in : m_inTuns) {
    if (!in->IsLoopTunnel())
      LOG_FAIL("replacement for throw call");
    ((LoopTunnel*)in)->StartFillingSizes();
  }
}

void RealLoop::ClearTunnelCaches() {
  for (auto tun : m_inTuns) {
    tun->ClearDataTypeCache();
  }
}

void RealLoop::FillTunnelSizes()
{
  SplitBase *control = GetControl();
  if (!control) {
    LOG_FAIL("replacement for throw call");
  }
  bool upToDate = true;
  TunVecIter iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    LoopTunnel *tun = (LoopTunnel*)(*iter);
#if TWOD
    if (!tun->m_msizes)
#else
      if (!tun->m_sizes)
#endif
	{
	  upToDate = false;
	  break;
	}
  }
  if (upToDate)
    return;
    
  ClearTunnelCaches();
  StartFillingTunnels();

  vector<int> numIters;
#if DOBLIS
  control->BuildSizes(true, numIters, NumGroupsInComm(m_comm));
#else
  control->BuildSizes(true, numIters, 1);
#endif
#if DODM
      ((LoopTunnel*)control)->UpdateLocalSizes();
#endif

  for (auto inTun : m_inTuns) {
    if (inTun != control) {
#if DOBLIS
      ((LoopTunnel*)inTun)->BuildSizes(false, numIters,  NumGroupsInComm(m_comm));
#else
      ((LoopTunnel*)inTun)->BuildSizes(false, numIters, 1);
#endif
#if DODM
      ((LoopTunnel*)inTun)->UpdateLocalSizes();
#endif
    }
  }
}

void RealLoop::BuildDataTypeCache()
{
  FillTunnelSizes();
  IntLoop<RealPSet>::BuildDataTypeCache();
}


void HardDeleteNode(Node *node)
{
  node->m_poss->RemoveFromGraphNodes(node);
  NodeConnVecIter iter = node->m_children.begin();
  for(; iter != node->m_children.end(); ++iter)
    delete *iter;
  node->m_children.clear();
  iter = node->m_inputs.begin();
  for(; iter != node->m_inputs.end(); ++iter)
    delete *iter;
  node->m_inputs.clear();
}

void RealLoop::TryToDeleteLoopTunnelSetAndCleanUp(LoopTunnel *tun)
{
  if (!m_shadows.empty())
    LOG_FAIL("replacement for throw call");
  if (tun->GetNodeClass() != LoopTunnel::GetClass()) {
    cout << "only handling LoopTunnels specifically, not splits\n";
    LOG_FAIL("replacement for throw call");
  }
  if (tun->m_tunType != POSSTUNIN)
    LOG_FAIL("replacement for throw call");
  
  LoopTunnel *setTunIn = (LoopTunnel*)(tun->Input(0));
  
  //If there is an poss that uses the LoopTunnel's input,
  // then we shouldn't remove it from all
  //Check all children of the set input (poss inputs)
  // by checking all of their children - if all children's
  // children are TunOuts, then we're good.
  NodeConnVecIter childIter = setTunIn->m_children.begin();
  for(; childIter != setTunIn->m_children.end(); ++childIter) {
    Node *child = (*childIter)->m_n;
    NodeConnVecIter childIter2 = child->m_children.begin();
    for(; childIter2 != child->m_children.end(); ++childIter2) {
      Node *childChild = (*childIter2)->m_n;
      if (!childChild->IsTunnel(POSSTUNOUT))
        return;
    }
  }
  
  //At this point, we've checked, and there is no reason
  // to keep this set of loop tunnels around
  unsigned int tunNum = UINT_MAX;
  for(unsigned int i = 0; i < m_inTuns.size(); ++i) {
    if (m_inTuns[i] == setTunIn) {
      tunNum = i;
      break;
    }
  }
  if (tunNum == UINT_MAX)
    LOG_FAIL("replacement for throw call");
  
  Node *setTunOut = m_outTuns[tunNum];
  
  if (setTunOut->m_children.size() > 0)
    LOG_FAIL("replacement for throw call");
  
  NodeConnVecIter outIter = setTunOut->m_inputs.begin();
  for(; outIter != setTunOut->m_inputs.end(); ++outIter) {
    HardDeleteNode((*outIter)->m_n);
  }
  
  m_outTuns.erase(m_outTuns.begin()+tunNum);
  HardDeleteNode(setTunOut);
  
  childIter = setTunIn->m_children.begin();
  for(; childIter != setTunIn->m_children.end(); ++childIter) {
    Node *child = (*childIter)->m_n;
    HardDeleteNode(child);
    delete *childIter;
  }
  setTunIn->m_children.clear();
  
  m_inTuns.erase(m_inTuns.begin()+tunNum);
  
  setTunIn->m_poss->DeleteChildAndCleanUp(setTunIn, true);
}

#if DOBLIS
void RealLoop::Parallelize(Comm comm)
{
  if (NumGroupsInComm(comm) <= 1)
    LOG_FAIL("replacement for throw call");
  if (RemoveParallelization(comm))
    LOG_FAIL("replacement for throw call");
  m_comm = comm;
  m_hasProped = false;
  if (!HasIndepIters()) {
    bool found = false;
    //If this code changes, reflace in OnlyParallelizedOnNonIndependentData
    NodeVecIter iter = m_inTuns.begin();
    for(; iter != m_inTuns.end(); ++iter) {
      LoopTunnel *tun = (LoopTunnel*)(*iter);
      if (!tun->IndepIters() && !tun->InputIsTemp()) {
        if (found)
          LOG_FAIL("replacement for throw call");
        found = true;
        unsigned int numOut = tun->NumOutputs();
        if (tun->IsSplit()) {
	  if (tun->GetNodeClass() != SplitSingleIter::GetClass())
	    LOG_FAIL("replacement for throw call");
	  --numOut;
	}
        int i;
        for(i = 0; i < (int)(tun->m_children.size()); ++i) {
          Node *possTun = (tun->m_children[i])->m_n;
          Poss *poss = possTun->m_poss;
          NodeSet nodeSet;
          for (unsigned int j = 0; j < numOut; ++j) {
            AddUsersOfLiveOutput(possTun, j, nodeSet);
          }
          if (!nodeSet.size())
            LOG_FAIL("replacement for throw call");
          poss->FillClique(nodeSet);
	  LOG_FAIL("replacement for throw call");
#if 0
          CritSect *crit = (CritSect*)(poss->FormSetForClique(nodeSet, true));
          if (crit->RemoveParallelization(CORECOMM)) {
            //This critical section is around some hierarchy of PSets
            // from which parallel code cannot be removed without getting
            // rid of all code
            RemoveAndDeletePoss(poss, true);
            --i;
          }
          if (HasParallelCode(crit->m_posses[0])) {
            LOG_FAIL("replacement for throw call");
          }
#endif
        }
      }
    }
  }
  
  ClearDataTypeCache();
  
  //If we're parallelizing a loop that is on a poss
  // that just got duplicated as part of a transformation,
  // then that duplicated poss doesn't have its size cache.
  //We need to form it and this loop's size cache (which will be
  // different thanks to paralellization.
  if (m_ownerPoss)
    m_ownerPoss->BuildDataTypeCache();
}

bool ContainsOnlyParallelization(PSet *set)
{
  if (!set->IsReal())
    return ContainsOnlyParallelization(set->GetReal());
  if (set->IsLoop()) {
    Loop *loop = (Loop*)set;
    if (loop->IsParallel()) {
      return true;
    }
  }
  
  PossMMapIter iter = set->m_posses.begin();
  for(; iter != set->m_posses.end(); ++iter) {
    Poss *poss = (*iter).second;
    bool foundParSet = false;
    PSetVecIter setIter = poss->m_sets.begin();
    for(; !foundParSet && setIter != poss->m_sets.end(); ++setIter) {
      PSet *set = *setIter;
      if (ContainsOnlyParallelization(set))
        foundParSet = true;
    }
    if (!foundParSet) {
      bool foundParNode = false;
      NodeVecIter nodeIter = poss->m_possNodes.begin();
      for(; !foundParNode && nodeIter != poss->m_possNodes.end(); ++nodeIter) {
        if ((*nodeIter)->IsParallel()) {
          foundParNode = true;
        }
      }
      if (!foundParNode)
        return false;
    }
  }
  return true;
}

bool ContainsOnlyParallelization(const NodeSet &set)
{
  PSetSet setSet;
  NodeSetIter iter = set.begin();
  for(; iter != set.end(); ++iter) {
    Node *node = *iter;
    if (node->IsParallel())
      return true;
    if (node->IsTunnel(SETTUNIN)) {
      Tunnel *tun = (Tunnel*)node;
      if (setSet.find(tun->m_pset) == setSet.end()) {
        setSet.insert(tun->m_pset);
        if (ContainsOnlyParallelization(tun->m_pset))
          return true;
      }
    }
  }
  return false;
}

bool RealLoop::HasIndepIters() const
{
  NodeVecConstIter iter = m_inTuns.begin();
  for (; iter != m_inTuns.end(); ++iter) {
    const LoopTunnel *in = (LoopTunnel*)(*iter);
    if (!in->IndepIters()) {
      return false;
    }
  }
  return true;
}
#endif

#if TWOD
 void RealLoop::SetDimName(DimName dim)
 {
   m_functionality += (char)(48+dim);
   m_dim = dim;
 }
#endif



ShadowPSet* RealLoop::GetNewShadow()
{
  ShadowLoop *shadow = new ShadowLoop;
  shadow->m_realPSet = this;
  m_shadows.push_back(shadow);
  return shadow;
}

 BasePSet* RealLoop::GetNewInst() 
{
  RealLoop *real = new RealLoop(m_type);
  real->m_bsSize = m_bsSize;
  real->m_label = m_label;
#if DOBLIS
  real->m_comm = m_comm;
#endif
  return (BasePSet*)real;
}

 
 RealLoop::~RealLoop()
 {
#if PRINTTRACKING
   cout << "migrating in real loop " << this << endl;
#endif
   Migrate();
#if PRINTTRACKING
   cout << "done migrating in real loop " << this << endl;
#endif
 }
 
Cost RealLoop::Prop()
{
  const Poss *poss = m_posses.begin()->second;
  for(auto in : poss->m_inTuns) {
    if (in->GetNodeClass() == SplitSingleIter::GetClass()) {
      const SplitSingleIter *split = (SplitSingleIter*)in;
      for (unsigned int i = 0; i < (split->NumOutputs()-1); ++i) {
	string name = split->GetNameStr(i);
	for(auto in2 : m_inTuns) {
	  //Spliting a variable and one of the partitions has
	  // the same name as another input
	  //This can happen with nested loops where an input
	  // is split twice on different loops
	  if (!in2->IsSplit() && name == in2->GetInputNameStr(0)) {
	    m_ownerPoss->PrintTransVecUp();
	    m_posses.begin()->second->PrintTransVecUp();
	    LOG_FAIL("replacement for throw call");
	  }
	}
      }
    }
    else {
      for (auto in2 : m_inTuns) {
	if (in != in2 && in2->GetNodeClass() != SplitSingleIter::GetClass()) {
	  if (in->Input(0) == in2->Input(0))
	    LOG_FAIL("replacement for throw call");
	}
      }
    }
  }
  
  return IntLoop<RealPSet>::Prop();
}
