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

#include "critSect.h"
#include "loopSupport.h"

void CritSect::PrintCurrPoss(IndStream &out, unsigned int &graphNum)
{
  if (!m_ownerPoss->m_pset->IsLoop())
    throw;
  Loop *loop = (Loop*)(m_ownerPoss->m_pset);
  Comm comm = loop->m_comm;
  if (comm == CORECOMM)
    throw;
  *out << "Critical section with communicator " << CommToStr(comm) << "; need correct output code\n";
  out.Indent();
  *out << "GetMutex(" << CommToStr(comm) << ");\n";
  
  PSet::PrintCurrPoss(out, graphNum);
  
  out.Indent();
  *out << "ReleaseMutex(" << CommToStr(comm) << ");\n";
}

bool HasParallelCode(Poss *poss)
{
  PSetVecIter setIter = poss->m_sets.begin();
  for(; setIter != poss->m_sets.end(); ++setIter) {
    PSet *set = *setIter;
    if (set->IsLoop()) {
      Loop *loop = (Loop*)set;
      if (loop->IsParallel()) {
        return true;
      }
    }
    PossVecIter possIter = set->m_posses.begin();
    for(; possIter != set->m_posses.end(); ++possIter) {
      if (HasParallelCode(*possIter))
        return true;
    }
  }
  
  NodeVecIter nodeIter = poss->m_possNodes.begin();
  for(; nodeIter != poss->m_possNodes.end(); ++nodeIter) {
    Node *node = *nodeIter;
    if (node->IsParallel()){
      return true;
    }
  }
  
  return false;
}

void CritSect::BuildSizeCache()
{
  NodeVecIter iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    (*iter)->BuildSizeCacheRecursive();
  }
  PSet::BuildSizeCache();
}

void CritSect::SanityCheck()
{
  //  throw;
  PSet::SanityCheck();
  if (!m_ownerPoss->m_pset->IsLoop())
    throw;
  PossVecIter iter = m_posses.begin();
  for( ; iter != m_posses.end(); ++iter) {
    Poss *poss = *iter;
    if (HasParallelCode(poss)) {
      PSet *set = m_ownerPoss->m_pset;
      if (set->IsLoop()) {
        Loop *loop = (Loop*)set;
        cout << "parallelism: " << CommToStr(loop->m_comm);
      }
      cout.flush();
      poss->ForcePrint();
      throw;
    }
  }
  NodeVecIter iter2 = m_inTuns.begin();
  for(; iter2 != m_inTuns.end(); ++iter2) {
    if ((*iter2)->GetNodeClass() !=
        CritSectTunnel::GetClass()) {
      cout << (*iter2)->GetNodeClass() << endl;
      throw;
    }
  }
  iter2 = m_outTuns.begin();
  for(; iter2 != m_outTuns.end(); ++iter2) {
    if ((*iter2)->GetNodeClass() !=
        CritSectTunnel::GetClass())
      throw;
  }
}

CritSectTunnel::CritSectTunnel()
{
  m_msizes = NULL;
  m_nsizes = NULL;
  m_mlsizes = NULL;
  m_nlsizes = NULL;
}

CritSectTunnel::CritSectTunnel(PossTunType type) 
 : PossTunnel(type) 
{
  m_msizes = NULL;
  m_nsizes = NULL;
  m_mlsizes = NULL;
  m_nlsizes = NULL;
}

CritSectTunnel::~CritSectTunnel()
{
  if (m_msizes) {
    delete m_msizes;
    m_msizes = NULL;
    delete m_mlsizes;
    m_mlsizes = NULL;
    delete m_nsizes;
    m_nsizes = NULL;
    delete m_nlsizes;
    m_nlsizes = NULL;
  }
}

const Sizes* CritSectTunnel::GetM(unsigned int num) const
{
  if (num > 0)
    throw;
  switch(m_tunType)
  {
    case (POSSTUNOUT):
      return GetInputM(0);
    case (SETTUNOUT):
      return m_msizes;
    case (POSSTUNIN):
      return GetInputM(0);
    case (SETTUNIN):
      return m_msizes;
    default:
      throw;
  }
}


const Sizes* CritSectTunnel::GetN(unsigned int num) const
{
  if (num > 0)
    throw;
  switch(m_tunType)
  {
    case (POSSTUNOUT):
      return GetInputN(0);
    case (SETTUNOUT):
      return m_nsizes;
    case (POSSTUNIN):
      return GetInputN(0);
    case (SETTUNIN):
      return m_nsizes;
    default:
      throw;
  }
}

const Sizes* CritSectTunnel::LocalM(unsigned int num) const
{
  if (num > 0)
    throw;
  switch(m_tunType)
  {
    case (POSSTUNOUT):
      return InputLocalM(0);
    case (SETTUNOUT):
      return m_mlsizes;
    case (POSSTUNIN):
      return InputLocalM(0);
    case (SETTUNIN):
      return m_mlsizes;
    default:
      throw;
  }
}

const Sizes* CritSectTunnel::LocalN(unsigned int num) const
{
  if (num > 0)
    throw;
  switch(m_tunType)
  {
    case (POSSTUNOUT):
      return InputLocalN(0);
    case (SETTUNOUT):
      return m_nlsizes;
    case (POSSTUNIN):
      return InputLocalN(0);
    case (SETTUNIN):
      return m_nlsizes;
    default:
      throw;
  }
}

void CritSectTunnel::BuildSizeCache()
{
  if (m_msizes)
    return;
  if (m_tunType == POSSTUNIN || m_tunType == POSSTUNOUT)
    return;
  m_msizes = new Sizes(*GetInputM(0));
  m_mlsizes = new Sizes(*InputLocalM(0));
  m_nsizes = new Sizes(*GetInputN(0));
  m_nlsizes = new Sizes(*InputLocalN(0));
  if (m_tunType == SETTUNIN) {
    m_msizes->SetParFactor(1);
    m_nsizes->SetParFactor(1);
    m_mlsizes->SetParFactor(1);
    m_nlsizes->SetParFactor(1);
  }
  else {
    if (!m_poss->m_pset->IsLoop())
      throw;
    Loop *loop = (Loop*)(m_poss->m_pset);
    unsigned int p = NumGroupsInComm(loop->m_comm);
    m_msizes->AddParFactor(p);
    m_nsizes->AddParFactor(p);
    m_mlsizes->AddParFactor(p);
    m_nlsizes->AddParFactor(p);
  }
}

void CritSectTunnel::ClearSizeCache()
{
  if (!m_msizes)
    return;
  delete m_msizes;
  m_msizes = NULL;
  delete m_mlsizes;
  m_mlsizes = NULL;
  delete m_nsizes;
  m_nsizes = NULL;
  delete m_nlsizes;
  m_nlsizes = NULL;
}

void CritSectTunnel::SanityCheck()
{
  PossTunnel::SanityCheck();
  if (m_tunType == SETTUNIN || m_tunType == SETTUNOUT) {
    if (!m_pset->IsCritSect())
      throw;
  }
  else {
    if (!m_poss->m_pset->IsCritSect())
      throw;
  }
}
