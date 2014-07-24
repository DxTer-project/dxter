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
#include "shadowPSet.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "loopSupport.h"

extern unsigned int M_phase;

ShadowPSet::ShadowPSet()
  : m_realPSet(NULL)
{
}

ShadowPSet::~ShadowPSet()
{
  m_realPSet->RemoveShadow(this);
  m_realPSet = NULL;
}

bool ShadowPSet::operator==(const BasePSet &rhs) const
{
  return *m_realPSet == rhs;
}

void ShadowPSet::Prop()
{
  if(m_hasProped)
    return;

  if (!m_realPSet)
    throw;

  if (!m_realPSet->IsReal())
    throw;

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

  m_realPSet->Prop();

  m_hasProped = true;
}

GraphNum ShadowPSet::TotalCount() const
{
  return m_realPSet->TotalCount();
}


void ShadowPSet::BuildDataTypeCache()
{

}

void ShadowPSet::ClearDataTypeCache()
{

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

const string& ShadowPSet::GetFunctionalityString() const
{
  return m_realPSet->GetFunctionalityString();
}


BasePSet* ShadowPSet::GetNewShadow()
{
  return m_realPSet->GetNewShadow();
}
