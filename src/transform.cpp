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



#include "transform.h"

void MultiTrans::AddTrans(SingleTrans *trans)
{
  if (m_trans.empty())
    m_isRef = trans->IsRef();
  else {
    if (m_isRef != trans->IsRef())
      throw;
  }
  m_trans.push_back(trans);
}

TransConstVec* MultiTrans::GetApplicableTrans(const Poss *poss, const Node *node) const
{
  TransConstVec *applicable = new TransConstVec;
  CostVec costs;
  TransConstVecIter iter = m_trans.begin();
  for(; iter != m_trans.end(); ++iter) {
    if (!(*iter)->IsSingle())
      throw;
    const SingleTrans *trans = (SingleTrans*)(*iter);
    if (trans->CanApply(poss, node) && trans->WorthApplying(node)) {
      applicable->push_back(trans);
      Cost cost = trans->RHSCostEstimate(node);
      costs.push_back(cost);
    }    
  }
  
  if (costs.size() != applicable->size())
    throw;

  if (applicable->size() < MAXNUMBEROFREFINEMENTS)
    return applicable;

  
  //sort
  list<const Transformation*> tList;
  list<Cost> cList;
  iter = applicable->begin();
  CostVecIter costIter = costs.begin();
  for(; iter != applicable->end(); ++iter, ++costIter) {
    const Transformation *trans = *iter;
    Cost cost = *costIter;
    bool added = false;
    list<const Transformation*>::iterator tListIter = tList.begin();
    list<Cost>::iterator cListIter = cList.begin();
    for(; tListIter != tList.end() && !added; ++tListIter, ++cListIter) {
      if (cost < *cListIter) {
	added = true;
	tList.insert(tListIter, trans);
	cList.insert(cListIter, cost);
      }
    }
    if (!added && tList.size() < MAXNUMBEROFREFINEMENTS) {
      tList.push_back(trans);
      cList.push_back(cost);
      added = true;
    }
    if (added && tList.size() > MAXNUMBEROFREFINEMENTS) {
      tListIter = tList.begin();
      advance(tListIter,MAXNUMBEROFREFINEMENTS);
      tList.erase(tListIter,tList.end());
      cListIter = cList.begin();
      advance(cListIter,MAXNUMBEROFREFINEMENTS);
      cList.erase(cListIter,cList.end());
    }
  }

  applicable->clear();
  list<const Transformation*>::iterator tListIter = tList.begin();
  for(; tListIter != tList.end(); ++tListIter)
    applicable->push_back(*tListIter);

  return applicable;
}

int MultiTrans::CanApply(const Poss *poss, const Node *node, void **cache) const
{
  TransConstVec *vec = GetApplicableTrans(poss, node);
  *cache = vec;
  return vec->size();
}

void MultiTrans::CleanCache(void **cache) const
{
  delete (TransConstVec*)(*cache);
}

void MultiTrans::Apply(Poss *poss, int num, Node *node, void **cache) const
{
  TransConstVec *vec = (TransConstVec*)(*cache);
  if ((unsigned int)num >= vec->size())
    throw;
  const Transformation *trans = (*vec)[num];
  if (!trans->IsSingle())
    throw;
  ((SingleTrans*)trans)->Apply(poss, node);
}
