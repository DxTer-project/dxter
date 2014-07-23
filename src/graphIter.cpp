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



#include "graphIter.h"

GraphIter::GraphIter(Poss *poss)
{
  m_poss = poss;
  m_setIters = new PossMMapIter[poss->m_sets.size()];
  m_subIters = new GraphIterPtr[poss->m_sets.size()];
  for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
    m_setIters[i] = m_poss->m_sets[i]->m_posses.begin();
    m_subIters[i] = new GraphIter(m_setIters[i]->second);
  }
}

GraphIter::GraphIter(const GraphIter &iter) 
{
  m_poss = NULL;
  *this = iter;
}

GraphIter::~GraphIter()
{
  for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
    delete m_subIters[i];
  }
  delete [] m_setIters;
  delete [] m_subIters;
  m_poss = NULL;
}

GraphIter& GraphIter::operator=(const GraphIter &rhs)
{
  if (m_poss) {
    for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
      delete m_subIters[i];
    }
    delete [] m_setIters;
    delete [] m_subIters;
  }
  m_poss = rhs.m_poss;
  m_setIters = new PossMMapIter[m_poss->m_sets.size()];
  m_subIters = new GraphIterPtr[m_poss->m_sets.size()];
  for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
    m_setIters[i] = rhs.m_setIters[i];
    m_subIters[i] = new GraphIter(rhs.m_subIters[i]);
  }
}

bool GraphIter::Increment()
{
  for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
    bool ret = m_subIters[i]->Increment();
    if (!ret) 
      return false;
    else {
      PSet *set = m_poss->m_sets[i];
      set->m_currHasPrinted = false;
      delete m_subIters[i];
      ++(m_setIters[i]);
      if (m_setIters[i] == set->m_posses.end()) {
	m_setIters[i] = set->m_posses.begin();
	m_subIters[i] = new GraphIter(m_setIters[i]->second());
      }
      else {
	m_subIters[i] = new GraphIter(m_setIters[i]->second());
	return false;
      }
    }
  }
  return true;
}


void GraphIter::GetCurrTransVec(TransVec &transVec)
{
  transVec.insert(transVec.end(),
		  m_poss->m_transVec.begin(),m_poss->m_transVec.end());
  for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
    m_subIters[i]->GetCurrTransVec(transVec);
  }
}

void GraphIter::AddCurrPossVars(VarSet &set) const
{
  NodeVecConstIter iter = m_poss->m_possNodes.begin();
  for(; iter != m_poss->m_possNodes.end(); ++iter) {
    (*iter)->AddVariables(set);
  }
  for (unsigned int i = 0; i < m_poss->m_sets.size(); ++i) {
    m_subIters[i]->GetCurrTransVec(transVec);
  }
}
