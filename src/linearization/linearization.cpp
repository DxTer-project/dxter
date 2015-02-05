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

#include "linearization.h"
#include "clearLinElem.h"
#include "base.h"
#include "setLinElem.h"
#include "realPSet.h"

Linearization::~Linearization()
{
  Clear();
}

void Linearization::Clear()
{
  for(auto clear : m_clears) {
    delete clear;
  }
  m_clears.clear();
  m_order.clear();
  m_cost = -1;
}


//GetCostNoRecursion
void Linearization::InsertVecClearing(const StrSet &stillLive, const StrSet &alwaysLive)
{
  if (!m_clears.empty())
    return;
  StrSet live;
  for(int i = m_order.size()-1; i >= 0; --i) {
    LinElem *elem = m_order[i];
    if (!elem->IsClear()) {
      if (elem->IsSet()) {
	/*
	  If input and used later, do not clear
	  If input, not output and not used later, clear within
	  If input and output but no children and not used later,
	       clear here
	  If input and output but has children or used later,
	       don't clear
	*/
	BasePSet *set = ((SetLinElem*)elem)->m_set;
	//If the variable is output, then it's live within

	StrSet liveHere = live;
	StrVec needClears;
	for (auto output : set->m_outTuns) {
	  string outName = output->GetNameStr(0); 
	  if (stillLive.find(outName) == stillLive.end() &&
	      alwaysLive.find(outName) == alwaysLive.end()) 
	    {
	      //erase outputs.  
	      //if they're still live, then will be input, too, and
	      // re-added below
	      live.erase(outName);
	      
	      if (live.find(outName) == live.end())
		{
		  needClears.push_back(outName);
		}
	    }

	  if (set->IsLoop()) {
	    Node *possTunOut = output->GetRealTunnel()->Input(0);
	    for(auto inToOut : possTunOut->m_inputs) {
	      liveHere.insert(inToOut->m_n->GetNameStr(inToOut->m_num));
	    }
	  }
	}

	for (auto input : set->m_inTuns) {
	  string inName = input->GetInputNameStr(0);
	  if (stillLive.find(inName) == stillLive.end()
	      && alwaysLive.find(inName) == alwaysLive.end())
	    {
	      live.insert(inName);
	    }
	  else
	    liveHere.insert(inName);
	}
	
	
	elem->CacheLiveVars(liveHere);
	
	for(string clear : needClears) {
	  ClearLinElem *newClear = new ClearLinElem(clear);
	  m_order.insert(m_order.begin()+i+1, newClear);
	  m_clears.push_back(newClear);
	}
      }
      else {
	StrSet goingLive;
	Node *node = ((NodeLinElem*)elem)->m_node;
	for (auto inConn : node->m_inputs) {
	  string inName = inConn->m_n->GetNameStr(inConn->m_num);
	  if (stillLive.find(inName) == stillLive.end()
	      && alwaysLive.find(inName) == alwaysLive.end())  {
	    goingLive.insert(inName);
	    if (live.find(inName) == live.end()) {
	      ClearLinElem *clear = new ClearLinElem(inName);
	      m_order.insert(m_order.begin()+i+1, clear);
	      m_clears.push_back(clear);
	    }
	  }
	}
	for (unsigned int i = 0; i < node->NumOutputs(); ++i) {
	  live.erase(node->GetNameStr(i));
	}
	live.insert(goingLive.begin(), goingLive.end());
      }
    }
  }
}


//Reflect in EnforceMemConstraint
Cost Linearization::GetCostNoRecursion(const StrSet &stillLive, const StrSet &alwaysLive)
{
  if (m_cost <= 0) {
    m_cost = 0;
    Cost currCost = 0;
    VarCostMap map;
    LinElemVecIter iter = m_order.begin();
    for( ; iter != m_order.end(); ++iter) {
      LinElem *elem = *iter;

      VarCostMap newVars = elem->NewVarsAndCosts();
      map.insert(newVars.begin(), newVars.end());
      for(auto &var : map) {
#if !DOTENSORS
	throw;
	/*
	  for domains other than tensors, we could be
	  in a loop where one variable is increasing
	  in size while another is decreasing.  If this
	  is the case, then the summ of their 
	  max size across iterations is about double
	  what their actual max size is
	*/
#endif
	currCost += var.second;
      }

      if (currCost > m_cost) {
	m_cost = currCost;
      }

      StrVec possiblyDyingVars = elem->PossiblyDyingVars();
      for(int i = 0; i < possiblyDyingVars.size(); ++i) {
	string var = possiblyDyingVars[i];
	if (stillLive.find(var) != stillLive.end()
	    || alwaysLive.find(var) != alwaysLive.end()) 
	  {
	    possiblyDyingVars.erase(possiblyDyingVars.begin()+i);
	    --i;
	  }
	else {
	  LinElemVecIter iter2 = iter;
	  ++iter2;
	  for(; iter2 != m_order.end(); ++iter2) {
	    if ((*iter2)->UsesInputVar(var)) {
	      possiblyDyingVars.erase(possiblyDyingVars.begin()+i);
	      --i;
	      break;
	    }
	  }
	}
      }
      //at this point all vars in possiblyDyingVars are actually dying
      for (auto var : possiblyDyingVars) {
	VarCostMapIter find = map.find(var);
	if (find != map.end()) {
	  m_cost -= find->second;
	  map.erase(find);
	}
      }
    }
  }
  return m_cost;
}

void Linearization::operator=(const Linearization &rhs)
{
  Clear();
  if (rhs.m_clears.empty()) {
    m_order = rhs.m_order;
  }
  else {
    m_order.reserve(rhs.m_order.size());
    m_clears.reserve(rhs.m_clears.size());
    for(auto elem : rhs.m_order) {
      if (elem->IsClear()) {
	ClearLinElem *clear = new ClearLinElem(*((ClearLinElem*)elem));
	m_order.push_back(clear);
	m_clears.push_back(clear);
      }
      else {
	m_order.push_back(elem);
      }
    }
  }
  m_cost = rhs.m_cost;
}

void Linearization::Print(IndStream &out)
{
  for(auto elem : m_order)
    elem->Print(out);
}

bool Linearization::LiveAfter(unsigned int loc, const string &name, const StrSet &alwaysLive) const
{
  if (alwaysLive.find(name) != alwaysLive.end())
    return true;
  for(; loc < m_order.size(); ++loc) {
    if (m_order[loc]->UsesInputVar(name))
      return true;
  }
  return false;
}

//Reflect in GetCostNoRecursion
bool Linearization::EnforceMemConstraint(Cost maxMem, const StrSet &stillLive, const StrSet &alwaysLive)
{
  Cost currCost = 0;
  StrSet live;
  for(int i = m_order.size() - 1; i >= 0; --i) {
    LinElem *elem = m_order[i];
    if (!elem->IsClear()) {
      if (elem->IsSet() && (((SetLinElem*)elem)->m_set)->IsReal()) {
	/*
	  If input and used later, do not clear
	  If input, not output and not used later, clear within
	  If input and output but no children and not used later,
	  clear here
	  If input and output but has children or used later,
	  don't clear
	*/
	RealPSet *set = (RealPSet*)(((SetLinElem*)elem)->m_set);
	//If the variable is output, then it's live within

	Cost costDiff = 0;

	StrSet liveHere = live;
	for (auto output : set->m_outTuns) {
	  string outName = output->GetNameStr(0); 
	  if (stillLive.find(outName) == stillLive.end() &&
	      alwaysLive.find(outName) == alwaysLive.end()) 
	    {
	      live.erase(outName);
#if DOTENSORS
	      costDiff -= ((DLANode*)output)->MaxNumberOfLocalElements(0);
#else
	/*
	  for domains other than tensors, we could be
	  in a loop where one variable is increasing
	  in size while another is decreasing.  If this
	  is the case, then the summ of their 
	  max size across iterations is about double
	  what their actual max size is
	*/
	      throw;
	      costDiff -=((DLANode*)output)->MaxNumberOfElements(0);
#endif
	    }
	  
	  if (set->IsLoop()) {
	    Node *possTunOut = output->GetRealTunnel()->Input(0);
	    for(auto inToOut : possTunOut->m_inputs) {
	      liveHere.insert(inToOut->m_n->GetNameStr(inToOut->m_num));
	    }
	  }
	}
	
	for (auto input : set->m_inTuns) {
	  string inName = input->GetInputNameStr(0);
	  if (stillLive.find(inName) == stillLive.end()
	      && alwaysLive.find(inName) == alwaysLive.end()
	      && live.find(inName) == live.end())
	    {
	      live.insert(inName);
#if DOTENSORS
	      costDiff += ((DLANode*)(input->Input(0)))->MaxNumberOfLocalElements(input->InputConnNum(0));
#else
	/*
	  for domains other than tensors, we could be
	  in a loop where one variable is increasing
	  in size while another is decreasing.  If this
	  is the case, then the summ of their 
	  max size across iterations is about double
	  what their actual max size is
	*/
	      throw;
	      costDiff += ((DLANode*)(input->Input(0)))->MaxNumberOfElements(input->InputConnNum(0));
#endif
	    }
	  else
	    liveHere.insert(inName);
	}
	
	if (set->EnforceMemConstraint(maxMem-currCost, liveHere))
	  return true;
	
	currCost += costDiff;
      }
      else {
	Node *node = ((NodeLinElem*)elem)->m_node;
	Cost costDiff = 0;
	for (unsigned int i = 0; i < node->NumOutputs(); ++i) {
	  string name = node->GetNameStr(i);
	  if (stillLive.find(name) == stillLive.end()
	       && alwaysLive.find(name) == alwaysLive.end())
            {
              live.insert(name);
#if DOTENSORS
              costDiff -= ((DLANode*)node)->MaxNumberOfLocalElements(i);
#else
              throw;
              costDiff -= ((DLANode*)node)->MaxNumberOfElements(i);
#endif
	    }
	}
        for(auto inputConn : node->m_inputs) {
	  string inName = inputConn->m_n->GetInputNameStr(inputConn->m_num);
	  if (stillLive.find(inName) == stillLive.end()
	      && alwaysLive.find(inName) == alwaysLive.end()
	      && live.find(inName) == live.end())
	    {
	      live.insert(inName);
#if DOTENSORS
	      costDiff += ((DLANode*)(inputConn->m_n))->MaxNumberOfLocalElements(inputConn->m_num);
#else
	/*
	  for domains other than tensors, we could be
	  in a loop where one variable is increasing
	  in size while another is decreasing.  If this
	  is the case, then the summ of their 
	  max size across iterations is about double
	  what their actual max size is
	*/
	      throw;
	      costDiff += ((DLANode*)(inputConn->m_n))->MaxNumberOfElements(inputConn->m_num);
#endif
	    }
	}
	  currCost += costDiff;
	  
      if (currCost >= maxMem)
	throw;
	}
      }
      }
      
  return false;
}


