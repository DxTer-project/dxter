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

#include "linearizer.h"
#include "poss.h"
#include "nodeLinElem.h"
#include "setLinElem.h"
#include "BasePSet.h"
#include "node.h"

void AddAndRecurse(Linearization &curr, LinElemVec &readyToAdd, LinElem *currAdd, Linearization &opt) ;
void RecursivelyFindOpt(Linearization &curr, const LinElemVec &readyToAdd, Linearization &opt) ;

Linearizer::Linearizer(const Poss *poss)
{
  PtrToLinElemMap map;
  for(auto node : poss->m_possNodes) {
    FindOrAdd(node, map);
  }
  for(auto set : poss->m_sets) {
    FindOrAdd(set, map);
  }
}

Linearizer::~Linearizer()
{
  for(auto elem : m_elems)
    delete elem;
}

LinElem* Linearizer::FindOrAdd(Node *node, PtrToLinElemMap &map)
{
  PtrToLinElemMapIter find = map.find(node);
  if (find != map.end())
    return find->second;

  if (node->IsTunnel()) {
    const Tunnel *tun = (Tunnel*)node;
    switch (tun->m_tunType)
      {
      case (POSSTUNIN):
      case (POSSTUNOUT):
	return NULL;
      case (SETTUNIN):
      case (SETTUNOUT):
	return FindOrAdd(tun->m_pset, map);
      default:
	throw;
      }
  }


  NodeLinElem *elem = new NodeLinElem(node);
  map[node] = elem;
  m_elems.push_back(elem);

  for(auto inputConn : node->m_inputs) {
    LinElem *inElem = FindOrAdd(inputConn->m_n, map);
    if (inElem)
      elem->AddInputIfUnique(inElem);
  }

  if (node->NumOutputs() > 1)
    throw;
  
  LinElem *overwriter = NULL;
  LinElemSet others;
  for(auto childConn : node->m_children) {
    LinElem *outElem = FindOrAdd(childConn->m_n, map);
    if (outElem) {
      elem->AddChildIfUnique(outElem);
      if (childConn->m_n->Overwrites(node, 0)) {
	if (overwriter) {
	  if (overwriter != outElem) {
	    cout << "two children of set overwrite output\n";
	    throw;
	  }
	}
	else {
	  overwriter = outElem;
	}
      }
      else {
	others.insert(outElem);
      }
    }
  }
  
  return elem;
}

LinElem* Linearizer::FindOrAdd(BasePSet *set, PtrToLinElemMap &map)
{
  PtrToLinElemMapIter find = map.find(set);
  if (find != map.end())
    return find->second;
  
  SetLinElem *elem = new SetLinElem(set);
  map[set] = elem;
  m_elems.push_back(elem);

  for(auto inTun : set->m_inTuns) {
    for(auto inputConn : inTun->m_inputs) {
      LinElem *inElem = FindOrAdd(inputConn->m_n, map);
      if (inElem)
	elem->AddInputIfUnique(inElem);
    }
  }

  for(auto outTun : set->m_outTuns) {
    LinElem *overwriter = NULL;
    LinElemSet others;
    for(auto childConn : outTun->m_children) {
      LinElem *outElem = FindOrAdd(childConn->m_n, map);
      if (outElem) {
	elem->AddChildIfUnique(outElem); 
	if (childConn->m_n->Overwrites(outTun, childConn->m_num)) {
	  if (overwriter) {
	    if (overwriter != outElem) {
	      cout << "two children of set overwrite output\n";
	      throw;
	    }
	  }
	  else {
	    overwriter = outElem;
	  }
	}
	else {
	  others.insert(outElem);
	}
      }
    }
    if (overwriter) {
      LinElemSetIter iter = others.find(overwriter);
      if (iter != others.end()) {
	others.erase(iter);
      }
      overwriter->m_preds.insert(overwriter->m_preds.begin(),
				 others.begin(),
				 others.end());
      for(auto other : others) {
	other->m_succs.push_back(overwriter);
      }
    }
  }
  
  return elem;
}

void Linearizer::FindOptimalLinearization()
{
  ClearCurrLinearization();
  
  m_lin.m_cost = 0;

  LinElemVec readyToAdd;

  for(auto elem : m_elems) {
    if (elem->CanAddToLinearOrder()) {
      readyToAdd.push_back(elem);
    }
  }

  if (readyToAdd.empty())
    throw;

  Linearization curr;

  RecursivelyFindOpt(curr, readyToAdd, m_lin);

  if (!curr.m_order.empty())
    throw;

  if (m_lin.GetCost() <= 0)
    throw;

  if (m_lin.m_order.size() != m_elems.size())
    throw;
}

void RecursivelyFindOpt(Linearization &curr, const LinElemVec &readyToAdd, Linearization &opt) 
{
  if (readyToAdd.empty()) {
    if (opt.m_cost == 0) {
      opt = curr;
    }
    else {
      if (opt.m_order.size() != curr.m_order.size())
	throw;
      if (opt.GetCost() > curr.GetCost())
	opt = curr;
    }
  }
  else {
    for(unsigned int i = 0; i < readyToAdd.size(); ++i) {
      LinElemVec newReadyToAdd = readyToAdd;
      newReadyToAdd.erase(newReadyToAdd.begin()+i);
      LinElem *currAdd = readyToAdd[i];

      AddAndRecurse(curr, newReadyToAdd, currAdd, opt); 
    }
  }
}

void AddAndRecurse(Linearization &curr, LinElemVec &readyToAdd, LinElem *currAdd, Linearization &opt) 
{
  curr.m_order.push_back(currAdd);
  currAdd->SetAdded();
  if (currAdd->m_children.size() > 1 || !currAdd->m_succs.empty()) {
    for(auto child : currAdd->m_children) {
      if (child->CanAddToLinearOrder()) {
	readyToAdd.push_back(child);
      }
    }
    for (auto succ : currAdd->m_succs) {
      if (succ->CanAddToLinearOrder()) {
	readyToAdd.push_back(succ);
      }
    }
    RecursivelyFindOpt(curr, readyToAdd, opt);
  }
  else if (currAdd->m_children.size() == 1) {
    LinElem *child = currAdd->m_children[0];
    if (child->CanAddToLinearOrder())
      AddAndRecurse(curr, readyToAdd, child, opt);
    else
      RecursivelyFindOpt(curr, readyToAdd, opt);
  }
  else {
    RecursivelyFindOpt(curr, readyToAdd, opt);
  }
  currAdd->ClearAdded();
  curr.m_order.pop_back();
}


void Linearizer::FindAnyLinearization()
{
  ClearCurrLinearization();

  LinElemVec readyToAdd;

  for(auto elem : m_elems) {
    if (elem->CanAddToLinearOrder()) {
      readyToAdd.push_back(elem);
    }
  }

  while (!readyToAdd.empty()) {
    LinElem *elem = readyToAdd.back();
    readyToAdd.pop_back();
    m_lin.m_order.push_back(elem);
    elem->SetAdded();
    for(auto succ : elem->m_succs) {
      if (succ->CanAddToLinearOrder())
	readyToAdd.push_back(succ);
    }
    for(auto child : elem->m_children) {
      if (child->CanAddToLinearOrder())
	readyToAdd.push_back(child);
    }
  }

  if (m_lin.m_order.size() != m_elems.size()) {
    PrintConnections();
    throw;
  }
}

void Linearizer::ClearCurrLinearization()
{
  m_lin.Clear();
  for(auto elem : m_elems) {
    elem->ClearCache();
    elem->ClearAdded();
  }
}

bool Linearizer::HasCurrLinearization() const
{
  if (m_lin.m_cost < 0) {
    return false;
  }
  else {
    if (m_lin.m_order.size() != m_elems.size())
      throw;
    return true;
  }
}

void Linearizer::InsertVecClearing(const StrSet &stillLive)
{
  if (!HasCurrLinearization())
    throw;
  m_lin.InsertVecClearing(stillLive);
}

void Linearizer::PrintConnections() const
{
  for(auto elem : m_elems) {
    cout << elem << ":\n";
    for(auto in : elem->m_inputs) {
      cout << "\tin:\t" << in << endl;
    }
    for(auto child : elem->m_children) {
      cout << "\tout:\t" << child << endl;
    }
  }
}
