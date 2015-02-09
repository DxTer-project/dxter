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
#include "basePSet.h"
#include "node.h"
#include "helperNodes.h"

bool FoundInVec(const LinElem *elem, const LinElemVec &vec)
{
  for(auto t : vec)
    if (t == elem)
      return true;
  return false;
}

Linearizer::Linearizer(const Poss *poss)
{
  m_alwaysLiveCost = 0;
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
  
  if (node->GetNodeClass() == OutputNode::GetClass()) {
#if DOTENSORS
    m_alwaysLiveCost += ((DLANode*)(node->Input(0)))->MaxNumberOfLocalElements(node->InputConnNum(0));
#else
    m_alwaysLiveCost += ((DLANode*)(node->Input(0)))->MaxNumberOfElements(node->InputConnNum(0));
#endif
    m_alwaysLive.insert(node->GetInputNameStr(0));
  }
  
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
      else if (childConn->m_n->IsDataDependencyOfInput()) {
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
      if (other->m_succ) {
        cout << "already has succesor\n";
        throw;
      }
      else
        other->m_succ = overwriter;
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
        if (other->m_succ) {
          cout << "has succs\n";
          throw;
        }
        else
          other->m_succ = overwriter;
      }
    }
  }
  
  return elem;
}

void Linearizer::FindOptimalLinearization(const StrSet &stillLive)
{
  ClearCurrLinearization();
  
  m_lin.m_cost = -1;
  
  LinElemSet readyToAdd;
  
  Linearization curr;
  
  for(auto elem : m_elems) {
    if (elem->IsNode()) {
      NodeLinElem *nodeElem = (NodeLinElem*)elem;
      if (nodeElem->m_node->GetNodeClass() == InputNode::GetClass()) {
        curr.m_order.push_back(elem);
        elem->SetAdded();
      }
    }
  }
  
  //output nodes don't print any code and for tensors are often just connected to inputs, so
  // put them in the order now to reduce computation complexity of finding all topological sorts
  for(auto elem : m_elems) {
    if (elem->IsNode()) {
      NodeLinElem *nodeElem = (NodeLinElem*)elem;
      if (nodeElem->m_node->GetNodeClass() == OutputNode::GetClass()) {
        if (nodeElem->CanAddToLinearOrder())
        {
          curr.m_order.push_back(elem);
          elem->SetAdded();
        }
      }
    }
  }
  
#if 0
  int clumpCount = 0;
  int count = 0;
  cout << "\\\\  Scheduling the following ( " << m_elems.size() << " )\n";
  for(auto elem : m_elems) {
    if (!elem->HasAdded()) {
      bool clump = elem->ShouldClump();
      ++count;
      if (clump)
        ++clumpCount;
      if (elem->IsNode()) {
        cout << "\\\\" << ((NodeLinElem*)elem)->m_node->GetNodeClass() << (clump ? " will clump\n" : "\n");
      }
      else if (elem->IsSet())
        cout << "\\\\set " << ((SetLinElem*)elem)->m_set->GetFunctionalityString() << (clump ? " will clump\n" : "\n");
    }
  }
  cout << "\\\\" << count << " to schedule with " << clumpCount << " clumped\n";
#endif
  
  for(auto elem : m_elems) {
    if (elem->CanAddToLinearOrder()) {
      if (!elem->ShouldClump() ||
          !elem->OtherInputInClumpIsAlreadyRead(readyToAdd))
      {
        readyToAdd.insert(elem);
      }
    }
  }
  
  if (readyToAdd.empty())
    throw;
  
  //  cout << "\n\n\nfinding optimal with " << m_elems.size() << endl;
  //  cout.flush();
  
  RecursivelyFindOpt(curr, readyToAdd, m_lin, stillLive);
  
  for(auto elem : curr.m_order) {
    if (elem->IsNode()) {
      NodeLinElem *nodeElem = (NodeLinElem*)elem;
      if (nodeElem->m_node->GetNodeClass() != InputNode::GetClass() &&
          nodeElem->m_node->GetNodeClass() != OutputNode::GetClass())
        throw;
    }
    else
      throw;
  }
  
  if (m_lin.m_order.size() != m_elems.size()) {
    cout << m_lin.m_order.size() << endl;
    cout << m_elems.size() << endl;
    throw;
  }
  
  if (m_lin.GetCostNoRecursion(stillLive, m_alwaysLive) < 0) {
    m_lin.m_cost = -1;
    cout << m_lin.GetCostNoRecursion(stillLive, m_alwaysLive) << endl;
    throw;
  }
}

void Linearizer::RecursivelyFindOpt(Linearization &curr, const LinElemSet &readyToAdd, Linearization &opt, const StrSet &stillLive) const
{
  if (readyToAdd.empty()) {
    if (opt.m_cost < 0) {
      if (curr.m_order.size() != m_elems.size()) {
        cout << curr.m_order.size() << endl;
        cout << m_elems.size() << endl;
        LinElemSet tmp;
        tmp.insert(m_elems.begin(), m_elems.end());
        for(auto n : curr.m_order)
          tmp.erase(n);
        for(auto n : tmp) {
          cout << n->CanAddToLinearOrder() << endl;
          if (n->IsNode()) {
            cout << ((NodeLinElem*)n)->m_node->GetType() << endl;
          }
          else if (n->IsSet()) {
            cout <<((SetLinElem*)n)->m_set->GetFunctionalityString() << endl;
          }
        }
        
        throw;
      }
      opt = curr;
      opt.GetCostNoRecursion(stillLive, m_alwaysLive);
    }
    else {
      if (curr.m_order.size() != m_elems.size()) {
        cout << curr.m_order.size() << endl;
        cout << m_elems.size() << endl;
        LinElemSet tmp;
        tmp.insert(m_elems.begin(), m_elems.end());
        for(auto n : curr.m_order)
          tmp.erase(n);
        for(auto n : tmp) {
          cout << n->CanAddToLinearOrder() << endl;
          if (n->IsNode()) {
            cout << ((NodeLinElem*)n)->m_node->GetType() << endl;
          }
          else if (n->IsSet()) {
            cout <<((SetLinElem*)n)->m_set->GetFunctionalityString() << endl;
          }
        }
        
        throw;
        
      }
      if (opt.GetCostNoRecursion(stillLive, m_alwaysLive) > curr.GetCostNoRecursion(stillLive, m_alwaysLive)) {
        opt = curr;
      }
    }
  }
  else {
    for(auto elem : readyToAdd) {
      if (!elem->CreatesNewVars()) {
        LinElemSet tmp = readyToAdd;
        tmp.erase(elem);
        AddAndRecurse(curr, tmp, elem, opt, stillLive, false);
        return;
      }
    }
    
    LinElemSet skip;
    
    LinElemSetIter iter = readyToAdd.begin();
    for(; iter != readyToAdd.end(); ++iter) {
      LinElem *currAdd = *iter;
      if (skip.find(currAdd) != skip.end())
        continue;
      LinElemSet newReadyToAdd = readyToAdd;
      newReadyToAdd.erase(currAdd);
      
      //matching AddAndRecurse
      if (currAdd->m_children.size() == 1 && !currAdd->m_succ) {
        LinElem *child = currAdd->m_children[0];
        LinElemSetIter iter2 = iter;
        ++iter2;
        for(; iter2 != readyToAdd.end(); ++iter2) {
          LinElem *other = *iter2;
          //matching in AddAndRecurse
          if (other->m_children.size() == 1 && !other->m_succ) {
            if (other->m_children[0] == child) {
              skip.insert(other);
            }
          }
        }
      }
      
      AddAndRecurse(curr, newReadyToAdd, currAdd, opt, stillLive, true);
    }
  }
}

void Linearizer::AddClumpingInputsToReadyList(LinElemSet &readyToAdd, LinElem *child) const
{
  for (auto input : child->m_inputs) {
    if (input->CanAddToLinearOrder() && input->ShouldClump()) {
      readyToAdd.insert(input);
      break;
    }
  }
}

bool Linearizer::AddUp(Linearization &curr, LinElemSet &readyToAdd,
                       LinElem *currAdd, LinElem *ignore,
                       unsigned int &countOfAdded) const
{
  if (currAdd->HasAdded())
    throw;
  
  bool keepGoing = true;
  while (keepGoing && !currAdd->CanAddToLinearOrder()) {
    keepGoing = false;
    for (auto input : currAdd->m_inputs) {
      if (!input->HasAdded()) {
        keepGoing |= AddUp(curr, readyToAdd, input, currAdd, countOfAdded);
      }
    }
    for (auto pred : currAdd->m_preds) {
      if (!pred->HasAdded())
        keepGoing |= AddUp(curr, readyToAdd,
                           pred, currAdd,
                           countOfAdded);
    }
  }
  
  if (currAdd->CanAddToLinearOrder()) {
    readyToAdd.erase(currAdd);
    ++countOfAdded;
    curr.m_order.push_back(currAdd);
    currAdd->SetAdded();
    
    if (currAdd->m_succ) {
      if (currAdd->m_succ->CanAddToLinearOrder()) {
        readyToAdd.insert(currAdd->m_succ);
      }
      else {
        AddClumpingInputsToReadyList(readyToAdd, currAdd->m_succ);
      }
    }
    for(auto child : currAdd->m_children) {
      if (ignore != child)  {
        if (child->CanAddToLinearOrder()) {
          readyToAdd.insert(child);
        }
        else {
          AddClumpingInputsToReadyList(readyToAdd, child);
        }
      }
    }
    return true;
  }
  else {
    return false;
  }
}

void Linearizer::AddAndRecurse(Linearization &curr, LinElemSet &readyToAdd, LinElem *currAdd,
                               Linearization &opt, const StrSet &stillLive, bool addSingleChildImmediately) const
{
  if (!currAdd->CanAddToLinearOrder())
    throw;
  curr.m_order.push_back(currAdd);
  currAdd->SetAdded();
  
  //matching RecursivelyFindOpt...
  if (currAdd->m_children.size() > 1 || (currAdd->m_children.size() == 1 && (!addSingleChildImmediately && currAdd->m_children[0]->CreatesNewVars()))
      || currAdd->m_succ)
  {
    //Go through children and add those that should be added immediately
    int printedImmediately = 0;
    for(auto child : currAdd->m_children) {
      if (child->CanAddToLinearOrder()) {
        bool done = false;
        if (child->IsNode()) {
          NodeLinElem *childNode = (NodeLinElem*)child;
          Node *node = childNode->m_node;
          if (node->GetNodeClass() == OutputNode::GetClass()) {
            curr.m_order.push_back(child);
            child->SetAdded();
            ++printedImmediately;
            done = true;
            if (child->m_succ)
              throw;
          }
        }
        if (!done)
          readyToAdd.insert(child);
      }
      else {
        //If that was the last input to (e.g.) a set
        // other than temps, then now the temps can
        // add, followed by the set
        AddClumpingInputsToReadyList(readyToAdd, child);
      }
    }
    if (currAdd->m_succ
        && currAdd->m_succ->CanAddToLinearOrder()
        && readyToAdd.find(currAdd->m_succ) == readyToAdd.end())
    {
      LinElem *succ = currAdd->m_succ;
      bool skip = false;
      if (succ->ShouldClump()) {
        for (auto input : succ->m_children[0]->m_inputs) {
          if (input->ShouldClump()
              && (readyToAdd.find(input) != readyToAdd.end()))
          {
            skip = true;
            break;
          }
        }
      }
      if (!skip)
        readyToAdd.insert(currAdd->m_succ);
    }
    RecursivelyFindOpt(curr, readyToAdd, opt, stillLive);
    while (printedImmediately > 0) {
      LinElem *elem = curr.m_order.back();
      elem->ClearAdded();
      curr.m_order.pop_back();
      --printedImmediately;
    }
  }
  else if (currAdd->m_children.size() == 1) {
    LinElem *child = currAdd->m_children[0];
    bool keepGoing = true;
    unsigned int countOfAdded = 0;
    while (addSingleChildImmediately && keepGoing) {
      keepGoing = false;
      //this is a tempvar node
      // if the child is a set, print all of its other input temp var nodes and it
      // if it's not a set, try to print it
      for (auto input : child->m_inputs) {
        if (!input->HasAdded()
            && input != currAdd)
        {
          keepGoing |= AddUp(curr, readyToAdd, input, child, countOfAdded);
        }
      }
    }
    if (child->CanAddToLinearOrder()) {
      AddAndRecurse(curr, readyToAdd, child, opt, stillLive, addSingleChildImmediately);
    }
    else
      RecursivelyFindOpt(curr, readyToAdd, opt, stillLive);
    for(; countOfAdded > 0; --countOfAdded) {
      LinElem *elem = curr.m_order.back();
      curr.m_order.pop_back();
      elem->ClearAdded();
    }
  }
  else {
    //don't have children, so continue working with readyToAdd list
    RecursivelyFindOpt(curr, readyToAdd, opt, stillLive);
  }
  currAdd->ClearAdded();
  curr.m_order.pop_back();
}


void Linearizer::FindAnyLinearization()
{
  ClearCurrLinearization();
  
  LinElemSet readyToAdd;
  
  for(auto elem : m_elems) {
    if (elem->CanAddToLinearOrder()) {
      readyToAdd.insert(elem);
    }
  }
  
  while (!readyToAdd.empty()) {
    LinElem *elem = *(readyToAdd.begin());
    readyToAdd.erase(readyToAdd.begin());
    m_lin.m_order.push_back(elem);
    elem->SetAdded();
    if (elem->m_succ && elem->m_succ->CanAddToLinearOrder())
      readyToAdd.insert(elem->m_succ);
    for(auto child : elem->m_children) {
      if (child->CanAddToLinearOrder())
        readyToAdd.insert(child);
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
  m_lin.InsertVecClearing(stillLive, m_alwaysLive);
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
