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



#include "transform.h"
#include "poss.h"
#include "pset.h"
#include "loopSupport.h"
#include <sstream>
#include "elemRedist.h"
#include "tensorRedist.h"
#include <iomanip>
#include "twoSidedTrxm.h"
#include "pack.h"
#include "critSect.h"

#define CHECKFORLOOPS

#define ALWAYSFUSE !DOLLDLA

StrSet Poss::M_fusedSets;

unsigned int Poss::M_count = 1;

Poss::Poss()
{
  m_num = M_count;
  m_parent = 0;
  ++M_count;
  m_fullyExpanded = false;
  m_pset = NULL;
  m_hashValid = false;
  m_isSane = true;
}

Poss::Poss(PossTunnel *tunn)
{
  if (tunn->m_tunType != POSSTUNOUT) {
    cout << "tunn->m_tunType != POSSTUNOUT\n";
    throw;
  }
  m_num = M_count;
  m_parent = 0;
  ++M_count;
  AddUp(m_possNodes, tunn, true, false);
  m_fullyExpanded = false;
  m_pset = NULL;
  m_hashValid = false;
  m_isSane = true;
}

Poss::Poss(Node *node, bool goUp)
{
  m_num = M_count;
  m_parent = 0;
  ++M_count;
  m_fullyExpanded = false;
  m_pset = NULL;
  m_hashValid = false;
  m_isSane = true;
  
  if (!goUp) {
    node->SetPoss(this);
    m_possNodes.push_back(node);
    
    for (unsigned int i = 0; i < node->NumOutputs(); ++i) {
      PossTunnel *out = new PossTunnel(POSSTUNOUT);
      out->AddInput(node, i);
      out->SetPoss(this);
      m_possNodes.push_back(out);
      m_outTuns.push_back(out);
    }
    
    
    for (unsigned int i = 0; i < node->m_inputs.size(); ++i) {
      NodeConn *conn = node->InputConn(i);
      PossTunnel *in = new PossTunnel(POSSTUNIN);
      in->AddInput(conn->m_n, conn->m_num);
      conn->m_n->RemoveChild(node, conn->m_num);
      conn->m_n = in;
      conn->m_num = 0;
      in->AddChild(node, 0);
      
      in->SetPoss(this);
      m_possNodes.push_back(in);
      m_inTuns.push_back(in);
    }
  }
  else {
    for (unsigned int i = 0; i < node->NumOutputs(); ++i) {
      PossTunnel *out = new PossTunnel(POSSTUNOUT);
      out->AddInput(node, i);
      AddUp(m_possNodes, out, true, false);
    }
  }
}


Poss::Poss(int numArgs, ...)
{
  NodeVec nodes;
  va_list listPointer;
  va_start (listPointer, numArgs);
  for(int i = 0; i < numArgs; i++) {
    Node *node = va_arg(listPointer, Node* );
    nodes.push_back(node);
  }
  InitHelper(nodes, true, false);
}

Poss::Poss(const NodeVec &nodes, bool outTuns, bool disconnectFromOwner)
{
  InitHelper(nodes, outTuns, disconnectFromOwner);
}


void Poss::MarkInsane(bool wrongPhase)
{
  if (!wrongPhase) {
    cout << "!!!!!!!wrongPhase but is insane\n";
    this->PrintTransVec();
    throw;
  }
  m_isSane=false;
}

void Poss::InitHelper(const NodeVec &vec, bool outTuns, bool disconnectFromOwner)
{
  m_num = M_count;
  m_parent = 0;
  ++M_count;
  m_fullyExpanded = false;
  m_pset = NULL;
  m_hashValid = false;
  m_isSane = true;
  
  
  NodeVecConstIter iter = vec.begin();
  for(; iter != vec.end(); ++iter) {
    Node *node = *iter;
    if (outTuns) {
      if (!node->IsPossTunnel()) {
        PossTunnel *out = new PossTunnel(POSSTUNOUT);
        out->AddInput(node);
        AddUp(m_possNodes, out, true, disconnectFromOwner);
      }
      else
        AddUp(m_possNodes, node, true, disconnectFromOwner);
    }
    else {
      AddUp(m_possNodes, node, false, disconnectFromOwner);
    }
  }
  
}

Poss::~Poss()
{
  {
    NodeVecIter iter = m_possNodes.begin();
    for( ; iter != m_possNodes.end(); ++iter) {
      Node *node = *iter;
      delete node;
    }
    m_possNodes.clear();
  }
  PSetVecIter iter = m_sets.begin();
  for(; iter != m_sets.end(); ++iter)
    delete *iter;
}

void Poss::ForcePrint()
{
  ClearBeforeProp();
  Prop();
  bool hasPrinted = true;
  ClearNodesPrinted();
  NodeVecConstIter iter = m_inTuns.begin();
  for( ; iter != m_inTuns.end(); ++iter)
    (*iter)->SetPrinted();
  while(hasPrinted) {
    hasPrinted = false;
    NodeVecConstIter nodeIter = m_possNodes.begin();
    for( ; nodeIter != m_possNodes.end(); ++nodeIter) {
      if (!(*nodeIter)->HasPrinted()) {
#ifdef DOELEM
        IndStream out(&cout,ELEMSTREAM);
#elif DOSQM
        IndStream out(&cout,BLISSTREAM);
#endif
        (*nodeIter)->Print(out, 1);
        hasPrinted |= (*nodeIter)->HasPrinted();
      }
    }
  }
}

bool Poss::CanPrint() const
{
  if (m_hasPrinted)
    return false;
  NodeVecConstIter iter = m_inTuns.begin();
  for( ; iter != m_inTuns.end(); ++iter) {
    if (!(*iter)->CanPrintCode()) {
      return false;
    }
  }
  return true;
}

void Poss::Duplicate(const Poss *orig, NodeMap &map, bool possMerging)
{
  /**********
	     Any changes to this function should be reflected
	     in PSet::InlinePoss
  **********/
  Poss *poss = (Poss*)orig;
  m_parent = poss->m_num;
  //If this changes, update PSet::InlinePoss
  NodeVecConstIter iter = orig->m_possNodes.begin();
  for( ; iter != orig->m_possNodes.end(); ++iter) {
    Node *oldNode = *iter;
    Node *newNode = oldNode->GetNewInst();
    newNode->Duplicate(oldNode,false, possMerging);
    if (newNode->m_inputs.size() != oldNode->m_inputs.size())
      throw;
    m_possNodes.push_back(newNode);
    newNode->m_poss = this;
    map[oldNode] = newNode;
  }
  PSetVecIter setIter = poss->m_sets.begin();
  for(; setIter != poss->m_sets.end(); ++setIter) {
    PSet *newSet = (*setIter)->GetNewInst();
    newSet->Duplicate(*setIter, map,possMerging);
    AddPSet(newSet, true);
  }
  setIter = m_sets.begin();
  for(; setIter != m_sets.end(); ++setIter) {
    (*setIter)->PatchAfterDuplicate(map);
  }
  iter = poss->m_inTuns.begin();
  for(; iter != poss->m_inTuns.end(); ++iter)  {
    Node *node = map[*iter];
    if (!node) {
      cout << "!node in dup\n";
      throw;
    }
    m_inTuns.push_back(node);
    if (!(*iter)->m_poss) {
      cout << "!(*iter)->m_poss for " << *iter << endl;
      throw;
    }
  }
  iter = poss->m_outTuns.begin();
  for(; iter != poss->m_outTuns.end(); ++iter) {
    Node *node = map[*iter];
    if (!(*iter)->m_poss) {
      cout << "!(*iter)->m_poss\n";
      throw;
    }
    if (!node) {
      cout << "!node in dup\n";
      throw;
    }
    m_outTuns.push_back(node);
  }
  m_transVec = poss->m_transVec;
}

bool Poss::operator==(Poss &rhs)
{
  Poss &poss = (Poss&)rhs;
  if (poss.GetHash() != GetHash()) {
    return false;
  }
  if (m_possNodes.size() != poss.m_possNodes.size()
      || m_inTuns.size() != poss.m_inTuns.size()
      || m_outTuns.size() != poss.m_outTuns.size()
      || m_sets.size() != poss.m_sets.size())
    return false;
  for(unsigned int i = 0; i < m_outTuns.size(); ++i)
    if (*OutTun(i) != *poss.OutTun(i))
      return false;
  return true;
}

void Poss::AddNode(Node *node)
{
  InvalidateHash();
  if (node->m_poss) {
    cout << "node already has a poss\n";
    throw;
  }
  node->m_poss = this;
#if 0
  NodeVecIter iter = m_possNodes.begin();
  for(; iter != m_possNodes.end(); ++iter) {
    if (*iter == node)  {
      cout << "node already in my list\n";
      throw;
    }
  }
#endif
  m_possNodes.push_back(node);
}

void Poss::AddNodes(int numNodes, ...)
{
  va_list listPointer;
  va_start (listPointer, numNodes);
  for(int i = 0; i < numNodes; i++) {
    Node *node = va_arg(listPointer, Node*);
    AddNode(node);
  }
}

void Poss::AddPSet(PSet *pset, bool expectToBeNew)
{
  bool isNew = AddElemToVec(m_sets,pset,false);
  if (!isNew && expectToBeNew) {
    cout << "didn't add\n";
    throw;
  }
  else {
    if (pset->m_ownerPoss && pset->m_ownerPoss != this) {
      cout << "already has an owner!\n";
      throw;
    }
    pset->m_ownerPoss = this;
  }
  InvalidateHash();
}

void Poss::DeleteNode(Node *node)
{
  InvalidateHash();
  NodeVecIter iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter) {
    if (*iter == node) {
      m_inTuns.erase(iter);
      break;
    }
  }
  iter = m_outTuns.begin();
  for(; iter != m_outTuns.end(); ++iter){
    if (*iter == node) {
      m_outTuns.erase(iter);
      break;
    }
  }
  iter = m_possNodes.begin();
  unsigned int i = 0;
  for( ; iter != m_possNodes.end(); ++iter) {
    if (*iter == node) {
      m_possNodes.erase(iter);
      delete node;
      for(iter = m_possNodes.begin(); iter != m_possNodes.end(); ++iter) {
        if (*iter == node) {
          cout << "found second instance of node!!!\n";
        }
      }
      return;
    }
    ++i;
  }
  printf("node %s %p not found in %p\n", node->GetType().c_str(), node, this);
  throw;
}

void Poss::DeleteChildAndCleanUp(Node *output,
                                 bool goThroughTunnels, bool handleTunnelsAsNormalNodes)
{
  InvalidateHash();
  bool found = false;
  NodeVecConstIter tmp = m_possNodes.begin();
  for(; !found && tmp != m_possNodes.end(); ++tmp)
    if (*tmp == output) {
      found = true;
    }
  if (!found) {
    cout << "DeleteChildAndCleanUp called on wrong poss!\n";
    throw;
  }
  if (output->IsLoopTunnel() && !handleTunnelsAsNormalNodes) {
    if (!output->IsPossTunnel(SETTUNIN)) {
      if (!goThroughTunnels) {
        cout << "calling DeleteChildAndCleanUp on PossTunnel!\n";
        throw;
      }
      else if (output->GetNodeClass() == LoopTunnel::GetClass()){
        ((Loop*)((LoopTunnel*)output)->m_pset)->TryToDeleteLoopTunnelSetAndCleanUp((LoopTunnel*)output);
      }
      else
        throw;
    }
  }
  
  NodeConnVecIter iter = output->m_inputs.begin();
  for( ; iter != output->m_inputs.end(); ++iter) {
    Node *input = (*iter)->m_n;
    if ((input->m_poss != output->m_poss)
        && !input->IsPossTunnel()
        && !output->IsPossTunnel())
      {
	throw;
	cout << "input->m_poss != output->m_poss\n";
      }
    input->RemoveChild(output,(*iter)->m_num);
    if(input->m_children.empty()) {
      input->m_poss->DeleteChildAndCleanUp(input);
    }
    delete *iter;
  }
  output->m_inputs.clear();
  DeleteNode(output);
}




void Poss::AddUp(NodeVec &vec, Node *node, bool start, bool disconnectFromOwner)
{
  InvalidateHash();
  if (node->IsPossTunnel(POSSTUNIN) && !start && !node->m_poss) {
    if (node->m_poss && node->m_poss != this) {
      cout << "node already on a poss\n";
      throw;
    }
    if (AddElemToVec(vec, node, false)) {
      m_inTuns.push_back(node);
      if (node->m_poss && node->m_poss != this) {
        if (disconnectFromOwner) {
          node->m_poss->RemoveFromGraphNodes(node);
        }
        else
          throw;
      }
      node->SetPoss(this);
    }
  }
  else if (AddElemToVec(vec, node, false)) {
    if (node->IsPossTunnel(SETTUNOUT)) {
      //Found a PSet to be added to this poss
      // add the output tunnels of the set, add the
      // set to my set list, and add the input set
      // tunnels and everything preceeding them
      PSet *pset = ((PossTunnel*)node)->m_pset;
      if (pset->m_ownerPoss && pset->m_ownerPoss != this && disconnectFromOwner)
        pset->m_ownerPoss->RemoveFromSets(pset);
      AddPSet(pset, false);
      
      NodeVecIter iter = pset->m_outTuns.begin();
      for(; iter != pset->m_outTuns.end(); ++iter) {
        Node *out = *iter;
        if (out->m_poss && out->m_poss != this) {
          if (disconnectFromOwner) {
            Poss *poss = out->m_poss;
            poss->RemoveFromGraphNodes(out);
          }
          else
            throw;
        }
        out->SetPoss(this);
      }
      
      iter = pset->m_inTuns.begin();
      for(; iter != pset->m_inTuns.end(); ++iter) {
        AddUp(vec, *iter, false, disconnectFromOwner);
      }
      iter = pset->m_outTuns.begin();
      for(; iter != pset->m_outTuns.end(); ++iter) {
        if (AddElemToVec(vec, *iter, false)) {
          if (node->m_poss && node->m_poss != this) {
            if (disconnectFromOwner) {
              node->m_poss->RemoveFromGraphNodes(node);
            }
            else
              throw;
          }
          (*iter)->SetPoss(this);
        }
      }
      return;
    }
    else {
      if (node->IsPossTunnel(POSSTUNOUT) && start) {
        m_outTuns.push_back(node);
      }
      if (node->m_poss && node->m_poss != this) {
        if (disconnectFromOwner) {
          node->m_poss->RemoveFromGraphNodes(node);
        }
        else
          throw;
      }
      node->SetPoss(this);
      NodeConnVecIter iter = node->m_inputs.begin();
      for( ; iter != node->m_inputs.end(); ++iter) {
        AddUp(vec, (*iter)->m_n, false, disconnectFromOwner);
      }
    }
  }
}

void Poss::AddLoop(Loop *loop)
{
  NodeVecIter iter = loop->m_outTuns.begin();
  for(; iter != loop->m_outTuns.end(); ++iter) {
    AddUp(m_possNodes, *iter, true, false);
  }
}

bool Poss::Simplify(const TransMap &simplifiers)
{
  bool didSomething = false;
  for(unsigned int nodeIdx = 0; nodeIdx < m_possNodes.size(); ++nodeIdx) {
    Node *node = m_possNodes[nodeIdx];
    TransMapConstIter iter = simplifiers.find(node->GetNodeClass());
    if (iter != simplifiers.end()) {
      TransVecConstIter transIter = (*iter).second->begin();
      for(; transIter != (*iter).second->end(); ++transIter) {
        const Transformation *trans = *transIter;
        if (!trans->IsSingle())
          throw;
        if (((SingleTrans*)trans)->CanApply(node)) {
          didSomething = true;
          InvalidateHash();
          ((SingleTrans*)trans)->Apply(node);
          m_transVec.push_back(const_cast<Transformation*>(trans));
          nodeIdx = -1;
          break;
        }
      }
    }
  }
  return didSomething;
}

void Poss::ClearNodesPrinted()
{
  NodeVecConstIter nodeIter = m_possNodes.begin();
  for( ; nodeIter != m_possNodes.end(); ++nodeIter)
    (*nodeIter)->ClearPrinted();
}


void Poss::PatchAfterDuplicate(NodeMap &map)
{
  NodeVecIter iter = m_possNodes.begin();
  for( ; iter != m_possNodes.end(); ++iter) {
    (*iter)->PatchAfterDuplicate(map);
  }
}


void Poss::Prop()
{
  NodeVecIter nodeIter = m_possNodes.begin();
  for( ; nodeIter != m_possNodes.end(); ++nodeIter) {
    Node *node = *nodeIter;
    node->Prop();
    if (node->m_poss != this) {
      cout << node->GetType() << endl;
      cout << node->GetNameStr(0) << endl;
      cout << node->GetNodeClass() << endl;
      throw;
    }
  }
  PSetVecIter setIter = m_sets.begin();
  for( ; setIter != m_sets.end(); ++setIter) {
    PSet *set = *setIter;
    set->Prop();
  }

  if (m_inTuns.size() != m_pset->m_inTuns.size()) {
    cout << m_inTuns.size() << endl;
    cout << m_pset->m_inTuns.size() << endl;
    cout << this << " " << m_pset << endl;
    throw;
  }
  
  unsigned int i = 0;
  nodeIter = m_inTuns.begin();
  for(; nodeIter != m_inTuns.end(); ++nodeIter, ++i) {
    if (m_pset->InTun(i) != (*nodeIter)->Input(0))
      throw;
    if (!m_pset->InTun(i)->InChildren(*nodeIter,0)) {
      cout << "!m_pset->m_inputs[i]->InChildren(*iter,0)\n";
      cout << "poss tunnel " << *nodeIter << " on poss " << this << " on set " << m_pset << endl;
      cout << "expected on set tunnel " << m_pset->InTun(i) << endl;
      cout << "i = " << i << endl;
      for(unsigned int j = 0; j < m_pset->InTun(0)->m_children.size(); ++j) {
        cout << m_pset->InTun(i)->m_children[j]->m_n << " is " << j << " child\n";
        if ( m_pset->InTun(i)->m_children[j]->m_num)
          cout << "!!!!!" << m_pset->InTun(i)->m_children[j]->m_num << " is " << j << " num\n";
      }
      throw;
    }
  }
}

void Poss::Cull(Phase phase)
{
  if (!m_isSane)
    throw;
  NodeVecIter nodeIter = m_possNodes.begin();
  for( ; nodeIter != m_possNodes.end(); ++nodeIter) {
    (*nodeIter)->Cull(phase);
  }
  PSetVecIter possIter = m_sets.begin();
  for( ; possIter != m_sets.end(); ++possIter) {
    (*possIter)->Cull(phase);
  }
}

void Poss::PrintTransVec()
{
  TransVecConstIter iter = m_transVec.begin();
  for(; iter != m_transVec.end(); ++iter)
    cout << (*iter)->GetType() << endl;
  if (m_sets.size()) {
    //    cout << "not supported for hierarchical posses\n";
    //    throw;
  }
}


void Poss::RemoveConnectionToSet()
{
  InvalidateHash();
  for (unsigned int i = 0; i < m_inTuns.size(); ++i) {
    Node *tun = InTun(i);
    for (unsigned int j = 0; j < tun->m_inputs.size(); ++j) {
      delete tun->InputConn(j);
    }
    InTun(i)->m_inputs.clear();
  }
  for (unsigned int i = 0; i < m_outTuns.size(); ++i) {
    Node *tun = OutTun(i);
    for (unsigned int j = 0; j < tun->m_children.size(); ++j) {
      delete tun->m_children[j];
    }
    OutTun(i)->m_children.clear();
  }
}

void Poss::ExpandTunnels()
{
  InvalidateHash();
  if (m_pset) {
    cout << "m_pset is set\n";
    throw;
  }
  for (unsigned int i = 0; i < m_inTuns.size(); ++i) {
    Node *tunIn = InTun(i);
    while (tunIn->m_inputs.size() > 1 && tunIn->GetNodeClass() != Split::GetClass()) {
      cout << "splitting in\n";
      PossTunnel *tun = (PossTunnel*)(tunIn->GetNewInst());
      tun->m_tunType = POSSTUNIN;
      unsigned int inNum = tunIn->m_inputs.size()-1;
      NodeConn *conn = tunIn->InputConn(inNum);
      tun->AddInput(conn->m_n,conn->m_num);
      conn->m_n->RemoveChild(tunIn,conn->m_num);
      AddNode(tun);
      m_inTuns.insert(m_inTuns.begin()+i,tun);
      delete conn;
      tunIn->m_inputs.erase(tunIn->m_inputs.begin()+inNum);
      tunIn->RedirectChildren(inNum, tun, 0);
    }
  }
  for (unsigned int i = 0; i < m_outTuns.size(); ++i) {
    Node *tunOut = OutTun(i);
    while (tunOut->m_inputs.size() > 1 && !tunOut->IsLoopTunnel()) {
      if (tunOut->m_children.size()) {
        cout << "tunOut has children\n";
        throw;
      }
      cout << "splitting out\n";
      PossTunnel *tun = (PossTunnel*)(tunOut->GetNewInst());
      tun->m_tunType = POSSTUNOUT;
      unsigned int inNum = tunOut->m_inputs.size()-1;
      NodeConn *conn = tunOut->InputConn(inNum);
      tun->AddInput(conn->m_n,conn->m_num);
      conn->m_n->RemoveChild(tunOut,conn->m_num);
      AddNode(tun);
      m_outTuns.insert(m_outTuns.begin()+i,tun);
      delete conn;
      tunOut->m_inputs.erase(tunOut->m_inputs.begin()+inNum);
      
    }
  }
}

bool Poss::MergePosses(PossMMap &newPosses,const TransMap &simplifiers, CullFunction cullFunc)
{
  /*
    First, recurse through my PSet's then check if
    any of them can be merged if the recursion doesn't change
    anything
  */
  InvalidateHash();
  bool didMerge = false;
  PSetVecIter iter = m_sets.begin();
  for(; iter != m_sets.end(); ++iter) {
    PSet *pset = *iter;
    if (pset->m_ownerPoss != this) {
      cout << "bad owner\n";
      throw;
    }
    didMerge |= pset->MergePosses(simplifiers, cullFunc);
  }
  if (!didMerge) {
    //Didn't make any changes in the recursion
    //First see if there are any loops to fuse
    // (check that the fusion hasn't already been done)
    //If so, duplicate this and do the fusion
    //If not, then merge non-loop posses
    if (m_sets.size() >= 2) {
      for (unsigned int left = 0; left < m_sets.size()-1; ++left) {
        //check if the poss set is transparent poss (i.e. not a loop)
        PSet *leftSet = m_sets[left];
        if (leftSet->IsLoop()) {
          for (unsigned int right = left + 1; right < m_sets.size(); ++right) {
            PSet *rightSet = m_sets[right];
            if (rightSet->IsLoop()) {
              if (!((Loop*)(m_sets[left]))->WorthFusing((Loop*)(m_sets[right]))) {
                continue;
              }
              if (HasFused((Loop*)(m_sets[left]),(Loop*)(m_sets[right]))) {
                //		cout << "has fused\n\n";
                continue;
              }
              if (m_sets[left]->CanMerge(m_sets[right])) {
                NodeMap nodeMap;
                NodeVecIter tunIter = m_pset->m_inTuns.begin();
                for(; tunIter != m_pset->m_inTuns.end(); ++tunIter) {
                  nodeMap[*tunIter] = *tunIter;
                }
                tunIter = m_pset->m_outTuns.begin();
                for(; tunIter != m_pset->m_outTuns.end(); ++tunIter) {
                  nodeMap[*tunIter] = *tunIter;
                }
#if ALWAYSFUSE
                FuseLoops(left,right,simplifiers,cullFunc);
                didMerge = true;
#else
                Poss *newPoss = new Poss;
                newPoss->Duplicate(this,nodeMap,true);
                newPoss->PatchAfterDuplicate(nodeMap);
                //cout << "fusing " << left << " and " << right << endl;
                //throw;
                newPoss->FuseLoops(left,right,simplifiers,cullFunc);
                SetFused((Loop*)(m_sets[left]),(Loop*)(m_sets[right]));
                if (AddPossToMMap(newPosses,newPoss,newPoss->GetHash()))
                  didMerge = true;
                else
                  delete newPoss;
                break;
#endif
              }
            }
          }
          if (didMerge)
            break;
        }
      }
    }
  }
  if (!didMerge) {
    if (m_sets.size() >= 2) {
      for (unsigned int left = 0; left < m_sets.size()-1 && !didMerge; ++left) {
        //check if the poss set is transparent poss (i.e. not a loop)
        if (m_sets[left]->IsTransparent()) {
          for (unsigned int right = left + 1; right < m_sets.size(); ++right) {
            if (m_sets[right]->IsTransparent()) {
              if (m_sets[left]->CanMerge(m_sets[right])) {
                MergePosses(left,right,simplifiers,cullFunc);
                didMerge = true;
                break;
              }
            }
          }
        }
      }
    }
  }
  
  if (didMerge)
    m_fullyExpanded = false;
  
  return didMerge;
}

void Poss::MergePosses(unsigned int left, unsigned int right, const TransMap &simplifiers, CullFunction cullFunc)
{
  InvalidateHash();
  if (left >= m_sets.size() || right >= m_sets.size()) {
    cout << "left or right is too big\n";
    throw;
  }
  if (left >= right) {
    cout << "left >= right\n";
    throw;
  }
  PSet *leftSet = m_sets[left];
  PSet *rightSet = m_sets[right];
  leftSet->Cull(cullFunc);
  rightSet->Cull(cullFunc);
  if (leftSet->m_ownerPoss != this ||
      rightSet->m_ownerPoss != this) {
    cout << "Bad owner\n";
    throw;
  }
  
  NodeMap map;
  PSet *newSet = leftSet->GetNewInst();

  newSet->m_functionality = leftSet->m_functionality+rightSet->m_functionality;
  
  m_sets.erase(m_sets.begin()+left);
  if (*(m_sets.begin()+right-1) != rightSet) {
    cout << "error;";
    throw;
  }
  m_sets.erase(m_sets.begin()+right-1);
  
  newSet->m_ownerPoss = this;
  m_sets.insert(m_sets.begin()+left,newSet);
  
  NodeVecConstIter iter;
  
  //Create output set tunnels on the new set to match
  // those on the old sets
  iter  = leftSet->m_outTuns.begin();
  for (; iter != leftSet->m_outTuns.end(); ++iter) {
    PossTunnel *tun = (PossTunnel*)((*iter)->GetNewInst());
    tun->Duplicate(*iter,true,true);
    newSet->m_outTuns.push_back(tun);
    tun->m_pset = newSet;
    map[*iter] = tun;
    RemoveFromGraphNodes(*iter);
  }
  
  iter  = rightSet->m_outTuns.begin();
  for (; iter != rightSet->m_outTuns.end(); ++iter) {
    PossTunnel *tun = (PossTunnel*)((*iter)->GetNewInst());
    tun->Duplicate(*iter,true,true);
    newSet->m_outTuns.push_back(tun);
    tun->m_pset = newSet;
    map[*iter] = tun;
    RemoveFromGraphNodes(*iter);
  }
  
  //Create input set tunnels from left set
  iter = leftSet->m_inTuns.begin();
  for (; iter != leftSet->m_inTuns.end(); ++iter) {
    PossTunnel *tun = (PossTunnel*)((*iter)->GetNewInst());
    tun->Duplicate(*iter,true,true);
    newSet->m_inTuns.push_back(tun);
    tun->m_pset =  newSet;
    map[*iter] = tun;
    for (unsigned int i = 0; i < (*iter)->m_inputs.size(); ++i) {
      Node *input = (*iter)->Input(i);
      if (map[input]) {
        NodeConn *conn = new NodeConn(map[input],(*iter)->InputConnNum(i));
        tun->m_inputs.push_back(conn);
      }
      else {
        tun->AddInput((*iter)->Input(i),(*iter)->InputConnNum(i));
      }
    }
    RemoveFromGraphNodes(*iter);
  }
  
  //Create input set tunnels from right set
  iter  = rightSet->m_inTuns.begin();
  for (; iter != rightSet->m_inTuns.end(); ++iter) {
    PossTunnel *tun = (PossTunnel*)((*iter)->GetNewInst());
    tun->Duplicate(*iter,true,true);
    newSet->m_inTuns.push_back(tun);
    tun->m_pset = newSet;
    map[*iter] = tun;
    for (unsigned int i = 0; i < (*iter)->m_inputs.size(); ++i) {
      Node *input = (*iter)->Input(i);
      if (map[input]) {
        NodeConn *conn = new NodeConn(map[input],(*iter)->InputConnNum(i));
        tun->m_inputs.push_back(conn);
      }
      else {
        tun->AddInput((*iter)->Input(i),(*iter)->InputConnNum(i));
      }
    }
    RemoveFromGraphNodes(*iter);
  }
  
  for (PossMMapIter leftIter = leftSet->m_posses.begin(); leftIter != leftSet->m_posses.end(); ++leftIter) {
    for (PossMMapIter rightIter = rightSet->m_posses.begin(); rightIter != rightSet->m_posses.end(); ++rightIter) {
      Poss *newLeft = new Poss;
      Poss *newRight = new Poss;
      newLeft->Duplicate((*leftIter).second,map,true);
      newRight->Duplicate((*rightIter).second,map,true);
      
      for(unsigned int i = 0; i < newLeft->m_inTuns.size(); ++i) {
        map[newLeft->m_inTuns[i]->Input(0)]->AddChild(newLeft->m_inTuns[i],0);
      }
      for(unsigned int i = 0; i < newRight->m_inTuns.size(); ++i) {
        map[newRight->m_inTuns[i]->Input(0)]->AddChild(newRight->m_inTuns[i],0);
      }
      
      for(unsigned int i = 0; i < newLeft->m_outTuns.size(); ++i)
        map[newLeft->m_outTuns[i]->m_children[0]->m_n]->m_inputs.push_back(new NodeConn(newLeft->m_outTuns[i],0));
      
      for(unsigned int i = 0; i < newRight->m_outTuns.size(); ++i)
        map[newRight->m_outTuns[i]->m_children[0]->m_n]->m_inputs.push_back(new NodeConn(newRight->m_outTuns[i],0));
      
      while(newRight->m_sets.size()) {
        PSet *set = newRight->m_sets[0];
        set->m_ownerPoss = newLeft;
        newLeft->m_sets.push_back(set);
        newRight->m_sets.erase(newRight->m_sets.begin());
      }
      
      NodeVecIter nodeIter = newRight->m_possNodes.begin();
      for(; nodeIter != newRight->m_possNodes.end(); ++nodeIter) {
        (*nodeIter)->m_poss = newLeft;
        newLeft->m_possNodes.push_back(*nodeIter);
      }
      newLeft->m_inTuns.insert(newLeft->m_inTuns.end(),newRight->m_inTuns.begin(),newRight->m_inTuns.end());
      newLeft->m_outTuns.insert(newLeft->m_outTuns.end(),newRight->m_outTuns.begin(),newRight->m_outTuns.end());
      newLeft->m_transVec.insert(newLeft->m_transVec.end(),newRight->m_transVec.begin(),newRight->m_transVec.end());
      newLeft->PatchAfterDuplicate(map);
      newRight->m_possNodes.clear();
      newRight->m_inTuns.clear();
      newRight->m_outTuns.clear();
      
      if (newRight->m_sets.size()) {
        cout << "are the sets getting copied?\n";
        throw;
      }
      
      newSet->m_posses.insert(PossMMapPair(newLeft->GetHash(),newLeft));
      newLeft->m_pset = newSet;
      delete newRight;
    }
  }
  
  NodeVec newInputTunnelsToFix;
  
  iter  = leftSet->m_outTuns.begin();
  for (; iter != leftSet->m_outTuns.end(); ++iter) {
    PossTunnel *tun = (PossTunnel*)(map[*iter]);
    for(unsigned int i = 0; i < (*iter)->m_children.size(); ++i) {
      Node *child = (*iter)->m_children[i]->m_n;
      if (map[child] != NULL)  {
        if (!map[child]->IsPossTunnel() ||
            ((PossTunnel*)map[child])->m_tunType != SETTUNIN) {
          cout << "mapped child not of expected type\n";
          throw;
        }
        newInputTunnelsToFix.push_back(map[child]);
      }
      else {
        child->ChangeInput1Way(*iter, (*iter)->m_children[i]->m_num, tun, (*iter)->m_children[i]->m_num);
      }
      delete (*iter)->m_children[i];
      (*iter)->m_children.erase((*iter)->m_children.begin()+i);
      --i;
    }
    if (!(*iter)->m_children.empty()) {
      cout << "!(*iter)->m_children.empty() 1\n";
      cout << (*iter)->m_children.size() << endl;
      throw;
    }
  }
  
  
  iter  = rightSet->m_outTuns.begin();
  for (; iter != rightSet->m_outTuns.end(); ++iter) {
    PossTunnel *tun = (PossTunnel*)(map[*iter]);
    for(unsigned int i = 0; i < (*iter)->m_children.size(); ++i) {
      Node *child = (*iter)->m_children[i]->m_n;
      if (map[child] != NULL)  {
        if (!map[child]->IsPossTunnel() ||
            ((PossTunnel*)map[child])->m_tunType != SETTUNIN) {
          cout << "mapped child not of expected type\n";
          throw;
        }
        newInputTunnelsToFix.push_back(map[child]);
      }
      else {
        child->ChangeInput1Way(*iter, (*iter)->m_children[i]->m_num, tun, (*iter)->m_children[i]->m_num);
      }
      delete (*iter)->m_children[i];
      (*iter)->m_children.erase((*iter)->m_children.begin()+i);
      --i;
    }
    if (!(*iter)->m_children.empty()) {
      cout << "!(*iter)->m_children.empty() 1\n";
      cout << (*iter)->m_children.size() << endl;
      throw;
    }
  }
  
  
  NodeVecIter setTunInIter = newInputTunnelsToFix.begin();
  for(; setTunInIter != newInputTunnelsToFix.end(); ++setTunInIter) {
    Node *newSetInput = *setTunInIter;
    unsigned int inputInputNum = 0;
    NodeConnVecIter inputInputConIter = newSetInput->m_inputs.begin();
    NodeSet set;
    for(unsigned int i = 0; i < newSetInput->m_inputs.size(); ++i) {
      Node *newSetOutput = NULL;
      if (newSetInput->Input(i)->IsPossTunnel() && ((PossTunnel*)((*inputInputConIter)->m_n))->m_pset == newSet)
        newSetOutput = newSetInput->Input(i);
      unsigned int outputTunnelOutputNum = newSetInput->InputConnNum(i);
      if (newSetOutput) {
        if (!newSetOutput->IsPossTunnel()) {
          cout << "!newSetOutput->IsPossTunnel()\n";
          throw;
        }
        if (((PossTunnel*)newSetOutput)->m_tunType != SETTUNOUT) {
          cout << "((PossTunnel*)newSetOutput)->m_tunType != SETTUNOUT\n";
          throw;
        }
        NodeConnVecIter newSetInputChildIter = newSetInput->m_children.begin();
        NodeConnVecIter newSetOutputInputIter = newSetOutput->m_inputs.begin();
        //	jth input to the children should be wired to outputTunnelOutputNum input to input
        for(; newSetInputChildIter != newSetInput->m_children.end()
	      && newSetOutputInputIter != newSetOutput->m_inputs.end();
            ++newSetInputChildIter, ++newSetOutputInputIter)
	  {
	    Node *possInput = (*newSetInputChildIter)->m_n;
	    Node *possOutput = (*newSetOutputInputIter)->m_n;
	    Node *inputToPossOutput = possOutput->m_inputs[outputTunnelOutputNum]->m_n;
	    unsigned int inputToPossOutputNum = possOutput->m_inputs[outputTunnelOutputNum]->m_num;
	    if (possInput->m_poss != possOutput->m_poss) {
	      cout << "(possInput->m_poss != possOutput->m_poss)\n";
	      cout << possInput->m_poss << " != " << possOutput->m_poss << endl;
	      throw;
	    }
	    for(unsigned int j = 0; j < possInput->m_children.size(); ++j) {
	      if (possInput->m_children[j]->m_num == inputInputNum) {
		Node *userOfInput = possInput->m_children[j]->m_n;
		if (userOfInput->m_poss != inputToPossOutput->m_poss) {
		  cout << "userOfInput->m_poss != inputToPossOutput->m_poss\n";
		  throw;
		}
		userOfInput->ChangeInput1Way(possInput, inputInputNum, inputToPossOutput, inputToPossOutputNum);
		possInput->m_children.erase(possInput->m_children.begin()+j);
		--j;
	      }
	      else if (possInput->m_children[j]->m_num > inputInputNum) {
		set.insert(possInput->m_children[j]->m_n);
		--(possInput->m_children[j]->m_num);
	      }
	    }
	    NodeSetIter setIter = set.begin();
	    for(; setIter != set.end(); ++setIter) {
	      NodeConnVecIter inputsIter = (*setIter)->m_inputs.begin();
	      for(; inputsIter != (*setIter)->m_inputs.end(); ++inputsIter) {
		if ((*inputsIter)->m_n == possInput) {
		  if ((*inputsIter)->m_num == inputInputNum) {
		    cout << "(*inputsIter)->m_num == inputInputNum\n";
		    throw;
		  }
		  else if ((*inputsIter)->m_num > inputInputNum) {
		    --((*inputsIter)->m_num);
		  }
		}
	      }
	    }
	  }
        if (newSetInputChildIter != newSetInput->m_children.end()
            || newSetOutputInputIter != newSetOutput->m_inputs.end())
	  {
	    cout << "unbalanced posses\n";
	    throw;
	  }
        newSetInput->m_inputs.erase(newSetInput->m_inputs.begin()+i);
        --i;
        --inputInputNum;
      }
    }
  }
  
  for(unsigned int i = 0; i < leftSet->m_inTuns.size(); ++i) {
    Node *node = leftSet->m_inTuns[i];
    for(unsigned int j = 0; j < node->m_inputs.size(); ++j) {
      NodeConn *conn = node->m_inputs[j];
      if (!map[conn->m_n])
        conn->m_n->RemoveChild(node,conn->m_num);
      delete conn;
    }
    node->m_inputs.clear();
  }
  
  for(unsigned int i = 0; i < rightSet->m_inTuns.size(); ++i) {
    Node *node = rightSet->m_inTuns[i];
    for(unsigned int j = 0; j < node->m_inputs.size(); ++j) {
      NodeConn *conn = node->m_inputs[j];
      if (!map[conn->m_n])
        conn->m_n->RemoveChild(node,conn->m_num);
      delete conn;
    }
    node->m_inputs.clear();
  }
  delete leftSet;
  delete rightSet;
  
  newSet->CombineAndRemoveTunnels();
  
  NodeVecIter tunIter = newSet->m_inTuns.begin();
  for(; tunIter != newSet->m_inTuns.end(); ++tunIter) {
    AddNode(*tunIter);
  }
  
  tunIter = newSet->m_outTuns.begin();
  for(; tunIter != newSet->m_outTuns.end(); ++tunIter) {
    AddNode(*tunIter);
  }
  
  newSet->Simplify(simplifiers);
  newSet->BuildDataTypeCache();
}

bool AddNodesDown(Node *edgeStart, unsigned int childNum, NodeVec &outputTuns, NodeSet &possNodes)
{
  bool ret = false;
  NodeConn *conn = edgeStart->m_children[childNum];
  Node *child = conn->m_n;
  if (possNodes.find(child) != possNodes.end())
    return false;
  if (
#if DODM
      (child->GetNodeClass() == RedistNode::GetClass()) || 
#endif
      child->IsPossTunnel()) {
    if (child->IsPossTunnel(POSSTUNIN)) {
      cout << "Whoa!\n";
      throw;
    }
    //The child will be grouped with another poss set
    PossTunnel *tun = new PossTunnel(POSSTUNOUT);
    possNodes.insert(tun);
    outputTuns.push_back(tun);
    child->ChangeInput1Way(edgeStart, conn->m_num, tun, 0);
    conn->SetNode(tun);
    tun->m_inputs.push_back(new NodeConn(edgeStart,conn->m_num));
  }
  else {
    unsigned int i;
    for(i = 0; i < child->m_inputs.size(); ++i) {
      if (possNodes.find(child->Input(i)) == possNodes.end())
        return false;
    }
    ret = true;
    possNodes.insert(child);
    child->m_poss->RemoveFromGraphNodes(child);
    child->m_poss = NULL;
    for(unsigned int i = 0; i < child->m_children.size(); ++i)
      ret |= AddNodesDown(child, i, outputTuns, possNodes);
  }
  return ret;
}


void AddTunnelDown(Node *edgeStart, unsigned int childNum, NodeVec &outputTuns, NodeSet &possNodes)
{
  NodeConn *conn = edgeStart->m_children[childNum];
  Node *child = conn->m_n;
  if (possNodes.find(child) != possNodes.end())
    return;
  PossTunnel *tun = new PossTunnel(POSSTUNOUT);
  possNodes.insert(tun);
  outputTuns.push_back(tun);
  child->ChangeInput1Way(edgeStart, conn->m_num, tun, 0);
  conn->SetNode(tun);
  tun->m_inputs.push_back(new NodeConn(edgeStart,conn->m_num));
}

void AddPossTunnels(Node *node, Node *ignore, NodeVec &outputTuns, NodeSet &possNodes)
{
  if (possNodes.find(node) != possNodes.end())
    return;
  possNodes.insert(node);
  node->m_poss->RemoveFromGraphNodes(node);
  node->m_poss = NULL;
  NodeConnVecIter iter = node->m_inputs.begin();
  for(; iter != node->m_inputs.end(); ++iter) {
    NodeConn *conn = *iter;
    if (conn->m_n != ignore) {
#if DODM
      if (conn->m_n->GetNodeClass() == RedistNode::GetClass()) {
        AddPossTunnels(conn->m_n, node, outputTuns, possNodes);
      }
      else 
#endif
	{
	  PossTunnel *tun = new PossTunnel(POSSTUNIN);
	  possNodes.insert(tun);
	  conn->m_n->RemoveChild(node, conn->m_num);
	  tun->AddInput(conn->m_n, conn->m_num);
	  conn->SetNode(tun);
	  conn->m_num=0;
	  tun->AddChild(node, 0);
	}
    }
  }
#if DODM
  unsigned int i = 0;
  for (; i < node->m_children.size(); ++i) {
    NodeConn *conn = node->m_children[i];
    if (conn->m_n != ignore) {
      if (conn->m_n->GetNodeClass() == RedistNode::GetClass()) {
        AddPossTunnels(conn->m_n, node, outputTuns, possNodes);
      }
    }
  }
#endif
}

bool FoundLoop(Node *node, NodeVec &queue)
{
  NodeVecIter checkIter = queue.begin();
  for(; checkIter != queue.end(); ++checkIter) {
    if (*checkIter == node) {
      cout << "recursion on node " << node << " " << node->GetType() << endl;
      for(; checkIter != queue.end(); ++checkIter) {
        cout << "On queue: " << (*checkIter) << " " << (*checkIter)->GetType() << endl;
      }
      return true;
    }
  }
  if (node->IsPossTunnel(POSSTUNIN)) {
    return false;
  }
  queue.push_back(node);
  if (node->IsPossTunnel(SETTUNOUT)) {
    const PossTunnel *tunOut = (PossTunnel*)node;
    const PSet *foundSet = tunOut->m_pset;
    NodeVecConstIter iter = foundSet->m_inTuns.begin();
    for(; iter != foundSet->m_inTuns.end(); ++iter) {
      if (FoundLoop(*iter, queue)) {
        queue.pop_back();
        return true;
      }
    }
  }
  else {
    NodeConnVecConstIter iter = node->m_inputs.begin();
    for (; iter != node->m_inputs.end(); ++iter) {
      if (FoundLoop((*iter)->m_n, queue)) {
        queue.pop_back();
        return true;
      }
    }
  }
  queue.pop_back();
  return false;
}

PSet* Poss::FormSubPSet(NodeVec &outputTuns, bool isCritSect)
{
  Poss *newPoss = new Poss(outputTuns, true, true);
  PSet *set;
  if (isCritSect)
    throw;
#if 0
  if (isCritSect)
    set = new CritSect(newPoss);
  else
#endif
    set = new PSet(newPoss);
  
  AddPSet(set, true);
  
  for (unsigned int j = 0; j < set->m_inTuns.size(); ++j) {
    Node *tun = set->m_inTuns[j];
    if (!AddElemToVec(m_possNodes, tun, false))
      throw;
    tun->SetPoss(this);
  }
  
  for (unsigned int j = 0; j < set->m_outTuns.size(); ++j) {
    Node *tun = set->m_outTuns[j];
    if (!AddElemToVec(m_possNodes, tun, false))
      throw;
    tun->SetPoss(this);
  }
  
  return set;
}

void AddUsersOfLiveOutput(Node *node, unsigned int connNum, NodeSet &set)
{
  NodeConnVecIter iter = node->m_children.begin();
  for(; iter != node->m_children.end(); ++iter) {
    NodeConn *conn = *iter;
    if (conn->m_num == connNum) {
      Node *child = conn->m_n;
      if (child->IsPossTunnel(SETTUNIN)) {
        if (set.insert(child).second) {
          PossTunnel *tun = (PossTunnel*)child;
          PSet *pset = tun->m_pset;
          NodeVecIter iter2 = pset->m_inTuns.begin();
          for(; iter2 != pset->m_inTuns.end(); ++iter2)
            set.insert(*iter2);
          iter2 = pset->m_outTuns.begin();
          for(; iter2 != pset->m_outTuns.end(); ++iter2) {
            set.insert(*iter2);
          }
        }
        if (child->IsLoopTunnel()) {
          AddUsersOfLiveOutput(((LoopTunnel*)child)->GetMatchingOutTun(), 0, set);
        }
        else {
          //For normal PSets, we can't just find the matching output tunnel
          // to know along which edge the variable stays live
          //Here, we're trying to figure that out
          Node *currNode = child->Child(0);
          unsigned int numIn = 0;
          unsigned int numOut = 0;
          Node *nextNode = NULL;
          bool found = true;
          while(found) {
            found = false;
            for(unsigned int j = 0; j < currNode->m_children.size(); ++j) {
              NodeConn *childConn = currNode->m_children[j];
              if (childConn->m_num == numIn) {
                unsigned int tempOut;
                if (childConn->m_n->KeepsInputVarLive(currNode, numIn, tempOut)) {
                  if (found) {
                    if (!childConn->m_n->IsPossTunnel(POSSTUNOUT)) {
                      childConn->m_n->m_poss->PrintSetConnections();
                      throw;
                    }
                  }
                  else {
                    found = true;
                    numOut = tempOut;
                    nextNode = childConn->m_n;
                  }
                }
              }
            }
            if (found) {
              if (!nextNode)
                throw;
              currNode = nextNode;
              numIn = numOut;
              if (currNode->IsPossTunnel(POSSTUNOUT))
                found = false;
              if (currNode->IsPossTunnel(SETTUNIN)) {
                if (currNode->IsLoopTunnel()) {
                  nextNode = ((LoopTunnel*)currNode)->GetMatchingOutTun();
                  if (!nextNode)
                    throw;
                  currNode = nextNode;
                  numIn = 0;
                }
                else if (currNode->GetNodeClass() == PossTunnel::GetClass()) {
                  nextNode = currNode->Child(0);
                  if (!nextNode)
                    throw;
                  currNode = nextNode;
                  numIn = 0;
                }
              }
            }
          }
          if (currNode->IsPossTunnel(POSSTUNOUT)) {
            currNode = currNode->Child(0);
            AddUsersOfLiveOutput(currNode, 0, set);
          }
        }
      }
      else {
        if (child->IsDataDependencyOfInput() && !child->IsPossTunnel(POSSTUNOUT)) {
          set.insert(child);
          unsigned int num;
          if (conn->m_n->KeepsInputVarLive(node, connNum, num))
            AddUsersOfLiveOutput(conn->m_n, num, set);
        }
      }
    }
  }
}

bool CheckPath(Node *node, NodeVec &vec, NodeSet &set)
{
  bool ret = false;
  NodeConnVecIter iter = node->m_children.begin();
  for(; iter != node->m_children.end(); ++iter) {
    NodeConn *conn = *iter;
    Node *child = conn->m_n;
    if (child->IsPossTunnel(POSSTUNOUT)) {
      return false;
    }
    else if (set.find(child) != set.end()) {
      NodeVecIter iter2 = vec.begin();
      for(; iter2 != vec.end(); ++iter2) {
        set.insert(*iter2);
        return true;
      }
    }
    else {
      NodeVec vec2 = vec;
      vec2.push_back(child);
      ret |= CheckPath(child, vec2, set);
    }
  }
  return ret;
}

void Poss::FillClique(NodeSet &set)
{
  NodeSetIter iter = set.begin();
  NodeVec vec;
  for(; iter != set.end(); ++iter) {
    if (CheckPath(*iter, vec, set)) {
      FillClique(set);
      return;
    }
  }
}

PSet* Poss::FormSetForClique(NodeSet &set, bool isCritSect)
{
  NodeVec outputTuns;
  NodeSetIter iter = set.begin();
  for(; iter != set.end(); ++iter) {
    Node *node = *iter;
    if (!node->IsPossTunnel(SETTUNIN) && !node->IsPossTunnel(POSSTUNOUT)) {
      for(int i = 0; i < (int)(node->m_children.size()); ++i) {
        if (set.find(node->Child(i)) == set.end()) {
          unsigned int connNum = node->ChildConnNum(i);
          PossTunnel *tun;
	  if (isCritSect)
	    throw;
#if 0
	  tun = new CritSectTunnel(POSSTUNOUT);
#endif
	  else
	    tun = new PossTunnel(POSSTUNOUT);
          outputTuns.push_back(tun);
          set.insert(tun);
          for(int j = i; j < (int)(node->m_children.size()); ++j) {
            if (set.find(node->Child(j)) == set.end()) {
              node->RedirectChild(j, tun, 0);
              --j;
            }
          }
          tun->AddInput(node, connNum);
          i=-1;
        }
      }
    }
    if (!node->IsPossTunnel(SETTUNOUT) && !node->IsPossTunnel(POSSTUNIN)) {
      for(int i = 0; i < (int)(node->m_inputs.size()); ++i) {
        Node *input = node->Input(i);
        unsigned int connNum = node->InputConnNum(i);
        PossTunnel *tun = NULL;
        if (set.find(input) == set.end()) {
	  if (isCritSect)
	    throw;
#if 0
	  tun = new CritSectTunnel(POSSTUNIN);
#endif
	  else
	    tun = new PossTunnel(POSSTUNIN);
          //          cout << "creating input for child " << node->GetNodeClass() << endl;
          set.insert(tun);
          node->ChangeInput2Way(input, connNum, tun, 0);
          tun->AddInput(input, connNum);
          --i;
          for(int j = 0; j < (int)(input->m_children.size()); ++j) {
            NodeConn *conn = input->m_children[j];
            if (conn->m_num == connNum && conn->m_n != tun) {
              if (set.find(conn->m_n) != set.end()
                  && !conn->m_n->IsPossTunnel(POSSTUNIN))
		{
		  conn->m_n->ChangeInput2Way(input, connNum, tun, 0);
		}
            }
          }
        }
      }
    }
  }
  /*
    iter = set.begin();
    for(; iter != set.end(); ++iter) {
    cout << "node " << *iter << " " << (*iter)->GetNodeClass() << endl;
    cout << "  first owned by " << (*iter)->m_poss << endl;
    }
  */
  //   cout << outputTuns.size() << " outputs\n";
  //   NodeVecIter iter3 = outputTuns.begin();
  //   for(; iter3 != outputTuns.end(); ++iter3) {
  //     cout << "output node " << *iter3 << " " << (*iter3)->GetNodeClass() << endl;
  //   }
  PSet *pset = FormSubPSet(outputTuns, isCritSect);
  //   NodeVecIter iter2 = m_possNodes.begin();
  //   for(; iter2 != m_possNodes.end(); ++iter2) {
  //     if (set.find(*iter2) != set.end()) {
  //       Node *node = *iter2;
  //       cout << node <<" found in poss\n";
  //       cout <<"  " << (*iter2)->GetNodeClass() << endl;
  //     }
  //   }
  //   iter = set.begin();
  //   for(; iter != set.end(); ++iter) {
  //     cout << "node " << *iter << " " << (*iter)->GetNodeClass() << endl;
  //     cout << "  now owned by " << (*iter)->m_poss << endl;
  //   }
  return pset;
}


void Poss::FormSets(unsigned int phase)
{
#if DOSOPHASE
  if (phase == SOPHASE) {
    InvalidateHash();
    
    //It's important to recurse first; otherwise,
    // we'd form sets, recurse into them, and form similar sets again
    
    for(unsigned int i=0; i < (unsigned int)(m_sets.size()); ++i) {
      (m_sets[i])->FormSets(phase);
    }
    
    if (m_pset->m_isTopLevel)
      return;
  }
#endif //DOSOPHASE
#if DOSUMSCATTERTENSORPHASE
  if (phase == SUMSCATTERTENSORPHASE) {
   
    //       Be careful to prevent the following
    //       B = A
    //       C = B
    //       D = Op1 (C)
    //       E = F
    //       F = Op2 (B, E)
   
    //       Where E=F is in a different poss then everything else, forming a loop
    //       because of F = Op2(...)
   
   
#ifdef CHECKFORLOOPS
    NodeVec vec;
    for(unsigned int i = 0; i < m_possNodes.size(); ++i) {
      if (FoundLoop(m_possNodes[i],vec)) {
	cout << "Found loop 1\n";
	cout.flush();
	throw;
      }
      if (!vec.empty()) {
	cout << "vec not empty\n";
	throw;
      }
    }
#endif
   
    InvalidateHash();
   
   
    //It's important to recurse first; otherwise,
    // we'd form a bunch of sets of just redist nodes,
    // recurse into it and form a set of those redist nodes
    // ad infinitum
   
    int i;
    for(i=0; i < (int)(m_sets.size()); ++i) {
      (m_sets[i])->FormSets(phase);
    }
    
    if (m_pset->m_isTopLevel)
      return;
    
    for(i=0; i < (int)(m_sets.size()); ++i) {
      if (m_sets[i]->IsLoop()) {
	m_sets[i]->FormSetAround();
	i = -1;
      }
    }
   
   
    for(i=0; i < (int)(m_possNodes.size()); ++i) {
      Node *node = m_possNodes[i];
      if (!node->IsPossTunnel() && node->GetNodeClass() == RedistNode::GetClass()) {
	//Found a node that isn't a poss tunnel
	//Let's form a new set!
	NodeSet possNodes;
	NodeVec outputTuns;
	AddPossTunnels(node, NULL, outputTuns, possNodes);
	
#if DOSUMSCATTERTENSORPHASE
       	if (false && phase == SUMSCATTERTENSORPHASE) {
	  bool newNode = true;
	  do {
	    newNode = false;
	    NodeSetIter nodeIter = possNodes.begin();
	    for(; nodeIter != possNodes.end(); ++nodeIter) {
	      Node *currNode = *nodeIter;
	      if (!currNode->IsPossTunnel()) {
		bool ret = false;
		for (unsigned int j = 0; j < currNode->m_children.size(); ++j) {
		  ret |= AddNodesDown(currNode, j, outputTuns, possNodes);
		}
		if (ret) {
		  nodeIter = possNodes.begin();
		  newNode = true;
		}
	      }
	    }
	  } while (newNode);
	}
#else
	throw;
#endif
   
	NodeSetIter nodeIter = possNodes.begin();
	for(; nodeIter != possNodes.end(); ++nodeIter) {
	  Node *currNode = *nodeIter;
	  if (!currNode->IsPossTunnel()) {
	    for (unsigned int j = 0; j < currNode->m_children.size(); ++j) {
	      AddTunnelDown(currNode, j, outputTuns, possNodes);
	    }
	  }
	}
   
	Poss *newPoss = new Poss(outputTuns, true, true);
	PSet *set = new PSet(newPoss);
   
	AddPSet(set, true);
   
	for (unsigned int j = 0; j < set->m_inTuns.size(); ++j) {
	  Node *tun = set->m_inTuns[j];
	  if (!AddElemToVec(m_possNodes, tun, false))
	    throw;
	  tun->SetPoss(this);
	}
   
	for (unsigned int j = 0; j < set->m_outTuns.size(); ++j) {
	  Node *tun = set->m_outTuns[j];
	  if (!AddElemToVec(m_possNodes, tun, false))
	    throw;
	  tun->SetPoss(this);
	}
   
	i = -1;
      }
    }
   
    
//     for(i=0; i < (int)(m_possNodes.size()); ++i) {
//       Node *node = m_possNodes[i];
//       if (!node->IsPossTunnel()) {
// 	//Found a node that isn't part of a poss set
// 	//Let's form a new set!
   
// 	RemoveFromGraphNodes(node);
// 	node->m_poss = NULL;
   
// 	NodeVec outputTuns;
   
// 	NodeConnVecIter iter = node->m_inputs.begin();
// 	for(; iter != node->m_inputs.end(); ++iter) {
// 	  NodeConn *conn = *iter;
// 	  PossTunnel *tun = new PossTunnel(POSSTUNIN);
// 	  conn->m_n->RemoveChild(node, conn->m_num);
// 	  tun->AddInput(conn->m_n, conn->m_num);
// 	  conn->SetNode(tun);
// 	  conn->m_num=0;
// 	  tun->AddChild(node, 0);
// 	}
   
// 	iter = node->m_children.begin();
// 	for(; iter != node->m_children.end(); ++iter) {
// 	  NodeConn *conn = *iter;
// 	  Node *child = conn->m_n;
// 	  PossTunnel *tun = new PossTunnel(POSSTUNOUT);
// 	  outputTuns.push_back(tun);
// 	  child->ChangeInput1Way(node, conn->m_num, tun, 0);
// 	  conn->SetNode(tun);
// 	  tun->m_inputs.push_back(new NodeConn(node,conn->m_num));
// 	}
   
   
// 	Poss *newPoss = new Poss(outputTuns, true, true);
// 	PSet *set = new PSet(newPoss);
   
// 	AddPSet(set, true);
   
// 	for (unsigned int j = 0; j < set->m_inTuns.size(); ++j) {
// 	  Node *tun = set->m_inTuns[j];
// 	  if (!AddElemToVec(m_possNodes, tun, false))
// 	    throw;
// 	  tun->SetPoss(this);
// 	}
   
// 	for (unsigned int j = 0; j < set->m_outTuns.size(); ++j) {
// 	  Node *tun = set->m_outTuns[j];
// 	  if (!AddElemToVec(m_possNodes, tun, false))
// 	    throw;
// 	  tun->SetPoss(this);
// 	}
   
// 	i = -1;
//       }
//     }
   
#ifdef CHECKFORLOOPS
    for(unsigned int k = 0; k < m_outTuns.size(); ++k) {
      if (FoundLoop(m_outTuns[k],vec)) {
	cout << "Found loop 3\n";
	cout.flush();
	throw;
      }
      if (!vec.empty()) {
	cout << "vec not empty\n";
	throw;
      }
    }
#endif
    }
#endif
}

void Poss::FuseLoops(unsigned int left, unsigned int right, const TransMap &simplifiers, CullFunction cullFunc)
{
  InvalidateHash();
  
  unsigned int size = m_sets.size();
  if (left >= size || right >= size) {
    cout << "left or right is too big\n";
    throw;
  }
  if (left >= right) {
    cout << "left >= right\n";
    throw;
  }
  
  
  if (!m_sets[left]->IsLoop() || !m_sets[right]->IsLoop())
    throw;
  Loop *leftSet = (Loop*)(m_sets[left]);
  Loop *rightSet = (Loop*)(m_sets[right]);
  leftSet->Cull(cullFunc);
  rightSet->Cull(cullFunc);
  if (!leftSet->Size() || !rightSet->Size())
    throw;
  
  
  NodeMap tunMap;
  Loop *newSet = (Loop*)leftSet->GetNewInst();
  newSet->m_functionality = leftSet->m_functionality + rightSet->m_functionality;
#if TWOD
  if (leftSet->GetDimName() == rightSet->GetDimName())
    newSet->SetDimName(leftSet->GetDimName());
#endif
  newSet->m_bsSize = (((Loop*)leftSet)->m_bsSize);
  newSet->m_label.clear();
  newSet->m_label.insert(leftSet->m_label.begin(),leftSet->m_label.end());
  newSet->m_label.insert(rightSet->m_label.begin(),rightSet->m_label.end());
  
  m_sets.erase(m_sets.begin()+left);
  m_sets.erase(m_sets.begin()+right-1);
  
  newSet->m_ownerPoss = this;
  m_sets.insert(m_sets.begin()+left,newSet);
  
  NodeVecConstIter iter;
  
  //Create output set tunnels on the new set to match
  // those on the old sets
  iter  = leftSet->m_outTuns.begin();
  for (; iter != leftSet->m_outTuns.end(); ++iter) {
    PossTunnel *tun = (PossTunnel*)((*iter)->GetNewInst());
    tun->Duplicate(*iter,true,true);
    newSet->m_outTuns.push_back(tun);
    tun->m_pset = newSet;
    tunMap[*iter] = tun;
    RemoveFromGraphNodes(*iter);
  }
  
  iter  = rightSet->m_outTuns.begin();
  for (; iter != rightSet->m_outTuns.end(); ++iter) {
    PossTunnel *tun = (PossTunnel*)((*iter)->GetNewInst());
    tun->Duplicate(*iter,true,true);
    newSet->m_outTuns.push_back(tun);
    tun->m_pset = newSet;
    tunMap[*iter] = tun;
    RemoveFromGraphNodes(*iter);
  }
  
  //Create input set tunnels from left set
  iter = leftSet->m_inTuns.begin();
  for (; iter != leftSet->m_inTuns.end(); ++iter) {
    PossTunnel *tun = (PossTunnel*)((*iter)->GetNewInst());
    tun->Duplicate(*iter,true,true);
    newSet->m_inTuns.push_back(tun);
    tun->m_pset =  newSet;
    tunMap[*iter] = tun;
    for (unsigned int i = 0; i < (*iter)->m_inputs.size(); ++i) {
      Node *input = (*iter)->Input(i);
      if (tunMap[input]) {
        tun->AddInput(tunMap[input],(*iter)->InputConnNum(i));
      }
      else {
        tun->AddInput((*iter)->Input(i),(*iter)->InputConnNum(i));
      }
    }
    RemoveFromGraphNodes(*iter);
  }
  
  //Create input set tunnels from right set
  iter  = rightSet->m_inTuns.begin();
  for (; iter != rightSet->m_inTuns.end(); ++iter) {
    PossTunnel *tun = (PossTunnel*)((*iter)->GetNewInst());
    tun->Duplicate(*iter,true,true);
    newSet->m_inTuns.push_back(tun);
    tun->m_pset = newSet;
    tunMap[*iter] = tun;
    for (unsigned int i = 0; i < (*iter)->m_inputs.size(); ++i) {
      Node *input = (*iter)->Input(i);
      if (tunMap[input]) {
        tun->AddInput(tunMap[input],(*iter)->InputConnNum(i));
      }
      else {
        tun->AddInput((*iter)->Input(i),(*iter)->InputConnNum(i));
      }
    }
    RemoveFromGraphNodes(*iter);
  }
  
  for (PossMMapIter leftIter = leftSet->m_posses.begin(); leftIter != leftSet->m_posses.end(); ++leftIter) {
    for (PossMMapIter rightIter = rightSet->m_posses.begin(); rightIter != rightSet->m_posses.end(); ++rightIter) {
      Poss *newLeft = new Poss;
      Poss *newRight = new Poss;
      NodeMap map = tunMap;
      newLeft->Duplicate((*leftIter).second,map,true);
      newRight->Duplicate((*rightIter).second,map,true);
      
      for(unsigned int i = 0; i < newLeft->m_inTuns.size(); ++i) {
        map[newLeft->m_inTuns[i]->Input(0)]->AddChild(newLeft->m_inTuns[i],0);
      }
      for(unsigned int i = 0; i < newRight->m_inTuns.size(); ++i) {
        map[newRight->m_inTuns[i]->Input(0)]->AddChild(newRight->m_inTuns[i],0);
      }
      
      for(unsigned int i = 0; i < newLeft->m_outTuns.size(); ++i)
        map[newLeft->m_outTuns[i]->m_children[0]->m_n]->m_inputs.push_back(new NodeConn(newLeft->m_outTuns[i],0));
      
      for(unsigned int i = 0; i < newRight->m_outTuns.size(); ++i)
        map[newRight->m_outTuns[i]->m_children[0]->m_n]->m_inputs.push_back(new NodeConn(newRight->m_outTuns[i],0));
      
      while(newRight->m_sets.size()) {
        PSet *set = newRight->m_sets[0];
        set->m_ownerPoss = newLeft;
        newLeft->m_sets.push_back(set);
        newRight->m_sets.erase(newRight->m_sets.begin());
      }
      
      NodeVecIter nodeIter = newRight->m_possNodes.begin();
      for(; nodeIter != newRight->m_possNodes.end(); ++nodeIter) {
        (*nodeIter)->m_poss = newLeft;
        newLeft->m_possNodes.push_back(*nodeIter);
      }
      newLeft->m_inTuns.insert(newLeft->m_inTuns.end(),newRight->m_inTuns.begin(),newRight->m_inTuns.end());
      newLeft->m_outTuns.insert(newLeft->m_outTuns.end(),newRight->m_outTuns.begin(),newRight->m_outTuns.end());
      newLeft->m_transVec.insert(newLeft->m_transVec.end(),newRight->m_transVec.begin(),newRight->m_transVec.end());
      //      newLeft->m_transVec.push_back("Fused loops");
      newLeft->PatchAfterDuplicate(map);
      newRight->m_possNodes.clear();
      newRight->m_inTuns.clear();
      newRight->m_outTuns.clear();
      
      newSet->m_posses.insert(PossMMapPair(newLeft->GetHash(),newLeft));     
      newLeft->m_pset = newSet;
      delete newRight;
    }
  }
  
  NodeVec newInputTunnelsToFix;
  
  iter  = leftSet->m_outTuns.begin();
  for (; iter != leftSet->m_outTuns.end(); ++iter) {
    PossTunnel *tun = (PossTunnel*)(tunMap[*iter]);
    for(unsigned int i = 0; i < (*iter)->m_children.size(); ++i) {
      Node *child = (*iter)->m_children[i]->m_n;
      if (tunMap[child] != NULL)  {
        if (!tunMap[child]->IsLoopTunnel() ||
            ((PossTunnel*)tunMap[child])->m_tunType != SETTUNIN) {
          cout << "tunMapped child not of expected type\n";
          throw;
        }
        newInputTunnelsToFix.push_back(tunMap[child]);
      }
      else {
        child->ChangeInput1Way(*iter, (*iter)->m_children[i]->m_num, tun, (*iter)->m_children[i]->m_num);
      }
      delete (*iter)->m_children[i];
      (*iter)->m_children.erase((*iter)->m_children.begin()+i);
      --i;
    }
    if (!(*iter)->m_children.empty()) {
      cout << "!(*iter)->m_children.empty() 1\n";
      cout << (*iter)->m_children.size() << endl;
      throw;
    }
  }
  
  
  iter  = rightSet->m_outTuns.begin();
  for (; iter != rightSet->m_outTuns.end(); ++iter) {
    PossTunnel *tun = (PossTunnel*)(tunMap[*iter]);
    for(unsigned int i = 0; i < (*iter)->m_children.size(); ++i) {
      Node *child = (*iter)->m_children[i]->m_n;
      if (tunMap[child] != NULL)  {
        if (!tunMap[child]->IsPossTunnel() ||
            ((PossTunnel*)tunMap[child])->m_tunType != SETTUNIN) {
          cout << "tunMapped child not of expected type\n";
          throw;
        }
        newInputTunnelsToFix.push_back(tunMap[child]);
      }
      else {
        child->ChangeInput1Way(*iter, (*iter)->m_children[i]->m_num, tun, (*iter)->m_children[i]->m_num);
      }
      delete (*iter)->m_children[i];
      (*iter)->m_children.erase((*iter)->m_children.begin()+i);
      --i;
    }
    if (!(*iter)->m_children.empty()) {
      cout << "!(*iter)->m_children.empty() 1\n";
      cout << (*iter)->m_children.size() << endl;
      throw;
    }
  }
  
  NodeVecIter setTunInIter = newInputTunnelsToFix.begin();
  for(; setTunInIter != newInputTunnelsToFix.end(); ++setTunInIter) {
    if (!(*setTunInIter)->IsLoopTunnel())
      throw;
    LoopTunnel *rightSetInput = (LoopTunnel*)(*setTunInIter);
    LoopTunnel *leftSetOutput = NULL;
    if (rightSetInput->m_inputs.size() != 1)
      throw;
    if (rightSetInput->Input(0)->IsPossTunnel() && ((PossTunnel*)rightSetInput->Input(0))->m_pset == newSet) {
      LoopTunnel *newOutputToUse = ((LoopTunnel*)rightSetInput)->GetMatchingOutTun();
      if (rightSetInput->InputConnNum(0) > 0)
        throw;
      leftSetOutput = (LoopTunnel*)(rightSetInput->Input(0));
      if (!leftSetOutput->IsLoopTunnel()) {
        cout << "!leftSetOutput->IsLoopTunnel()\n";
        throw;
      }
      if (((PossTunnel*)leftSetOutput)->m_tunType != SETTUNOUT) {
        cout << "((PossTunnel*)leftSetOutput)->m_tunType != SETTUNOUT\n";
        throw;
      }
      ClassType rightInputType = rightSetInput->GetNodeClass();
      ClassType leftOutputType = leftSetOutput->GetNodeClass();
      if ((rightInputType == Split::GetClass() && leftOutputType != Combine::GetClass())
          || (leftOutputType == Combine::GetClass() && rightInputType != Split::GetClass()))
	{
	  LoopTunnel *inTun = (LoopTunnel*)leftSetOutput->GetMatchingInTun();
	  rightSetInput->ChangeInput1Way(leftSetOutput,0,inTun->Input(0),inTun->InputConnNum(0));
	  leftSetOutput->RemoveChild(rightSetInput,0);
	  leftSetOutput->CopyTunnelInfo(rightSetInput);
	  leftSetOutput->GetMatchingInTun()->CopyTunnelInfo(rightSetInput);
	}
      else {
        NodeConnVecIter rightSetInputChildIter = rightSetInput->m_children.begin();
        NodeConnVecIter leftSetOutputInputIter = leftSetOutput->m_inputs.begin();
        //	jth input to the children should be wired to outputTunnelOutputNum input to input
        for(; rightSetInputChildIter != rightSetInput->m_children.end()
	      && leftSetOutputInputIter != leftSetOutput->m_inputs.end();
            ++rightSetInputChildIter, ++leftSetOutputInputIter)
	  {
	    Node *rightPossInput = (*rightSetInputChildIter)->m_n;
	    Node *leftPossOutput = (*leftSetOutputInputIter)->m_n;
	    Poss *poss = (Poss*)rightPossInput->m_poss;
	    if (rightPossInput->m_poss != leftPossOutput->m_poss) {
	      cout << "(rightPossInput->m_poss != leftPossOutput->m_poss)\n";
	      cout << rightPossInput->m_poss << " != " << leftPossOutput->m_poss << endl;
	      throw;
	    }
	    while (rightPossInput->m_children.size()) {
	      NodeConn *childConn = rightPossInput->m_children[0];
	      Node *inputToPossOutput = leftPossOutput->Input(childConn->m_num);
	      unsigned int inputToPossOutputNum = leftPossOutput->InputConnNum(childConn->m_num);
	      Node *userOfInput = childConn->m_n;
	      if (userOfInput->m_poss != inputToPossOutput->m_poss) {
		cout << "userOfInput->m_poss != inputToPossOutput->m_poss\n";
		throw;
	      }
	      userOfInput->ChangeInput1Way(rightPossInput, childConn->m_num, inputToPossOutput, inputToPossOutputNum);
	      delete childConn;
	      rightPossInput->m_children.erase(rightPossInput->m_children.begin());
	    }
	    delete leftPossOutput->m_children[0];
	    leftPossOutput->m_children.erase(leftPossOutput->m_children.begin());
	    poss->DeleteChildAndCleanUp(leftPossOutput, false, true);
	    if (leftPossOutput->m_children.size() || leftPossOutput->m_inputs.size())
	      throw;
          
	    if (rightPossInput->m_children.size())
	      throw;
	    poss->DeleteNode(rightPossInput);
	  }
        
        if (rightSetInputChildIter != rightSetInput->m_children.end()
            || leftSetOutputInputIter != leftSetOutput->m_inputs.end())
	  {
	    cout << "unbalanced posses\n";
	    throw;
	  }
        leftSetOutput->RemoveChild(rightSetInput,0);
        if (leftSetOutput->m_children.size())
          leftSetOutput->RedirectChildren(newOutputToUse);
        delete leftSetOutput;
        newSet->RemoveOutTun(leftSetOutput);
        newSet->RemoveInTun(rightSetInput);
        rightSetInput->m_inputs.erase(rightSetInput->m_inputs.begin());
        if (rightSetInput->m_inputs.size()) {
          //	rightSetInput->m_inputs.erase(rightSetInput->m_inputs.begin());
          throw;
        }
        delete rightSetInput;
      }
    }
  }
  
  for(unsigned int i = 0; i < leftSet->m_inTuns.size(); ++i) {
    Node *node = leftSet->m_inTuns[i];
    for(unsigned int j = 0; j < node->m_inputs.size(); ++j) {
      NodeConn *conn = node->m_inputs[j];
      if (!tunMap[conn->m_n])
        conn->m_n->RemoveChild(node,conn->m_num);
      delete conn;
    }
    node->m_inputs.clear();
  }
  
  for(unsigned int i = 0; i < rightSet->m_inTuns.size(); ++i) {
    Node *node = rightSet->m_inTuns[i];
    for(unsigned int j = 0; j < node->m_inputs.size(); ++j) {
      NodeConn *conn = node->m_inputs[j];
      if (!tunMap[conn->m_n])
        conn->m_n->RemoveChild(node,conn->m_num);
      delete conn;
    }
    node->m_inputs.clear();
  }
  
  delete leftSet;
  delete rightSet;
  
  newSet->CombineAndRemoveTunnels();
  
  NodeVecIter tunIter = newSet->m_inTuns.begin();
  bool foundControl = false;
  for(; tunIter != newSet->m_inTuns.end(); ++tunIter) {
    Node *tun = *tunIter;
    AddNode(tun);
    if (tun->GetNodeClass() == Split::GetClass()) {
      Split *split = (Split*)tun;
      if (split->m_isControlTun) {
        if (foundControl)
          split->m_isControlTun = false;
        else
          foundControl = true;
      }
    }
  }
  if (!foundControl)
    throw;
  
  tunIter = newSet->m_outTuns.begin();
  for(; tunIter != newSet->m_outTuns.end(); ++tunIter) {
    AddNode(*tunIter);
  }
  
  for(unsigned int i = 0; i < newSet->m_inTuns.size(); ++i) {
    Node *tun = newSet->m_inTuns[i];
    if (tun->GetNodeClass() == Split::GetClass()) {
      for(unsigned int j = i+1; j < newSet->m_inTuns.size(); ++j) {
        Node *tun2 = newSet->m_inTuns[j];
        if (tun2->GetNodeClass() == Split::GetClass()) {
          NodeConn *conn1 = tun->m_inputs[0];
          NodeConn *conn2 = tun2->m_inputs[0];
          if (conn1->m_n == conn2->m_n && conn1->m_num == conn2->m_num) {
            Split *split1 = (Split*)tun;
            Split *split2 = (Split*)tun2;
#if TWOD
            if (split1->m_dir == split2->m_dir)
#else
	      if (split1->m_partDim == split2->m_partDim)
#endif
		throw;
	      else {
		split1->SetAddDir();
		split2->SetAddDir();
	      }
          }
        }
      }
    }
  }
  
  newSet->Simplify(simplifiers);
  //adding below.  Seems necessary
  newSet->BuildDataTypeCache();
  BuildDataTypeCache();
}

void Poss::ClearBeforeProp()
{
  m_isSane = true;
  NodeVecIter nodeIter = m_possNodes.begin();
  for( ; nodeIter != m_possNodes.end(); ++nodeIter) {
    (*nodeIter)->ClearBeforeProp();
  }
  PSetVecIter iter = m_sets.begin();
  for(; iter != m_sets.end(); ++iter)
    (*iter)->ClearBeforeProp();
}

void Poss::ClearFullyExpanded()
{
  PSetVecIter iter = m_sets.begin();
  for(; iter != m_sets.end(); ++iter)
    (*iter)->ClearFullyExpanded();
  m_fullyExpanded = false;
}

Cost Poss::EvalCurr(TransConstVec &transList)
{
  unsigned int numPSets = m_sets.size();
  Cost tot = 0;
  
  TransVecIter iter2 = m_transVec.begin();
  for(; iter2 != m_transVec.end(); ++iter2) {
    transList.push_back(*iter2);
  }
  
  NodeVecConstIter iter = m_possNodes.begin();
  for( ; iter != m_possNodes.end(); ++iter) {
    DLANode *node =  (DLANode*)(*iter);
    double cost = node->GetCost();
    tot += cost;
  }
  for(unsigned int i = 0; i < numPSets; ++i) {
    tot += m_sets[i]->EvalCurrPoss(transList);
  }
  
  return tot;
}

Cost Poss::EvalAndSetBest()
{
  unsigned int numPSets = m_sets.size();
  Cost tot = 0;
  
  NodeVecConstIter iter = m_possNodes.begin();
  for( ; iter != m_possNodes.end(); ++iter) {
    DLANode *node =  (DLANode*)(*iter);
    double cost = node->GetCost();
    tot += cost;
  }

  for(unsigned int i = 0; i < numPSets; ++i) {
    tot += m_sets[i]->EvalAndSetBest();
  }
  
  return tot;
}

void Poss::Print(IndStream &out, unsigned int &graphNum)
{
  unsigned int numPSets = m_sets.size();
  
  NodeVecConstIter nodeIter = m_inTuns.begin();
  for(; nodeIter != m_inTuns.end(); ++nodeIter) {
    (*nodeIter)->Print(out, graphNum);
    /*
      if ((*nodeIter)->GetNodeClass() == Split::GetClass()) {
      const Sizes *m = ((DLANode*)(*nodeIter))->GetM(1);
      const Sizes *n = ((DLANode*)(*nodeIter))->GetN(1);
      m->Print();
      n->Print();
      }
    */
    if (!(*nodeIter)->HasPrinted()) {
      cout << "tunnel input " << (*nodeIter)->GetType()
	   << "hasn't printed even though he should have\n";
    }
  }
  out.Indent();
  *out << "//------------------------------------//\n" << endl;
  bool hasPrinted = true;
  while(hasPrinted) {
    hasPrinted = false;
    NodeVecConstIter nodeIter = m_possNodes.begin();
    for( ; nodeIter != m_possNodes.end(); ++nodeIter) {
      Node *node = *nodeIter;
      //Don't print the poss out tunnels until the end
      // so the repartitioning code all goes after the loop body
      if (!node->HasPrinted()
          && !node->IsPossTunnel(POSSTUNOUT)
          && node->CanPrintCode())
	{
	  (*nodeIter)->Print(out, graphNum);
	  hasPrinted |= (*nodeIter)->HasPrinted();
	}
    }
    for(unsigned int i = 0; i < numPSets; ++i) {
      if (m_sets[i]->GetCurrPoss()->CanPrint()) {
        m_sets[i]->PrintCurrPoss(out, graphNum);
        hasPrinted = true;
      }
    }
  }
  
  
  *out << endl;
  out.Indent();
  *out << "//------------------------------------//" << endl;
  
  nodeIter = m_outTuns.begin();
  for(; nodeIter != m_outTuns.end(); ++nodeIter) {
    (*nodeIter)->Print(out, graphNum);
    (*nodeIter)->SetPrinted();
  }
  
  nodeIter = m_pset->m_outTuns.begin();
  for(; nodeIter != m_pset->m_outTuns.end(); ++nodeIter) {
    (*nodeIter)->Print(out, graphNum);
    (*nodeIter)->SetPrinted();
  }
  *out << endl;
  
  bool bad = false;
  
  nodeIter = m_inTuns.begin();
  for(; nodeIter != m_inTuns.end(); ++nodeIter) {
    if (!(*nodeIter)->HasPrinted()) {
      cout << (*nodeIter)->GetType() << " hasn't printed\n";
      bad = true;
    }
  }
  
  for(unsigned int i = 0; i < numPSets; ++i) {
    if (!m_sets[i]->GetCurrPoss()->m_hasPrinted) {
      cout << "set " << m_sets[i]->GetCurrPoss() << " not printed\n";
      for(unsigned int j = 0; j < m_sets[i]->m_inTuns.size(); ++j) {
        Node *tun = m_sets[i]->m_inTuns[j];
        cout << "in tun " << tun << endl;
        if (!tun->CanPrintCode())
          tun->CanPrintCode();
      }
      m_sets[i]->GetCurrPoss()->CanPrint();
      m_sets[i]->GetCurrPoss()->ForcePrint();
      bad = true;
    }
  }
  
  nodeIter = m_possNodes.begin();
  for(; nodeIter != m_possNodes.end(); ++nodeIter) {
    if (!(*nodeIter)->HasPrinted()) {
      cout << (*nodeIter)->GetType() << " " << *nodeIter << " hasn't printed\n";
      cout << "Inputs are\n";
      (*nodeIter)->PrintInputs();
      cout << "Is it possible that the node is read only but ReadOnly doesn't return true?\n\n";
      bad = true;
      //      PrintSetConnections();
    }
  }
  
  if (bad) {
    cout << this << " is bad\n";
    cout << "contains " << m_sets.size() << " posses\n";
    throw;
  }
  
  m_hasPrinted = true;
}

void Poss::EvalRoot(IndStream &out, unsigned int &graphNum, unsigned int whichGraph, unsigned int &optGraphs, double &optCosts)
{
  bool keepGoing = true;
  
  ClearCurrPoss();
  
  while (keepGoing) {
    if (whichGraph <= 0 || whichGraph == graphNum) {
      TransConstVec transList;
      Cost tot = EvalCurr(transList);
#ifdef MATLAB
      *out << "cost(" << graphNum << ") = "
	   << setprecision(15) << tot << ";\n";
#else
      *out << "cost[" << graphNum << ",1] = "
	   << setprecision(15) << tot << ";\n";
#endif
      
#ifdef MATLAB
      *out << "refs(" << graphNum << ") = {[ ";
      TransConstVecIter iter = transList.begin();
      for(; iter != transList.end(); ++iter) {
        *out << (size_t)(*iter) << " ";
      }
      *out << "]};\n";
      if (!(graphNum % 1000)) {
        *out << "'loaded " << graphNum << "'\n";
      }
#endif
      if (optCosts <= 0 || tot < optCosts) {
        optCosts = tot;
        optGraphs = graphNum;
      }
    }
    
    ++graphNum;
    keepGoing = !IncrementCurrPoss();
  }
}

void Poss::PrintRoot(IndStream &out, unsigned int &graphNum, unsigned int whichGraph)
{
  bool keepGoing = true;
  unsigned int numPSets = m_sets.size();
  
  ClearCurrPoss();
  
  while (keepGoing) {
    if (whichGraph <= 0 || whichGraph == graphNum) {
      ClearPrinted();
      
      *out << "/*** Algorithm " << graphNum << " ***" << endl;
      *out << "\tUnique Num: " << m_num << endl;
      *out << "\tChild of: " << m_parent << endl;
      *out << "\tResult of transformations:" << endl;
      TransVec transVec;
      GetCurrTransVec(transVec);
      TransVecConstIter transIter = transVec.begin();
      for( ; transIter != transVec.end(); ++transIter)
        *out << "\t" << (*transIter)->GetType() << endl;
      *out << "*****************************************/" << endl;

      VarSet set;
      AddCurrPossVars(set);
      VarSetIter varIter = set.begin();
#if DOTENSORS
      out.Indent();
      *out << "ObjShape tempShape;\n";
#endif
      for(; varIter != set.end(); ++varIter) {
	(*varIter).PrintDecl(out);
      }
      
      if (m_pset && m_pset->IsLoop()
          && ((Loop*)m_pset)->GetType() == BLISLOOP)
	{
	  string loopLevel = out.LoopLevel(1);
	  string idx = "idx" + loopLevel;
	  string dimLen = "dimLen" + loopLevel;
	  string bs = "bs" + loopLevel;
	  out.Indent();
	  *out << "dim_t " << idx << ", " << dimLen << ", " << bs << ";\n";
	}
      
      
      //This actualy sets some stuff so it can print
      if (!CanPrint()) {
        cout << "couldn't print\n";
        throw;
      }
      
      NodeVecConstIter nodeIter = m_inTuns.begin();
      for(; nodeIter != m_inTuns.end(); ++nodeIter) {
        (*nodeIter)->Print(out, graphNum);
      }
      bool hasPrinted = true;
      while(hasPrinted) {
        hasPrinted = false;
        NodeVecConstIter nodeIter = m_possNodes.begin();
        for( ; nodeIter != m_possNodes.end(); ++nodeIter) {
          //Don't bring the poss out tunnels until the end
          // so the repartitioning code all goes after the loop body
          if (!(*nodeIter)->HasPrinted() && !(*nodeIter)->IsPossTunnel(POSSTUNOUT)) {
            (*nodeIter)->Print(out, graphNum);
            hasPrinted |= (*nodeIter)->HasPrinted();
          }
        }
        for(unsigned int i = 0; i < numPSets; ++i) {
          if (m_sets[i]->GetCurrPoss()->CanPrint()) {
            m_sets[i]->PrintCurrPoss(out, graphNum);
            hasPrinted = true;
          }
        }
      }
      
      nodeIter = m_outTuns.begin();
      for(; nodeIter != m_outTuns.end(); ++nodeIter) {
        (*nodeIter)->Print(out, graphNum);
        (*nodeIter)->SetPrinted();
      }
      
      nodeIter = m_pset->m_outTuns.begin();
      for(; nodeIter != m_pset->m_outTuns.end(); ++nodeIter) {
        (*nodeIter)->Print(out, graphNum);
        (*nodeIter)->SetPrinted();
      }
      *out << endl;
      
      out.Indent();
      *out << "/*****************************************/" << endl;
    }
    
    ++graphNum;
    keepGoing = !IncrementCurrPoss();
  }
  
  m_hasPrinted = true;
}

void Poss::PrintCurrRoot(IndStream &out)
{
  unsigned int numPSets = m_sets.size();
  
  ClearPrinted();
  
  *out << "\tUnique Num: " << m_num << endl;
  *out << "\tChild of: " << m_parent << endl;
  *out << "\tResult of transformations:" << endl;
  TransVec transVec;
  GetCurrTransVec(transVec);
  TransVecConstIter transIter = transVec.begin();
  for( ; transIter != transVec.end(); ++transIter)
    *out << "\t" << (*transIter)->GetType() << endl;
  *out << "*****************************************" << endl;

  VarSet set;
  AddCurrPossVars(set);
  VarSetIter varIter = set.begin();
#if DOTENSORS
  out.Indent();
  *out << "ObjShape tempShape;\n";
#endif
  for(; varIter != set.end(); ++varIter) {
    (*varIter).PrintDecl(out);
  }

  
  if (m_pset && m_pset->IsLoop()
      && ((Loop*)m_pset)->GetType() == BLISLOOP)
    {
      string loopLevel = out.LoopLevel(1);
      string idx = "idx" + loopLevel;
      string dimLen = "dimLen" + loopLevel;
      string bs = "bs" + loopLevel;
      out.Indent();
      *out << "dim_t " << idx << ", " << dimLen << ", " << bs << ";\n";
    }
      
      
  //This actualy sets some stuff so it can print
  if (!CanPrint()) {
    cout << "couldn't print\n";
    throw;
  }

  unsigned int graphNum = -1;
      
  NodeVecConstIter nodeIter = m_inTuns.begin();
  for(; nodeIter != m_inTuns.end(); ++nodeIter) {
    (*nodeIter)->Print(out, graphNum);
  }
  bool hasPrinted = true;
  while(hasPrinted) {
    hasPrinted = false;
    NodeVecConstIter nodeIter = m_possNodes.begin();
    for( ; nodeIter != m_possNodes.end(); ++nodeIter) {
      //Don't print the poss out tunnels until the end
      // so the repartitioning code all goes after the loop body
      if (!(*nodeIter)->HasPrinted() && !(*nodeIter)->IsPossTunnel(POSSTUNOUT)) {
	(*nodeIter)->Print(out, graphNum);
	hasPrinted |= (*nodeIter)->HasPrinted();
      }
    }
    for(unsigned int i = 0; i < numPSets; ++i) {
      if (m_sets[i]->GetCurrPoss()->CanPrint()) {
	m_sets[i]->PrintCurrPoss(out, graphNum);
	hasPrinted = true;
      }
    }
  }
      
  nodeIter = m_outTuns.begin();
  for(; nodeIter != m_outTuns.end(); ++nodeIter) {
    (*nodeIter)->Print(out, graphNum);
    (*nodeIter)->SetPrinted();
  }
  
  nodeIter = m_pset->m_outTuns.begin();
  for(; nodeIter != m_pset->m_outTuns.end(); ++nodeIter) {
    (*nodeIter)->Print(out, graphNum);
    (*nodeIter)->SetPrinted();
  }
  *out << endl;
  
  out.Indent();
  *out << "*****************************************" << endl;
}

void Poss::ClearPrinted()
{
  NodeVecConstIter nodeIter = m_possNodes.begin();
  for( ; nodeIter != m_possNodes.end(); ++nodeIter)
    (*nodeIter)->ClearPrinted();
  PSetVecIter iter = m_sets.begin();
  for(; iter != m_sets.end(); ++iter)
    (*iter)->ClearPrinted();
  m_hasPrinted = false;
}

unsigned int Poss::TotalCount() const
{
  unsigned int tot = 1;
  PSetVecConstIter iter = m_sets.begin();
  for(; iter != m_sets.end(); ++iter)
    tot *= (*iter)->TotalCount();
  return tot;
}

void Poss::GetTransVec(TransVec &transVec) const
{
  transVec.insert(transVec.end(),m_transVec.begin(), m_transVec.end());
  if (m_sets.size()) {
    cout << "not supported for hierarchical posses\n";
    throw;
  }
}

void Poss::GetCurrTransVec(TransVec &transVec) const
{
  transVec.insert(transVec.end(),m_transVec.begin(), m_transVec.end());
  for(unsigned int i = 0; i < m_sets.size(); ++i) {
    m_sets[i]->GetCurrTransVec(transVec);
  }
}

bool Poss::TakeIter(const TransMap &transMap, const TransMap &simplifiers,
                    PossMMap &newPosses)
{
  bool didSomething = false;
  PSetVecIter iter = m_sets.begin();
  for (; iter != m_sets.end(); ++iter) {
    didSomething |= (*iter)->TakeIter(transMap, simplifiers);
  }
  if (!didSomething) {
    NodeMap setTunnels;
    NodeVecIter iter2 = m_pset->m_inTuns.begin();
    for(; iter2 != m_pset->m_inTuns.end(); ++iter2)
      setTunnels[*iter2] = *iter2;
    iter2 = m_pset->m_outTuns.begin();
    for(; iter2 != m_pset->m_outTuns.end(); ++iter2)
      setTunnels[*iter2] = *iter2;
    for(unsigned int nodeIdx = 0; nodeIdx < m_possNodes.size(); ++nodeIdx) {
      Node *node = m_possNodes[nodeIdx];
      TransMapConstIter transMapIter = transMap.find(node->GetNodeClass());
      if (transMapIter != transMap.end()) {
        TransVecConstIter transIter = transMapIter->second->begin();
        for(; transIter != transMapIter->second->end(); ++transIter) {
          const Transformation *trans = *transIter;
	  if (trans->IsSingle()) {
	    if (node->HasApplied(trans)) {
	      //            cout << "skipping because " << trans->GetType() << " has been applied\n";
	      continue;
	    }
	    const SingleTrans *single = (SingleTrans*)trans;
	    if (single->CanApply(node)) {
	      node->Applied(single);
	      if (single->IsRef())
		node->SetHasRefined();
	      Poss *newPoss = new Poss;
	      NodeMap nodeMap = setTunnels;
	      newPoss->Duplicate(this,nodeMap,false);
	      newPoss->PatchAfterDuplicate(nodeMap);
	      Node *newNode = nodeMap[node];
	      single->Apply(newNode);
	      newPoss->m_transVec.push_back(const_cast<Transformation*>(trans));
	      newPoss->Simplify(simplifiers);
	      //newPoss->BuildDataTypeCache();
	      if(!AddPossToMMap(newPosses,newPoss,newPoss->GetHash())) {
		delete newPoss;
	      }
	      else {
		didSomething = true;
	      }
	    }
	  }
	  else if (trans->IsVarRef()) {
	    const VarTrans *var = (VarTrans*)trans;
	    void *cache = NULL;
	    if (node->HasApplied(var))
	      continue;
	    int count = var->CanApply(node, &cache);
	    if (count > 0) {
	      node->Applied(var);
	      
	      if (trans->IsRef())
		node->SetHasRefined();
	      
	      for (int i = 0; i < count; ++ i) {
		//Need to check if this transformationhas been applied since this graph
		// could have been loaded from disk, and the transformation could have been
		// applied previously, but a different MultiTrans was marked as applied on the node
		const Transformation *marking = (var->IsMultiRef() ? 
						 ((MultiTrans*)var)->GetTrans(&cache,i) :
						 var );
	      
		Poss *newPoss = new Poss;
		NodeMap nodeMap = setTunnels;
		newPoss->Duplicate(this,nodeMap,false);
		newPoss->PatchAfterDuplicate(nodeMap);
		Node *newNode = nodeMap[node];
		var->Apply(i, newNode, &cache);
		newPoss->m_transVec.push_back(const_cast<Transformation*>(marking));
		newPoss->Simplify(simplifiers);
		//newPoss->BuildDataTypeCache();
		if(!AddPossToMMap(newPosses,newPoss,newPoss->GetHash())) {
		  delete newPoss;
		}
		else {
		  didSomething = true;
		}
	      }
	    }
	  }
	  else {
	    throw;
	  }
	}
      }
    }
    m_fullyExpanded = true;
  }
  if (!didSomething)
    m_fullyExpanded = true;
  else {
    InvalidateHash();
  }
  return didSomething;
}

bool Poss::GlobalSimplification(const TransMap &globalSimplifiers, const TransMap &simplifiers)
{
  if (!globalSimplifiers.empty()) {
    //recurse first since a global simplifier could affect this poss's nodes
    // (that also means this can't be easily parallelized)
    bool didSomething = false;
    PSetVecIter iter = m_sets.begin();
    for(; iter != m_sets.end(); ++iter) {
      didSomething |= (*iter)->GlobalSimplification(globalSimplifiers, simplifiers);
    }
    didSomething |= Simplify(globalSimplifiers);
    didSomething |= Simplify(simplifiers);
    return didSomething;
  }
  else
    return false;
}

string GetFusedString(const IntSet *set)
{
  std::stringstream str;
  
  IntSetConstIter iter = set->begin();
  for(; iter != set->end(); ++iter) {
    str << *iter << ";";
  }
  
  return str.str();
}

bool Poss::HasFused(const Loop *left, const Loop *right) const
{
  IntSet fusedSet(left->m_label.begin(),left->m_label.end());
  fusedSet.insert(right->m_label.begin(),right->m_label.end());
  string str = GetFusedString(&fusedSet);
  return M_fusedSets.find(str) != M_fusedSets.end();
}

void Poss::SetFused(const Loop *left, const Loop *right)
{
  IntSet fusedSet(left->m_label.begin(),left->m_label.end());
  fusedSet.insert(right->m_label.begin(),right->m_label.end());
  string str = GetFusedString(&fusedSet);
  M_fusedSets.insert(str);
}


void Poss::ClearCurrPoss()
{
  PSetVecIter iter = m_sets.begin();
  for(; iter != m_sets.end(); ++iter) {
    (*iter)->ClearCurrPoss();
  }
}

bool Poss::IncrementCurrPoss()
{
  unsigned int i;
  for (i = 0; i < m_sets.size(); ++i) {
    if (!m_sets[i]->IncrementCurrPoss())
      return false;
  }
  return true;
}

string Poss::GetFunctionalityString() const
{
  string str;
  NodeVecConstIter iter = m_outTuns.begin();
  for(; iter != m_outTuns.end(); ++iter) {
    const Node *node = *iter;
    str += node->GetFunctionalityString();
  }
  return str;
}

size_t Poss::Hash(const string &str)
{
  static std::hash<std::string> hasher;
  size_t tmp = 0;
  for (int i = 0; i < str.size(); i+=10) {
    tmp += hasher(str.substr(i,10));
  }
  return tmp;
}

size_t Poss::GetHash()
{
  if (m_hashValid)
    return m_hash;
  else {
    m_hash = Hash(GetFunctionalityString());
    m_hashValid = true;
    return m_hash;
  }
}

void Poss::RemoveFromGraphNodes(Node *node)
{
  InvalidateHash();
  NodeVecIter iter = m_possNodes.begin();
  for(; iter != m_possNodes.end(); ++iter) {
    if (*iter == node) {
      m_possNodes.erase(iter);
      node->m_poss = NULL;
      return;
    }
  }
  cout << "node not found\n";
  throw;
}

void Poss::RemoveFromSets(PSet *set)
{
  InvalidateHash();
  PSetVecIter iter = m_sets.begin();
  for(; iter != m_sets.end(); ++iter) {
    if (*iter == set) {
      m_sets.erase(iter);
      set->m_ownerPoss = NULL;
      return;
    }
  }
  cout << "set not found\n";
  throw;
}

void Poss::PrintNodeAddresses() const
{
  cout << "Nodes on " << this << endl;
  for(unsigned int i = 0; i < m_possNodes.size(); ++i)
    cout << m_possNodes[i]->GetNodeClass() << " " << m_possNodes[i] << endl;
}

void Poss::Flatten(ofstream &out) const
{
  unsigned int size;
  WRITE(START);
  WRITE(m_pset);
  size = m_transVec.size();
  WRITE(size);
  TransVecConstIter iter2 = m_transVec.begin();
  for(; iter2 != m_transVec.end(); ++iter2)
    WRITE(*iter2);
  size = m_sets.size();
  WRITE(size);
  PSetVecConstIter iter4 = m_sets.begin();
  for(; iter4 != m_sets.end(); ++iter4) {
    bool isLoop = (*iter4)->IsLoop();
    WRITE(isLoop);
    WRITE(*iter4);
  }
  FullyFlatten(m_possNodes, out);
  size = m_inTuns.size();
  WRITE(size);
  NodeVecConstIter iter = m_inTuns.begin();
  for(; iter != m_inTuns.end(); ++iter)
    WRITE(*iter);
  size = m_outTuns.size();
  WRITE(size);
  iter = m_outTuns.begin();
  for(; iter != m_outTuns.end(); ++iter)
    WRITE(*iter);
  iter4 = m_sets.begin();
  for(; iter4 != m_sets.end(); ++iter4) {
    (*iter4)->Flatten(out);
  }
  WRITE(END);
}

void Poss::FlattenStatic(ofstream &out)
{
  unsigned int size = M_fusedSets.size();
  WRITE(size);
  StrSetConstIter iter = M_fusedSets.begin();
  for(; iter != M_fusedSets.end(); ++iter)
    out << *iter << endl;
}


void Poss::Unflatten(ifstream &in, SaveInfo &info)
{
  char tmp;
  READ(tmp);
  if (tmp != START)
    throw;
  unsigned int size;
  READ(m_pset);
  Swap(&m_pset,info.psetMap);
  READ(size);
  for(unsigned int i = 0; i < size; ++i) {
    Transformation *trans;
    READ(trans);
    Swap(&trans,info.transMap);
    m_transVec.push_back(trans);
  }
  READ(size);
  for(unsigned int i = 0; i < size; ++i) {
    bool isLoop;
    READ(isLoop);
    PSet *newSet;
    if (isLoop)
      newSet = new Loop;
    else
      newSet = new PSet;
    PSet *oldSet;
    READ(oldSet);
    (*(info.psetMap))[oldSet] = newSet;
    m_sets.push_back(newSet);
  }
  FullyUnflatten(m_possNodes, in, info);
  READ(size);
  for(unsigned int i = 0; i < size; ++i) {
    Node *tun;
    READ(tun);
    Swap(&tun,info.nodeMap);
    m_inTuns.push_back(tun);
  }
  READ(size);
  for(unsigned int i = 0; i < size; ++i) {
    Node *tun;
    READ(tun);
    Swap(&tun,info.nodeMap);
    m_outTuns.push_back(tun);
  }
  PSetVecIter iter = m_sets.begin();
  for(; iter != m_sets.end(); ++iter) {
    (*iter)->Unflatten(in,info);
  }
  READ(tmp);
  if (tmp != END)
    throw;
  PatchAfterDuplicate(*(info.nodeMap));
}

void Poss::UnflattenStatic(ifstream &in)
{
  unsigned int size;
  READ(size);
  for (unsigned int i = 0; i < size; ++i) {
    string str;
    getline(in, str);
    M_fusedSets.insert(str);
  }
}

void Poss::BuildDataTypeCache()
{
  PSetVecIter setIter;
#if DOBLIS
  setIter = m_sets.begin();
  for(; setIter != m_sets.end(); ++setIter) {
    if ((*setIter)->IsCritSect())
      (*setIter)->BuildDataTypeCache();
  }
#endif
  NodeVecIter iter1 = m_possNodes.begin();
  for(; iter1 != m_possNodes.end(); ++iter1) {
    Node *node = *iter1;
    node->m_flags &= ~BUILDFLAG;
  }
  iter1 = m_possNodes.begin();
  for(; iter1 != m_possNodes.end(); ++iter1) {
    Node *node = *iter1;
    node->BuildDataTypeCacheRecursive();
  }
  setIter = m_sets.begin();
  for(; setIter != m_sets.end(); ++setIter) {
#if DOBLIS
    if (!(*setIter)->IsCritSect())
#endif
      (*setIter)->BuildDataTypeCache();
  }
}

void Poss::ClearDataTypeCache()
{
  NodeVecIter iter1 = m_possNodes.begin();
  for(; iter1 != m_possNodes.end(); ++iter1)
    (*iter1)->ClearDataTypeCache();
  PSetVecIter iter2 = m_sets.begin();
  for(; iter2 != m_sets.end(); ++iter2)
    (*iter2)->ClearDataTypeCache();
}

void PrintSetOrNodeInputs(Node *node)
{
  NodeConnVecIter iter = node->m_inputs.begin();
  for(; iter != node->m_inputs.end(); ++iter) {
    Node *input = (*iter)->m_n;
    if (input->IsPossTunnel(SETTUNOUT)) {
      cout << "\tInput: Set " << ((PossTunnel*)input)->m_pset << endl;
    }
    else if (input->IsPossTunnel(POSSTUNIN)) {
      cout << "\tInput: PossTunIn " << input << " ("
      << (*iter)->m_num << ") " << input->GetNodeClass() << endl;
    }
    else if (input->IsPossTunnel())
      throw;
    else {
      cout << "\tInput: Node " << input << " ("
      << (*iter)->m_num << ") " << input->GetNodeClass() << endl;
    }
  }
}

void PrintSetOrNodeChildren(Node *node)
{
  NodeConnVecIter iter = node->m_children.begin();
  for(; iter != node->m_children.end(); ++iter) {
    Node *child = (*iter)->m_n;
    if (child->IsPossTunnel(POSSTUNOUT))
      continue;
    else if (child->IsPossTunnel(SETTUNIN)) {
      cout << "\tChild: Set " << ((PossTunnel*)child)->m_pset
      << " (" << (*iter)->m_num << ")" << endl;
      cout << "\t\t" << node->GetNameStr((*iter)->m_num) << endl;
    }
    else if (child->IsPossTunnel())
      throw;
    else {
      cout << "\tChild: Node " << child << " ("
      << (*iter)->m_num << ") " << child->GetNodeClass() << endl;
      cout << "\t\t" << node->GetNameStr((*iter)->m_num) << endl;
    }
  }
}

void Poss::PrintSetConnections()
{
  PSetSet psetSet;
  NodeSet nodeSet;
  bool printed = true;
  while (printed) {
    printed = false;
    for(unsigned int i = 0; i < m_possNodes.size(); ++i) {
      Node *node = m_possNodes[i];
      if (nodeSet.find(node) == nodeSet.end()) {
        if (node->IsPossTunnel()) {
          PossTunnel *tun = (PossTunnel*)node;
          switch(tun->m_tunType)
          {
            case (POSSTUNIN):
              cout << "POSSTUNIN " << node << " " << node->GetNodeClass() << endl;
              PrintSetOrNodeChildren(tun);
              nodeSet.insert(node);
              break;
            case (POSSTUNOUT):
            case (SETTUNOUT):
              nodeSet.insert(node);
              break;
            case (SETTUNIN):
            {
              PSet *set = tun->m_pset;
              if(psetSet.find(set) != psetSet.end())
                throw;
              psetSet.insert(set);
              cout << "Set " << set << "\t" << set->GetFunctionalityString() << endl;
              if (set->IsLoop())
                cout << "\tIs Loop\n";
              else
                cout << "\tIs not Loop\n";
              NodeVecIter iter = set->m_inTuns.begin();
              for(; iter != set->m_inTuns.end(); ++iter) {
                nodeSet.insert(*iter);
                PrintSetOrNodeInputs(*iter);
              }
              iter = set->m_outTuns.begin();
              for(; iter != set->m_outTuns.end(); ++iter) {
                nodeSet.insert(*iter);
                PrintSetOrNodeChildren(*iter);
              }
            }
              break;
            default:
              throw;
          }
          
        }
        else {
          cout << "Node " << node << " " << node->GetNodeClass() << endl;
          PrintSetOrNodeInputs(node);
          PrintSetOrNodeChildren(node);
          nodeSet.insert(node);
          
        }
      }
    }
  }
  if (psetSet.size() != m_sets.size())
    throw;
  if (nodeSet.size() != m_possNodes.size())
    throw;
}

void Poss::AddCurrPossVars(VarSet &set) const
{
  PSetVecConstIter iter = m_sets.begin();
  for(; iter != m_sets.end(); ++iter) {
    (*iter)->AddCurrPossVars(set);
  }

  NodeVecConstIter iter2 = m_possNodes.begin();
  for(; iter2 != m_possNodes.end(); ++iter2) {
    (*iter2)->AddVariables(set);
  }
}
