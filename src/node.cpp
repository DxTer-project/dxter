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
#include "elemRedist.h"
#include "loop.h"
#include "gemm.h"
#include "trxm.h"
#include "chol.h"
#include "node.h"
#include "blas.h"
#include "split.h"
#include "twoSidedTrxm.h"
#include "tensorSumScatter.h"

//#define PRINTCOSTS

static unsigned int currNum = 0;

Node::Node()
  :m_flags(0), m_poss(NULL) 
{
  m_num = currNum;
  ++currNum;
#ifdef TRACKORIG
  m_orig = NULL;
#endif
}

#include "helperNodes.h"

Node::~Node()
{
  NodeConnVecIter iter = m_inputs.begin();
  for(; iter != m_inputs.end(); ++iter)
    delete *iter;
  m_inputs.clear();
  
  iter = m_children.begin();
  for(; iter != m_children.end(); ++iter)
    delete *iter;
  m_children.clear();
  //  memset(this,0,sizeof(Node));
}

void Node::Cull(Phase phase)
{
  if (MaxPhase() < phase) {
    if (!HasRefined()) {
#if DOELEM
      if (GetNodeClass() == Gemm::GetClass()) {
        if (((Gemm*)this)->GetLayer() == SMLAYER)
          return;
      }
      if (GetNodeClass() == Axpy::GetClass()) {
        if (((Axpy*)this)->GetLayer() == DMLAYER)
          return;
      }
      if (GetNodeClass() == Trxm::GetClass()) {
        if (((Trxm*)this)->m_invert && (((Trxm*)this)->GetLayer() == SMLAYER))
          return;
      }
      if (GetNodeClass() == Chol::GetClass()) {
        if (((Chol*)this)->GetLayer() == SMLAYER)
          return;
      }
#endif
      cout << "skipping culling for invalid node type " << GetType() << endl;
      cout << MaxPhase() << " is MaxPhase\n";
#if DOTENSORS
      if (GetNodeClass() == SumScatterUpdateNode::GetClass()) {
	SumScatterUpdateNode *sum = (SumScatterUpdateNode*)this;
	cout << "sumscatter " << sum->InputDataType(0).m_dist.str() 
	     << " -> " << sum->DataType(0).m_dist.str() << endl;
      }

#endif
      for (unsigned int i = 0; i < m_inputs.size(); ++i) {
        DLANode *in = (DLANode*)Input(i);
        cout << "Input " << i 
#if DOELEM
	     << " " << DistTypeToStr(((DLANode*)this)->InputDataType(i).m_dist) << endl;
#else
	<<endl;
#endif
        ClassType type = in->GetNodeClass();
#if DOELEM
        if (type == RedistNode::GetClass()) {
          cout << DistTypeToStr(in->InputDataType(0).m_dist) << endl;
        }
#endif
      }
      m_poss->PrintTransVec();
      throw;
    }
    else
      m_poss->MarkInsane(true);
  }
}

void Node::AddChild(Node *node, unsigned int num)
{
  NodeConn *conn = new NodeConn(node,num);
  m_children.push_back(conn);
}

void Node::RemoveChild(Node *node, unsigned int num)
{
  NodeConnVecIter iter = m_children.begin();
  for( ; iter != m_children.end(); ++iter) {
    if ((*iter)->m_n == node && (*iter)->m_num == num) {
      delete *iter;
      m_children.erase(iter);
      return;
    }
  }
  printf("in RemoveChild child %s %p not found on %s %p\n", node->GetType().c_str(), node, GetType().c_str(), this);
  cout << "only has " << m_children.size() << " children\n";
  cout << "child on " << node->m_poss << " parent on " << node->m_poss << endl;
  iter = m_children.begin();
  for (; iter != m_children.end(); ++iter) {
    cout << "child: " << (*iter)->m_n->GetNodeClass() << endl;
  }
  fflush(stdout);
  throw;
}

void Node::RemoveInput(Node *node, unsigned int num)
{
  /*****
   Notice that this invalidates/makes insane the node
   Only to be used while picking apart a poss!
   *****/
  NodeConnVecIter iter = m_inputs.begin();
  for( ; iter != m_inputs.end(); ++iter) {
    if ((*iter)->m_n == node && (*iter)->m_num == num) {
      delete *iter;
      m_inputs.erase(iter);
      return;
    }
  }
  printf("in RemoveInput input %s %p not found on %s %p\n", node->GetType().c_str(), node, GetType().c_str(), this);
  cout << "Has " << m_inputs.size() << " inputs\n";
  cout << "Input(0) = " << Input(0) << " " << Input(0)->GetType() << endl;
  
  throw;
}

bool Node::operator==(const Node &rhs) const
{
  if (GetType() != rhs.GetType())
    return false;
  if (m_inputs.size() != rhs.m_inputs.size())
    return false;
  //  if (m_children.size() != rhs.m_children.size())
  //    return false;
  NodeConnVecConstIter iter1, iter2;
  iter1 = m_inputs.begin();
  iter2 = rhs.m_inputs.begin();
  for( ; iter1 != m_inputs.end(); ++iter1, ++iter2) {
    Node *node1 = (*iter1)->m_n;
    Node *node2 = (*iter2)->m_n;
    if ((*iter1)->m_num != (*iter2)->m_num)
      return false;
    if (node1->IsPossTunnel()) {
      if (!node2->IsPossTunnel())
        return false;
      if (((PossTunnel*)node1)->m_tunType != ((PossTunnel*)node2)->m_tunType)
        return false;
      if (node1->IsPossTunnel(SETTUNOUT)) {
        PossTunnel *tun1 = (PossTunnel*)node1;
        PossTunnel *tun2 = (PossTunnel*)node2;
        PSet *set1 = tun1->m_pset;
        PSet *set2 = tun2->m_pset;
        if (!(*set1 == *set2))
          return false;
        if (set1->m_inTuns.size() != set2->m_inTuns.size()) {
          return false;
        }
        for(unsigned int i = 0; i < set1->m_inTuns.size(); ++i) {
          if (! (*(set1->m_inTuns[i]) == *(set2->m_inTuns[i])))
            return false;
        }
      }
    }
    else if (node2->IsPossTunnel())
      return false;
    else if (!(**iter1 == **iter2))
      return false;
  }
  return true;
}

void Node::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  //  cout << "duplicated " << orig << endl;
  
#ifdef TRACKORIG
  if (orig->m_orig)
    m_orig = orig->m_orig;
  else
    m_orig = orig;
#endif
  
  if (!shallow)
    m_poss = orig->m_poss;
  if (possMerging)
    m_applications = orig->m_applications;
  m_inverseOps = orig->m_inverseOps;

  m_flags = orig->m_flags;
  
  //Don't duplicate this multiple times through double inheritance
  if (!shallow && m_inputs.empty()) {
    //need to go through children, create new nodeconn, and copy
    NodeConnVecConstIter iter = orig->m_inputs.begin();
    for(; iter != orig->m_inputs.end(); ++iter)
      m_inputs.push_back(new NodeConn(*iter));
  }
  if (!shallow && m_children.empty()) {
    NodeConnVecConstIter iter = orig->m_children.begin();
    for(; iter != orig->m_children.end(); ++iter)
      m_children.push_back(new NodeConn(*iter));
  }
}

void Node::PatchAfterDuplicate(NodeMap &map)
{
  for(unsigned int i = 0; i < m_inputs.size(); ++i) {
    if (!Input(i)) {
      cout<<"!m_inputs[i]->m_n\n";
      throw;
    }
    if (!map[Input(i)]) {
      cout << "map[m_inputs[i]->m_n] size = " << map.size() << endl;
      cout << "Input's poss = " << Input(i)->m_poss << endl;
      cout << "my poss = " << m_poss << endl;
      printf("didn't find input %s, %p for %s, %p\n", Input(i)->GetType().c_str(), Input(i), GetType().c_str(), this);
      throw;
    }
    m_inputs[i]->SetNode(map[Input(i)]);
  }
  for(unsigned int i = 0; i < m_children.size(); ++i) {
    Node *child = m_children[i]->m_n;
    Node *mapped = map[child];
    if (!mapped) {
      cout << "this " << this << " " << this->GetNodeClass() << endl;
      cout << "child " << child << " " << child->GetNodeClass() << endl;
      cout << "map[m_children[i]->m_n] size = " << map.size() << endl;
      printf("didn't find child %p %s, %p %s\n", m_children[i]->m_n, m_children[i]->m_n->GetType().c_str(), this, GetType().c_str());
      cout << "m_poss = " << m_poss << endl;
      cout << "child's poss = " << m_children[i]->m_n->m_poss << endl;
      throw;
    }
    m_children[i]->SetNode(mapped);
  }
}

void Node::Print(IndStream &out, unsigned int graphNum)
{
  if (!HasPrinted()) {
    if(CanPrintCode()) {
      SetPrinted();
      PrintCode(out);
#ifdef PRINTCOSTS
      if (((DLANode*)this)->m_cost) {
	(*out).precision(5);
	*out << scientific << "cost of previous node: " << ((DLANode*)this)->m_cost << endl;
      }
#endif
      
      if (!IsPossTunnel()) {
        //If there are redist nodes as my child, go ahead
        // and print them to make less likely a weird situation
        // that arrises from DAGs turned into code
        // where you can have something like
        // opA(A10_D_VR_STAR)
        // opB(A10_D_STAR_VR)
        // A10 = A10_D_VR_STAR
        // A10 = A10_D_STAR_VR
        // A10_D_MC_STAR = A10
        // opC(A10_D_MC_STAR)
        // the 1st, 3rd, 5th, and 6th lines are a stream
        // of ops and the 2nd and 4th are, where the
        // 4th should be the final update of A10
        //I think this will work often, but not always
        NodeConnVecConstIter childIter;
#if DOELEM
	childIter = m_children.begin();
        for( ; childIter != m_children.end(); ++childIter) {
          Node *child = (*childIter)->m_n;
          if (child->GetNodeClass() == RedistNode::GetClass())
            child->Print(out, graphNum);
        }
#endif
        //Try to print the read-only users of this first
        childIter = m_children.begin();
        for( ; childIter != m_children.end(); ++childIter) {
          Node *child = (*childIter)->m_n;
          if (!child->IsPossTunnel() &&
              (child->Input(child->m_inputs.size()-1) != this))
          {
            child->Print(out, graphNum);
          }
        }
      }
    }
  }
}

bool Node::CanPrintCode() const
{
  NodeConnVecConstIter iter = m_inputs.begin();
  for( ; iter != m_inputs.end(); ++iter) {
    const NodeConn *conn = *iter;
    const Node *input = conn->m_n;
    unsigned int num = conn->m_num;
    if (!input->HasPrinted()) {
      return false;
    }
    //If input is SETTUNIN, then this is a POSTUNIN
    // there will be other children of input that 
    // also overwrite just like this.
    else if (!input->IsPossTunnel(SETTUNIN) && Overwrites(input, num))
    {
      NodeConnVecConstIter childIter = input->m_children.begin();
      for(; childIter != input->m_children.end(); ++childIter) {
        const NodeConn *childConn = *childIter;
        const Node *child = childConn->m_n;
        unsigned int childNum = childConn->m_num;
        if (childNum == num && child != this && !child->IsPossTunnel(POSSTUNOUT)) {
          if (child->Overwrites(input, num)) {
            cout << "two children of " << input->GetNodeClass() << " overwrite output " << childNum << "\n";
            cout << GetNodeClass() << endl;
            cout << child->GetNodeClass() << endl;
            throw;
          }
          else {
            if (!child->HasPrinted()) {
              if (IsPossTunnel(SETTUNIN) && child->IsPossTunnel(SETTUNIN)) {
                const PossTunnel *tun1 = (PossTunnel*)this;
                const PossTunnel *tun2 = (PossTunnel*)child;
                if (tun1->m_pset != tun2->m_pset)
                  return false;
              }
              else
                return false;
            }
            else if (child->IsPossTunnel(SETTUNIN)) {
              const PossTunnel *tun2 = (PossTunnel*)child;
              const PSet *set = tun2->m_pset;
              if (this->IsPossTunnel(SETTUNIN)) {
                if (set == ((PossTunnel*)this)->m_pset)
                  continue;
              }
              if (!set->GetCurrPoss()->m_hasPrinted)
                return false;
            }
          }
        }
      }
    }
  }
  
  return true;
}

void Node::AddInput(Node *node)
{
  if (!node) {
    cout << "!node\n";
    cout.flush();
  }
  if (node->NumOutputs() == 1)
    AddInput(node, 0);
  else {
    cout << "bad call to AddInput on " << node->GetType() << "\n";
    throw;
  }
}

void Node::AddInput(Node *node, unsigned int num)
{
  //  cout << "adding input " << node << " to " << this << endl;
  NodeConn *conn = new NodeConn(node, num);
  m_inputs.push_back(conn);
  node->AddChild(this, num);
}

void Node::AddInputs(int numArgs, ...)
{
  va_list listPointer;
  va_start (listPointer, numArgs);
  for(int i = 0; i < numArgs; i+=2) {
    Node *node = va_arg(listPointer, Node* );
    int num = va_arg(listPointer, int);
    AddInput(node, num);
  }
}

void Node::ChangeInput1Way(Node *oldInput, unsigned int oldNum, Node *newInput, unsigned int newNum)
{
  //Reflect in 2Way
  if (oldInput == newInput && oldNum == newNum)
    throw;
  unsigned int i = 0;
  NodeConnVecIter iter = m_inputs.begin();
  for ( ; iter != m_inputs.end(); ++iter) {
    if ((*iter)->m_n == oldInput && (*iter)->m_num == oldNum) {
      m_inputs[i]->SetNode(newInput);
      m_inputs[i]->m_num = newNum;
      newInput->AddChild(this, newNum);
      return;
    }
    ++i;
  }
  cout << "Didn't find oldInput " << oldInput->GetType() << " " << oldInput << " to " <<
  GetType() << " " << this << "!\n";
  for (unsigned int i = 0; i < m_inputs.size(); ++i)
    cout << "input " << i << " is " << m_inputs[i]->m_n->GetType() << " " <<  m_inputs[i] << " on " << m_inputs[i]->m_n->m_poss << endl;
  throw;
}

void Node::ChangeInput2Way(Node *oldInput, unsigned int oldNum, Node *newInput, unsigned int newNum)
{
  //Reflect in 1Way
  if (oldInput == newInput && oldNum == newNum)
    throw;
  unsigned int i = 0;
  NodeConnVecIter iter = m_inputs.begin();
  for ( ; iter != m_inputs.end(); ++iter) {
    if ((*iter)->m_n == oldInput && (*iter)->m_num == oldNum) {
      NodeConn *conn = m_inputs[i];
      conn->m_n->RemoveChild(this, conn->m_num);
      conn->SetNode(newInput);
      conn->m_num = newNum;
      newInput->AddChild(this, newNum);
      return;
    }
    ++i;
  }
  cout << "Didn't find oldInput " << oldInput->GetType() << " " << oldInput << " to " <<
  GetType() << " " << this << "!\n";
  for (unsigned int i = 0; i < m_inputs.size(); ++i)
    cout << "input " << i << " is " << m_inputs[i]->m_n->GetType() << " " <<  m_inputs[i] << " on " << m_inputs[i]->m_n->m_poss << endl;
  throw;
}

void Node::RedirectChildren(Node *newInput)
{
  if (newInput->NumOutputs() != 1) {
    cout << "Bad call to RedirectChildren 1\n";
    throw;
  }
  else
    RedirectChildren(newInput, 0);
}

void Node::RedirectChildren(Node *newInput, unsigned int newNum)
{
  if (NumOutputs() != 1) {
    cout << "Bad call to RedirectChildren 2\n";
    throw;
  }
  else
    RedirectChildren(0, newInput, newNum);
}

void Node::RedirectChildren(unsigned int oldNum, Node *newInput, unsigned int newNum)
{
  bool found = false;
  if (!m_children.size()) {
    cout << "no children for " << this << endl;
    throw;
  }
  for(unsigned int i = 0; i < m_children.size(); ++i) {
    NodeConn *output = m_children[i];
    if (output->m_num == oldNum) {
      //      cout << "on " << GetType() << " redirecting output->m_n " << output->m_n->GetType() << endl;
      RedirectChild(i, newInput, newNum);
      --i;
      found = true;
    }
  }
  if (!found) {
    cout << "!found\n";
    throw;
  }
}

void Node::RedirectAllChildren(Node *newInput)
{
  while (m_children.size()) {
    RedirectChild(0, newInput, m_children[0]->m_num);
  }
}

void Node::RedirectChild(unsigned int childNum, Node *newInput, unsigned int newNum)
{
  NodeConn *conn = m_children[childNum];
  conn->m_n->ChangeInput1Way(this, conn->m_num, newInput, newNum);
  delete conn;
  m_children.erase(m_children.begin() + childNum);
}

bool Node::HasApplied(const Transformation *trans) const
{
  return (m_applications.find(trans) != m_applications.end())
  || (m_inverseOps.find(trans) != m_inverseOps.end());
}

bool Node::Applied(const Transformation *trans)
{
  m_applications.insert(trans);
  return true;
}

void Node::CheckConnections()
{
  NodeConnVecIter iter1 = m_inputs.begin();
  for( ; iter1 != m_inputs.end(); ++iter1 ) {
    //check that all inputs are either in the same poss or i'm a posstunnel
    if (!IsPossTunnel(SETTUNOUT) && !IsPossTunnel(POSSTUNIN)) {
      if ((*iter1)->m_n->m_poss != m_poss) {
        cout << "Input to " << GetType() << " isn't in the same poss (" << (*iter1)->m_n->GetType() << ")\n";
        throw;
      }
    }
    
    if ((*iter1)->m_n->IsPossTunnel(SETTUNOUT)) {
      PossTunnel *tun = (PossTunnel*)((*iter1)->m_n);
      PSetVecIter iter = m_poss->m_sets.begin();
      bool found = false;
      for(; iter != m_poss->m_sets.end() && !found; ++iter) {
        if (*iter == tun->m_pset)
          found = true;
      }
      if (!found)
        throw;
    }
    
    
    if ((*iter1)->m_n->m_poss) {
      NodeVecConstIter nodeIter = (*iter1)->m_n->m_poss->m_possNodes.begin();
      bool found = false;
      for( ; !found && nodeIter != (*iter1)->m_n->m_poss->m_possNodes.end(); ++nodeIter) {
        if (*nodeIter == (*iter1)->m_n)
          found = true;
      }
      if (!found) {
        cout << "Input " << (*iter1)->m_n->GetType() << " not in Alg's list\n";
        throw;
      }
    }
    
    if (!(*iter1)->m_n->InChildren(this, (*iter1)->m_num)) {
      cout << "Input doesn't know about child\n";
      cout << "Input is " << (*iter1)->m_n << " " << (*iter1)->m_num << endl;
      cout << "Child is " << GetType() << " " << this << endl;
      throw;
    }
  }
  
  NodeConnVecConstIter iter2 = m_children.begin();
  for( ; iter2 != m_children.end(); ++iter2 ) {
    if ((*iter2)->m_num >= NumOutputs()) {
      cout << "m_children->m_num >= NumOutputs()\n";
      cout << "m_children->m_num = " << (*iter2)->m_num << "   NumOutputs() = " << NumOutputs() << endl;
      throw;
    }
    if ((*iter2)->m_n->m_poss) {
      NodeVecConstIter nodeIter = (*iter2)->m_n->m_poss->m_possNodes.begin();
      bool found = false;
      for( ; !found && nodeIter != (*iter2)->m_n->m_poss->m_possNodes.end(); ++nodeIter) {
        if (*nodeIter == (*iter2)->m_n)
          found = true;
      }
      if (!found) {
        cout << "Output not in Alg's list\n";
        throw;
      }
    }
  }
  
  {
    NodeConnVecConstIter iter = m_children.begin();
    for(; iter != m_children.end(); ++iter) {
      Node *child = (*iter)->m_n;
      if (!child->InInputs(this,(*iter)->m_num)) {
        cout << "Didn't find " << this << " " << GetType() << " as input to child " << child << endl;
        cout << "Didn't find " << this << " " << GetType() << " as input to child " << child << " " << child->GetType() << " input list\n";
        throw;
      }
    }
  }
  
  if (!IsPossTunnel()
      || ((PossTunnel*)this)->m_tunType == POSSTUNIN
      || ((PossTunnel*)this)->m_tunType == POSSTUNOUT)
  {
    bool found = false;
    NodeVecConstIter iter = m_poss->m_possNodes.begin();
    for (; !found && iter != m_poss->m_possNodes.end(); ++iter)
      if (*iter == this)
        found = true;
    if (!found) {
      cout << "Poss doesn't know about " << GetType() << endl;
      throw;
    }
  }
}

NodeType Node::GetType() const
{
  cout << GetNodeClass() << endl;
  throw;
}

Name Node::GetInputName(unsigned int num) const
{
  if (num >= m_inputs.size()) {
    cout << "num too big " << GetType() << "\n";
  }
  return Input(num)->GetName(InputConnNum(num));
}

void Node::AddToPoss(Poss *poss)
{
  m_poss = poss;
  poss->m_possNodes.push_back(this);
}

Node* Node::Input(unsigned int num) const
{
  if (num >= m_inputs.size()) {
    cout << "num bad size on " << GetType() << "\n";
    cout << num << " >= " << m_inputs.size() << endl;
    throw;
  }
  if (!(m_inputs[num]->m_n)) {
    cout << "(!(m_inputs[num]->m_n))\n";
    throw;
  }
  NodeConn *conn = m_inputs[num];
  return conn->m_n;
}

NodeConn* Node::InputConn(unsigned int num) const
{
  if (num >= m_inputs.size()) {
    cout << "num bad size on " << GetType() << "\n";
    cout << "num " << num << " >= " << m_inputs.size() << endl;
    throw;
  }
  if (!(m_inputs[num]->m_n)) {
    cout << "(!(m_inputs[num]->m_n))\n";
    throw;
  }
  return m_inputs[num];
}

unsigned int Node::InputConnNum(unsigned int num) const
{
  if (m_inputs.size() <= num)
    throw;
  return m_inputs[num]->m_num;
}

Node* Node::Child(unsigned int num) const
{
  return m_children[num]->m_n;
}

unsigned int Node::ChildConnNum(unsigned int num) const
{
  return m_children[num]->m_num;
}

unsigned int Node::NumChildrenOfOutput(unsigned int num) const
{
  unsigned int count = 0;
  NodeConnVecConstIter iter = m_children.begin();
  for(; iter != m_children.end(); ++iter)
    if ((*iter)->m_num == num)
      ++count;
  return count;
}

bool Node::InChildren(Node *node, unsigned int num) const
{
  NodeConnVecConstIter iter = m_children.begin();
  for(; iter != m_children.end(); ++iter) {
    if ((*iter)->m_n == node && (*iter)->m_num == num)
      return true;
  }
  return false;
}

bool Node::InInputs(Node *node, unsigned int num) const
{
  NodeConnVecConstIter iter = m_inputs.begin();
  for(; iter != m_inputs.end(); ++iter)
    if ((*iter)->m_n == node && (*iter)->m_num == num)
      return true;
  return false;
}

void Node::PrintChildren()
{
  for (unsigned int i = 0; i < m_children.size(); ++i)
    cout << "Child " << i << " " << Child(i)->GetType() << " " << Child(i) << endl;
}

void Node::PrintInputs()
{
  for (unsigned int i = 0; i < m_inputs.size(); ++i) {
    const Node *node = Input(i);
    cout << "Input " << i << " " << node->GetType() << " " << node << endl;
  }
}

void Node::Flatten(ofstream &out) const
{
  unsigned int size = m_applications.size();
  WRITE(size);
  TransSetConstIter iter = m_applications.begin();
  for(; iter != m_applications.end(); ++iter)
    WRITE(*iter);
  size = m_inverseOps.size();
  WRITE(size);
  iter = m_inverseOps.begin();
  for(; iter != m_inverseOps.end(); ++iter)
    WRITE(*iter);
  size = m_inputs.size();
  WRITE(size);
  NodeConnVecConstIter iter2 = m_inputs.begin();
  for(; iter2 != m_inputs.end(); ++iter2) {
    (*iter2)->Flatten(out);
  }
  size = m_children.size();
  WRITE(size);
  iter2 = m_children.begin();
  for(; iter2 != m_children.end(); ++iter2) {
    (*iter2)->Flatten(out);
  }
  WRITE(m_poss);
  FlattenCore(out);
}


void Node::Unflatten(ifstream &in, SaveInfo &info)
{
  unsigned int size;
  READ(size);
  for(unsigned int i = 0; i < size; ++i) {
    Transformation *trans;
    READ(trans);
    Swap(&trans, info.transMap);
    m_applications.insert(trans);
  }
  READ(size);
  for(unsigned int i = 0; i < size; ++i) {
    Transformation *trans;
    READ(trans);
    Swap(&trans, info.transMap);
    m_inverseOps.insert(trans);
  }
  READ(size);
  for(unsigned int i = 0; i < size; ++i) {
    NodeConn *conn = new NodeConn;
    conn->Unflatten(in);
    m_inputs.push_back(conn);
  }
  READ(size);
  for(unsigned int i = 0; i < size; ++i) {
    NodeConn *conn = new NodeConn;
    conn->Unflatten(in);
    m_children.push_back(conn);
  }
  READ(m_poss);
  Swap(&m_poss,info.possMap);
  UnflattenCore(in,info);
}


void FullyFlatten(const NodeVec &vec, ofstream &out)
{
  unsigned int size = vec.size();
  WRITE(size);
  NodeVecConstIter iter = vec.begin();
  for(; iter != vec.end(); ++iter) {
    WRITE(*iter);
    out << (*iter)->GetNodeClass() << endl;
  }
  iter = vec.begin();
  for(; iter != vec.end(); ++iter) {
    (*iter)->Flatten(out);
    WRITE(END);
  }
}

void FullyUnflatten(NodeVec &vec, ifstream &in, SaveInfo &info)
{
  unsigned int size;
  READ(size);
  for(unsigned int i = 0; i < size; ++i) {
    Node *node;
    READ(node);
    string className;
    getline(in, className);
    Node *newNode = Universe::GetBlankClassInst(className);
    (*(info.nodeMap))[node] = newNode;
    vec.push_back(newNode);
  }
  char tmp;
  NodeVecIter iter = vec.begin();
  for(; iter != vec.end(); ++iter) {
    (*iter)->Unflatten(in, info);
    READ(tmp);
    if (tmp != END)
      throw;
  }
}


void Node::BuildDataTypeCacheRecursive()
{
  if (!(m_flags & BUILDFLAG)) {
    NodeConnVecIter iter = m_inputs.begin();
    for(; iter != m_inputs.end(); ++iter) {
      Node *node = (*iter)->m_n;
      if (node->m_poss == m_poss)
	{
	  node->BuildDataTypeCacheRecursive();
	}
    }
    m_flags |= BUILDFLAG;
    BuildDataTypeCache();
  }
}

#if DOBLIS
bool Node::InCriticalSection() const
{
  Poss *poss = m_poss;
  while (poss && poss->m_pset) {
    PSet *set = poss->m_pset;
    if (set->IsCritSect()) {
      return true;
    }
    poss = set->m_ownerPoss;
  }
  return false;
}

Comm Node::WithinParallelism() const
{
  const Poss *poss = m_poss;
  while (poss && poss->m_pset) {
    const PSet *set = poss->m_pset;
    if (set->IsCritSect()) {
      return CORECOMM;
    }
#if DOBLIS
    else if (set->IsLoop()) {
      const Loop *loop = (Loop*)set;
      if (loop->IsParallel())
	return loop->m_comm;
    }
#endif
    poss = set->m_ownerPoss;
  }
  return CORECOMM;
}
#endif //DOBLIS

void Node::AddVariables(VarSet &set) const
{
#if DOTENSORS
  for(unsigned int i = 0; i < NumOutputs(); ++i) {
    Name name = GetName(i);
    Var var(name);
    set.insert(var);
    Var var2(name.m_type);
    set.insert(var2);
  }
#endif
}

const DataTypeInfo& Node::InputDataType(unsigned int num) const
{
  const NodeConn *conn = m_inputs[num];
  return conn->m_n->DataType(conn->m_num);
}
