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

#include "layers.h"
#if DOBLIS||DOLLDLA

#include "transpose.h"
#include "loopSupport.h"
#include "helperNodes.h"


Transpose::Transpose(Trans trans, bool objectTrans)
: m_trans(trans), m_objTrans(objectTrans)
{
}

void Transpose::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<1,1>::Duplicate(orig, shallow, possMerging);
  const Transpose *trans = (Transpose*)orig;
  m_trans = trans->m_trans;
  m_objTrans = trans->m_objTrans;
}

void Transpose::Flatten(ofstream &out) const
{
  DLAOp<1,1>::Flatten(out);
  WRITE(m_trans);
  WRITE(m_objTrans);
}

void Transpose::Unflatten(ifstream &in, SaveInfo &info)
{
  DLAOp<1,1>::Unflatten(in, info);
  READ(m_trans);
  READ(m_objTrans);
}


const Sizes* Transpose::GetM(unsigned int num) const
{
  if (m_trans == NORMAL || m_trans == CONJ)
    return GetInputM(0);
  else
    return GetInputN(0);
}

const Sizes* Transpose::GetN(unsigned int num) const
{
  if (m_trans == NORMAL || m_trans == CONJ)
    return GetInputN(0);
  else
    return GetInputM(0);
}

const Sizes* Transpose::LocalM(unsigned int num) const
{
  if (m_trans == NORMAL || m_trans == CONJ)
    return InputLocalM(0);
  else
    return InputLocalN(0);
}

const Sizes* Transpose::LocalN(unsigned int num) const
{
  if (m_trans == NORMAL || m_trans == CONJ)
    return InputLocalN(0);
  else
    return InputLocalM(0);
}

void Transpose::PrintCode(IndStream &out)
{
  string inputName = GetInputName(0).str();
  out.Indent();
  if (!m_objTrans) {
#if DOSMPPHASE
    bool barrierBack = false;
    Node *in = this->Input(0);
    unsigned int num = this->InputConnNum(0);
    while (true) {
      if (in->GetNodeClass() == InputNode::GetClass()) {
        barrierBack = true;
        break;
      }
      else if (in->IsLoopTunnel())
        break;
      else if (in->IsPossTunnel(POSSTUNIN)) {
        num = in->InputConnNum(0);
        in = in->Input(0);
      }
      else if (in->IsPossTunnel(SETTUNIN)) {
        num = in->InputConnNum(0);
        in = in->Input(0);
      }
      else {
        break;
      }
    }
    if (!barrierBack) {
      if (Child(0)->GetNodeClass() == OutputNode::GetClass()) {
        *out << "th_barrier( GlobalComm );\n";
        out.Indent();
        *out << "if (th_am_root(GlobalComm))\n";
        out.Indent(1);
      }
    }
    if (barrierBack) {
      *out << "if (th_am_root(GlobalComm))\n";
      out.Indent(1);
    }
    
#endif //DOSMPPHASE
    if (m_trans == TRANS) {
      *out << "bli_obj_induce_trans( ";
    }
    else if (m_trans == CONJTRANS)
      *out << "bli_obj_induce_conjtrans( ";
    else if (m_trans == CONJ)
      *out << "bli_obj_induce_conj( ";
    else
      throw;
    *out << inputName << " );\n";
#if DOSMPPHASE
    if (barrierBack) {
      out.Indent();
      *out << "th_barrier( GlobalComm );\n";
    }
#endif //DOSMPPHASE
  }
  else {
    string name = GetNameStr(0);
    *out << "obj_t " << name << ";\n";
    out.Indent();
    if (m_trans == TRANS || m_trans == CONJTRANS)
      *out << "bli_obj_alias_with_trans( ";
    else
      *out << "bli_obj_alias_with_conj( ";
    if (m_trans == TRANS)
      *out << "BLIS_TRANSPOSE";
    else if (m_trans == CONJTRANS)
      *out << "BLIS_CONJ_TRANSPOSE";
    else if (m_trans == CONJ)
      *out << "BLIS_CONJUGATE";
    else
      throw;
    *out << ", " << inputName << ", "
    << name << ");\n";
  }
}

bool Transpose::Overwrites(const Node *input, unsigned int num) const
{
  if (m_objTrans)
    return false;
  else {
    const NodeConn *conn = m_inputs[0];
    return conn->m_n == input && conn->m_num == num;
  }
}

Name Transpose::GetName(unsigned int num) const
{
  if (num > 0)
    throw;
  if (m_objTrans) {
    Name name = GetInputName(0);
    if (m_trans == TRANS)
      name.m_name += "T";
    else if (m_trans == CONJ)
      name.m_name += "C";
    if (m_trans == CONJTRANS)
      name.m_name += "H";
    return name;
  }
  else
    return GetInputName(0);
}

void Transpose::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<1,1>::Prop();
    m_cost = ZERO;
  }
}

Transpose* AddTranspose(Trans transVal, bool objTrans, Node *input, unsigned int num, bool addToPoss)
{
  Transpose *trans = new Transpose(transVal, objTrans);
  trans->AddInput(input, num);
  if (addToPoss)
    input->m_poss->AddNode(trans);
  return trans;
}

bool CombineTranspose::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Transpose::GetClass())
    return false;
  const Transpose *trans = (Transpose*)node;
  NodeConnVecConstIter iter = trans->m_inputs.begin();
  for(; iter != trans->m_inputs.end(); ++iter) {
    const Node *input = (*iter)->m_n;
    unsigned int inputNum = (*iter)->m_num;
    NodeConnVecConstIter iter2 = input->m_children.end();
    for(; iter2 != input->m_children.end(); ++iter2) {
      if ((*iter2)->m_n == node)
        continue;
      else if ((*iter2)->m_num == inputNum) {
        if ((*iter2)->m_n->GetNodeClass() == Transpose::GetClass()) {
          if (((Transpose*)((*iter2)->m_n))->m_trans == trans->m_trans)
            return true;
        }
      }
    }
  }
  return false;
}

void CombineTranspose::Apply(Node *node) const
{
  Transpose *trans = (Transpose*)node;
  NodeConnVecConstIter iter = trans->m_inputs.begin();
  for(; iter != trans->m_inputs.end(); ++iter) {
    Node *input = (*iter)->m_n;
    unsigned int inputNum = (*iter)->m_num;
    NodeConnVecConstIter iter2 = input->m_children.end();
    for(; iter2 != input->m_children.end(); ++iter2) {
      if ((*iter2)->m_n == node)
        continue;
      else if ((*iter2)->m_num == inputNum) {
        if ((*iter2)->m_n->GetNodeClass() == Transpose::GetClass()) {
          if (((Transpose*)((*iter2)->m_n))->m_trans == trans->m_trans) {
            (*iter2)->m_n->RedirectChildren(trans, 0);
            (*iter2)->m_n->m_poss->DeleteChildAndCleanUp((*iter2)->m_n);
            return;
          }
        }
      }
    }
  }
  throw;
}


Transpose* InsertTranspose(Trans trans, bool objTrans,
                           Node *node, unsigned int inNum, bool addToPoss)
{
  Node *input = node->Input(inNum);
  unsigned int num = node->InputConnNum(inNum);
  Transpose *newTrans = AddTranspose(trans, objTrans, input, num, addToPoss);
  node->ChangeInput2Way(input, num, newTrans, 0);
  return newTrans;
}

#endif
