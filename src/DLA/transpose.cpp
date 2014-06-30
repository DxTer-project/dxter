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
#include "var.h"

#if DOBLIS
Transpose::Transpose(Trans trans, bool objectTrans)
: m_trans(trans), m_objTrans(objectTrans)
{
}
#elif DOLLDLA
Transpose::Transpose(Trans trans)
: m_trans(trans)
{
}
#endif

void Transpose::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<1,1>::Duplicate(orig, shallow, possMerging);
  const Transpose *trans = (Transpose*)orig;
  m_trans = trans->m_trans;
#if DOBLIS
  m_objTrans = trans->m_objTrans;
#endif
}

void Transpose::Flatten(ofstream &out) const
{
  DLAOp<1,1>::Flatten(out);
  WRITE(m_trans);
#if DOBLIS
  WRITE(m_objTrans);
#endif
}

void Transpose::Unflatten(ifstream &in, SaveInfo &info)
{
  DLAOp<1,1>::Unflatten(in, info);
  READ(m_trans);
#if DOBLIS
  READ(m_objTrans);
#endif
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

#if DODM
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
#endif

void Transpose::PrintCode(IndStream &out)
{
  string inputName = GetInputName(0).str();

  out.Indent();

#if DOBLIS
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
#elif DOLLDLA
  string inName = GetInputNameStr(0);
  *out << LLDLATransVarName(inName,m_trans) << " = " << inName << ";\n";
#endif
}

bool Transpose::Overwrites(const Node *input, unsigned int num) const
{
#if DOBLIS
  if (m_objTrans)
    return false;
  else 
    {
    const NodeConn *conn = m_inputs[0];
    return conn->m_n == input && conn->m_num == num;
    }
#elif DOLLDLA
  return false;
#endif
}

bool Transpose::KeepsInputVarLive(Node *input, unsigned int numInArg, unsigned int &numOutArg) const
{
#if DOBLIS
  if (Overwrites(input, numInArg)) {
    numOutArg = 0;
    return true;
  }
  else
    return false;
#elif DOLLDLA
  return false;
#endif
}

Name Transpose::GetName(unsigned int num) const
{
  if (num > 0)
    throw;
#if DOBLIS
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
#elif DOLLDLA
  Name in = GetInputName(0);
  in.m_name = LLDLATransVarName(in.m_name, m_trans);
  return in;
#endif
}

#if DOLLDLA
void Transpose::AddVariables(VarSet &set) const
{
  DLAOp<1,1>::AddVariables(set);
  Var var(GetInputNameStr(0), m_trans);
  set.insert(var);
}
#endif

void Transpose::Prop()
{
  if (!IsValidCost(m_cost)) {
    if (m_trans == NORMAL)
      throw;
    DLAOp<1,1>::Prop();
    m_cost = ZERO;
  }
}

#if DOLLDLA
const DataTypeInfo& Transpose::DataType(unsigned int num) const
{
  return m_info;
}

void Transpose::ClearDataTypeCache()
{
  m_info.m_rowStride = BADSTRIDE;
}

void Transpose::BuildDataTypeCache()
{
  const DataTypeInfo &in = InputDataType(0);
  if (m_trans == CONJ) {
    m_info = in;
  }
  else {
    m_info.m_rowStride = in.m_colStride;
    m_info.m_colStride = in.m_rowStride;
    m_info.m_numRowsVar = in.m_numColsVar;
    m_info.m_numColsVar = in.m_numRowsVar;
    m_info.m_rowStrideVar = in.m_colStrideVar;
    m_info.m_colStrideVar = in.m_rowStrideVar;
  }
}
#endif


Transpose* AddTranspose(Trans transVal,
#if DOBLIS
			bool objTrans,
#endif
			Node *input, unsigned int num, bool addToPoss)
{
#if DOBLIS
  Transpose *trans = new Transpose(transVal, objTrans);
#elif DOLLDLA
  Transpose *trans = new Transpose(transVal);
#endif
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



Transpose* InsertTranspose(Trans trans, 
#if DOBLIS
			   bool objTrans,
#endif
                           Node *node, unsigned int inNum, bool addToPoss)
{
  Node *input = node->Input(inNum);
  unsigned int num = node->InputConnNum(inNum);
  Transpose *newTrans = AddTranspose(trans, 
#if DOBLIS
				     objTrans, 
#endif
				     input, num, addToPoss);
  node->ChangeInput2Way(input, num, newTrans, 0);
  return newTrans;
}



#endif
