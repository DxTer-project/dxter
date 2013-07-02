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

#define PORTIONPARALLELIZABLE .8

#include "blis.h"
#include "loopSupport.h"

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

bool CombineTranspose::CanApply(const Poss *poss, const Node *node) const
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

void CombineTranspose::Apply(Poss *poss, Node *node) const
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

GetUpToDiag::GetUpToDiag(Tri tri, PartDir dir)
{
  m_tri = tri;
  m_sizes = NULL;
  m_lsizes = NULL;
  if (dir != PARTRIGHT && dir != PARTDOWN)
    throw;
  m_dir = dir;
}

unsigned int GetUpToDiag::NumOutputs() const
{
  if (m_inputs.size() == 3)
    return 2;
  else if (m_inputs.size() == 4)
    return 3;
  else
    throw;
}

const Sizes* GetUpToDiag::GetM(unsigned int num) const
{
  if (m_dir == PARTDOWN)
    return GetInputM(num+1);
  else
    return m_sizes;
}

const Sizes* GetUpToDiag::GetN(unsigned int num) const
{
  if (m_dir == PARTDOWN)
    return m_sizes;
  else
    return GetInputN(num+1);
}

const Sizes* GetUpToDiag::LocalM(unsigned int num) const
{
  if (m_dir == PARTDOWN)
    return InputLocalM(num+1);
  else
    return m_lsizes;
}

const Sizes* GetUpToDiag::LocalN(unsigned int num) const
{
  if (m_dir == PARTDOWN)
    return m_lsizes;
  else
    return InputLocalN(num+1);
}

void GetUpToDiag::ClearSizeCache()
{
  if (m_sizes) {
    delete m_sizes;
    m_sizes = NULL;
    delete m_lsizes;
    m_lsizes = NULL;
  }
}

void GetUpToDiag::BuildSizeCache()
{
  if (m_sizes)
    return;
  else {
    if (m_dir == PARTDOWN) {
      m_sizes = new Sizes;
      m_sizes->PairwiseSum(*GetInputM(0), *GetInputM(1));
      m_lsizes = new Sizes;
      m_lsizes->PairwiseSum(*InputLocalM(0), *InputLocalM(1));
    }
    else if (m_dir == PARTRIGHT) {
      m_sizes = new Sizes;
      m_sizes->PairwiseSum(*GetInputN(0), *GetInputN(1));
      m_lsizes = new Sizes;
      m_lsizes->PairwiseSum(*InputLocalN(0), *InputLocalN(1));
    }
    else
      throw;
  }
}

Name GetUpToDiag::GetName(unsigned int num) const
{
  Name name = GetInputName(num+1);
  if (m_dir == PARTDOWN) {
    if (m_tri == LOWER) {
      name.m_name += "_L";
    }
    else if (m_tri == UPPER) {
      name.m_name += "_R";
    }
    else
      throw;
  }
  else if (m_dir == PARTRIGHT) {
    if (m_tri == LOWER) {
      name.m_name += "_B";
    }
    else if (m_tri == UPPER) {
      name.m_name += "_T";
    }
    else
      throw;
  }
  else
    throw;
  return name;
}

void GetUpToDiag::Prop()
{
  if (!IsValidCost(m_cost)) {
    const unsigned int numIn = m_inputs.size();
    if (numIn != 3 && numIn != 4)
      throw;
    for (unsigned int i = 0; i < numIn; ++i)
      Input(i)->Prop();
    
    m_cost = ZERO;
  }
}

void GetUpToDiag::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  GetUpToDiag *diag = (GetUpToDiag*)orig;
  m_tri = diag->m_tri;
  m_dir = diag->m_dir;
}

void GetUpToDiag::PrintCode(IndStream &out)
{
  char triChar;
  bool lower = m_tri == LOWER;
  if (m_dir == PARTDOWN) {
    if (lower)
      triChar = 'L';
    else
      triChar = 'R';
  }
  else {
    if (lower)
      triChar = 'B';
    else
      triChar = 'T';
  }
  
  out.Indent();
  if (m_inputs.size() == 3)
    *out << "obj_t " << GetNameStr(0) << ", " << GetNameStr(1) << ";\n";
  else
    *out << "obj_t " << GetNameStr(0) << ", " << GetNameStr(1) << ", " << GetNameStr(2) << ";\n";
  out.Indent();
  *out << "dim_t off" << triChar << " = ";
  if (m_dir == PARTDOWN) {
    if (lower) {
      *out << "0;\n";
      out.Indent();
      *out << "dim_t n" << triChar << " = bli_min( bli_obj_width_after_trans( "
      << GetInputNameStr( 1 ) << " ), \n";
      out.Indent(2);
      *out << "bli_obj_diag_offset_after_trans( " << GetInputNameStr( 1 ) << " ) + "
      << "bs" << out.LoopLevel() << " );\n";
    }
    else {
      *out << "bli_max( 0, bli_obj_diag_offset_after_trans( " << GetInputNameStr( 1 )
      << " ) );\n";
      out.Indent();
      *out << "dim_t n" << triChar << " = bli_obj_width_after_trans( "
      << GetInputNameStr( 1 ) << " ) - off" << triChar << ";\n";
    }
  }
  else if (m_dir == PARTRIGHT) {
    if (lower) {
      *out << "bli_max( 0, -bli_obj_diag_offset_after_trans( " << GetInputNameStr( 1 )
      << " ) );\n";
      out.Indent();
      *out << "dim_t m" << triChar << " = bli_obj_length_after_trans( "
      << GetInputNameStr( 1 ) << " ) - off" << triChar << ";\n";
    }
    else {
      *out << "0;\n";
      out.Indent();
      *out << "dim_t m" << triChar << " = bli_min( bli_obj_length_after_trans( "
      << GetInputNameStr( 1 ) << " ), \n";
      out.Indent(2);
      *out << "-bli_obj_diag_offset_after_trans( " << GetInputNameStr( 1 ) << " ) + "
      << "bs" << out.LoopLevel() << " );\n";
    }
  }
  
  out.Indent();
  if (m_dir == PARTDOWN)
    *out << "bli_acquire_mpart_l2r( BLIS_SUBPART1,\n";
  else if (m_dir == PARTRIGHT)
    *out << "bli_acquire_mpart_t2b( BLIS_SUBPART1,\n";
  out.Indent(2);
  *out << "off" << triChar << ", "
  << (m_dir == PARTDOWN ? "n" : "m") << triChar << ", &"
  << GetInputNameStr(1) << ", &" << GetNameStr(0) << " );\n";
  //  out.Indent(2);
  out.Indent();
  if (m_dir == PARTDOWN)
    *out << "bli_acquire_mpart_l2r( BLIS_SUBPART1,\n";
  else if (m_dir == PARTRIGHT)
    *out << "bli_acquire_mpart_t2b( BLIS_SUBPART1,\n";
  out.Indent(2);
  *out << "off" << triChar << ", "
  << (m_dir == PARTDOWN ? "n" : "m") << triChar << ", &"
  << GetInputNameStr(2) << ", &" << GetNameStr(1) << " );\n";
  if (m_inputs.size() == 4) {
    out.Indent();
    if (m_dir == PARTDOWN)
      *out << "bli_acquire_mpart_l2r( BLIS_SUBPART1,\n";
    else if (m_dir == PARTRIGHT)
      *out << "bli_acquire_mpart_t2b( BLIS_SUBPART1,\n";
    out.Indent(2);
    *out << "off" << triChar << ", "
    << (m_dir == PARTDOWN ? "n" : "m") << triChar << ", &"
    << GetInputNameStr(3) << ", &" << GetNameStr(2) << " );\n";
  }
}

void GetUpToDiag::Flatten(ofstream &out) const
{
  DLANode::Flatten(out);
  WRITE(m_tri);
  WRITE(m_dir);
}

void GetUpToDiag::Unflatten(ifstream &in, SaveInfo &info)
{
  DLANode::Unflatten(in, info);
  READ(m_tri);
  READ(m_dir);
}

void CombineDiag::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();
    m_cost = ZERO;
  }
}

SetObjProps::SetObjProps(Tri tri, Diag diag, TriStruct triStruct)
{
  m_tri = tri;
  m_diag = diag;
  m_struct = triStruct;
}

void SetObjProps::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<1,1>::Prop();
    m_cost = ZERO;
  }
}

void SetObjProps::PrintCode(IndStream &out)
{
  if (m_struct != GEN) {
    out.Indent();
    *out << "bli_obj_set_struc( ";
    switch(m_struct)
    {
      case (HERM):
        *out << "BLIS_HERMITIAN";
        break;
      case (SYMM):
        *out << "BLIS_SYMMETRIC";
        break;
      case (TRI):
        *out << "BLIS_TRIANGULAR";
        break;
      default:
        throw;
    }
    *out << ", " << GetNameStr(0) << " );\n";
  }
  if (m_tri != NOTTRI) {
    out.Indent();
    *out << "bli_obj_set_uplo( ";
    if (m_tri == LOWER)
      *out << "BLIS_LOWER";
    else
      *out << "BLIS_UPPER";
    *out << ", " << GetNameStr(0) << " );\n";
  }
  if (m_diag == UNIT) {
    out.Indent();
    *out << "bli_obj_set_diag( BLIS_UNIT_DIAG, "
    << GetNameStr(0) << " );\n";
  }
}

void SetObjProps::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<1,1>::Duplicate(orig, shallow, possMerging);
  const SetObjProps *props = (SetObjProps*)orig;
  m_tri = props->m_tri;
  m_struct = props->m_struct;
}

void SetObjProps::Flatten(ofstream &out) const
{
  DLAOp<1,1>::Flatten(out);
  WRITE(m_tri);
  WRITE(m_struct);
}

void SetObjProps::Unflatten(ifstream &in, SaveInfo &info)
{
  DLAOp<1,1>::Unflatten(in,info);
  READ(m_tri);
  READ(m_struct);
}

void Copy::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();
    m_cost = (PSIRVAL+PSIWVAL) * LocalM(0)->SumProds11(*LocalN(0));
  }
}

void Copy::PrintCode(IndStream &out)
{
  out.Indent();
  *out << "bli_copym( &" << GetInputNameStr(0) << ", &"
  << GetInputNameStr(1) << " );\n";
}

bool ParallelizeMDim::CanApply(const Poss *poss, const Node *node) const
{
  const PackBuff *buff = (PackBuff*)node;
  if (buff->m_packMat == PACKBPANEL) {
    if (buff->m_children.size() != 1) {
      cout << "more than one child of B panel packbuff!\n";
      throw;
    }
    if (buff->Child(0)->GetNodeClass() != Pack::GetClass()) {
      cout << "child of B panel packbuff isn't Pack\n";
      throw;
    }
    const Pack *pack = (Pack*)buff->Child(0);
    if (!LegalParallelizationNestingUp(pack, m_comm))
      return false;
    if (pack->m_children.size() < 1)
      throw;
    bool found = false;
    NodeConnVecConstIter iter = pack->m_children.begin();
    for(; iter != pack->m_children.end(); ++iter) {
      const Node *child = (*iter)->m_n;
      if (!child->IsLoopTunnel()) {
        if (child->IsPossTunnel()) {
          PossTunnel *tun = (PossTunnel*)child;
          if (!tun->m_pset->IsCritSect())
            throw;
          else
            return false;
        }
        else
          throw;
      }
      if (!child->IsPossTunnel(SETTUNIN)) {
        cout << "Child of pack is loop tunnel but not SETTUNIN\n";
        throw;
      }
      const LoopTunnel *tun = (LoopTunnel*)child;
      
      const Loop *loop = (Loop*)(tun->m_pset);
      if (!loop->HasIndepIters()) {
	continue;
      }
      else
	found = true;
      if (loop->m_dim != DIMM)
        throw;
      if (loop->m_comm == m_comm)
        return false;
      else if (loop->m_comm != CORECOMM)
        throw;
      if (!LegalParallelizationNestingDown(loop, m_comm))
        return false;
      /*
      const Split *control = loop->GetControl();
      const unsigned int numExecs = control->NumberOfLoopExecs();
      unsigned int parFactor = NumGroupsInComm(m_comm);
      int numParallelizable = 0;
      for(unsigned int i = 0; i < numExecs; ++i) {
        unsigned int numIters = control->NumIters(i);
        if (numIters >= parFactor)
          ++numParallelizable;
      }
      if ((((double)numParallelizable) / numExecs) < PORTIONPARALLELIZABLE)
        return false;
      */
    }
    return true;
  }
  return false;
}

void ParallelizeMDim::Apply(Poss *poss, Node *node) const
{
  PackBuff *buff = (PackBuff*)node;
  buff->Parallelize(m_comm);
  Pack *pack = (Pack*)buff->Child(0);
  pack->Parallelize(m_comm);
  NodeConnVecConstIter iter = pack->m_children.begin();
  for(; iter != pack->m_children.end(); ++iter) {
    Node *child = (*iter)->m_n;
    LoopTunnel *tun = (LoopTunnel*)child;
    Loop *loop = (Loop*)(tun->m_pset);
    if (!loop->HasIndepIters()) {
      continue;
    }
    else {
      loop->Parallelize(m_comm);
    }
  }
}


bool ParallelizeInnerNDim::CanApply(const Poss *poss, const Node *node) const
{
  const PackBuff *buff = (PackBuff*)node;
  if (buff->m_packMat == PACKABLOCK) {
    if (buff->m_children.size() != 1) {
      cout << "more than one child of A block packbuff!\n";
      throw;
    }
    if (buff->Child(0)->GetNodeClass() != Pack::GetClass()) {
      cout << "child of A block packbuff isn't Pack\n";
      throw;
    }
    const Pack *pack = (Pack*)buff->Child(0);
    if (!pack->m_children.size())
      throw;
    if (!LegalParallelizationNestingUp(pack, m_comm))
      return false;
    NodeConnVecConstIter iter = pack->m_children.begin();
    for(; iter != pack->m_children.end(); ++iter) {
      const Node *child = (*iter)->m_n;
      if (!child->IsDLA())
        throw;
      const DLANode *dla = (DLANode*)child;
      if (dla->IsBLISParallelizable()) {
        return true;
      }
    }
  }
  return false;
}

void ParallelizeInnerNDim::Apply(Poss *poss, Node *node) const
{
  PackBuff *buff = (PackBuff*)node;
  buff->Parallelize(m_comm);
  Pack *pack = (Pack*)buff->Child(0);
  pack->Parallelize(m_comm);
  NodeConnVecConstIter iter = pack->m_children.begin();
  for(; iter != pack->m_children.end(); ++iter) {
    DLANode *child = (DLANode*)((*iter)->m_n);
    if (child->IsBLISParallelizable())
      child->Parallelize(m_comm);
  }
}

bool ParallelizeOuterNDim::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != LoopTunnel::GetClass())
    throw;
  const LoopTunnel *tun = (LoopTunnel*)node;
  if (tun->m_tunType != SETTUNIN)
    return false;
  if (!LegalParallelizationNestingUp(tun, m_comm))
    return false;
  const Loop *loop = (Loop*)(tun->m_pset);
  if (loop->m_comm == m_comm)
    return false;
  if (loop->m_dim != DIMN)
    return false;
  if (loop->m_comm == m_comm)
    return false;
  if (!LegalParallelizationNestingDown(loop, m_comm))
    return false;
  if (loop->m_comm != CORECOMM) {
    //Need to handle multiple par factors on loop
    throw;
  }
  if (!loop->HasIndepIters())
    throw;
  
  /*  const Split *control = loop->GetControl();
   const unsigned int numExecs = control->NumberOfLoopExecs();
   unsigned int parFactor = NumGroupsInComm(m_comm);
   int numParallelizable = 0;
   for(unsigned int i = 0; i < numExecs; ++i) {
   unsigned int numIters = control->NumIters(i);
   if (numIters >= parFactor)
   ++numParallelizable;
   }
   if ((((double)numParallelizable) / numExecs) < PORTIONPARALLELIZABLE)
   return false;
   */
  return true;
}

void ParallelizeOuterNDim::Apply(Poss *poss, Node *node) const
{
  LoopTunnel *tun = (LoopTunnel*)node;
  Loop *loop = (Loop*)(tun->m_pset);
  loop->Parallelize(m_comm);
}


bool ParallelizeK::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != LoopTunnel::GetClass())
    throw;
  const LoopTunnel *tun = (LoopTunnel*)node;
  if (tun->m_tunType != SETTUNIN)
    return false;
  if (!LegalParallelizationNestingUp(tun, m_comm))
    return false;
  const Loop *loop = (Loop*)(tun->m_pset);
  if (loop->m_comm == m_comm)
    return false;
  if (loop->m_dim != DIMK)
    return false;
  if (!LegalParallelizationNestingDown(loop, m_comm))
    return false;
  if (loop->m_comm != CORECOMM) {
    //Need to handle multiple par factors on loop
    throw;
  }
  if (loop->HasIndepIters())
    throw;
  if (loop->OnlyParallelizedOnNonIndependentData())
    return false;
  else {
    return true;
  }
}

void ParallelizeK::Apply(Poss *poss, Node *node) const
{
  LoopTunnel *tun = (LoopTunnel*)node;
  Loop *loop = (Loop*)(tun->m_pset);
  loop->Parallelize(m_comm);
}

bool LegalParallelizationNestingUp(const Node *node, Comm comm)
{
  Poss *poss = node->m_poss;
  while(poss) {
    PSet *pset = poss->m_pset;
    if (!pset)
      break;
    if (pset->IsLoop()) {
      return CommAllowedWithin(((Loop*)pset)->m_comm, comm);
    }
    else if (pset->IsCritSect()) {
      return false;
    }
    poss = pset->m_ownerPoss;
  }
  return true;
}

bool LegalParallelizationNestingDown(const PSet *pset, Comm comm)
{
  PossVecConstIter iter = pset->m_posses.begin();
  bool foundGood = false;
  for(; !foundGood && iter != pset->m_posses.end(); ++iter) {
    const Poss *poss = *iter;
    bool foundBad = false;
    PSetVecConstIter iter2 = poss->m_sets.begin();
    for(; !foundBad && iter2 != poss->m_sets.end(); ++iter2) {
      const PSet *pset = *iter2;
      if (pset->IsLoop()) {
        if (!CommAllowedWithin(comm,((Loop*)pset)->m_comm))
          foundBad = true;
      }
      else if (!pset->IsCritSect()) {
        if (!LegalParallelizationNestingDown(pset, comm))
          foundBad = true;
      }
    }
    if (!foundBad)
      foundGood = true;
  }
  return foundGood;
}
