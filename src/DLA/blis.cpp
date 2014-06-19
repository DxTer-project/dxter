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
#if DOBLIS

#define PORTIONPARALLELIZABLE .8

#include "blis.h"
#include "loopSupport.h"
#include "helperNodes.h"

GetUpToDiag::GetUpToDiag(Tri tri, PartDir dir)
{
  m_tri = tri;
  m_sizes = NULL;
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
void GetUpToDiag::ClearDataTypeCache()
{
  if (m_sizes) {
    delete m_sizes;
    m_sizes = NULL;
  }
}

void GetUpToDiag::BuildDataTypeCache()
{
  if (m_sizes)
    return;
  else {
    if (m_dir == PARTDOWN) {
      m_sizes = new Sizes;
      m_sizes->PairwiseSum(*GetInputM(0), *GetInputM(1));
    }
    else if (m_dir == PARTRIGHT) {
      m_sizes = new Sizes;
      m_sizes->PairwiseSum(*GetInputN(0), *GetInputN(1));
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
    m_cost = (PSIRVAL+PSIWVAL) * GetM(0)->SumProds11(*GetN(0));
  }
}

void Copy::PrintCode(IndStream &out)
{
  out.Indent();
  *out << "bli_copym( &" << GetInputNameStr(0) << ", &"
  << GetInputNameStr(1) << " );\n";
}

bool ParallelizeMDim::CanApply(const Node *node) const
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
    const Pack *pack = (Pack*)(buff->Child(0));
    if (!LegalParallelizationNestingUp(pack, m_comm)) {
      return false;
    }
    if (buff->m_comm != CORECOMM || pack->m_comm != CORECOMM)
      return false;
    if (pack->m_children.size() < 1)
      throw;
    return true;
  }
  return false;
}

void ParallelizeMDim::Apply(Node *node) const
{
  PackBuff *buff = (PackBuff*)node;
  buff->Parallelize(m_comm);
  Pack *pack = (Pack*)buff->Child(0);
  pack->Parallelize(m_comm);
  NodeConnVecIter iter = pack->m_children.begin();
  for(; iter != pack->m_children.end(); ++iter) {
    Node *child = (*iter)->m_n;
    LoopTunnel *tun = (LoopTunnel*)child;
    Loop *loop = (Loop*)(tun->m_pset);
    if (loop->HasIndepIters()
        && LegalParallelizationNestingDown(loop, m_comm)
        && loop->m_comm == CORECOMM)
    {
      loop->Parallelize(m_comm);
    }
    //For Trsm, the packed B panel goes into the
    // Trsm loop, which doesn't have indep iters
    //That panel is updated in the loop and then
    // input to the Gemm loop, which does have
    // indep iters
    tun = tun->GetMatchingOutTun();
    NodeConnVecIter iter2 = tun->m_children.begin();
    for(; iter2 != tun->m_children.end(); ++iter2) {
      Node *child2 = (*iter2)->m_n;
      bool skip = false;
      while (!skip && !child2->IsLoopTunnel()) {
        if (child2->IsPossTunnel(SETTUNIN) ||
            child2->IsPossTunnel(POSSTUNIN))
          child2 = child2->m_children[0]->m_n;
        else
          skip = true;
      }
      if (child2->IsLoopTunnel()) {
        loop = (Loop*)(((LoopTunnel*)child2)->m_pset);
        if (loop->HasIndepIters()
            && LegalParallelizationNestingDown(loop, m_comm)
            && loop->m_comm == CORECOMM)
        {
          loop->Parallelize(m_comm);
        }
      }
    }
  }
}


bool ParallelizeInnerNDim::CanApply(const Node *node) const
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
    if (pack->InCriticalSection())
      return false;
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
      if (dla->IsBLISParallelizable() && !dla->IsParallel()) {
        return true;
      }
    }
  }
  return false;
}

void ParallelizeInnerNDim::Apply(Node *node) const
{
  PackBuff *buff = (PackBuff*)node;
  buff->Parallelize(m_comm);
  Pack *pack = (Pack*)buff->Child(0);
  pack->Parallelize(m_comm);
  NodeConnVecConstIter iter = pack->m_children.begin();
  for(; iter != pack->m_children.end(); ++iter) {
    DLANode *child = (DLANode*)((*iter)->m_n);
    if (child->IsBLISParallelizable()) {
      child->Parallelize(m_comm);
    }
  }
}

bool ParallelizeOuterNDim::CanApply(const Node *node) const
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
  if (loop->GetDimName() != DIMN)
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

void ParallelizeOuterNDim::Apply(Node *node) const
{
  LoopTunnel *tun = (LoopTunnel*)node;
  Loop *loop = (Loop*)(tun->m_pset);
  loop->Parallelize(m_comm);
}


bool ParallelizeK::CanApply(const Node *node) const
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
  if (loop->GetDimName() != DIMK)
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

void ParallelizeK::Apply(Node *node) const
{
  LoopTunnel *tun = (LoopTunnel*)node;
  Loop *loop = (Loop*)(tun->m_pset);
  loop->Parallelize(m_comm);
}

bool LegalParallelizationNestingUp(const Node *node, Comm comm)
{
  Poss *poss = node->m_poss;
#if DOSM
  bool foundProcComm = comm != L2COMM;
#endif
  while(poss) {
    PSet *pset = poss->m_pset;
    if (!pset)
      break;
    if (pset->IsLoop()) {
      Loop *loop = (Loop*)pset;
      if (loop->IsParallel()) {
        if (!CommAllowedWithin(loop->m_comm, comm))
          return false;
#if DOSM
        else if (!foundProcComm)
          foundProcComm = loop->m_comm == PROCCOMM;
#endif
      }
    }
    else if (pset->IsCritSect()) {
      return false;
    }
    poss = pset->m_ownerPoss;
  }
#if DOSM
  return foundProcComm;
#else
  return true;
#endif
}

bool LegalParallelizationNestingDown(const PSet *pset, Comm comm)
{
  PossMMapConstIter iter = pset->m_posses.begin();
  bool foundGood = false;
  for(; !foundGood && iter != pset->m_posses.end(); ++iter) {
    const Poss *poss = (*iter).second;
    bool foundBad = false;
    PSetVecConstIter iter2 = poss->m_sets.begin();
    for(; !foundBad && iter2 != poss->m_sets.end(); ++iter2) {
      const PSet *pset = *iter2;
      if (pset->IsLoop()) {
        if (!CommAllowedWithin(comm,((Loop*)pset)->m_comm))
          foundBad = true;
      }
      if (!pset->IsCritSect()) {
        if (!LegalParallelizationNestingDown(pset, comm))
          foundBad = true;
      }
    }
    NodeVecConstIter iter3 = poss->m_possNodes.begin();
    for (; !foundBad && iter3 != poss->m_possNodes.end(); ++iter3) {
      const Node *node = *iter3;
      if (node->IsParallel())
        if (!CommAllowedWithin(comm, node->ParallelComm()))
          foundBad = true;
    }
    if (!foundBad)
      foundGood = true;
  }
  return foundGood;
}
/*
 bool IncreaseParallelizedLoop::CanApply(const Node *node) const
 {
 if (node->GetNodeClass() != Split::GetClass())
 throw;
 const Split *split = (Split*)node;
 if (!split->IsPossTunnel(SETTUNIN))
 return false;
 const Loop *loop = split->GetMyLoop();
 if (loop->m_comm == CORECOMM)
 return false;
 Comm comm = node->WithinParallelism();
 if (comm == CORECOMM)
 #if NUMPROCS > 1
 return loop->m_comm != ALLPROCCOMM;
 #elif NUML2PERPROC > 1
 return loop->m_comm != PROCCOMM;
 #else
 throw;
 #endif
 return !IsImmediateSubComm(comm, loop->m_comm);
 }
 
 void IncreaseParallelizedLoop::Apply(Node *node) const
 {
 Split *split = (Split*)node;
 Loop *loop = split->GetMyLoop();
 Comm comm = node->WithinParallelism();
 if (comm == CORECOMM)
 #if NUMPROCS > 1
 {
 loop->m_comm = ALLL2COMM;
 loop->ReplaceAllComms(L2COMM,L2COMMSUBALLL2);
 }
 #else
 throw;
 #endif
 else
 loop->m_comm = GetSubComm(comm);
 }
 */

bool FoundBarrier(const Node *node, unsigned int input, Comm comm)
{
  NodeConnVecConstIter iter = node->m_inputs.begin();
  for(; iter != node->m_inputs.end(); ++iter) {
    const NodeConn *conn = *iter;
    const Node *in = conn->m_n;
    //Assume this was handles outside of the pset
    if (in->IsPossTunnel(POSSTUNIN))
      continue;
    else if (in->IsPossTunnel(SETTUNOUT)) {
      return false;
    }
    else {
      Comm bar = node->HasBarrier();
      if (bar != CORECOMM &&
          !CommAllowedWithin(bar, comm))
        return false;
    }
  }
  return true;
}

Cost AdditionalCostForBringingIntoL2(Node *node, unsigned int num, Size numAElems, Comm comm)
{
#if DOSM
  DLANode *input = (DLANode*)(node->Input(num));
  if (input->GetNodeClass() != Pack::GetClass())
    throw;
  Pack *pack = (Pack*)input;
  if (pack->m_comm != comm)
    throw;
#if 1
  switch (comm) 
    {
    case (ALLPROCCOMM):
    case (ALLL2COMM):
#if NUMPROCS>1
      return numAElems * PSIRVAL;
#elif NUML2PERPROC>1
      return numAElems * PSIRVAL;
#else
      break;
#endif
    case (PROCCOMM):
#if NUML2PERPROC>1
      return numAElems * PSIRVAL;
#else
      return 0;
#endif	    
	  case (L2COMM):
	  case (L2COMMSUBALLL2):
	  case (CORECOMM):
	    return 0;
	  default:
	    throw;
	  }
#endif
#else
  throw;
#endif
}

Cost AdditionalCostOfBarrier(Comm comm, unsigned int numHits)
{
  if (comm != CORECOMM)
    return 20 * numHits * NumCoresInComm(comm);
  else
    return 0;
}
#endif
