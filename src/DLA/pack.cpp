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



#include "pack.h"
#include "loopSupport.h"
#include "trxm.h"

string PackTypeToStr(PackType type)
{
  switch(type)
  {
    case (PACKROWPANS):
      return "BLIS_PACKED_ROW_PANELS";
    case (PACKCOLPANS):
      return "BLIS_PACKED_COL_PANELS";
    default:
      throw;
  }
}

string PackSizeToStr(PackSize size)
{
  switch(size)
  {
    case (USEMRSIZE):
      return "gemm_mr";
    case (USENRSIZE):
      return "gemm_nr";
    case (USEKRSIZE):
      return "gemm_kr";
    default:
      throw;
  }
}
/*
 string PackSizeForExtToStr(PackSize size)
 {
 switch(size)
 {
 case (USEMRSIZE):
 return "gemm_extmr";
 case (USENRSIZE):
 return "gemm_extnr";
 case (USEKRSIZE):
 return "gemm_extkr";
 default:
 throw;
 }
 }
 */

Pack::Pack(PackType pack, unsigned int var,
           bool scaleAlpha, bool densify, bool invertDiag,
           bool revUpper, bool revLower)
: m_pack(pack),
m_var(var),
m_scaleAlpha(scaleAlpha),
m_densify(densify),
m_invertDiag(invertDiag),
m_revUpper(revUpper),
  m_revLower(revLower),
  m_comm(CORECOMM)
{
  if (scaleAlpha) {
    cout << "need to scale by alpha but alpha is not a parameter yet!!!\n";
    throw;
  }
  SetLayer(S3LAYER);
}

void Pack::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  const Pack *pack = (Pack*)orig;
  m_pack = pack->m_pack;
  m_var = pack->m_var;
  m_scaleAlpha = pack->m_scaleAlpha;
  m_densify = pack->m_densify;
  m_invertDiag = pack->m_invertDiag;
  m_revUpper = pack->m_revUpper;
  m_revLower = pack->m_revLower;
  m_comm = pack->m_comm;
}

void Pack::FlattenCore(ofstream &out) const
{
  DLANode::FlattenCore(out);
  WRITE(m_pack);
  WRITE(m_var);
  WRITE(m_scaleAlpha);
  WRITE(m_densify);
  WRITE(m_invertDiag);
  WRITE(m_revUpper);
  WRITE(m_revLower);
  WRITE(m_comm);
}

void Pack::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLANode::UnflattenCore(in, info);
  READ(m_pack);
  READ(m_var);
  READ(m_scaleAlpha);
  READ(m_densify);
  READ(m_invertDiag);
  READ(m_revUpper);
  READ(m_revLower);
  READ(m_comm);
}

NodeType Pack::GetType() const
{
  return "Pack " + PackTypeToStr(m_pack);
}

DistType Pack::GetDistType(unsigned int num) const
{
  return InputDistType(1);
}

Phase Pack::MaxPhase() const
{
  return NUMPHASES;
}

void Pack::PrintCode(IndStream &out)
{
  if (m_scaleAlpha)
    throw;
  out.Indent();
  if (m_var == 2)
    *out << "bli_packm_blk_var2";
  else if (m_var == 3)
    *out << "bli_packm_blk_var3";
  else
    throw;
  if (m_comm == CORECOMM) {
    *out << "( ";
  }
  else {
    *out << "_par( ";
  }
  *out << "&BLIS_ONE, &"
       << GetInputName(0).str() << ", &"
       << GetInputName(1).str();
  if (m_comm != CORECOMM) {
    *out << ", " << CommToStr(m_comm);
  }
  *out << " );\n";
}

Name Pack::GetName(unsigned int num) const
{
  if (num > 0)
    throw;
  return GetInputName(1);
}

void Pack::Prop()
{
  if (m_inputs.size() != 2)
    throw;
  if (!IsValidCost(m_cost)) {
    Input(0)->Prop();
    Input(1)->Prop();
    const Sizes *size1 = InputLocalM(0);
    const Sizes *size2 = InputLocalN(0);
    m_cost = (size1->SumProds11(*size2) * (PSIWVAL + PSIRVAL)) / NumCoresInComm(m_comm);
  }
}

void Pack::SanityCheck()
{
  if (m_inputs.size() != 2) {
    cout << "Error on " << GetNodeClass() << endl;
    cout << "Has " << m_inputs.size() << " inputs\n";
    throw;
  }
  DLANode::SanityCheck();
}

void PackBuff::UpdateChildrensInnerMultiple(PackSize size)
{
  for(unsigned int i = 0; i < m_children.size(); ++i) {
    Node *pack = Child(i);
    if (pack->GetNodeClass() != Pack::GetClass()) {
      throw;
    }
    NodeConnVecIter packChildIter = pack->m_children.begin();
    for( ; packChildIter != pack->m_children.end(); ++packChildIter) {
      DLANode *packChild = (DLANode*)((*packChildIter)->m_n);
      if (packChild->IsPossTunnel()) {
        if (!packChild->IsPossTunnel(SETTUNIN))
          throw;
        NodeConnVecIter possTunInIter = packChild->m_children.begin();
        for(; possTunInIter != packChild->m_children.end(); ++possTunInIter) {
          DLANode *possTunIn = (DLANode*)((*possTunInIter)->m_n);
          if (!possTunIn->IsPossTunnel(POSSTUNIN))
            throw;
          NodeConnVecIter childIter = possTunIn->m_children.begin();
          for(; childIter != possTunIn->m_children.end(); ++childIter) {
            DLANode *child = (DLANode*)((*childIter)->m_n);
            if (!child->IsPossTunnel(POSSTUNOUT))
              child->UpdateInnerPackingMultiple(size);
          }
        }
      }
      else
        packChild->UpdateInnerPackingMultiple(size);
    }
  }
}

void PackBuff::PrintCode(IndStream &out)
{
  string name = GetNameStr(0);
  string inputName = GetInputName(0).str();

  if (m_triStruct != GEN) {
    out.Indent();
    
    *out << "bli_obj_set_struc( ";
    switch(m_triStruct) {
      case (TRI):
        *out << "BLIS_TRIANGULAR";
        break;
      case (HERM):
        *out << "BLIS_HERMITIAN";
        break;
      case (SYMM):
        *out << "BLIS_SYMMETRIC";
        break;
      default:
        throw;
    }
    *out << ", " << inputName << " );\n";
  }
  
  if (m_tri != NOTTRI) {
    out.Indent();
    if (m_tri == LOWER) {
      *out << "bli_obj_set_uplo( BLIS_LOWER, " << inputName << " );\n";
    }
    else if (m_tri == UPPER) {
      *out << "bli_obj_set_uplo( BLIS_UPPER, " << inputName << " );\n";
    }
    else
      throw;
  }
  
  if (m_diag == UNIT) {
    out.Indent();
    *out << "bli_obj_set_diag( BLIS_UNIT_DIAG, "
    << inputName << " );\n";
  }
  
  unsigned int indentOffset = 0;
  if (m_comm != CORECOMM) {
    out.Indent();
    *out << "th_barrier( " << CommToStr(m_comm) << " );\n";
    out.Indent();
    *out << "if (th_am_root(" << CommToStr(m_comm) << ")) {\n";
    indentOffset = 1;
  }
  
  out.Indent(indentOffset);
  *out << "bli_packm_init_pack( ";
  if (m_densify)
    *out << "TRUE, ";
  else
    *out << "FALSE, ";
  if (m_invertDiag)
    *out << "BLIS_INVERT_DIAG, ";
  else
    *out << "BLIS_NO_INVERT_DIAG, ";
  *out << PackTypeToStr(m_pack) << ", \n";
  out.Indent(indentOffset+2);
  if (m_revUpper)
    *out << "BLIS_PACK_REV_IF_UPPER, ";
  else
    *out << "BLIS_PACK_FWD_IF_UPPER, ";
  if (m_revLower)
    *out << "BLIS_PACK_REV_IF_LOWER, ";
  else
    *out << "BLIS_PACK_FWD_IF_LOWER, ";
  *out << endl;
  out.Indent(indentOffset+2);
  if (m_packMat == PACKABLOCK)
    *out << "BLIS_BUFFER_FOR_A_BLOCK,";
  else if (m_packMat == PACKBPANEL)
    *out << "BLIS_BUFFER_FOR_B_PANEL,";
  else
    throw;
  *out << endl;
  out.Indent(indentOffset+2);
  *out << PackSizeToStr(m_m) << ", "
  //       << PackSizeForExtToStr(m_m) << ", "
  << PackSizeToStr(m_n) << ", ";
  //       << PackSizeForExtToStr(m_n) << ", ";
  *out << endl;
  out.Indent(indentOffset+2);
  *out << "&" << Input(0)->GetName(InputConnNum(0)).str() << ", &"
  << name << " );\n";

  if (m_comm != CORECOMM) {
    out.Indent();
    *out << "}\n";
    out.Indent();
    *out << "th_broadcast_without_second_barrier(" << CommToStr(m_comm)
	 << ", 0, (void*)(&" << name << "), sizeof(" << name << "));\n";
  }
}

void PackBuff::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<1,1>::Prop();
    m_cost = 0;

    const Sizes *size1 = InputLocalM(0);
    const Sizes *size2 = InputLocalN(0);
    if (m_packMat == PACKABLOCK) {
      if (!(*size1 <= BLIS_MC_BS))
	throw;
      if (!(*size2 <= BLIS_KC_BS))
	throw;
    }
    else if (m_packMat == PACKBPANEL) {
      if (!(*size1 <= BLIS_KC_BS))
	throw;
      if (!(*size2 <= BLIS_NC_BS))
	throw;
    }
    else
      throw;
  }
}

void PackBuff::SanityCheck()
{
  if (m_inputs.size() != 1) {
    cout << "Error on " << GetNodeClass() << endl;
    cout << "Has " << m_inputs.size() << " inputs\n";
    throw;
  }
  DLANode::SanityCheck();
}

PackBuff::PackBuff(string name,
                   PackType pack,
                   PackMat mat,
                   Tri tri, Diag diag, TriStruct triStruct,
                   bool densify, bool invertDiag, bool revUpper, bool revLower,
                   PackSize mSize, PackSize nSize)
  : m_comm(CORECOMM)
{
  m_name.m_type = UNKNOWN;
  if (name.find("_packed") == string::npos) {
    m_name.m_name = name + "_packed";
  }
  else
    m_name.m_name = name;
  m_pack = pack;
  m_packMat = mat;
  m_tri = tri;
  m_diag = diag;
  m_triStruct = triStruct;
  m_densify = densify;
  m_invertDiag = invertDiag;
  m_revUpper = revUpper;
  m_revLower = revLower;
  m_m = mSize;
  m_n = nSize;
}

Name PackBuff::GetName(unsigned int num) const
{
  return m_name;
}

void PackBuff::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  const PackBuff *buff = (PackBuff*)orig;
  m_pack = buff->m_pack;
  m_packMat = buff->m_packMat;
  m_tri = buff->m_tri;
  m_diag = buff->m_diag;
  m_triStruct = buff->m_triStruct;
  m_densify = buff->m_densify;
  m_invertDiag = buff->m_invertDiag;
  m_revUpper = buff->m_revUpper;
  m_revLower = buff->m_revLower;
  m_m = buff->m_m;
  m_n = buff->m_n;
  m_name = buff->m_name;
  m_comm = buff->m_comm;
  DLAOp<1,1>::Duplicate(orig, shallow, possMerging);
}

void PackBuff::FlattenCore(ofstream &out) const
{
  DLAOp<1,1>::FlattenCore(out);
  out << m_name.m_name << endl;
  WRITE(m_pack);
  WRITE(m_packMat);
  WRITE(m_tri);
  WRITE(m_diag);
  WRITE(m_triStruct);
  WRITE(m_m);
  WRITE(m_n);
  WRITE(m_densify);
  WRITE(m_invertDiag);
  WRITE(m_revUpper);
  WRITE(m_revLower);
  WRITE(m_comm);
}
void PackBuff::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<1,1>::UnflattenCore(in, info);
  getline(in, m_name.m_name);
  m_name.m_type = UNKNOWN;
  READ(m_pack);
  READ(m_packMat);
  READ(m_tri);
  READ(m_diag);
  READ(m_triStruct);
  READ(m_m);
  READ(m_n);
  READ(m_densify);
  READ(m_invertDiag);
  READ(m_revUpper);
  READ(m_revLower);
  READ(m_comm);
}

bool LoopInvariantPackBuffMotion::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != PackBuff::GetClass())
    return false;
  if (node->Input(0)->GetNodeClass() != LoopTunnel::GetClass())
    return false;
  const Node *input = node->Input(0);
  const Node *inputInput = input->Input(0);
  if (!input->IsPossTunnel(POSSTUNIN) || !inputInput->IsPossTunnel(SETTUNIN)) {
    return false;
  }
  const LoopTunnel *tun = (LoopTunnel*)(node->Input(0));
  return tun->IsConst();
}

void LoopInvariantPackBuffMotion::Apply(Poss *poss, Node *node) const
{
  LoopTunnel *possTunIn = (LoopTunnel*)(node->Input(0));
  LoopTunnel *setTunIn = (LoopTunnel*)(possTunIn->Input(0));
  PackBuff *buff = (PackBuff*)node;
  Loop *loop = (Loop*)(poss->m_pset);
  Poss *owningPoss = loop->m_ownerPoss;
  PackBuff *newBuff = new PackBuff(buff->m_name.m_name,
                                   buff->m_pack,
                                   buff->m_packMat,
                                   buff->m_tri, buff->m_diag, buff->m_triStruct,
                                   buff->m_densify, buff->m_invertDiag,
                                   buff->m_revUpper, buff->m_revLower,
                                   buff->m_m, buff->m_n);
  owningPoss->AddNode(newBuff);
  newBuff->AddInput(setTunIn->Input(0), setTunIn->InputConnNum(0));
  
  Node *newTun = loop->CreateNewLoopTunnels(newBuff, 0, poss, FULLUP);
  
  buff->RedirectChildren(newTun, 0);
  poss->DeleteChildAndCleanUp(buff, true);
}

bool LoopInvariantPackMotion::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != Pack::GetClass())
    throw;
  const Node *input = node->Input(0);
  if (input->GetNodeClass() != LoopTunnel::GetClass() || !input->IsPossTunnel(POSSTUNIN))
    return false;
  const Node *inputInput = input->Input(0);
  if (inputInput->GetNodeClass() != LoopTunnel::GetClass() || !inputInput->IsPossTunnel(SETTUNIN))
    throw;
  return true;
}

void LoopInvariantPackMotion::Apply(Poss *poss, Node *node) const
{
  if (!poss->m_pset->IsLoop())
    throw;
  Pack *oldPack = (Pack*)node;
  Loop *loop = (Loop*)(poss->m_pset);
  Poss *loopOwner = loop->m_ownerPoss;
  LoopTunnel *oldPossTunIn0 = (LoopTunnel*)(node->Input(0));
  LoopTunnel *oldSetTunIn0 = (LoopTunnel*)(oldPossTunIn0->Input(0));
  LoopTunnel *oldPossTunIn1 = (LoopTunnel*)(node->Input(1));
  LoopTunnel *oldSetTunIn1 = (LoopTunnel*)(oldPossTunIn1->Input(0));
  
  UpStat stat = FULLUP;
  
  NodeConnVecIter iter = oldPack->m_children.begin();
  for(; iter != oldPack->m_children.end(); ++iter) {
    Node *node = (*iter)->m_n;
    if (node->GetNodeClass() == Trxm::GetClass()) {
      const Trxm *trxm = (Trxm*)node;
      if (trxm->m_invert && trxm->GetLayer() == S3LAYER) {
        stat = PARTUP;
        break;
      }
    }
  }
  
  Pack *newPack = new Pack(oldPack->m_pack, oldPack->m_var,
                           oldPack->m_scaleAlpha, oldPack->m_densify,
                           oldPack->m_invertDiag,
                           oldPack->m_revUpper, oldPack->m_revLower);
  newPack->AddInput(oldSetTunIn0->Input(0),oldSetTunIn0->InputConnNum(0));
  newPack->AddInput(oldSetTunIn1->Input(0),oldSetTunIn1->InputConnNum(0));
  loopOwner->AddNode(newPack);
  
  LoopTunnel *newPossTunIn = loop->CreateNewLoopTunnels(newPack,
                                                        0, poss, stat);
  if (!newPossTunIn)
    throw;
  
  node->RedirectChildren(newPossTunIn,0);
  
  node->m_poss->DeleteChildAndCleanUp(node, true);
}

bool CombinePacking::CanApply(const Poss *poss, const Node *node) const
{
  const Pack *pack = (Pack*)node;
  if (node->GetNodeClass() != Pack::GetClass())
    throw;
  const Node *par = node->Input(0);
  unsigned int num = node->InputConnNum(0);
  NodeConnVecConstIter iter = par->m_children.begin();
  for(; iter != par->m_children.end(); ++iter) {
    const NodeConn *conn = *iter;
    const Node *child = conn->m_n;
    if (child != node) {
      if (conn->m_num == num) {
        if (child->GetNodeClass() == Pack::GetClass()) {
          const Pack *pack2 = (Pack*)(child);
          if (pack->m_pack != pack2->m_pack)
            continue;
          if (pack->m_var != pack2->m_var)
            continue;
          if (pack->m_scaleAlpha != pack2->m_scaleAlpha)
            continue;
          if (pack->m_densify != pack2->m_densify)
            continue;
          if (pack->m_invertDiag != pack2->m_invertDiag)
            continue;
          if (pack->m_revUpper != pack2->m_revUpper)
            continue;
          if (pack->m_revLower != pack2->m_revLower)
            continue;
          return true;
        }
      }
    }
  }
  return false;
}

void CombinePacking::Apply(Poss *poss, Node *node) const
{
  Pack *pack = (Pack*)node;
  Node *par = node->Input(0);
  unsigned int num = node->InputConnNum(0);
  NodeConnVecIter iter = par->m_children.begin();
  for(; iter != par->m_children.end(); ++iter) {
    NodeConn *conn = *iter;
    if (conn->m_n != node) {
      if (conn->m_num == num) {
        if (conn->m_n->GetNodeClass() == Pack::GetClass()) {
          Pack *pack2 = (Pack*)(conn->m_n);
          if (pack->m_pack != pack2->m_pack)
            continue;
          if (pack->m_var != pack2->m_var)
            continue;
          if (pack->m_scaleAlpha != pack2->m_scaleAlpha)
            continue;
          if (pack->m_densify != pack2->m_densify)
            continue;
          if (pack->m_invertDiag != pack2->m_invertDiag)
            continue;
          if (pack->m_revUpper != pack2->m_revUpper)
            continue;
          if (pack->m_revLower != pack2->m_revLower)
            continue;
          pack2->RedirectChildren(pack, 0);
          pack2->m_poss->DeleteChildAndCleanUp(pack2, true);
          return;
        }
      }
    }
  }
  throw;
}

bool FoundNodeUp(const Node *node, const Node *find)
{
  if (node == find)
    return true;
  NodeConnVecConstIter iter = node->m_inputs.begin();
  for(; iter != node->m_inputs.end(); ++iter) {
    if (node->IsPossTunnel(POSSTUNIN))
      return false;
    else if (node->IsPossTunnel(SETTUNOUT)) {
      PSet *set = ((PossTunnel*)node)->m_pset;
      NodeVecConstIter iter2 = set->m_inTuns.begin();
      for(; iter2 != set->m_inTuns.end(); ++iter2) {
        if (FoundNodeUp(*iter2,find))
          return true;
      }
    }
    else {
      if (FoundNodeUp((*iter)->m_n, find))
        return true;
    }
  }
  return false;
}

bool CombinePackBuff::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != PackBuff::GetClass())
    return false;
  const PackBuff *packBuff = (PackBuff*)node;
  NodeVecConstIter iter = poss->m_possNodes.begin();
  for(; iter != poss->m_possNodes.end(); ++iter) {
    const Node *possNode = *iter;
    if (node != possNode) {
      if (possNode->GetNodeClass() == PackBuff::GetClass()) {
        const PackBuff *packBuff2 = (PackBuff*)(possNode);
        if (packBuff->m_name.str() != packBuff2->m_name.str())
          continue;
        if (packBuff->m_pack != packBuff2->m_pack
            || packBuff->m_packMat != packBuff2->m_packMat
            || packBuff->m_tri != packBuff2->m_tri
            || packBuff->m_diag != packBuff2->m_diag
            || packBuff->m_triStruct != packBuff2->m_triStruct
            || packBuff->m_densify != packBuff2->m_densify
            || packBuff->m_invertDiag != packBuff2->m_invertDiag
            || packBuff->m_revUpper != packBuff2->m_revUpper
            || packBuff->m_revLower != packBuff2->m_revLower)
        {
          cout << "Two packBuff's with same name"
          << " and different parameters\n";
          cout << packBuff->m_name.str() << endl;
          throw;
        }
        else if (packBuff->m_m == packBuff2->m_m
                 && packBuff->m_n == packBuff2->m_n)
        {
          if (FoundNodeUp(packBuff,packBuff2))
            return false;
          else
            return true;
        }
      }
    }
  }
  return false;
}

void CombinePackBuff::Apply(Poss *poss, Node *node) const
{
  PackBuff *packBuff = (PackBuff*)node;
  NodeVecConstIter iter = poss->m_possNodes.begin();
  for(; iter != poss->m_possNodes.end(); ++iter) {
    Node *possNode = *iter;
    if (node != possNode) {
      if (possNode->GetNodeClass() == PackBuff::GetClass()) {
        PackBuff *packBuff2 = (PackBuff*)(possNode);
        if (packBuff->m_name.str() != packBuff2->m_name.str())
          continue;
        if (packBuff->m_pack != packBuff2->m_pack
            || packBuff->m_packMat != packBuff2->m_packMat
            || packBuff->m_tri != packBuff2->m_tri
            || packBuff->m_diag != packBuff2->m_diag
            || packBuff->m_triStruct != packBuff2->m_triStruct
            || packBuff->m_densify != packBuff2->m_densify
            || packBuff->m_invertDiag != packBuff2->m_invertDiag
            || packBuff->m_revUpper != packBuff2->m_revUpper
            || packBuff->m_revLower != packBuff2->m_revLower)
        {
          cout << "Two packBuff's with same name"
          << " and different sizes\n";
          cout << packBuff->m_name.str() << endl;
          throw;
        }
        else if (packBuff->m_m == packBuff2->m_m
                 && packBuff->m_n == packBuff2->m_n)
        {
          /*
           cout << "packBuff " << packBuff << endl;
           PrintSetOrNodeInputs(packBuff);
           PrintSetOrNodeChildren(packBuff);
           cout << "packBuff2 " << packBuff2 << endl;
           PrintSetOrNodeInputs(packBuff2);
           PrintSetOrNodeChildren(packBuff2);
           */
          
          packBuff2->RedirectChildren(packBuff, 0);
          
          /*
           cout << "packBuff after\n";
           PrintSetOrNodeInputs(packBuff);
           PrintSetOrNodeChildren(packBuff);
           */
          
          packBuff2->m_poss->DeleteChildAndCleanUp(packBuff2, true);
          return;
        }
      }
    }
  }
  throw;
}

bool UnifyPackBuffParams::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != PackBuff::GetClass())
    return false;
  const PackBuff *packBuff = (PackBuff*)node;
  NodeVecConstIter iter = poss->m_possNodes.begin();
  for(; iter != poss->m_possNodes.end(); ++iter) {
    const Node *possNode = *iter;
    if (node != possNode) {
      if (possNode->GetNodeClass() == PackBuff::GetClass()) {
        const PackBuff *packBuff2 = (PackBuff*)(possNode);
        if (packBuff->m_name.str() != packBuff2->m_name.str())
          continue;
        if (packBuff->m_pack != packBuff2->m_pack
            || packBuff->m_tri != packBuff2->m_tri
            || packBuff->m_diag != packBuff2->m_diag
            || packBuff->m_triStruct != packBuff2->m_triStruct
            || packBuff->m_densify != packBuff2->m_densify
            || packBuff->m_invertDiag != packBuff2->m_invertDiag
            || packBuff->m_revUpper != packBuff2->m_revUpper
            || packBuff->m_revLower != packBuff2->m_revLower)
        {
          cout << "Two packBuff's with same name"
          << " and different parameters\n";
          cout << packBuff->m_name.str() << endl;
          throw;
        }
        else if(packBuff->m_m != packBuff2->m_m
                || packBuff->m_n != packBuff2->m_n) {
          if (packBuff->m_m == USEMRSIZE && packBuff->m_n == USENRSIZE
              && packBuff2->m_m == USEKRSIZE) {
            if (packBuff->m_pack != PACKCOLPANS)
              throw;
            return true;
          }
        }
      }
    }
  }
  return false;
}

void UnifyPackBuffParams::Apply(Poss *poss, Node *node) const
{
  PackBuff *packBuff = (PackBuff*)node;
  NodeVecIter iter = poss->m_possNodes.begin();
  for(; iter != poss->m_possNodes.end(); ++iter) {
    Node *possNode = *iter;
    if (node != possNode) {
      if (possNode->GetNodeClass() == PackBuff::GetClass()) {
        PackBuff *packBuff2 = (PackBuff*)(possNode);
        if (packBuff->m_name.str() != packBuff2->m_name.str())
          continue;
        if (packBuff->m_pack != packBuff2->m_pack
            || packBuff->m_tri != packBuff2->m_tri
            || packBuff->m_diag != packBuff2->m_diag
            || packBuff->m_triStruct != packBuff2->m_triStruct
            || packBuff->m_densify != packBuff2->m_densify
            || packBuff->m_invertDiag != packBuff2->m_invertDiag
            || packBuff->m_revUpper != packBuff2->m_revUpper
            || packBuff->m_revLower != packBuff2->m_revLower)
        {
          cout << "Two packBuff's with same name"
          << " and different parameters\n";
          cout << packBuff->m_name.str() << endl;
          throw;
        }
        else if(packBuff->m_m != packBuff2->m_m
                || packBuff->m_n != packBuff2->m_n) {
          if (packBuff->m_m == USEMRSIZE && packBuff->m_n == USENRSIZE
              && packBuff2->m_m == USEKRSIZE) {
            packBuff2->m_m = USEMRSIZE;
            packBuff2->UpdateChildrensInnerMultiple(USEMRSIZE);
            return;
          }
        }
      }
    }
  }
  throw;
}

bool ReuseTrsmPacking::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != Pack::GetClass())
    throw;
  Node *input = node->Input(0);
  if (!input->IsPossTunnel(SETTUNOUT))
    return false;
  PossTunnel *setTun = (PossTunnel*)input;
  if (setTun->m_inputs.size() > 1)
    return false;
  PossTunnel *possTun = (PossTunnel*)(setTun->Input(0));
  if (possTun->GetNodeClass() != Combine::GetClass())
    return false;
  //possTun should be a combine of X in some direction and X_1 should be attached to TrsmBP's
  Node *source = possTun->Input(1);
  if (source->GetNodeClass() != TrxmBP::GetClass()
      || !((TrxmBP*)source)->m_invert)
    return false;
  if (((DLANode*)source)->GetLayer() != m_layer)
    return false;
  unsigned int num = possTun->InputConnNum(1);
  if (num != 1)
    return false;
  else
    return true;
}

void ReuseTrsmPacking::Apply(Poss *poss, Node *node) const
{
  PossTunnel *setTun = (PossTunnel*)(node->Input(0));
  PossTunnel *possTun = (PossTunnel*)(setTun->Input(0));
  TrxmBP *trsm = (TrxmBP*)(possTun->Input(1));
  NodeConnVecIter iter = trsm->m_children.begin();
  for(; iter != trsm->m_children.end(); ++iter) {
    NodeConn *conn = *iter;
    if (conn->m_num == 0) {
      if (conn->m_n->IsPossTunnel(POSSTUNOUT)) {
        PossTunnel *possTunOut = (PossTunnel*)(conn->m_n);
        if (possTunOut->m_children.size() != 1)
          throw;
        if (!possTunOut->m_children[0]->m_n->IsPossTunnel(SETTUNOUT))
          throw;
        PossTunnel *setTunOut = (PossTunnel*)(possTunOut->m_children[0]->m_n);
        node->RedirectChildren(setTunOut,0);
        poss->DeleteChildAndCleanUp(node, false);
        return;
      }
    }
  }
  throw;
}

const Node* FindOtherPackBuffs(const Poss *poss, PackMat pack, const Node *ignore)
{
  NodeVecConstIter iter = poss->m_possNodes.begin();
  for(; iter != poss->m_possNodes.end(); ++iter) {
    const Node *node = *iter;
    if (node != ignore && node->GetNodeClass() == PackBuff::GetClass()) {
      const PackBuff *buff = (PackBuff*)node;
      if (buff->m_packMat == pack)
        return buff;
    }
  }
  PSetVecConstIter iter2 = poss->m_sets.begin();
  for(; iter2 != poss->m_sets.end(); ++iter2) {
    const PSet *set = *iter2;
    PossVecConstIter iter3 = set->m_posses.begin();
    for(; iter3 != set->m_posses.end(); ++iter3) {
      const Node *node = FindOtherPackBuffs(*iter3, pack, ignore);
      if (node != NULL)
        return node;
    }
  }
  return NULL;
}

#if DOSOPHASE
string RenamePackBuff::GetNewName(const PackBuff *buff) const
{
  switch (buff->m_packMat)
  {
    case (PACKABLOCK):
      return "packed_A_blk";
    case (PACKBPANEL):
      return "packed_B_pan";
      break;
    default:
      throw;
  }
}

bool RenamePackBuff::CanApply(const Poss *poss, const Node *node) const
{
  if (CurrPhase != SOPHASE)
    return false;
  if (node->GetNodeClass() != PackBuff::GetClass())
    throw;
  const PackBuff *buff = (PackBuff*)node;
  if (buff->m_name.m_name == GetNewName(buff))
    return false;
  const PackBuff *buff2 = (PackBuff*)FindOtherPackBuffs(poss, buff->m_packMat, node);
  if (buff2 != NULL) {
    throw;
  }
  else
    return true;
}

void RenamePackBuff::Apply(Poss *poss, Node *node) const
{
  PackBuff *buff = (PackBuff*)node;
  buff->m_name.m_name = GetNewName(buff);
}

#endif //DOSOPHASE
