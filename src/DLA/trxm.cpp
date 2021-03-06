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

#include "layers.h"
#if DOBLIS||DOELEM

#include "trxm.h"
#include "blas.h"
#include "elemRedist.h"
#include "string.h"
#include "helperNodes.h"
#include "loopSupport.h"
#include "TriInv.h"
#include "pack.h"
#include "blis.h"

using namespace std;


RealLoop* TrmmLoopLeftVar1(Node *Ain, ConnNum Anum,
                       Node *Bin, ConnNum Bnum,
                       Tri tri, Diag diag, Trans trans,
                       Coef coeff, Type type,
                       Layer layer);

RealLoop* TrmmLoopLeftVar2(Node *Ain, ConnNum Anum,
                       Node *Bin, ConnNum Bnum,
                       Tri tri, Diag diag, Trans trans,
                       Coef coeff, Type type,
                       Layer layer);

RealLoop* TrxmLoopLeftVar3(Node *Ain, ConnNum Anum,
                       Node *Bin, ConnNum Bnum,
                       bool invert,
                       Tri tri, Diag diag, Trans trans,
                       Coef coeff, Type type,
                       Layer layer);

RealLoop* TrmmLoopRightVar1(Node *Ain, ConnNum Anum,
                        Node *Bin, ConnNum Bnum,
                        Tri tri, Diag diag, Trans trans,
                        Coef coeff, Type type,
                        Layer layer);

RealLoop* TrmmLoopRightVar2(Node *Ain, ConnNum Anum,
                        Node *Bin, ConnNum Bnum,
                        Tri tri, Diag diag, Trans trans,
                        Coef coeff, Type type,
                        Layer layer);

RealLoop* TrxmLoopRightVar3(Node *Ain, ConnNum Anum,
                        Node *Bin, ConnNum Bnum,
                        bool invert,
                        Tri tri, Diag diag, Trans trans,
                        Coef coeff, Type type,
                        Layer layer);


RealLoop* TrsmLoopLeftVar1(Node *Ain, ConnNum Anum,
                       Node *Bin, ConnNum Bnum,
                       Tri tri, Diag diag, Trans trans,
                       Coef coeff, Layer layer, Type type);

RealLoop* TrsmLoopLeftVar2(Node *Ain, ConnNum Anum,
                       Node *Bin, ConnNum Bnum,
                       Tri tri, Diag diag, Trans trans,
                       Coef coeff, Layer layer, Type type);

RealLoop* TrsmLoopRightVar1(Node *Ain, ConnNum Anum,
                        Node *Bin, ConnNum Bnum,
                        Tri tri, Diag diag, Trans trans,
                        Coef coeff, Type type,
                        Layer layer);

RealLoop* TrsmLoopRightVar2(Node *Ain, ConnNum Anum,
                        Node *Bin, ConnNum Bnum,
                        Tri tri, Diag diag, Trans trans,
                        Coef coeff, Layer layer, Type type);

RealLoop* Trmm3LoopLeftVar2(Node *Ain, ConnNum Anum,
                        Node *Bin, ConnNum Bnum,
                        Node *Cin, ConnNum Cnum,
                        Tri tri, Diag diag, Trans trans,
                        Coef coeff, Coef beta, Type type,
                        Layer layer);

RealLoop* Trmm3LoopLeftVar3(Node *Ain, ConnNum Anum,
                        Node *Bin, ConnNum Bnum,
                        Node *Cin, ConnNum Cnum,
                        Tri tri, Diag diag, Trans trans,
                        Coef coeff, Coef beta, Type type, Layer layer);



#if DOELEM
bool IsDMTrxm(const Node *node)
{
  if (node->GetNodeClass() == Trxm::GetClass()) {
    const Trxm *trxm = (Trxm*)node;
    return trxm->GetLayer() == DMLAYER;
  }
  return false;
}

bool IsDMTrxm(const Node *node, bool invert)
{
  if (node->GetNodeClass() == Trxm::GetClass()) {
    const Trxm *trxm = (Trxm*)node;
    return trxm->GetLayer() == DMLAYER && trxm->m_invert == invert;
  }
  return false;
}
#endif

Trxm::Trxm(bool invert, Layer layer, Side side, Tri tri, Diag diag, Trans trans, Coef coeff, Type type)
: TrProps(invert, side, tri, diag, trans, coeff, type)
{
  SetLayer(layer);
}

void Trxm::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<2,1>::Duplicate(orig, shallow, possMerging);
  TrProps::Duplicate((Trxm*)orig);
}

void Trxm::FlattenCore(ofstream &out) const
{
  DLAOp<2,1>::FlattenCore(out);
  TrProps::FlattenCore(out);
}


void Trxm::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<2,1>::UnflattenCore(in,info);
  TrProps::UnflattenCore(in,info);
}


Phase Trxm::MaxPhase() const
{
  switch(GetLayer()) {
#if DODPPHASE
    case (ABSLAYER):
    case (DMLAYER):
      return DPPHASE;
    case (SMLAYER):
      return NUMPHASES;
#elif DOSR1PHASE
    case (ABSLAYER):
      return SR1PHASE;
    case (S1LAYER):
      return SR2PHASE;
    case (S2LAYER):
      return SR3PHASE;
#endif
    default:
      throw;

  }
}

bool Trxm::DoNotCullDP() const
{
#if DOELEM
  return GetLayer() == DMLAYER;
#else
  return false;
#endif
}

NodeType Trxm::GetType() const
{
  return (m_invert ? "Trsm" : "Trmm") + LayerNumToStr(GetLayer()) + " " +
  SideToStr(m_side) + " " + TriToStr(m_tri) + " " + TransToStr(m_trans) 
#if DOELEM
    + " " + DistTypeToStr(InputDataType(1).m_dist);
#else
  ;
#endif  
}


void Trxm::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();

    if (m_invert) {
#if DODM
      const Sizes *localM = InputLocalM(0);
      const Sizes *localN = InputLocalN(0);
#else
      const Sizes *localM = GetInputM(0);
      const Sizes *localN = GetInputN(0);
#endif
      if (*localM != *localN) {
	cout << "Trsm::Prop 1\n";
	throw;
      }
#if DODM
      if ((m_side == LEFT && *localN != *InputLocalM(0))
	  || (m_side == RIGHT && *localM != *InputLocalN(0)))
#else
      if ((m_side == LEFT && *localN != *GetInputM(0))
	  || (m_side == RIGHT && *localM != *GetInputN(0)))
#endif
	{
	  cout << "Trsm::Prop 2\n";
	  throw;
	}
    
      if (GetLayer() == ABSLAYER) {
#if DOELEM
	if (InputDataType(1).m_dist != D_MC_MR) {
	  cout << "input not D_MC_MR 7";
	  throw;
	}
#endif
      }
#if DOELEM
      else if (GetLayer() == DMLAYER) {
	if (m_inputs.size() != 2)
	  cout << "1 m_inputs.size() != 2\n";
	else {
	  if (InputDataType(0).m_dist != D_MC_MR)
	    cout << "input not D_MC_MR 1 from " << Input(0)->GetType() << endl;
	  if (InputDataType(1).m_dist != D_MC_MR)
	    cout << "input not D_MC_MR 2 from " << Input(0)->GetType() << endl;
	}
      }
      else if (GetLayer() == SMLAYER) {
	if (InputDataType(0).m_dist == D_VC_STAR && InputDataType(1).m_dist == D_VC_STAR) {
	  if (m_tri != LOWER || m_side != LEFT || m_trans != NORMAL)
	    throw;
	  if (m_inputs.size() != 2)
	    cout << "2 m_inputs.size() != 2\n";
	}
	else {
	  if (m_inputs.size() != 2)
	    cout << "2 m_inputs.size() != 2\n";
	  else {
	    if (InputDataType(0).m_dist != D_STAR_STAR)
	      cout << "input not D_STAR_STAR";
	  }
	}
      }
#endif
    }
#if DOELEM
    else {
      if (GetLayer() == DMLAYER) {
	if (m_inputs.size() != 2)
	  cout << "1 m_inputs.size() != 2\n";
	else {
	  if (InputDataType(0).m_dist != D_MC_MR)
	    cout << "input not D_MC_MR 1 from " << Input(0)->GetType() << endl;
	  if (InputDataType(1).m_dist != D_MC_MR)
	    cout << "input not D_MC_MR 2 from " << Input(0)->GetType() << endl;
	}
      }
      else if (GetLayer() == SMLAYER) {
	if (m_inputs.size() != 2) {
	  cout << "2 m_inputs.size() != 2\n";
	  throw;
	}
	else {
	  if (InputDataType(0).m_dist != D_STAR_STAR) {
	    cout << "input not D_STAR_STAR";
	    throw;
	  }
	  DistType type = InputDataType(1).m_dist;
	  if (type != DataType(0).m_dist)
	    m_poss->MarkInsane();
	}
      }
    }
#endif


    switch (GetLayer()) {
      case (ABSLAYER):
        m_cost = ZERO;
        break;
#if DOELEM
      case (DMLAYER):
        m_cost = ZERO;
        break;
      case (SMLAYER):
        m_cost = GetCost(SMLAYER, m_side, LocalM(0), LocalN(0));
        break;
#elif DOBLIS
      case (S1LAYER):
      case (S2LAYER):
      case (S3LAYER):
        m_cost = ZERO;
        break;
#else
    default:
      throw;
#endif
    }
  }
}

void Trxm::PrintCode(IndStream &out)
{
  out.Indent();
  if (m_invert) {
    switch(GetLayer()) {
      case (ABSLAYER):
        *out << "AbsTrsm( ";
        break;
#if DOELEM
      case (DMLAYER):
        *out << "DistTrsm( ";
        break;
      case (SMLAYER):
        if (InputDataType(0).m_dist == D_VC_STAR && InputDataType(1).m_dist == D_VC_STAR)
          *out << "internal::SpecialLocalTrsm( ";
        else
          *out << "internal::LocalTrsm( ";
        break;
#endif
      default:
        throw;
    }
  }
  else {
    switch (GetLayer()) {
      case (ABSLAYER):
        *out << "AbsTrmm( ";
        break;
#if DOELEM
      case (DMLAYER):
        *out << "DistTrmm( ";
        break;
      case (SMLAYER):
        *out << "internal::LocalTrmm( ";
        break;
#elif DOBLIS
      case (S1LAYER):
        *out << "BLISTrmmLimitedN( ";
        break;
      case (S2LAYER):
        *out << "TrmmRankKUpdate( ";
        break;
      case (S3LAYER):
        *out << "PackedTrmmGebp( ";
        break;
#endif
      default:
        throw;
    }
  }
  
  *out << SideToStr(m_side)
  << ", " << TriToStr(m_tri) << ", "
  << TransToStr(m_trans) << ", "
  << DiagToStr(m_diag) << ", \n" << out.Tabs(2);
  out << m_coeff;
  *out << ", "
  << GetInputName(0).str() << ","
  << "\n" << out.Tabs(2) << "" << GetInputName(1).str()
  << ");\n";
}

Cost Trxm::GetCost(Layer layer, Side side, const Sizes *localMs, const Sizes *localNs)
{
#if DOELEM
  if (layer == SMLAYER) {
    if (side == LEFT)
      return GAMMA * localMs->SumProds21(*localNs);
    else
      return GAMMA * localNs->SumProds21(*localMs);
  }
  else
#endif
    throw;
}

Trmm3::Trmm3(Layer layer, Side side, Tri tri, Diag diag, Trans trans,
             Coef coeff, Coef beta, Type type)
: TrProps(false, side, tri, diag, trans, coeff, type),
m_beta(beta)
{
  SetLayer(layer);
}

void Trmm3::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<3,1>::Duplicate(orig, shallow, possMerging);
  TrProps::Duplicate((Trmm3*)orig);
  const Trmm3 *trmm3 = (Trmm3*)orig;
  m_beta = trmm3->m_beta;
}

void Trmm3::FlattenCore(ofstream &out) const
{
  DLAOp<3,1>::FlattenCore(out);
  TrProps::FlattenCore(out);
  WRITE(m_beta);
}


void Trmm3::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<3,1>::UnflattenCore(in,info);
  TrProps::UnflattenCore(in,info);
  READ(m_beta);
}


Phase Trmm3::MaxPhase() const
{
  switch(GetLayer()) {
#if DOELEM
      throw;
    case (ABSLAYER):
    case (DMLAYER):
      return DPPHASE;
    case (SMLAYER):
      return NUMPHASES;
#elif DOBLIS
    case (ABSLAYER):
      return SR1PHASE;
    case (S1LAYER):
      return SR2PHASE;
    case (S2LAYER):
      return SR3PHASE;
#endif
    default:
      throw;
  }
}

#if DOELEM
bool Trmm3::DoNotCullDP() const
{
  return GetLayer() == DMLAYER;
}
#endif

NodeType Trmm3::GetType() const
{
  return "Trmm3" + LayerNumToStr(GetLayer()) + " " +
  SideToStr(m_side) + " " + TriToStr(m_tri) + " " + TransToStr(m_trans) 
#if DOELEM
    + " " + DistTypeToStr(InputDataType(1).m_dist);
#elif DOBLIS
  ;
#endif
}


void Trmm3::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3,1>::Prop();
    switch (GetLayer()) {
      case (ABSLAYER):
        m_cost = ZERO;
        break;
#if DOELEM
      case (DMLAYER):
      case (SMLAYER):
        throw;
        break;
#elif DOBLIS
      case (S1LAYER):
      case (S2LAYER):
      case (S3LAYER):
        m_cost = ZERO;
        break;
#else
    default:
      throw;
#endif
    }
  }
}

void Trmm3::PrintCode(IndStream &out)
{
  out.Indent();
  switch (GetLayer()) {
    case (ABSLAYER):
      *out << "AbsTrmm3( ";
      break;
#if DOELEM
    case (DMLAYER):
      *out << "DistTrmm3( ";
      break;
    case (SMLAYER):
      *out << "internal::LocalTrmm3( ";
      break;
#elif DOBLIS
    case (S1LAYER):
      *out << "BLISTrmm3LimitedN( ";
      break;
    case (S2LAYER):
      *out << "Trmm3RankKUpdate( ";
      break;
#endif
    default:
      throw;
  }
  
  *out << SideToStr(m_side)
  << ", " << TriToStr(m_tri) << ", "
  << TransToStr(m_trans) << ", "
  << DiagToStr(m_diag) << ", \n" << out.Tabs(2);
  out << m_coeff;
  *out << ", "
  << GetInputName(0).str() << ","
  << "\n" << out.Tabs(2) << "" << GetInputName(1).str()  << ","
  << "\n" << out.Tabs(2) << "" << GetInputName(2).str()
  << ");\n";
}


string TrxmLoopExp::GetType() const
{
  string str;
  switch(m_dim) {
    case(1):
      str = "Trxm Loop Exp var 1";
      break;
    case(2):
      str = "Trxm Loop Exp var 2";
      break;
    case(3):
      str = "Trxm Loop Exp var 3";
      break;
    default:
      throw;
  }
  return str  + (m_side == LEFT ? " left" : " right");
}

bool TrxmLoopExp::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Trxm::GetClass()) {
    const Trxm *trxm = (Trxm*)node;
    if (trxm->GetLayer() != m_fromLayer)
      return false;
    if (trxm->m_side != m_side)
      return false;
    if (m_dim == 1 && trxm->m_side == LEFT 
	&& (trxm->m_tri != LOWER || trxm->m_trans != NORMAL))
      return false;
    return true;
  }
  return false;
}

void TrxmLoopExp::Apply(Node *node) const
{
  Trxm *trxm = (Trxm*)node;
  RealLoop *loop;
  
  NodeConn *connA, *connB;
  connA = trxm->m_inputs[0];
  connB = trxm->m_inputs[1];
  
  
  if (trxm->m_invert) {
    switch(m_dim) {
      case(1):
        if (trxm->m_side == LEFT)
          loop = TrsmLoopLeftVar1(connA->m_n, connA->m_num,
                                  connB->m_n, connB->m_num,
                                  trxm->m_tri, trxm->m_diag, trxm->m_trans,
                                  trxm->m_coeff, m_toLayer,
                                  trxm->m_type);
        else {
          loop = TrsmLoopRightVar1(connA->m_n, connA->m_num,
                                   connB->m_n, connB->m_num,
                                   trxm->m_tri, trxm->m_diag, trxm->m_trans,
                                   trxm->m_coeff, trxm->m_type, m_toLayer);
        }
        break;
      case(2):
        if (trxm->m_side == LEFT)
          loop = TrsmLoopLeftVar2(connA->m_n, connA->m_num,
                                  connB->m_n, connB->m_num,
                                  trxm->m_tri, trxm->m_diag, trxm->m_trans,
                                  trxm->m_coeff, m_toLayer,
                                  trxm->m_type);
        else
          loop = TrsmLoopRightVar2(connA->m_n, connA->m_num,
                                   connB->m_n, connB->m_num,
                                   trxm->m_tri, trxm->m_diag, trxm->m_trans,
                                   trxm->m_coeff, m_toLayer,
                                   trxm->m_type);
        break;
        
      case(3):
        if (trxm->m_side == LEFT)
          loop = TrxmLoopLeftVar3(connA->m_n, connA->m_num,
                                  connB->m_n, connB->m_num,
                                  true,
                                  trxm->m_tri, trxm->m_diag, trxm->m_trans,
                                  trxm->m_coeff, trxm->m_type, m_toLayer);
        else
          loop = TrxmLoopRightVar3(connA->m_n, connA->m_num,
                                   connB->m_n, connB->m_num,
                                   true,
                                   trxm->m_tri, trxm->m_diag, trxm->m_trans,
                                   trxm->m_coeff, trxm->m_type, m_toLayer);
        break;
      default:
        throw;
    }
  }
  else {
    switch(m_dim) {
      case(1):
        if (trxm->m_side == LEFT)
          loop = TrmmLoopLeftVar1(connA->m_n, connA->m_num,
                                  connB->m_n, connB->m_num,
                                  trxm->m_tri, trxm->m_diag, trxm->m_trans,
                                  trxm->m_coeff, trxm->m_type,
                                  m_toLayer);
        else
          loop = TrmmLoopRightVar1(connA->m_n, connA->m_num,
                                   connB->m_n, connB->m_num,
                                   trxm->m_tri, trxm->m_diag, trxm->m_trans,
                                   trxm->m_coeff, trxm->m_type,
                                   m_toLayer);
        break;
      case(2):
        if (trxm->m_side == LEFT)
          loop = TrmmLoopLeftVar2(connA->m_n, connA->m_num,
                                  connB->m_n, connB->m_num,
                                  trxm->m_tri, trxm->m_diag, trxm->m_trans,
                                  trxm->m_coeff, trxm->m_type,
                                  m_toLayer);
        else
          loop = TrmmLoopRightVar2(connA->m_n, connA->m_num,
                                   connB->m_n, connB->m_num,
                                   trxm->m_tri, trxm->m_diag, trxm->m_trans,
                                   trxm->m_coeff, trxm->m_type,
                                   m_toLayer);
        break;
      case(3):
        if (trxm->m_side == LEFT)
          loop = TrxmLoopLeftVar3(connA->m_n, connA->m_num,
                                  connB->m_n, connB->m_num,
                                  false,
                                  trxm->m_tri, trxm->m_diag, trxm->m_trans,
                                  trxm->m_coeff, trxm->m_type,
                                  m_toLayer);
        else
          loop = TrxmLoopRightVar3(connA->m_n, connA->m_num,
                                   connB->m_n, connB->m_num,
                                   false,
                                   trxm->m_tri, trxm->m_diag, trxm->m_trans,
                                   trxm->m_coeff, trxm->m_type,
                                   m_toLayer);
        break;
      default:
        throw;
    }
  }
  
  node->m_poss->AddPSet(loop);
  
  node->RedirectChildren(loop->OutTun(1),0);
  node->m_poss->DeleteChildAndCleanUp(node);
}


string Trmm3LoopExp::GetType() const
{
  return "Trmm3 Loop Exp " + LayerNumToStr(m_fromLayer)
  + " -> " + LayerNumToStr(m_toLayer);
}

bool Trmm3LoopExp::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Trmm3::GetClass()) {
    const Trmm3 *trmm3 = (Trmm3*)node;
    if (trmm3->m_side == RIGHT)
      throw;
    return trmm3->GetLayer() == m_fromLayer;
  }
  return false;
}

void Trmm3LoopExp::Apply(Node *node) const
{
  Trmm3 *trmm3 = (Trmm3*)node;
  RealLoop *loop = NULL;
  
  NodeConn *connA, *connB, *connC;
  connA = trmm3->m_inputs[0];
  connB = trmm3->m_inputs[1];
  connC = trmm3->m_inputs[2];
  
  if (m_var == 2)
    loop = Trmm3LoopLeftVar2(connA->m_n, connA->m_num,
                             connB->m_n, connB->m_num,
                             connC->m_n, connC->m_num,
                             trmm3->m_tri, trmm3->m_diag, trmm3->m_trans,
                             trmm3->m_coeff, trmm3->m_beta, trmm3->m_type,
                             m_toLayer);
  else if (m_var == 3)
    loop = Trmm3LoopLeftVar3(connA->m_n, connA->m_num,
                             connB->m_n, connB->m_num,
                             connC->m_n, connC->m_num,
                             trmm3->m_tri, trmm3->m_diag, trmm3->m_trans,
                             trmm3->m_coeff, trmm3->m_beta, trmm3->m_type,
                             m_toLayer);
  
  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(2),0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

#if DOELEM
string DistTrxmToLocalTrxm::GetType() const
{
  return "Distributed Trxm to Local Trxm " + DistTypeToStr(m_leftType) + ", " + DistTypeToStr(m_rightType);
}

bool DistTrxmToLocalTrxm::CanApply(const Node *node) const
{
  return IsDMTrxm(node);
}

void DistTrxmToLocalTrxm::Apply(Node *node) const
{
  Trxm *trxm = (Trxm*)node;
  RedistNode *node1 = new RedistNode(D_STAR_STAR);
  RedistNode *node2 = new RedistNode(trxm->m_side == LEFT ? m_leftType: m_rightType);
  Trxm *node3 = new Trxm(trxm->m_invert, SMLAYER, trxm->m_side, trxm->m_tri, trxm->m_diag, trxm->m_trans, trxm->m_coeff, trxm->m_type);
  RedistNode *node4 = new RedistNode(D_MC_MR);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node2->AddInput(node->Input(1),node->InputConnNum(1));
  node3->AddInput(node1,0);
  node3->AddInput(node2,0);
  node4->AddInput(node3,0);
  node->m_poss->AddNodes(4, node1, node2, node3, node4);
  node->RedirectChildren(node4,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

Cost DistTrxmToLocalTrxm::RHSCostEstimate(const Node *node) const
{
  const Trxm *trxm = (Trxm*)node;
  const Sizes *A1 = trxm->GetInputM(0);
  const Sizes *A2 = trxm->GetInputN(0);
  const Sizes *B1 = trxm->GetInputM(1);
  const Sizes *B2 = trxm->GetInputN(1);
  Cost cost = RedistNode::GetCost(D_MC_MR, D_STAR_STAR, A1, A2);
  DistType dist = trxm->m_side == LEFT ? m_leftType: m_rightType;
  cost += RedistNode::GetCost(D_MC_MR, dist, B1, B2);
  Sizes localB1, localB2;
  GetLocalSizes(dist, B1, B2, localB1, localB2);
  cost += Trxm::GetCost(SMLAYER, trxm->m_side, &localB1, &localB2);
  cost += RedistNode::GetCost(dist, D_MC_MR, B1, B2);
  return cost;
}

bool DistTrsmToSpecialLocalTrsm::CanApply(const Node *node) const
{
  if (!IsDMTrxm(node, true))
    return false;
  Trxm *trxm = (Trxm*)node;
  if (trxm->m_side != LEFT || trxm->m_tri != LOWER || trxm->m_trans != NORMAL)
    return false;
  return true;
}

void DistTrsmToSpecialLocalTrsm::Apply(Node *node) const
{
  Trxm *trxm = (Trxm*)node;
  if (trxm->m_side != LEFT || trxm->m_tri != LOWER || trxm->m_trans != NORMAL)
    throw;
  RedistNode *node1 = new RedistNode(D_VC_STAR);
  RedistNode *node2 = new RedistNode(D_VC_STAR);
  Trxm *node3 = new Trxm(true, SMLAYER,
                         trxm->m_side, trxm->m_tri,
                         trxm->m_diag, trxm->m_trans,
                         trxm->m_coeff, trxm->m_type);
  RedistNode *node4 = new RedistNode(D_MC_MR);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node2->AddInput(node->Input(1),node->InputConnNum(1));
  node3->AddInput(node1,0);
  node3->AddInput(node2,0);
  node4->AddInput(node3,0);
  node->m_poss->AddNodes(4, node1, node2, node3, node4);
  node->RedirectChildren(node4,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

LocalTrmmAcc::LocalTrmmAcc(Side side, Tri tri, Diag diag, Trans trans, Coef coeff, Type type)
: m_side(side), m_tri(tri), m_diag(diag), m_trans(trans), m_coeff(coeff), m_type(type)
{
  
}

void LocalTrmmAcc::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<3,1>::Duplicate(orig, shallow, possMerging);
  const LocalTrmmAcc *trmm = (LocalTrmmAcc*)orig;
  m_side = trmm->m_side;
  m_tri = trmm->m_tri;
  m_trans = trmm->m_trans;
  m_coeff = trmm->m_coeff;
  m_type = trmm->m_type;
}

void LocalTrmmAcc::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3,1>::Prop();
    bool isTrans = IsTransType(DataType(0).m_dist);
    m_cost = GetCost(m_side, isTrans, LocalM(0), LocalN(0));
  }
}

Cost LocalTrmmAcc::GetCost(Side side, bool isTrans, const Sizes *localMs, const Sizes *localNs)
{
  bool isLeft = side == LEFT;
  if (isLeft) {
    if (isTrans) {
      return GAMMA * localNs->SumProds21(*localMs);
    }
    else {
      return GAMMA * localMs->SumProds21(*localNs);
    }
  }
  else {
    if (isTrans) {
      return GAMMA * localMs->SumProds21(*localMs);
    }
    else {
      return GAMMA * localNs->SumProds21(*localNs);
    }
  }
  
}

void LocalTrmmAcc::PrintCode(IndStream &out)
{
  out.Indent();
  *out << "internal::LocalTrmmAccumulate";
  
  if (m_side == LEFT)
    *out << "L";
  else
    *out << "R";
  
  if (m_tri == LOWER)
    *out << "L";
  else
    *out << "U";
  
  if (m_trans == NORMAL)
    *out << "N";
  else
    *out << "T";
  
  *out << "( " << TransToStr(m_trans) << ", NON_UNIT, ";
  out << m_coeff;
  *out << ", " << endl << out.Tabs(1);
  
  *out << GetInputName(0).str() << ", " << GetInputName(1).str()
  << ", " << GetInputName(2).str() << " );\n";
}

void LocalTrmmAcc::FlattenCore(ofstream &out) const
{
  DLAOp<3,1>::FlattenCore(out);
  WRITE(m_side);
  WRITE(m_tri);
  WRITE(m_trans);
  WRITE(m_coeff);
  WRITE(m_type);
}

void LocalTrmmAcc::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<3,1>::UnflattenCore(in, info);
  READ(m_side);
  READ(m_tri);
  READ(m_trans);
  READ(m_coeff);
  READ(m_type);
}

bool DistTrmmToLocalTrmmStatA::CanApply(const Node *node) const
{
  return IsDMTrxm(node, false);
}

void DistTrmmToLocalTrmmStatA::Apply(Node *node) const
{
  Trxm *trmm = (Trxm*)node;
  
  Node *Ain = trmm->Input(0);
  ConnNum Anum = trmm->InputConnNum(0);
  Node *Bin = trmm->Input(1);
  ConnNum Bnum = trmm->InputConnNum(1);
  
  
  if (trmm->m_side == LEFT) {
    if (trmm->m_trans == NORMAL) {
      RedistNode *redist = new RedistNode(D_VR_STAR);
      redist->AddInput(Bin, Bnum);
      
      RedistNode *redist2 = new RedistNode(D_STAR_MR_T);
      redist2->AddInput(redist, 0);
      
      TempVarNode *tmp1 = new TempVarNode(D_MC_STAR);
      tmp1->AddInput(Bin, Bnum);
      
      LocalTrmmAcc *acc = new LocalTrmmAcc(LEFT, trmm->m_tri, trmm->m_diag, NORMAL, trmm->m_coeff, trmm->m_type);
      acc->AddInputs(6,
                     Ain, Anum,
                     redist2, 0,
                     tmp1, 0);
      
      SumScatterFrom *from = new SumScatterFrom;
      from->AddInputs(4,
                      acc, 0,
                      Bin, Bnum);
      
      node->m_poss->AddNodes(5, redist, redist2, tmp1, acc, from);
      
      node->RedirectChildren(from, 0);
      node->m_poss->DeleteChildAndCleanUp(node);
    }
    else {
      RedistNode *redist = new RedistNode(D_MC_STAR);
      redist->AddInput(Bin, Bnum);
      
      TempVarNode *tmp1 = new TempVarNode(D_MR_STAR);
      tmp1->AddInput(Bin, Bnum);
      
      LocalTrmmAcc *acc = new LocalTrmmAcc(LEFT, trmm->m_tri, trmm->m_diag,
                                           trmm->m_trans, trmm->m_coeff, trmm->m_type);
      acc->AddInputs(6,
                     Ain, Anum,
                     redist, 0,
                     tmp1, 0);
      
      TempVarNode *tmp2 = new TempVarNode(D_MR_MC);
      tmp2->AddInput(Bin, Bnum);
      
      SumScatterFrom *from = new SumScatterFrom;
      from->AddInputs(4,
                      acc, 0,
                      tmp2, 0);
      
      RedistNode *redist2 = new RedistNode(D_MC_MR);
      redist2->AddInput(from, 0);
      
      node->m_poss->AddNodes(6, redist, tmp1, acc, tmp2, from, redist2);
      
      node->RedirectChildren(redist2, 0);
      node->m_poss->DeleteChildAndCleanUp(node);
    }
  }
  else if (trmm->m_side == RIGHT) {
    if (trmm->m_trans == NORMAL) {
      RedistNode *redist = new RedistNode(D_STAR_VC);
      redist->AddInput(Bin, Bnum);
      
      RedistNode *redist2 = new RedistNode(D_STAR_MC);
      redist2->AddInput(redist, 0);
      
      TempVarNode *tmp1 = new TempVarNode(D_MR_STAR_T);
      tmp1->AddInput(Bin, Bnum);
      
      LocalTrmmAcc *acc = new LocalTrmmAcc(RIGHT, trmm->m_tri, trmm->m_diag, NORMAL,
                                           trmm->m_coeff, trmm->m_type);
      acc->AddInputs(6,
                     Ain, Anum,
                     redist2, 0,
                     tmp1, 0);
      
      TempVarNode *tmp2 = new TempVarNode(D_MR_MC_T);
      tmp2->AddInput(Bin, Bnum);
      
      SumScatterFrom *from = new SumScatterFrom;
      from->AddInputs(4,
                      acc, 0,
                      tmp2, 0);
      
      RedistNode *redist3 = new RedistNode(D_MC_MR);
      redist3->AddInput(from, 0);
      
      node->m_poss->AddNodes(7, redist, redist2, tmp1, acc, tmp2, from, redist3);
      
      node->RedirectChildren(redist3, 0);
      node->m_poss->DeleteChildAndCleanUp(node);
    }
    else {
      bool trans = trmm->m_trans == TRANS;
      RedistNode *redist = new RedistNode(trans ? D_MR_STAR_T : D_MR_STAR_H);
      redist->AddInput(Bin, Bnum);
      
      TempVarNode *tmp1 = new TempVarNode(trans ? D_MC_STAR_T : D_MC_STAR_H);
      tmp1->AddInput(Bin, Bnum);
      
      LocalTrmmAcc *acc = new LocalTrmmAcc(RIGHT, trmm->m_tri, trmm->m_diag,
                                           trmm->m_trans, trmm->m_coeff, trmm->m_type);
      acc->AddInputs(6,
                     Ain, Anum,
                     redist, 0,
                     tmp1, 0);
      
      TempVarNode *tmp2 = new TempVarNode(trans ? D_MC_MR_T : D_MC_MR_H);
      tmp2->AddInput(Bin, Bnum);
      
      SumScatterFrom *from = new SumScatterFrom;
      from->AddInputs(4,
                      acc, 0,
                      tmp2, 0);
      
      RedistNode *redist3 = new RedistNode(trans ? D_MR_MC_T : D_MR_MC_H);
      redist3->AddInput(from, 0);
      
      RedistNode *redist4 = new RedistNode(D_MC_MR);
      redist4->AddInput(redist3, 0);
      
      node->m_poss->AddNodes(7, redist, tmp1, acc, tmp2, from, redist3, redist4);
      
      node->RedirectChildren(redist4, 0);
      node->m_poss->DeleteChildAndCleanUp(node);
    }
  }
  else
    throw;
}

Cost DistTrmmToLocalTrmmStatA::RHSCostEstimate(const Node *node) const
{
  const Trxm *trmm = (Trxm*)node;
  const Sizes *B1 = trmm->GetInputM(1);
  const Sizes *B2 = trmm->GetInputN(1);
  if(trmm->m_side == LEFT) {
    if (trmm->m_trans == NORMAL) {
      Cost cost = RedistNode::GetCost(D_MC_MR, D_VR_STAR, B1, B2);
      cost += RedistNode::GetCost(D_VR_STAR, D_STAR_MR_T, B1, B2);
      Sizes localB1, localB2;
      GetLocalSizes(D_MC_STAR, B1, B2, localB1, localB2);
      cost += LocalTrmmAcc::GetCost(trmm->m_side, false, &localB1, &localB2);
      cost += SumScatterFrom::GetCost(D_MC_MR, D_MC_STAR, trmm->InputLocalM(1), trmm->InputLocalN(1));
      return cost;
    }
    else {
      Cost cost = RedistNode::GetCost(D_MC_MR, D_MC_STAR, B1, B2);
      Sizes localB1, localB2;
      GetLocalSizes(D_MR_STAR, B1, B2, localB1, localB2);
      cost += LocalTrmmAcc::GetCost(trmm->m_side, false, &localB1, &localB2);
      localB1.ClearSizes();
      localB2.ClearSizes();
      GetLocalSizes(D_MR_MC, B1, B2, localB1, localB2);
      cost += SumScatterFrom::GetCost(D_MR_MC, D_MR_STAR, &localB1, &localB2);
      cost += RedistNode::GetCost(D_MR_MC, D_MC_MR, B1, B2);
      return cost;
    }
  }
  else if (trmm->m_side == RIGHT) {
    if (trmm->m_trans == NORMAL) {
      Cost cost = RedistNode::GetCost(D_MC_MR, D_STAR_VR, B1, B2);
      cost += RedistNode::GetCost(D_STAR_VR, D_STAR_MC, B1, B2);
      Sizes localB1, localB2;
      GetLocalSizes(D_MR_STAR_T, B1, B2, localB1, localB2);
      cost += LocalTrmmAcc::GetCost(trmm->m_side, true, &localB1, &localB2);
      localB1.ClearSizes();
      localB2.ClearSizes();
      GetLocalSizes(D_MR_MC_T, B1, B2, localB1, localB2);
      cost += SumScatterFrom::GetCost(D_MR_MC_T, D_MR_STAR_T, &localB1, &localB2);
      cost += RedistNode::GetCost(D_MR_MC_T, D_MC_MR, B1, B2);
      return cost;
    }
    else {
      bool trans = trmm->m_trans == TRANS;
      DistType type = trans ? D_MR_STAR_T : D_MR_STAR_H;
      Cost cost = RedistNode::GetCost(D_MC_MR, type, B1, B2);
      type = trans ? D_MC_STAR_T : D_MC_STAR_H;
      Sizes localB1, localB2;
      GetLocalSizes(type, B1, B2, localB1, localB2);
      cost += LocalTrmmAcc::GetCost(trmm->m_side, true, &localB1, &localB2);
      DistType type2 = trans ? D_MC_MR_T : D_MC_MR_H;
      localB1.ClearSizes();
      localB2.ClearSizes();
      GetLocalSizes(type2, B1, B2, localB1, localB2);
      cost += SumScatterFrom::GetCost(type2, type, &localB1, &localB2);
      type = trans ? D_MR_MC_T : D_MR_MC_H;
      cost += RedistNode::GetCost(type2, type, B1, B2);
      cost += RedistNode::GetCost(type, D_MC_MR, B1, B2);
      return cost;
    }
  }
  else
    throw;
}


string TrxmTrans::GetTransType() const
{
  return "Trxm";
}

bool TrxmTrans::CanApply(const Node *node) const
{
  //Any change should also be made to TrsmTrans
  if (node->GetNodeClass() != Trxm::GetClass())
    return false;
  const Trxm *trxm = (Trxm*)node;
  if (trxm->GetLayer() != SMLAYER)
    return false;
  if (trxm->m_trans == TRANS && m_trans == CONJTRANS)
    return false;
  else if (trxm->m_trans == CONJTRANS && m_trans != CONJTRANS)
    return false;
  DLANode *source = (DLANode*)node->Input(1);
  if (!source->CanTrans())
    return false;
  for(ConnNum i = 0; i < node->m_children.size(); ++i) {
    DLANode *child = (DLANode*)(node->Child(i));
    if (child->GetNodeClass() != RedistNode::GetClass()) {
      if (!child->CanTransposeInputs())
        return false;
    }
    else {
      RedistNode *redist = (RedistNode*)(node->Child(i));
      DistType inType = ((DLANode*)node)->DataType(0).m_dist;
      if (!CanTrans(inType, redist->m_info.m_dist, true))
        return false;
    }
  }
  return true;
}

void TrxmTrans::PreApply(Node *node) const
{
  //Any change should also be made to TrsmTrans
  Trxm *trxm = (Trxm*)node;
  if (trxm->m_side == LEFT)
    trxm->m_side = RIGHT;
  else
    trxm->m_side = LEFT;
  if (trxm->m_trans == NORMAL)
    trxm->m_trans = m_trans;
  else
    trxm->m_trans = NORMAL;
}

void TrxmTrans::PostApply(Node *node) const
{
  
}


bool DTrmmToTrsm::CanApply(const Node *node) const
{
  if (!IsDMTrxm(node,false))
    return false;
  DLANode *parent = (DLANode*)(node->Input(0));
  if (parent &&
      parent->GetNodeClass() == TriInv::GetClass() &&
      parent->GetLayer() == DMLAYER)
  {
    return true;
  }
  else
    return false;
}

void DTrmmToTrsm::Apply(Node *node) const
{
  if (!IsDMTrxm(node,false))
    throw;
  Trxm *trmm = (Trxm*)node;
  TriInv *triInv = NULL;
  DLANode *parent = trmm->FindNonRedistParent(0);
  triInv = (TriInv*)parent;
  
  Trxm *trsm = new Trxm(true, DMLAYER, trmm->m_side, trmm->m_tri, trmm->m_diag,
                        trmm->m_trans, trmm->m_coeff,trmm->m_type);
  trsm->AddInput(triInv->Input(0),triInv->InputConnNum(0));
  trsm->AddInput(trmm->Input(1), trmm->InputConnNum(1));
  
  trmm->m_poss->AddNode(trsm);
  
  trmm->RedirectChildren(trsm);
  
  trmm->m_poss->DeleteChildAndCleanUp(trmm);
}


bool LTrmmToTrsm::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Trxm::GetClass())
    return false;
  const Trxm *trmm = (Trxm*)node;
  if (trmm->m_invert)
    return false;
  if (trmm->GetLayer() != SMLAYER)
    return false;
  DLANode *parent = ((DLANode*)node)->FindNonRedistParent(0);
  if (parent &&
      parent->GetNodeClass() == TriInv::GetClass() &&
      parent->GetLayer() == SMLAYER)
    return true;
  else
    return false;
}

void LTrmmToTrsm::Apply(Node *node) const
{
  Trxm *trmm = (Trxm*)node;
  DLANode *parent = trmm->FindNonRedistParent(0);
  TriInv *triInv = (TriInv*)parent;
  
  Trxm *trsm = new Trxm(true, SMLAYER, trmm->m_side, trmm->m_tri, trmm->m_diag,
                        trmm->m_trans, trmm->m_coeff,trmm->m_type);
  trsm->AddInput(triInv->Input(0),triInv->InputConnNum(0));
  trsm->AddInput(trmm->Input(1), trmm->InputConnNum(1));
  
  trmm->m_poss->AddNode(trsm);
  
  trmm->RedirectChildren(trsm);
  trmm->m_poss->DeleteChildAndCleanUp(trmm);
}

#endif


RealLoop* TrmmLoopLeftVar1(Node *Ain, ConnNum Anum,
                       Node *Bin, ConnNum Bnum,
                       Tri tri, Diag diag, Trans trans,
                       Coef coeff, Type type,
                       Layer layer)
{
  bool rev = (tri == UPPER && trans != NORMAL) || (tri == LOWER && trans == NORMAL);
  
  SplitSingleIter *splitA = new SplitSingleIter(rev ? PARTDIAGBACK: PARTDIAG, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(rev ? PARTUPWARD : PARTDOWN, POSSTUNIN, true);
  splitB->AddInput(Bin, Bnum);
  if (!rev)
    splitB->SetUpStats(FULLUP, FULLUP,
                       NOTUP, NOTUP);
  else
    splitB->SetUpStats(NOTUP, NOTUP,
                       FULLUP, FULLUP);
  
  Node *trmm;
  trmm = new Trxm(false, layer, LEFT, tri, diag, trans, coeff, type);
  
  trmm->AddInputs(4,
                  splitA, 4,
                  splitB, 1);
  
  Node *gemm;
  gemm = new Gemm(layer, trans, NORMAL, coeff, COEFONE, type);
  
  if (tri == LOWER) {
    if (trans == NORMAL) {
      gemm->AddInputs(6,
                      splitA, 1,
                      splitB, 0,
                      trmm, 0);
    }
    else {
      gemm->AddInputs(6,
                      splitA, 5,
                      splitB, 2,
                      trmm, 0);
    }
  }
  else {
    if (trans == NORMAL) {
      gemm->AddInputs(6,
                      splitA, 7,
                      splitB, 2,
                      trmm, 0);
    }
    else {
      gemm->AddInputs(6,
                      splitA, 3,
                      splitB, 0,
                      trmm, 0);
    }
  }
  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB = splitB->CreateMatchingCombine(1,
                                                1, gemm, 0);
  
  Poss *loopPoss = new Poss(2, comA, comB);
  RealLoop *loop;
#if DOELEM
  if (layer == DMLAYER)
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  else
    throw;
#elif DOBLIS
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISKC);
#endif
  
  return loop;
}

RealLoop* TrmmLoopLeftVar2(Node *Ain, ConnNum Anum,
                       Node *Bin, ConnNum Bnum,
                       Tri tri, Diag diag, Trans trans,
                       Coef coeff, Type type,
                       Layer layer)
{
  bool rev = (tri == UPPER && trans != NORMAL) || (tri == LOWER && trans == NORMAL);
  
  SplitSingleIter *splitA = new SplitSingleIter(rev ? PARTDIAGBACK : PARTDIAG, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(rev ? PARTUPWARD : PARTDOWN, POSSTUNIN, true);
  splitB->AddInput(Bin, Bnum);
  if (rev)
    splitB->SetUpStats(NOTUP, NOTUP,
                       PARTUP, PARTUP);
  else
    splitB->SetUpStats(PARTUP, PARTUP,
                       NOTUP, NOTUP);
  
  Node *gemm;
  gemm = new Gemm(layer, trans, NORMAL, coeff, COEFONE, type);
  
  if (tri == LOWER) {
    if (trans == NORMAL) {
      gemm->AddInputs(6,
                      splitA, 5,
                      splitB, 1,
                      splitB, 2);
    }
    else {
      gemm->AddInputs(6,
                      splitA, 1,
                      splitB, 1,
                      splitB, 0);
    }
  }
  else {
    if (trans == NORMAL) {
      gemm->AddInputs(6,
                      splitA, 3,
                      splitB, 1,
                      splitB, 0);
    }
    else {
      gemm->AddInputs(6,
                      splitA, 7,
                      splitB, 1,
                      splitB, 2);
    }
  }
  
  Node *trmm;
  trmm = new Trxm(false, layer, LEFT, tri, diag, trans, coeff, type);
  
  trmm->AddInputs(4,
                  splitA, 4,
                  splitB, 1);
  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB;
  if (rev)
    comB = splitB->CreateMatchingCombine(2,
                                         1, trmm, 0,
                                         2, gemm, 0);
  else
    comB = splitB->CreateMatchingCombine(2,
                                         0, gemm, 0,
                                         1, trmm, 0);
  
  Poss *loopPoss = new Poss(2, comA, comB);
  RealLoop *loop;
#if DOELEM
  if (layer == DMLAYER)
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  else
    throw;
#elif DOBLIS
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISKC);
#endif
  
  return loop;
}


RealLoop* TrmmLoopRightVar1(Node *Ain, ConnNum Anum,
                        Node *Bin, ConnNum Bnum,
                        Tri tri, Diag diag, Trans trans,
                        Coef coeff, Type type,
                        Layer layer)
{
  bool rev = (tri == UPPER && trans == NORMAL) || (tri == LOWER && trans != NORMAL);
  
  SplitSingleIter *splitA = new SplitSingleIter(rev ? PARTDIAGBACK: PARTDIAG, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(rev ? PARTLEFT : PARTRIGHT, POSSTUNIN, true);
  splitB->AddInput(Bin, Bnum);
  if (!rev)
    splitB->SetUpStats(FULLUP, NOTUP,
                       FULLUP, NOTUP);
  else
    splitB->SetUpStats(NOTUP, FULLUP,
                       NOTUP, FULLUP);
  
  Node *trmm = new Trxm(false, layer, RIGHT, tri, diag, trans, coeff, type);
  trmm->AddInputs(4,
                  splitA, 4,
                  splitB, 1);
  
  Node *gemm;
  gemm = new Gemm(layer, NORMAL, trans, coeff, COEFONE, type);
  
  if (tri == LOWER) {
    if (trans == NORMAL) {
      gemm->AddInputs(6,
                      splitB, 2,
                      splitA, 5,
                      trmm, 0);
    }
    else {
      gemm->AddInputs(6,
                      splitB, 0,
                      splitA, 1,
                      trmm, 0);
    }
  }
  else {
    if (trans == NORMAL) {
      gemm->AddInputs(6,
                      splitB, 0,
                      splitA, 3,
                      trmm, 0);
    }
    else {
      gemm->AddInputs(6,
                      splitB, 2,
                      splitA, 7,
                      trmm, 0);
    }
  }
  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB = splitB->CreateMatchingCombine(1,
                                                1, gemm, 0);
  
  Poss *loopPoss = new Poss(2, comA, comB);
  RealLoop *loop;
#if DOELEM
  if (layer == DMLAYER)
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  else
    throw;
#elif DOBLIS
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISNC);
#endif
  
  return loop;
}


RealLoop* TrmmLoopRightVar2(Node *Ain, ConnNum Anum,
                        Node *Bin, ConnNum Bnum,
                        Tri tri, Diag diag, Trans trans,
                        Coef coeff, Type type,
                        Layer layer)
{
  bool rev = (tri == UPPER && trans == NORMAL) || (tri == LOWER && trans != NORMAL);
  
  SplitSingleIter *splitA = new SplitSingleIter(rev ? PARTDIAGBACK : PARTDIAG, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(rev ? PARTLEFT : PARTRIGHT, POSSTUNIN, true);
  splitB->AddInput(Bin, Bnum);
  if (rev)
    splitB->SetUpStats(NOTUP, PARTUP,
                       NOTUP, PARTUP);
  else
    splitB->SetUpStats(PARTUP, NOTUP,
                       PARTUP, NOTUP);
  
  
  Node *gemm;
  gemm = new Gemm(layer, NORMAL, trans, coeff, COEFONE, type);
  
  if (tri == LOWER) {
    if (trans == NORMAL) {
      gemm->AddInputs(6,
                      splitB, 1,
                      splitA, 1,
                      splitB, 0);
    }
    else {
      gemm->AddInputs(6,
                      splitB, 1,
                      splitA, 5,
                      splitB, 2);
    }
  }
  else {
    if (trans == NORMAL) {
      gemm->AddInputs(6,
                      splitB, 1,
                      splitA, 7,
                      splitB, 2);
    }
    else {
      gemm->AddInputs(6,
                      splitB, 1,
                      splitA, 3,
                      splitB, 0);
    }
  }
  
  Node *trmm;
  trmm = new Trxm(false, layer, RIGHT, tri, diag, trans, coeff, type);
  trmm->AddInputs(4,
                  splitA, 4,
                  splitB, 1);
  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB;
  if ((tri == LOWER && trans != NORMAL) || (tri == UPPER && trans == NORMAL))
    comB = splitB->CreateMatchingCombine(2,
                                         1, trmm, 0,
                                         2, gemm, 0);
  else
    comB = splitB->CreateMatchingCombine(2,
                                         0, gemm, 0,
                                         1, trmm, 0);
  
  Poss *loopPoss = new Poss(2, comA, comB);
  RealLoop *loop;
#if DOELEM
  if (layer == DMLAYER)
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  else
    throw;
#elif DOBLIS
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISKC);
#endif
  
  return loop;
}


RealLoop* TrxmLoopRightVar3(Node *Ain, ConnNum Anum,
                        Node *Bin, ConnNum Bnum,
                        bool invert,
                        Tri tri, Diag diag, Trans trans,
                        Coef coeff, Type type,
                        Layer layer)
{
  LoopTunnel *Atun = new LoopTunnel(POSSTUNIN);
  Atun->AddInput(Ain, Anum);
  Atun->SetAllStats(FULLUP);
  Atun->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(PARTDOWN, POSSTUNIN, true);
  splitB->AddInput(Bin, Bnum);
  splitB->SetUpStats(FULLUP, FULLUP,
                     NOTUP, NOTUP);
  splitB->SetIndepIters();
  
  Node *trxm;
  trxm = new Trxm(invert, layer, RIGHT, tri, diag, trans, coeff, type);
  trxm->AddInputs(4,
                  Atun, 0,
                  splitB, 1);
  
  
  LoopTunnel *AtunOut = new LoopTunnel(POSSTUNOUT);
  AtunOut->AddInput(Atun, 0);
  AtunOut->AddInput(Atun, 0);
  AtunOut->CopyTunnelInfo(Atun);
  
#if DOELEM
  CombineSingleIter *comB = splitB->CreateMatchingCombine(1,
                                                1, trxm, 0);
  
  Poss *loopPoss = new Poss(2, AtunOut, comB);
  RealLoop *loop;
  if (layer == DMLAYER)
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  else
    throw;
  return loop;
#else
  throw;
#endif
  //    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISKC);
  

}

RealLoop* TrsmLoopLeftVar1(Node *Ain, ConnNum Anum,
                       Node *Bin, ConnNum Bnum,
                       Tri tri, Diag diag, Trans trans,
                       Coef coeff, Layer layer,
                       Type type)
{
  if (tri != LOWER || trans != NORMAL)
    throw;
  
  SplitSingleIter *splitA = new SplitSingleIter(PARTDIAG, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  ScaleNode *scale = new ScaleNode(layer, coeff);
  
  scale->AddInput(Bin, Bnum);
  
  SplitSingleIter *splitB = new SplitSingleIter(PARTDOWN, POSSTUNIN, true);
  splitB->AddInput(scale, 0);
  splitB->SetUpStats(FULLUP, FULLUP,
                     NOTUP, NOTUP);
  
  
  Node *gemm;
  
  gemm = new Gemm(layer, NORMAL, NORMAL, COEFNEGONE, COEFONE, type);
  
  gemm->AddInputs(6,
                  splitA, 1,
                  splitB, 0,
                  splitB, 1);
  
  Node *trsm = NULL;
  trsm = new Trxm(true, layer, LEFT, tri, diag, trans, COEFONE, type);
  
  trsm->AddInputs(4,
                  splitA, 4,
                  gemm, 0);
  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB;
  comB = splitB->CreateMatchingCombine(1,
                                       1, trsm, 0);
  
  Poss *loopPoss = new Poss(2, comA, comB);
  RealLoop *loop;
#if DOELEM
  if (layer == DMLAYER)
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  else
    throw;
#elif DOBLIS
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISKC);
#endif
  
  return loop;
}


RealLoop* TrsmLoopRightVar1(Node *Ain, ConnNum Anum,
                        Node *Bin, ConnNum Bnum,
                        Tri tri, Diag diag, Trans trans,
                        Coef coeff, Type type,
                        Layer layer)
{
  bool rev = (tri == LOWER && trans == NORMAL) || (tri == UPPER && trans != NORMAL);
  
  SplitSingleIter *splitA = new SplitSingleIter(rev ? PARTDIAGBACK: PARTDIAG, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(rev ? PARTLEFT : PARTRIGHT, POSSTUNIN, true);
  splitB->AddInput(Bin, Bnum);
  if (!rev)
    splitB->SetUpStats(FULLUP, NOTUP,
                       FULLUP, NOTUP);
  else
    splitB->SetUpStats(NOTUP, FULLUP,
                       NOTUP, FULLUP);
  
  Node *gemm;
  gemm = new Gemm(layer, NORMAL, trans, COEFNEGONE, COEFONE, type);
  
  if (tri == LOWER) {
    if (trans == NORMAL) {
      gemm->AddInputs(6,
                      splitB, 2,
                      splitA, 5,
                      splitB, 1);
    }
    else {
      gemm->AddInputs(6,
                      splitB, 0,
                      splitA, 1,
                      splitB, 1);
    }
  }
  else {
    if (trans == NORMAL) {
      gemm->AddInputs(6,
                      splitB, 0,
                      splitA, 3,
                      splitB, 1);
    }
    else {
      gemm->AddInputs(6,
                      splitB, 2,
                      splitA, 7,
                      splitB, 1);
    }
  }
  
  
  Node *trmm = new Trxm(true, layer, RIGHT, tri, diag, trans, coeff, type);
  trmm->AddInputs(4,
                  splitA, 4,
                  gemm, 0);
  
  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB = splitB->CreateMatchingCombine(1,
                                                1, trmm, 0);
  
  Poss *loopPoss = new Poss(2, comA, comB);
  RealLoop *loop;
#if DOELEM
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  #else
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISNC);
#endif
  
  return loop;
}


RealLoop* TrsmLoopLeftVar2(Node *Ain, ConnNum Anum,
                       Node *Bin, ConnNum Bnum,
                       Tri tri, Diag diag, Trans trans,
                       Coef coeff, Layer layer,
                       Type type)
{
  bool rev = (tri == UPPER && trans == NORMAL) || (tri == LOWER && trans != NORMAL);
  
  SplitSingleIter *splitA = new SplitSingleIter(rev ? PARTDIAGBACK : PARTDIAG, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  ScaleNode *scale = new ScaleNode(layer, coeff);
  scale->AddInput(Bin, Bnum);
  
  SplitSingleIter *splitB = new SplitSingleIter(rev ? PARTUPWARD : PARTDOWN, POSSTUNIN, true);
  splitB->AddInput(scale, 0);
  if (rev)
    splitB->SetUpStats(PARTUP, PARTUP,
                       FULLUP, FULLUP);
  else
    splitB->SetUpStats(FULLUP, FULLUP,
                       PARTUP, PARTUP);
  
  
  Node *trsm = NULL;
  trsm = new Trxm(true, layer, LEFT, tri, diag, trans, COEFONE, type);
  
  trsm->AddInputs(4,
                  splitA, 4,
                  splitB, 1);
  
  
  Node *gemm;
  
  gemm = new Gemm(layer, trans, NORMAL, COEFNEGONE, COEFONE, type);
  
  
  if (tri == LOWER) {
    if (trans == NORMAL) {
      gemm->AddInputs(6,
                      splitA, 5,
                      trsm, 0,
                      splitB, 2);
    }
    else {
      gemm->AddInputs(6,
                      splitA, 1,
                      trsm, 0,
                      splitB, 0);
    }
  }
  else {
    if (trans == NORMAL) {
      gemm->AddInputs(6,
                      splitA, 3,
                      trsm, 0,
                      splitB, 0);
    }
    else {
      gemm->AddInputs(6,
                      splitA, 7,
                      trsm, 0,
                      splitB, 2);
    }
  }
  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB;
  if (rev)
    comB = splitB->CreateMatchingCombine(2,
                                         0, gemm, 0,
                                         1, trsm, 0);
  else
    comB = splitB->CreateMatchingCombine(2,
                                         1, trsm, 0,
                                         2, gemm, 0);
  
  Poss *loopPoss = new Poss(2, comA, comB);
  RealLoop *loop;
#if DOELEM
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  #else
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISKC);
#endif
  
  loop->SetDimName(DIMK);
  
  return loop;
}

RealLoop* TrxmLoopLeftVar3(Node *Ain, ConnNum Anum,
                       Node *Bin, ConnNum Bnum,
                       bool invert,
                       Tri tri, Diag diag, Trans trans,
                       Coef coeff, Type type, Layer layer)
{
  LoopTunnel *Atun = new LoopTunnel(POSSTUNIN);
  Atun->AddInput(Ain, Anum);
  Atun->SetAllStats(FULLUP);
  Atun->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(PARTRIGHT, POSSTUNIN, true);
  splitB->AddInput(Bin, Bnum);
  splitB->SetUpStats(FULLUP, NOTUP,
                     FULLUP, NOTUP);
  splitB->SetIndepIters();
  
  Node *trxm;
  trxm = new Trxm(invert, layer, LEFT, tri, diag, trans, coeff, type);
  trxm->AddInputs(4,
                  Atun, 0,
                  splitB, 1);
  
  
  LoopTunnel *AtunOut = new LoopTunnel(POSSTUNOUT);
  AtunOut->AddInput(Atun, 0);
  AtunOut->AddInput(Atun, 0);
  AtunOut->CopyTunnelInfo(Atun);
  
  CombineSingleIter *comB = splitB->CreateMatchingCombine(1,
                                                1, trxm, 0);
  
  Poss *loopPoss = new Poss(2, AtunOut, comB);
  RealLoop *loop;
#if DOELEM
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  #else
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISNC);
#endif
  
  loop->SetDimName(DIMN);
  
  return loop;
}

RealLoop* TrsmLoopRightVar2(Node *Ain, ConnNum Anum,
                        Node *Bin, ConnNum Bnum,
                        Tri tri, Diag diag, Trans trans,
                        Coef coeff, Layer layer,
                        Type type)
{
  bool rev = ((tri == UPPER && trans != NORMAL) || (tri == LOWER && trans == NORMAL));
  
  SplitSingleIter *splitA = new SplitSingleIter(rev ? PARTDIAGBACK : PARTDIAG, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  ScaleNode *scale = new ScaleNode(layer, coeff);
  
  scale->AddInput(Bin, Bnum);
  
  SplitSingleIter *splitB = new SplitSingleIter(rev ? PARTLEFT : PARTRIGHT, POSSTUNIN, true);
  splitB->AddInput(scale, 0);
  if (rev)
    splitB->SetUpStats(PARTUP, FULLUP,
                       PARTUP, FULLUP);
  else
    splitB->SetUpStats(FULLUP, PARTUP,
                       FULLUP, PARTUP);
  
  
  Node *trsm = NULL;
  trsm = new Trxm(true, layer, RIGHT, tri, diag, trans, COEFONE, type);
  
  trsm->AddInputs(4,
                  splitA, 4,
                  splitB, 1);
  
  
  Node *gemm;
  
  gemm = new Gemm(layer, NORMAL, trans, COEFNEGONE, COEFONE, type);
  
  gemm->AddInput(trsm, 0);
  
  if (tri == LOWER) {
    if (trans == NORMAL)
      gemm->AddInput(splitA, 1);
    else
      gemm->AddInput(splitA, 5);
  }
  else {
    if (trans == NORMAL)
      gemm->AddInput(splitA, 7);
    else
      gemm->AddInput(splitA, 3);
  }
  
  if (rev)
    gemm->AddInput(splitB, 0);
  else
    gemm->AddInput(splitB, 2);
  
  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB;
  if (rev)
    comB = splitB->CreateMatchingCombine(2,
                                         0, gemm, 0,
                                         1, trsm, 0);
  else
    comB = splitB->CreateMatchingCombine(2,
                                         1, trsm, 0,
                                         2, gemm, 0);
  
  Poss *loopPoss = new Poss(2, comA, comB);
  RealLoop *loop;
#if DOELEM
  if (layer==DMLAYER)
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  else
    throw;
#elif DOBLIS
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISKC);
#endif
  
  return loop;
}


#if DOBLIS
TrxmBP::TrxmBP(bool invert, Layer layer, Side side, Tri tri, Trans trans,
               Coef coeff, Coef beta, Type type)
: TrProps(invert, side, tri, NOTTRIDIAG, trans, coeff, type),
m_beta(beta),
m_comm(CORECOMM)
{
  SetLayer(layer);
}

bool TrxmBP::IsBLISParallelizable() const
{
  return GetLayer() == S3LAYER;
}

bool TrxmBP::IsParallel() const
{
  return m_comm != CORECOMM;
}

bool TrxmBP::RemoveParallelization()
{
  m_comm = CORECOMM;
  return false;
}

void TrxmBP::Parallelize(Comm comm)
{
  if (GetLayer() == S3LAYER)
    m_comm = comm;
  else
    throw;
}

void TrxmBP::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<3,2>::Duplicate(orig, shallow, possMerging);
  TrProps::Duplicate((TrxmBP*)orig);
  TrxmBP *bp = (TrxmBP*)orig;
  m_beta = bp->m_beta;
  m_comm = bp->m_comm;
}

void TrxmBP::FlattenCore(ofstream &out) const
{
  DLAOp<3,2>::FlattenCore(out);
  TrProps::FlattenCore(out);
  WRITE(m_beta);
  WRITE(m_comm);
}

void TrxmBP::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<3,2>::UnflattenCore(in,info);
  TrProps::UnflattenCore(in,info);
  READ(m_beta);
  READ(m_comm);
}



NodeType TrxmBP::GetType() const
{
  return "TrxmBP " + LayerNumToStr(GetLayer()) + " "
  + SideToStr(m_side) + " "
  + TriToStr(m_tri) + " "
  + TransToStr(m_trans)
  + " " + CommToStr(m_comm);
}


void TrxmBP::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3,2>::Prop();
    switch (GetLayer()) {
      case (S3LAYER):
      {
        const Sizes *ms0 = GetInputM(0);
        const Sizes *ns0 = GetInputN(0);
	//        const Sizes *ms1 = GetInputM(1);
        const Sizes *ns1 = GetInputN(1);
        //	const Sizes *ms2 = GetInputM(2);
        //	const Sizes *ns2 = GetInputN(2);
	m_cost = GAMMA * ms0->SumProds111(*ns0,*ns1) / NumCoresInComm(m_comm);
	if (NumCoresInComm(m_comm) > 1) {
	    Size numAElems = ms0->SumProds11(*ns0);
	    m_cost += AdditionalCostForBringingIntoL2(this, 0, numAElems, m_comm);
        }
        break;
      }
      default:
        throw;
    }
  }
}

void TrxmBP::PrintCode(IndStream &out)
{
  out.Indent();
  if (m_invert)
    *out << "bli_trsm_";
  else
    *out << "bli_trmm_";
  if (m_side == RIGHT)
    *out << "r";
  else
    *out << "l";
  if (m_tri == LOWER)
    *out << "l";
  else if (m_tri == UPPER)
    *out << "u";
  else
    throw;
  *out << "_ker_var2";
  if (m_comm == CORECOMM)
    *out << "( ";
  else
    *out << "_par( ";
  out << m_coeff ;
  *out << ", &"
       << GetInputName(0).str() << ", &" << GetInputName(1).str() << ", \n"
       << out.Tabs(2);
  out << m_beta;
  *out << ", &" << GetInputName(2).str() << ", (tr"
       << (m_invert ? "s" : "m") << "m_t*)NULL";
  if (m_comm != CORECOMM) {
#if DOSM||DOSQM
    *out << ", ";
    if (m_invert)
      *out << CommToStr(m_comm);
    else
      *out << CommToStr(GetSubComm(m_comm));
#endif
  }
  *out <<  ");\n";
}

bool BLISTrxmLoopExp::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Trxm::GetClass())
    return ((Trxm*)node)->GetLayer() == m_fromLayer;
  return false;
}

void BLISTrxmLoopExp::Apply(Node *node) const
{
  Trxm *trxm = (Trxm*)node;
  bool isTrsm = trxm->m_invert;
  
  Tri tri = trxm->m_tri;
  bool isLeft = trxm->m_side == LEFT;
  
  PartDir lhsDir;
  PartDir outputDir;
  if (isTrsm) {
    if (isLeft) {
      if (trxm->m_trans != NORMAL && trxm->m_trans != CONJ) {
        if (tri == UPPER) {
          lhsDir = PARTRIGHT;
          outputDir = PARTDOWN;
        }
        else {
          lhsDir = PARTLEFT;
          outputDir = PARTUPWARD;
        }
      }
      else {
        if (tri == LOWER) {
          lhsDir = PARTDOWN;
          outputDir = PARTDOWN;
        }
        else {
          lhsDir = PARTUPWARD;
          outputDir = PARTUPWARD;
        }
      }
    }
    else {
      if (trxm->m_trans != NORMAL && trxm->m_trans != CONJ) {
        if (tri == UPPER) {
          lhsDir = PARTUPWARD;
          outputDir = PARTUPWARD;
        }
        else {
          lhsDir = PARTDOWN;
          outputDir = PARTDOWN;
        }
      }
      else {
        if (tri == LOWER) {
          lhsDir = PARTUPWARD;
          outputDir = PARTUPWARD;
        }
        else {
          lhsDir = PARTDOWN;
          outputDir = PARTDOWN;
        }
      }
    }
  }
  else {
    if (isLeft) {
      if (trxm->m_trans != NORMAL && trxm->m_trans != CONJ)
        lhsDir = PARTRIGHT;
      else
        lhsDir = PARTDOWN;
      outputDir = PARTDOWN;
    }
    else {
      lhsDir = PARTDOWN;
      outputDir = PARTDOWN;
    }
  }
  
  SplitSingleIter *splitLHS = new SplitSingleIter(lhsDir, POSSTUNIN, true);
  if (isLeft) {
    splitLHS->AddInput(node->Input(0), node->InputConnNum(0));
  }
  else {
    splitLHS->AddInput(node->Input(1), node->InputConnNum(1));
  }
  splitLHS->SetAllStats(FULLUP);
  splitLHS->SetIndepIters();
  
  SplitSingleIter *splitOutput = new SplitSingleIter(outputDir, POSSTUNIN);
  splitOutput->AddInput(trxm->Input(1), trxm->InputConnNum(1));
  splitOutput->SetUpStats(FULLUP, FULLUP,
                          NOTUP, NOTUP);
  if (!isTrsm || !isLeft)
  splitOutput->SetIndepIters();
  //else
  //Not independent
  
  Node *rhsSrc;
  ConnNum rhsSrcNum;
  if (isLeft) {
    rhsSrc = node->Input(1);
    rhsSrcNum = node->InputConnNum(1);
  }
  else {
    rhsSrc = node->Input(0);
    rhsSrcNum = node->InputConnNum(0);
    if (trxm->m_trans != NORMAL) {
      rhsSrc = AddTranspose(trxm->m_trans, true, rhsSrc, rhsSrcNum, true);
      rhsSrcNum = 0;
    }
  }
  
  PackBuff *rhsBuff = NULL;
  Pack *rhsPack = NULL;
  if (isLeft) {
    rhsBuff = new PackBuff(node->Input(1)->GetName(node->InputConnNum(1)).m_name,
                           PACKCOLPANS, PACKBPANEL, NOTTRI, NOTTRIDIAG, GEN,
                           false, false, false, false,
                           USEMRSIZE, USENRSIZE);
    rhsPack = new Pack(PACKCOLPANS, 2, false, false, false, false, false);
  }
  else {
    //check invertDiag
    rhsBuff = new PackBuff(node->Input(1)->GetName(node->InputConnNum(1)).m_name,
                           PACKCOLPANS, PACKBPANEL, trxm->m_tri, trxm->m_diag, TRI,
                           true, isTrsm, false, trxm->m_invert,
                           isTrsm ? USEMRSIZE : USENRSIZE,
                           isTrsm ? USEMRSIZE : USENRSIZE);
    rhsPack = new Pack(PACKCOLPANS, 3, false,
                       true, isTrsm, false, isTrsm);
  }
  rhsBuff->AddInput(rhsSrc, rhsSrcNum);
  rhsPack->AddInput(rhsSrc, rhsSrcNum);
  rhsPack->AddInput(rhsBuff, 0);
  
  node->m_poss->AddNode(rhsBuff);
  node->m_poss->AddNode(rhsPack);
  
  LoopTunnel *packTun = new LoopTunnel(POSSTUNIN);
  packTun->AddInput(rhsPack);
  if (isLeft) {
    //kind of lying here....
    packTun->SetAllStats(isTrsm ? PARTUP : FULLUP);
  }
  else
    packTun->SetAllStats(FULLUP);
  if (!isTrsm || !isLeft)
    packTun->SetIndepIters();
  
  Node *lhsSrc = splitLHS;
  ConnNum lhsSrcNum = 1;
  
  if (isLeft && trxm->m_trans != NORMAL) {
    lhsSrc = AddTranspose(trxm->m_trans, true, lhsSrc, lhsSrcNum, false);
    lhsSrcNum = 0;
  }
  
  PackBuff *lhsBuff;
  Pack *lhsPack;
  
  if (isLeft) {
    lhsBuff = new PackBuff(node->Input(0)->GetName(node->InputConnNum(0)).m_name,
                           PACKROWPANS, PACKABLOCK, tri, trxm->m_diag, TRI,
                           true, isTrsm, isTrsm, false,
                           USEMRSIZE, USEMRSIZE);
    lhsPack = new Pack(PACKROWPANS, 3, false, true, isTrsm, isTrsm, false);
  }
  else {
    lhsBuff = new PackBuff(node->Input(0)->GetName(node->InputConnNum(0)).m_name,
                           PACKROWPANS, PACKABLOCK, NOTTRI, NOTTRIDIAG, GEN,
                           false, false, false, false,
                           isTrsm ? USENRSIZE : USEMRSIZE,
                           isTrsm ? USEMRSIZE : USENRSIZE );
    lhsPack = new Pack(PACKROWPANS, 2, false, false, false, false, false);
  }
  lhsBuff->AddInput(lhsSrc, lhsSrcNum);
  lhsPack->AddInput(lhsSrc, lhsSrcNum);
  lhsPack->AddInput(lhsBuff, 0);
  
  if (trxm->m_trans != NORMAL && trxm->m_trans != CONJ)
    tri = SwapTri(tri);
  
  TrxmBP *trbp = new TrxmBP(isTrsm, m_toLayer, trxm->m_side, tri,
                            trxm->m_trans != CONJ ? trxm->m_trans : NORMAL,
                            trxm->m_coeff, COEFZERO, trxm->m_type);
  trbp->AddInput(lhsPack, 0);
  trbp->AddInput(packTun, 0);
  trbp->AddInput(splitOutput, 1);
  
  CombineSingleIter *comA = splitLHS->CreateMatchingCombine(0);
  
  CombineSingleIter *comB = splitOutput->CreateMatchingCombine(1,
                                                     1, trbp, isTrsm ? 1 : 0);
  
  LoopTunnel *packTunOut = new LoopTunnel(POSSTUNOUT);
  if (isTrsm)
    packTunOut->AddInput(trbp, 0);
  else
    packTunOut->AddInput(packTun, 0);
  packTunOut->AddInput(packTun, 0);
  packTunOut->CopyTunnelInfo(packTun);
  
  Poss *loopPoss = new Poss(3, comA, comB, packTunOut);
  RealLoop *loop = new RealLoop(BLISLOOP, loopPoss, USEBLISMC);
  
  loop->SetDimName(DIMM);
  
  node->m_poss->AddPSet(loop);
  
  node->RedirectChildren(loop->OutTun(1), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}
#endif //DOBLIS

RealLoop* Trmm3LoopLeftVar2(Node *Ain, ConnNum Anum,
                        Node *Bin, ConnNum Bnum,
                        Node *Cin, ConnNum Cnum,
                        Tri tri, Diag diag, Trans trans,
                        Coef coeff, Coef beta, Type type,
                        Layer layer)
{
  bool rev = (tri == UPPER && trans != NORMAL) || (tri == LOWER && trans == NORMAL);
  
  SplitSingleIter *splitA = new SplitSingleIter(rev ? PARTDIAGBACK : PARTDIAG, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(rev ? PARTUPWARD : PARTDOWN, POSSTUNIN);
  splitB->AddInput(Bin, Bnum);
  splitB->SetAllStats(FULLUP);
  splitB->SetIndepIters();
  
  SplitSingleIter *splitC = new SplitSingleIter(rev ? PARTUPWARD : PARTDOWN, POSSTUNIN, true);
  splitC->AddInput(Cin, Cnum);
  if (rev)
    splitC->SetUpStats(NOTUP, NOTUP,
                       PARTUP, PARTUP);
  else
    splitC->SetUpStats(PARTUP, PARTUP,
                       NOTUP, NOTUP);
  
  
  Node *gemm;
  gemm = new Gemm(layer, trans, NORMAL, coeff, COEFONE, type);
  
  if (tri == LOWER) {
    if (trans == NORMAL) {
      gemm->AddInputs(6,
                      splitA, 5,
                      splitB, 1,
                      splitC, 2);
    }
    else {
      gemm->AddInputs(6,
                      splitA, 1,
                      splitB, 1,
                      splitC, 0);
    }
  }
  else {
    if (trans == NORMAL) {
      gemm->AddInputs(6,
                      splitA, 3,
                      splitB, 1,
                      splitC, 0);
    }
    else {
      gemm->AddInputs(6,
                      splitA, 7,
                      splitB, 1,
                      splitC, 2);
    }
  }
  
  Node *trmm3;
  trmm3 = new Trmm3(layer, LEFT, tri, diag, trans, coeff, beta, type);
  
  trmm3->AddInputs(6,
                   splitA, 4,
                   splitB, 1,
                   splitC, 1);
  
#if DOELEM
    throw;
    return NULL;
#else
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB = splitB->CreateMatchingCombine(0);
  
  CombineSingleIter *comC;
  if (rev)
    comC = splitC->CreateMatchingCombine(2,
                                         1, trmm3, 0,
                                         2, gemm, 0);
  else
    comC = splitC->CreateMatchingCombine(2,
                                         0, gemm, 0,
                                         1, trmm3, 0);
  

  RealLoop *loop;

  Poss *loopPoss = new Poss(3, comA, comB, comC);
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISKC);
  
  return loop;
#endif

}

RealLoop* Trmm3LoopLeftVar3(Node *Ain, ConnNum Anum,
                        Node *Bin, ConnNum Bnum,
                        Node *Cin, ConnNum Cnum,
                        Tri tri, Diag diag, Trans trans,
                        Coef coeff, Coef beta, Type type, Layer layer)
{
  LoopTunnel *Atun = new LoopTunnel(POSSTUNIN);
  Atun->AddInput(Ain, Anum);
  Atun->SetAllStats(FULLUP);
  Atun->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(PARTRIGHT, POSSTUNIN, true);
  splitB->AddInput(Bin, Bnum);
  splitB->SetUpStats(FULLUP, FULLUP,
                     FULLUP, FULLUP);
  splitB->SetIndepIters();
  
  SplitSingleIter *splitC = new SplitSingleIter(PARTRIGHT, POSSTUNIN, true);
  splitC->AddInput(Bin, Bnum);
  splitC->SetUpStats(FULLUP, NOTUP,
                     FULLUP, NOTUP);
  splitC->SetIndepIters();
  
  Node *trmm3;
  trmm3 = new Trmm3(layer, LEFT, tri, diag, trans, coeff, beta, type);
  trmm3->AddInputs(6,
                   Atun, 0,
                   splitB, 1,
                   splitC, 1);
  
  
  LoopTunnel *AtunOut = new LoopTunnel(POSSTUNOUT);
  AtunOut->AddInput(Atun, 0);
  AtunOut->AddInput(Atun, 0);
  AtunOut->CopyTunnelInfo(Atun);
  
  CombineSingleIter *comB = splitB->CreateMatchingCombine(0);
  
  CombineSingleIter *comC = splitC->CreateMatchingCombine(1,
                                                1, trmm3, 0);
  
  Poss *loopPoss = new Poss(3, AtunOut, comB, comC);
  RealLoop *loop;
#if DOELEM
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  #else
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISNC);
#endif
  
  return loop;
}

#if DOBLIS
template<>
string TrxmRightToLeft<Trxm>::GetType() const
{
  return "Trxm right to left";
}

template<class TrxmType>
bool TrxmRightToLeft<TrxmType>::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == TrxmType::GetClass()) {
    const TrxmType *trxm = (TrxmType*)node;
    if (trxm->GetLayer() != m_layer)
      return false;
    return trxm->m_side == RIGHT;
  }
  return false;
}

template<class TrxmType>
void TrxmRightToLeft<TrxmType>::Apply(Node *node) const
{
  TrxmType *trxm = (TrxmType*)node;
  trxm->m_side = LEFT;
  if (trxm->m_trans == NORMAL)
    trxm->m_trans = TRANS;
  else if (trxm->m_trans == TRANS)
    trxm->m_trans = NORMAL;
  else if (trxm->m_trans == CONJTRANS) {
    trxm->m_trans = CONJ;
  }
  //else don't do anything
  
  InsertTranspose(TRANS, false, trxm, 1, true);
  
  Transpose *newTrans = new Transpose(TRANS, false);
  node->m_poss->AddNode(newTrans);
  trxm->RedirectAllChildren(newTrans);
  newTrans->AddInput(trxm, 0);
}

string Trmm3RightToLeft::GetType() const
{
  return "Trmm3 right to left";
}

bool Trmm3RightToLeft::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Trmm3::GetClass()) {
    const Trmm3 *trmm3 = (Trmm3*)node;
    if (trmm3->GetLayer() != m_layer)
      return false;
    return trmm3->m_side == RIGHT;
  }
  return false;
}


void Trmm3RightToLeft::Apply(Node *node) const
{
  Trmm3 *trmm3 = (Trmm3*)node;
  trmm3->m_side = LEFT;
  if (trmm3->m_trans == NORMAL)
    trmm3->m_trans = TRANS;
  else if (trmm3->m_trans == TRANS)
    trmm3->m_trans = NORMAL;
  else if (trmm3->m_trans == CONJTRANS) {
    trmm3->m_trans = CONJ;
  }
  //else don't do anything
  
  InsertTranspose(TRANS, true, trmm3, 1, true);
  InsertTranspose(TRANS, false, trmm3, 2, true);
  
  Transpose *newTrans = new Transpose(TRANS, false);
  node->m_poss->AddNode(newTrans);
  trmm3->RedirectAllChildren(newTrans);
  newTrans->AddInput(trmm3, 0);
}
#endif

template<>
string TrxmLowerLayer<Trxm>::GetType() const
{
  return "Trxm lower layer"  + LayerNumToStr(m_fromLayer)
  + " to " + LayerNumToStr(m_toLayer);
}

template<>
string TrxmLowerLayer<Trmm3>::GetType() const
{
  return "Trmm3 lower layer"  + LayerNumToStr(m_fromLayer)
  + " to " + LayerNumToStr(m_toLayer);
}

template<class TrxmType>
bool TrxmLowerLayer<TrxmType>::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == TrxmType::GetClass()) {
    const TrxmType *trxm = (TrxmType*)node;
    if (trxm->GetLayer() != m_fromLayer)
      return false;
    if (trxm->m_side == LEFT) {
      if (m_dim == DIMK)
        return (*(trxm->GetInputM(0)) <= m_bs
                && *(trxm->GetInputN(0)) <= m_bs);
      else if (m_dim == DIMN)
        return (*(trxm->GetInputN(1)) <= m_bs);
      else
        throw;
    }
    else {
      if (m_dim == DIMK)
        return (*(trxm->GetInputN(1)) <= m_bs);
      else if (m_dim == DIMN) {
        if (trxm->m_trans == NORMAL)
          return (*(trxm->GetInputN(0)) <= m_bs);
        else
          return (*(trxm->GetInputM(0)) <= m_bs);
      }
      else
        throw;
    }
    
  }
  return false;
  
}

template<class TrxmType>
void TrxmLowerLayer<TrxmType>::Apply(Node *node) const
{
  TrxmType *trxm = (TrxmType*)node;
  trxm->SetLayer(m_toLayer);
}


#if DOBLIS
bool BLISTrmm3LoopExp::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Trmm3::GetClass())
    return ((Trmm3*)node)->GetLayer() == m_fromLayer;
  return false;
}

void BLISTrmm3LoopExp::Apply(Node *node) const
{
  Trmm3 *trmm3 = (Trmm3*)node;
  Tri tri = trmm3->m_tri;
  
  Trans transA, transB;
  if (trmm3->m_side == LEFT) {
    transA = trmm3->m_trans;
    transB = NORMAL;
  }
  else if (trmm3->m_side == RIGHT) {
    throw;
  }
  else
    throw;
  
  PartDir aDir;
  PartDir bDir;
  if (transA != NORMAL)
    aDir = PARTRIGHT;
  else
    aDir = PARTDOWN;
  if (transB != NORMAL)
    bDir = PARTRIGHT;
  else
    bDir = PARTDOWN;
  
  SplitSingleIter *splitA = new SplitSingleIter(aDir, POSSTUNIN, true);
  splitA->AddInput(trmm3->Input(0), trmm3->InputConnNum(0));
  splitA->SetAllStats(FULLUP);
  //BAM IndepIters
  
  SplitSingleIter *splitB = new SplitSingleIter(bDir, POSSTUNIN);
  splitB->AddInput(trmm3->Input(1), trmm3->InputConnNum(1));
  splitB->SetAllStats(FULLUP);
  
  SplitSingleIter *splitC = new SplitSingleIter(PARTDOWN, POSSTUNIN);
  splitC->AddInput(trmm3->Input(2), trmm3->InputConnNum(2));
  splitC->SetUpStats(FULLUP, FULLUP,
                     NOTUP, NOTUP);
  
  Node *bSrc = node->Input(1);
  ConnNum bSrcNum = node->InputConnNum(1);
  
  if (transB != NORMAL) {
    bSrc = AddTranspose(transB, false, bSrc, bSrcNum, true);
    bSrcNum = 0;
  }
  
  PackBuff *bBuff = new PackBuff(node->Input(1)->GetName(node->InputConnNum(1)).m_name,
                                 PACKCOLPANS, PACKBPANEL, NOTTRI, NOTTRIDIAG, GEN,
                                 false, false, false, false,
                                 USEMRSIZE, USENRSIZE);
  Pack *bPack = new Pack(PACKCOLPANS, 2, false, false, false, false, false);
  bBuff->AddInput(bSrc, bSrcNum);
  bPack->AddInput(bSrc, bSrcNum);
  bPack->AddInput(bBuff, 0);
  
  node->m_poss->AddNode(bBuff);
  node->m_poss->AddNode(bPack);
  
  LoopTunnel *packTun = new LoopTunnel(POSSTUNIN);
  packTun->AddInput(bPack);
  packTun->SetAllStats(FULLUP);
  
  Node *aSrc = splitA;
  ConnNum aSrcNum = 1;
  
  if (transA != NORMAL) {
    aSrc = AddTranspose(transA, true, aSrc, aSrcNum, false);
    aSrcNum = 0;
  }
  
  PackBuff *aBuff;
  Pack *aPack;
  
  aBuff = new PackBuff(node->Input(0)->GetName(node->InputConnNum(0)).m_name,
                       PACKROWPANS, PACKABLOCK, tri, trmm3->m_diag, TRI,
                       true, false, false, false,
                       USEMRSIZE, USEMRSIZE);
  aPack = new Pack(PACKROWPANS, 3, false, true, false, false, false);
  
  aBuff->AddInput(aSrc, aSrcNum);
  aPack->AddInput(aSrc, aSrcNum);
  aPack->AddInput(aBuff, 0);
  
  if (transA != NORMAL)
    tri = SwapTri(tri);
  
  TrxmBP *trbp = new TrxmBP(false,
                            m_toLayer, trmm3->m_side, tri,
                            trmm3->m_trans != CONJ ? trmm3->m_trans : NORMAL,
                            trmm3->m_coeff, trmm3->m_beta, trmm3->m_type);
  trbp->AddInput(aPack, 0);
  trbp->AddInput(packTun, 0);
  trbp->AddInput(splitC, 1);
  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB = splitB->CreateMatchingCombine(0);
  
  CombineSingleIter *comC = splitC->CreateMatchingCombine(1,
                                                1, trbp, 0);
  
  LoopTunnel *packTunOut = new LoopTunnel(POSSTUNOUT);
  packTunOut->AddInput(packTun, 0);
  packTunOut->AddInput(packTun, 0);
  packTunOut->CopyTunnelInfo(packTun);
  
  Poss *loopPoss = new Poss(4, comA, comB, packTunOut, comC);
  RealLoop *loop = new RealLoop(BLISLOOP, loopPoss, USEBLISMC);
  
  loop->SetDimName(DIMM);
  
  node->m_poss->AddPSet(loop);
  
  node->RedirectChildren(loop->OutTun(3), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}
#endif

#if DOBLIS
bool TrmmAxpytoTrxm3::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Trxm::GetClass())
    throw;
  if (((Trxm*)node)->m_invert)
    return false;
  if (((DLANode*)node)->GetLayer() != m_layer)
    return false;
  if (node->m_children.size() != 1)
    return false;
  DLANode *child = (DLANode*)(node->Child(0));
  if (child->GetNodeClass() != Axpy::GetClass())
    return false;
  return child->GetLayer() == m_layer;
}

void TrmmAxpytoTrxm3::Apply(Node *node) const
{
  Trxm *trmm = (Trxm*)node;
  Axpy *axpy = (Axpy*)(trmm->Child(0));
  Trmm3 *trmm3 = new Trmm3(trmm->m_layer,
                           trmm->m_side, trmm->m_tri, trmm->m_diag,
                           trmm->m_trans, trmm->m_coeff*axpy->m_coeff, COEFONE,
                           trmm->m_type);
  
  trmm3->AddInputs(6,
                   trmm->Input(0), trmm->InputConnNum(0),
                   trmm->Input(1), trmm->InputConnNum(1),
                   axpy->Input(1), axpy->InputConnNum(1));
  
  
  node->m_poss->AddNode(trmm3);
  axpy->RedirectChildren(trmm3);
  node->m_poss->DeleteChildAndCleanUp(axpy);
}


bool CopyTrmmtoTrxm3::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Trxm::GetClass())
    throw;
  if (((Trxm*)node)->m_invert)
    return false;
  if (((DLANode*)node)->GetLayer() != m_layer)
    return false;
  DLANode *parent = (DLANode*)(node->Input(1));
  if (parent->GetNodeClass() != Copy::GetClass())
    return false;
  return parent->GetLayer() == m_layer;
}

void CopyTrmmtoTrxm3::Apply(Node *node) const
{
  Trxm *trmm = (Trxm*)node;
  Copy *copy = (Copy*)(trmm->Input(1));
  Trmm3 *trmm3 = new Trmm3(trmm->m_layer,
                           trmm->m_side, trmm->m_tri, trmm->m_diag,
                           trmm->m_trans, trmm->m_coeff, COEFZERO,
                           trmm->m_type);
  
  trmm3->AddInputs(6,
                   trmm->Input(0), trmm->InputConnNum(0),
                   copy->Input(0), copy->InputConnNum(0),
                   copy->Input(1), copy->InputConnNum(1));
  
  node->m_poss->AddNode(trmm3);
  trmm->RedirectChildren(trmm3);
  node->m_poss->DeleteChildAndCleanUp(trmm);
}

template class TrxmRightToLeft<Trxm>;
#endif


template class TrxmLowerLayer<Trxm>;

template class TrxmLowerLayer<Trmm3>;

#endif
