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

#include "blas.h"
#include "her2k.h"
#include "elemRedist.h"
#include "stdlib.h"
#include "helperNodes.h"
#include "loopSupport.h"

Her2kProps::Her2kProps(Tri tri, Trans trans, Coef alpha, Coef beta, Type type)
:m_tri(tri), m_trans(trans), m_alpha(alpha), m_beta(beta), m_type(type)
{
  switch(m_trans) {
    case(NORMAL):
    case(TRANS):
    case(CONJTRANS):
      break;
    default:
      throw;
  }
  
}

void Her2kProps::DuplicateProps(const Her2kProps *her2k)
{
  m_tri = her2k->m_tri;
  m_trans = her2k->m_trans;
  m_alpha = her2k->m_alpha;
  m_beta = her2k->m_beta;
  m_type = her2k->m_type;
  
  switch(m_trans) {
    case(NORMAL):
    case(TRANS):
    case(CONJTRANS):
      break;
    default:
      throw;
  }
}

void Her2kProps::FlattenCore(ofstream &out) const
{
  WRITE(m_tri);
  WRITE(m_trans);
  WRITE(m_alpha);
  WRITE(m_beta);
  WRITE(m_type);
}


void Her2kProps::UnflattenCore(ifstream &in, SaveInfo &info)
{
  READ(m_tri);
  READ(m_trans);
  READ(m_alpha);
  READ(m_beta);
  READ(m_type);
}

Her2k::Her2k(Layer layer, Tri tri, Trans trans, Coef alpha, Coef beta, Type type)
: Her2kProps(tri, trans, alpha, beta, type)
{
  SetLayer(layer);
}

void Her2k::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3,1>::Prop();

#if DOELEM
    if (GetLayer() == ABSLAYER) {
      if (InputDataType(0).m_dist != D_MC_MR)
	throw;
      if (InputDataType(1).m_dist != D_MC_MR)
	throw;
    }
    else if (GetLayer() == DMLAYER) {
      if (InputDataType(2).m_dist != D_MC_MR)
	cout << "input not D_MC_MR";
    
    }
    else if (GetLayer() == SMLAYER) {
      DistType t0, t1;
      t0 = InputDataType(0).m_dist;
      t1 = InputDataType(1).m_dist;
    
      if (m_trans == NORMAL) {
	if (t0 == D_STAR_VR) {
	  if (t1 != D_VR_STAR)
	    throw;
	}
	else if (t0 == D_STAR_VC) {
	  if (t1 != D_VC_STAR)
	    throw;
	}
	else
	  throw;
      }
      else {
	if (t0 == D_VR_STAR) {
	  if (t1 != D_STAR_VR)
	    throw;
	}
	else if (t0 == D_VC_STAR) {
	  if (t1 != D_STAR_VC)
	    throw;
	}
      }
    
      if (InputDataType(2).m_dist != D_STAR_STAR)
	throw;
    }
#endif
  
  
    if (m_type == REAL) {
      if (m_trans != NORMAL && m_trans != TRANS)
	throw;
    }
    else {
      if (m_trans != NORMAL && m_trans != CONJTRANS)
	throw;
    }
    
    if (GetLayer() == ABSLAYER)
      m_cost = 0;
#if DOELEM
    else if (GetLayer() == DMLAYER)
      m_cost = 0;
    else if (GetLayer() == SMLAYER) {
      if (m_trans == NORMAL)
        m_cost = GetCost(SMLAYER, InputLocalM(2),InputLocalN(0));
      else
        m_cost = GetCost(SMLAYER, InputLocalM(2),InputLocalM(0));
    }
#endif
    else
      m_cost = ZERO;
  }
}



void Her2k::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  const Her2k *her2k = (Her2k*)orig;
  DLAOp<3,1>::Duplicate(orig, shallow, possMerging);
  Her2kProps::DuplicateProps(her2k);
}

NodeType Her2k::GetType() const
{
  string str = "Her2k " + LayerNumToStr(GetLayer()) + TransToStr(m_trans) + TriToStr(m_tri);
  return str;
}


Cost Her2k::GetCost(Layer layer, const Sizes *m, const Sizes *k)
{
#if DOELEM
  if (layer == SMLAYER)
    return 2 * GAMMA * m->SumProds21(*k);
  else
#endif
    throw;
}

void Her2k::FlattenCore(ofstream &out) const
{
  DLAOp<3,1>::FlattenCore(out);
  Her2kProps::FlattenCore(out);
}
void Her2k::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<3,1>::UnflattenCore(in, info);
  Her2kProps::UnflattenCore(in,info);
}

Phase Her2k::MaxPhase() const
{
#if DODPPHASE
  if (GetLayer() == ABSLAYER || GetLayer() == DMLAYER)
    return DPPHASE;
  else if (GetLayer() == SMLAYER)
    return NUMPHASES;
  else
    throw;
#elif DOBLIS
  switch (GetLayer()) {
    case (ABSLAYER):
      return SR1PHASE;
    case (S1LAYER):
      return SR2PHASE;
    case (S2LAYER):
      return SR3PHASE;
    default:
      throw;
  }
#else
  throw;
#endif
}

#if DOELEM
bool Her2k::ShouldCullDP() const
{
  if (GetLayer() == ABSLAYER || GetLayer() == DMLAYER)
    return true;
  else
    return false;
}
#endif


void Her2k::PrintCode(IndStream &out)
{
  out.Indent();
  if (GetLayer() == ABSLAYER) {
    if (m_type == REAL)
      *out << "AbsSyr2k( " ;
    else
      *out << "AbsHer2k( " ;
  }
#if DOELEM
  else if (GetLayer() == DMLAYER) {
    if (m_type == REAL)
      *out << "DistSyr2k( " ;
    else
      *out << "DistHer2k( " ;
  }
  else if (GetLayer() == SMLAYER) {
    if (m_type == REAL)
      *out << "Syr2k( " ;
    else
      *out << "Her2k( " ;
  }
#endif
  *out << TriToStr(m_tri) << ", " << TransToStr(m_trans)
  << ", \n" << out.Tabs(1);
  out << m_alpha;
  *out << ", "
  << GetInputName(0).str() << ".LocalMatrix(), "
  << GetInputName(1).str() << ".LocalMatrix(), \n" << out.Tabs(1) << "";
  out << m_beta;
  *out << ", "
  << GetInputName(2).str() << ".LocalMatrix());\n";
}

string Her2kLoopExp::GetType() const
{
  switch(m_var) {
    case(1):
      return "Her2k Loop Exp - var 1";
    case(2):
      return "Her2k Loop Exp - var 2";
    case(3):
      return "Her2k Loop Exp - var 3";
    case(4):
      return "Her2k Loop Exp - var 4";
    case(9):
      return "Her2k Loop Exp - var 9";
    default:
      throw;
  }
}

bool Her2kLoopExp::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Her2k::GetClass()
      && ((Her2k*)node)->GetLayer() == m_fromLayer)
    return true;
  return false;
}

void Her2kLoopExp::Apply(Node *node) const
{
  Her2k *her2k = (Her2k*)node;
  RealLoop *loop;
  
  NodeConn *connA, *connB, *connC;
  connA = her2k->m_inputs[0];
  connB = her2k->m_inputs[1];
  connC = her2k->m_inputs[2];
  
  switch(m_var) {
    case(1):
      loop = Her2kLoopVar1(connA->m_n, connA->m_num,
                           connB->m_n, connB->m_num,
                           connC->m_n, connC->m_num,
                           her2k->m_tri,
                           her2k->m_trans,
                           her2k->m_alpha, her2k->m_beta, her2k->m_type, m_toLayer);
      break;
      
    case(2):
      loop = Her2kLoopVar2(connA->m_n, connA->m_num,
                           connB->m_n, connB->m_num,
                           connC->m_n, connC->m_num,
                           her2k->m_tri,
                           her2k->m_trans,
                           her2k->m_alpha, her2k->m_beta, her2k->m_type, m_toLayer);
      break;
      
    case(3):
      loop = Her2kLoopVar3(connA->m_n, connA->m_num,
                           connB->m_n, connB->m_num,
                           connC->m_n, connC->m_num,
                           her2k->m_tri,
                           her2k->m_trans,
                           her2k->m_alpha, her2k->m_beta, her2k->m_type, m_toLayer);
      break;
      
    case(4):
      loop = Her2kLoopVar4(connA->m_n, connA->m_num,
                           connB->m_n, connB->m_num,
                           connC->m_n, connC->m_num,
                           her2k->m_tri,
                           her2k->m_trans,
                           her2k->m_alpha, her2k->m_beta, her2k->m_type, m_toLayer);
      break;
      
      
    case(9):
      loop = Her2kLoopVar9(connA->m_n, connA->m_num,
                           connB->m_n, connB->m_num,
                           connC->m_n, connC->m_num,
                           her2k->m_tri,
                           her2k->m_trans,
                           her2k->m_alpha, her2k->m_beta, her2k->m_type, m_toLayer);
      break;
      
    default:
      throw;
  }
  
  node->m_poss->AddPSet(loop);
  
  node->RedirectChildren(loop->OutTun(2),0);
  node->m_poss->DeleteChildAndCleanUp(node);
}


Tri2k::Tri2k(Layer layer, Tri tri, Trans trans, Coef alpha, Coef beta, Type type)
: Her2kProps(tri,trans,alpha,beta,type)
{
  SetLayer(layer);
}

NodeType Tri2k::GetType() const
{
  string str = "Tri2k " + LayerNumToStr(GetLayer()) + TriToStr(m_tri);
  return str;
}

Phase Tri2k::MaxPhase() const
{
#if DODPPHASE
  if (GetLayer() == ABSLAYER || GetLayer() == DMLAYER)
    return DPPHASE;
  else if (GetLayer() == SMLAYER)
    return NUMPHASES;
  else
    throw;
#elif DOBLIS
  switch (GetLayer()) {
    case (ABSLAYER):
      return SR1PHASE;
    case (S1LAYER):
      return SR2PHASE;
    case (S2LAYER):
      return SR3PHASE;
    default:
      throw;
  }
#else
  throw;
#endif
}

void Tri2k::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  const Tri2k *tri2k = (Tri2k*)orig;
  Her2kProps::DuplicateProps(tri2k);
  DLAOp<5,1>::Duplicate(orig,shallow,possMerging);
}


void Tri2k::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<5,1>::Prop();
    /*
#if DOELEM
    if (GetLayer() == SMLAYER) {
    
      if (m_type == REAL) {
	if (m_trans != NORMAL && m_trans != TRANS)
	  throw;
      }
      else {
	if (m_trans != NORMAL && m_trans != CONJTRANS)
	  throw;
      }
    
      if(m_type == REAL && m_trans == CONJTRANS)
	throw;
      else if (m_type == COMPLEX && m_trans == TRANS)
	throw;
    
      DistType t;
    
      if (m_trans == NORMAL) {
	t = InputDataType(0).m_dist;
	if ((t != D_MC_STAR) && (t != D_STAR_MC_T) && (t != D_STAR_MC_H))
	  throw;
      
	t = InputDataType(1).m_dist;
	if ((t != D_MR_STAR)
	    && (m_type == REAL && t != D_STAR_MR_T)
	    && (m_type == COMPLEX && t != D_STAR_MR_H))
	  throw;
      
	t = InputDataType(2).m_dist;
	if ((t != D_MC_STAR) && (t != D_STAR_MC_T) && (t != D_STAR_MC_H))
	  throw;
      
	t = InputDataType(3).m_dist;
	if ((t != D_MR_STAR)
	    && (m_type == REAL && t != D_STAR_MR_T)
	    && (m_type == COMPLEX && t != D_STAR_MR_H))
	  throw;
      }
      else {
	t = InputDataType(0).m_dist;
	if ((t != D_STAR_MC)
	    && (m_type == REAL && t != D_MC_STAR_T)
	    && (m_type == COMPLEX && t != D_MC_STAR_H))
	  throw;
      
	t = InputDataType(1).m_dist;
	if ((t != D_STAR_MR) && (t != D_MR_STAR_T))
	  throw;
      
	t = InputDataType(2).m_dist;
	if ((t != D_STAR_MC)
	    && (m_type == REAL && t != D_MC_STAR_T)
	    && (m_type == COMPLEX && t != D_MC_STAR_H))
	  throw;
      
	t = InputDataType(3).m_dist;
	if ((t != D_STAR_MR) && (t != D_MR_STAR_T))
	  throw;
      }
    
      if (InputDataType(4).m_dist != D_MC_MR)
	throw;
    }
    */
#if DOBLIS
    if (GetLayer() == S1LAYER || GetLayer() == S2LAYER || GetLayer() == S3LAYER) {
      if (*GetInputM(0) != *GetInputM(2))
	throw;
      if (*GetInputM(0) != *GetInputM(4))
	throw;
    
      if (*GetInputN(0) != *GetInputM(1))
	throw;
      if (*GetInputN(2) != *GetInputM(3))
	throw;
    
      if (*GetInputN(1) != *GetInputN(3))
	throw;
      if (*GetInputN(1) != *GetInputN(4))
	throw;
    }
#endif



#if DOELEM
    if (GetLayer() == SMLAYER) {
      const Sizes *localM = InputLocalM(4);
      const Sizes *localN = InputLocalN(4);
      const Sizes *others = (UpdateTrans(m_trans, InputDataType(0).m_dist) == NORMAL ?
                             InputLocalN(0) : InputLocalM(0));
      
      m_cost = GetCost(SMLAYER, localM, localN, others);
    }
#elif DOBLIS
    if (GetLayer() == S1LAYER || GetLayer() == S2LAYER || GetLayer() == S3LAYER) {
      const Sizes *localM = GetInputM(4);
      const Sizes *localN = GetInputN(4);
      const Sizes *others = (GetInputN(0));
      
      m_cost = GetCost(S1LAYER, localM, localN, others);
    }
#endif
  }
}

Cost Tri2k::GetCost(Layer layer, const Sizes *localDim1, const Sizes *localDim2, const Sizes *localDim3)
{
#if DOELEM
  if (layer == SMLAYER)
#elif DOBLIS
   if (layer == S1LAYER)
#endif
    return 2 * GAMMA * localDim1->SumProds111(*localDim2, *localDim3);
  else
    throw;
}

void Tri2k::PrintCode(IndStream &out)
{
#if DOELEM
  if (GetLayer() == SMLAYER) {
    out.Indent();
    *out << "LocalTrr2k"
    << "( " << TriToStr(m_tri) << ", \n" << out.Tabs(1) << "";
    
    Trans otherExpected;
    if (m_trans == NORMAL) {
      if (m_type == REAL)
        otherExpected = TRANS;
      else
        otherExpected = CONJTRANS;
    }
    else
      otherExpected = NORMAL;
    
    Trans upTrans;
    upTrans = UpdateTrans(m_trans, InputDataType(0).m_dist);
    if (upTrans != NORMAL)
      *out << TransToStr(upTrans) << ", \n" << out.Tabs(1) << "";
    
    upTrans = UpdateTrans(otherExpected, InputDataType(1).m_dist);
    if (upTrans != NORMAL)
      *out << TransToStr(upTrans) << ", \n" << out.Tabs(1) << "";
    
    upTrans = UpdateTrans(m_trans, InputDataType(2).m_dist);
    if (upTrans != NORMAL)
      *out << TransToStr(upTrans) << ", \n" << out.Tabs(1) << "";
    
    upTrans = UpdateTrans(otherExpected, InputDataType(3).m_dist);
    if (upTrans != NORMAL)
      *out << TransToStr(upTrans) << ", \n" << out.Tabs(1) << "";
    
    out << m_alpha;
    *out << ", "
    << GetInputName(0).str() << ", "
    << GetInputName(1).str() << ", \n" << out.Tabs(1) << ""
    << GetInputName(2).str() << ", "
    << GetInputName(3).str() << ", \n" << out.Tabs(1) << "";
    out << m_beta;
    *out << ", "
    << GetInputName(4).str() << ");\n";
  }
#elif DOBLIS
  if (GetLayer() == S1LAYER || GetLayer() == S2LAYER || GetLayer() == S3LAYER) {
    out.Indent();
    *out << "BLISTrr2k" << LayerNumToStr(GetLayer())
    << "( " << TriToStr(m_tri) << ", \n" << out.Tabs(1) << "";
    out << m_alpha;
    *out << ", "
    << GetInputName(0).str() << ", "
    << GetInputName(1).str() << ", \n" << out.Tabs(1) << ""
    << GetInputName(2).str() << ", "
    << GetInputName(3).str() << ", \n" << out.Tabs(1) << "";
    out << m_beta;
    *out << ", "
    << GetInputName(4).str() << ");\n";
    
  }
#endif
}

void Tri2k::FlattenCore(ofstream &out) const
{
  DLAOp<5,1>::FlattenCore(out);
  Her2kProps::FlattenCore(out);
}

void Tri2k::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<5,1>::UnflattenCore(in, info);
  Tri2k::UnflattenCore(in,info);
}

#if DOELEM
bool DistHer2kToLocalTri2k::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Her2k::GetClass()
      && ((Her2k*)node)->GetLayer() == DMLAYER)
  {
    return true;
  }
  else
    return false;
}

void DistHer2kToLocalTri2k::Apply(Node *node) const
{
  Her2k *orig = (Her2k*)node;
  bool normal = orig->m_trans == NORMAL;
  RedistNode *redist1 = new RedistNode(normal ? D_MC_STAR : D_STAR_MC);
  RedistNode *redist2 = new RedistNode(normal ? D_MR_STAR : D_STAR_MR);
  RedistNode *redist3 = new RedistNode(normal ? D_MC_STAR : D_STAR_MC);
  RedistNode *redist4 = new RedistNode(normal ? D_MR_STAR : D_STAR_MR);
  switch(orig->m_trans) {
    case(NORMAL):
    case(TRANS):
    case(CONJTRANS):
      break;
    default:
      throw;
  }
  Tri2k *tri2k = new Tri2k(SMLAYER, orig->m_tri, orig->m_trans,
                           orig->m_alpha, orig->m_beta,
                           orig->m_type);
  redist1->AddInput(node->Input(0),node->InputConnNum(0));
  redist2->AddInput(node->Input(1),node->InputConnNum(1));
  redist3->AddInput(node->Input(1),node->InputConnNum(1));
  redist4->AddInput(node->Input(0),node->InputConnNum(0));
  tri2k->AddInput(redist1,0);
  tri2k->AddInput(redist2,0);
  tri2k->AddInput(redist3,0);
  tri2k->AddInput(redist4,0);
  tri2k->AddInput(node->Input(2),node->InputConnNum(2));
  node->m_poss->AddNodes(4, redist1, redist2, redist3, redist4);
  node->m_poss->AddNode(tri2k);
  node->RedirectChildren(tri2k,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

Cost DistHer2kToLocalTri2k::RHSCostEstimate(const Node *node) const
{
  Her2k *orig = (Her2k*)node;
  bool normal = orig->m_trans == NORMAL;
  Cost cost = RedistNode::GetCost(D_MC_MR, normal ? D_MC_STAR : D_STAR_MC, orig->GetInputM(0), orig->GetInputN(0));
  cost += RedistNode::GetCost(D_MC_MR, normal ? D_MR_STAR : D_STAR_MR,orig->GetInputM(1), orig->GetInputN(1));
  cost += RedistNode::GetCost(D_MC_MR, normal ? D_MC_STAR : D_STAR_MC,orig->GetInputM(1), orig->GetInputN(1));
  cost += RedistNode::GetCost(D_MC_MR, normal ? D_MR_STAR : D_STAR_MR,orig->GetInputM(0), orig->GetInputN(0));
  
  Sizes localM, localN;
  GetLocalSizes(normal ? D_MC_STAR : D_STAR_MC, orig->GetInputM(0), orig->GetInputN(0), localM, localN);
  cost += Tri2k::GetCost(SMLAYER, orig->InputLocalM(2), orig->InputLocalN(2),
                         (UpdateTrans(orig->m_trans, normal ? D_MC_STAR : D_STAR_MC) == NORMAL ?
                          &localN : &localM));
  return cost;
}

string DistHer2kToLocalHer2kContrib::GetType() const
{
  return "Distributed Her2k to Local Her2k " + DistTypeToStr(m_AType)
  + " " + DistTypeToStr(m_BType);
}

bool DistHer2kToLocalHer2kContrib::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Her2k::GetClass()
      && ((Her2k*)node)->GetLayer() == DMLAYER)
  {
    return true;
  }
  else
    return false;
}

void DistHer2kToLocalHer2kContrib::Apply(Node *node) const
{
  Her2k *orig = (Her2k*)node;
  bool normal = orig->m_trans == NORMAL;
  RedistNode *redist1 = new RedistNode(normal ? m_AType : m_BType);
  RedistNode *redist2 = new RedistNode(normal ? m_BType : m_AType);
  TempVarNode *tmp = new TempVarNode(D_STAR_STAR);
  SumScatterNode *node3 = new SumScatterNode(orig->m_beta);
  
  redist1->AddInput(node->Input(0), node->InputConnNum(0));
  redist2->AddInput(node->Input(1), node->InputConnNum(1));
  tmp->AddInput(node->Input(2),node->InputConnNum(2));
  
  Her2k *her2k = new Her2k(SMLAYER, orig->m_tri, orig->m_trans,
                           orig->m_alpha, COEFZERO,
                           orig->m_type);
  her2k->AddInput(redist1,0);
  her2k->AddInput(redist2,0);
  her2k->AddInput(tmp,0);
  node3->AddInput(her2k);
  node3->AddInput(orig->Input(2),node->InputConnNum(2));
  node->m_poss->AddNodes(5, redist1, redist2, node3, tmp, her2k);
  node->RedirectChildren(node3,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

Cost DistHer2kToLocalHer2kContrib::RHSCostEstimate(const Node *node) const
{
  const Her2k *orig = (Her2k*)node;
  bool normal = orig->m_trans == NORMAL;
  const Sizes *A1 = orig->GetInputM(0);
  const Sizes *A2 = orig->GetInputN(0);
  const Sizes *B1 = orig->GetInputM(1);
  const Sizes *B2 = orig->GetInputN(1);
  const Sizes *C1 = orig->GetInputM(2);
  const Sizes *C2 = orig->GetInputN(2);
  Cost cost = RedistNode::GetCost(D_MC_MR, normal ? m_AType : m_BType, A1, A2);
  cost += RedistNode::GetCost(D_MC_MR, normal ? m_BType : m_AType, B1, B2);
  Sizes localA1, localA2;
  GetLocalSizes(normal ? m_AType : m_BType, A1, A2, localA1, localA2);
  cost += Her2k::GetCost(SMLAYER, C1, normal ? &localA2 : &localA1);
  Sizes localC1, localC2;
  GetLocalSizes(D_MC_MR, C1, C2, localC1, localC2);
  cost += SumScatterNode::GetCost(&localC1, &localC2, D_MC_MR, D_STAR_STAR);
  return cost;
}



string Tri2kTrans::GetTransType() const
{
  return "Tri2k";
}

bool Tri2kTrans::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Tri2k::GetClass()
      || ((Tri2k*)node)->GetLayer() != SMLAYER)
    return false;
  Tri2k *tri2k = (Tri2k*)node;
  
  if (m_trans == CONJTRANS && tri2k->m_type != COMPLEX)
    return false;
  
  if (m_argNum == 0) {
    if (tri2k->m_trans == CONJTRANS && m_trans == TRANS)
      return false;
    else
      return ((DLANode*)tri2k->Input(0))->CanTrans();
  }
  else if (m_argNum == 1) {
    if (tri2k->m_trans == NORMAL && m_trans == TRANS && tri2k->m_type == COMPLEX)
      return false;
    else
      return ((DLANode*)tri2k->Input(1))->CanTrans();
  }
  else if (m_argNum == 2) {
    if (tri2k->m_trans == CONJTRANS && m_trans == TRANS)
      return false;
    else
      return ((DLANode*)tri2k->Input(2))->CanTrans();
  }
  else if (m_argNum == 3) {
    if (tri2k->m_trans == NORMAL && m_trans == TRANS && tri2k->m_type == COMPLEX)
      return false;
    else
      return ((DLANode*)tri2k->Input(3))->CanTrans();
  }
  else
    throw;
}
#endif

RealLoop* Her2kLoopVar1(Node *Ain, ConnNum Anum,
                    Node *Bin, ConnNum Bnum,
                    Node *Cin, ConnNum Cnum,
                    Tri tri,
                    Trans trans,
                    Coef alpha, Coef beta, Type type,
                    Layer layer)
{
  SplitSingleIter *splitA = new SplitSingleIter(trans==NORMAL ? PARTDOWN : PARTRIGHT, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(trans==NORMAL ? PARTDOWN : PARTRIGHT, POSSTUNIN);
  splitB->AddInput(Bin, Bnum);
  splitB->SetAllStats(FULLUP);
  splitB->SetIndepIters();
  
  SplitSingleIter *splitC = new SplitSingleIter(PARTDIAG, POSSTUNIN, true);
  splitC->AddInput(Cin, Cnum);
  if (tri == LOWER)
    splitC->SetUpStats(FULLUP, FULLUP,
                       NOTUP, NOTUP);
  else
    splitC->SetUpStats(FULLUP, FULLUP,
                       FULLUP, NOTUP);
  splitC->SetIndepIters();
  
  Node *gemm1;
  if (trans == NORMAL)
    gemm1 = new Gemm(layer, NORMAL, type == REAL ? TRANS : CONJTRANS, alpha, beta, type);
  else
    gemm1 = new Gemm(layer, type == REAL ? TRANS : CONJTRANS, NORMAL, alpha, beta, type);
  
  if (tri == LOWER)
    gemm1->AddInputs(6,
                     splitA, 1,
                     splitB, 0,
                     splitC, 1);
  else
    gemm1->AddInputs(6,
                     splitA, 1,
                     splitB, 2,
                     splitC, 7);
	
  
  Node *gemm2;
  if (trans == NORMAL)
    gemm2 = new Gemm(layer, NORMAL, type == REAL ? TRANS : CONJTRANS, alpha, COEFONE, type);
  else
    gemm2 = new Gemm(layer, type == REAL ? TRANS : CONJTRANS, NORMAL, alpha, COEFONE, type);
  if (tri == LOWER)
    gemm2->AddInputs(6,
                     splitB, 1,
                     splitA, 0,
                     gemm1, 0);
  else
    gemm2->AddInputs(6,
                     splitB, 1,
                     splitA, 2,
                     gemm1, 0);
  
  Node *her2k;
  her2k = new Her2k(layer, tri, trans, alpha, COEFONE, type);
  
  her2k->AddInputs(6, splitA, 1,
                   splitB, 1,
                   splitC, 4);
  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB = splitB->CreateMatchingCombine(0);
  
  CombineSingleIter *comC;
  if (tri == LOWER)
    comC = splitC->CreateMatchingCombine(2,
                                         1, gemm2, 0,
                                         4, her2k, 0);
  else
    comC = splitC->CreateMatchingCombine(2,
                                         4, her2k, 0,
                                         7, gemm2, 0);
  
  Poss *loopPoss = new Poss(3, comA, comB, comC);
  RealLoop *loop;
#if DOELEM
  if (layer == DMLAYER)
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  else
    throw;
#elif DOBLIS
    loop = new RealLoop(BLISLOOP, loopPoss, BlisKC);
#endif
  
  return loop;
}


RealLoop* Her2kLoopVar2(Node *Ain, ConnNum Anum,
                    Node *Bin, ConnNum Bnum,
                    Node *Cin, ConnNum Cnum,
                    Tri tri,
                    Trans trans,
                    Coef alpha, Coef beta, Type type,
                    Layer layer)
{
  ScaleTrapNode *scal = new ScaleTrapNode(layer, LEFT, tri, beta);
  scal->AddInput(Cin, Cnum);
  
  SplitSingleIter *splitA = new SplitSingleIter(trans==NORMAL ? PARTDOWN : PARTRIGHT, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(trans==NORMAL ? PARTDOWN : PARTRIGHT, POSSTUNIN);
  splitB->AddInput(Bin, Bnum);
  splitB->SetAllStats(FULLUP);
  splitB->SetIndepIters();
  
  SplitSingleIter *splitC = new SplitSingleIter(PARTDIAG, POSSTUNIN, true);
  splitC->AddInput(scal, 0);
  if (tri == LOWER)
    splitC->SetUpStats(FULLUP, FULLUP,
                       PARTUP, NOTUP);
  else
    splitC->SetUpStats(FULLUP, PARTUP,
                       FULLUP, NOTUP);
  
  Node *gemm1;
  if (trans == NORMAL)
    gemm1 = new Gemm(layer, NORMAL, type == REAL ? TRANS : CONJTRANS, alpha, COEFONE, type);
  else
    gemm1 = new Gemm(layer, type == REAL ? TRANS : CONJTRANS, NORMAL, alpha, COEFONE, type);
  
  if (tri == LOWER)
    gemm1->AddInputs(6,
                     splitB, 1,
                     splitA, 0,
                     splitC, 1);
  else
    gemm1->AddInputs(6,
                     splitB, 0,
                     splitA, 1,
                     splitC, 3);
	
  
  Node *gemm2;
  if (trans == NORMAL)
    gemm2 = new Gemm(layer, NORMAL, type == REAL ? TRANS : CONJTRANS, alpha, COEFONE, type);
  else
    gemm2 = new Gemm(layer, type == REAL ? TRANS : CONJTRANS, NORMAL, alpha, COEFONE, type);
  
  if (tri == LOWER)
    gemm2->AddInputs(6,
                     splitA, 2,
                     splitB, 1,
                     splitC, 5);
  else
    gemm2->AddInputs(6,
                     splitA, 1,
                     splitB, 2,
                     splitC, 7);
  
  Node *her2k;
  her2k = new Her2k(layer, tri, trans, alpha, COEFONE, type);
  
  her2k->AddInputs(6, splitA, 1,
                   splitB, 1,
                   splitC, 4);
  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB = splitB->CreateMatchingCombine(0);
  
  CombineSingleIter *comC;
  if (tri == LOWER)
    comC = splitC->CreateMatchingCombine(3,
                                         1, gemm1, 0,
                                         4, her2k, 0,
                                         5, gemm2, 0);
  else
    comC = splitC->CreateMatchingCombine(3,
                                         3, gemm1, 0,
                                         4, her2k, 0,
                                         7, gemm2, 0);
  
  Poss *loopPoss = new Poss(3, comA, comB, comC);
  RealLoop *loop;
#if DOELEM
  if (layer == DMLAYER)
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  else
    throw;
#elif DOBLIS
    loop = new RealLoop(BLISLOOP, loopPoss, BlisKC);
#endif
  
  return loop;
}

RealLoop* Her2kLoopVar3(Node *Ain, ConnNum Anum,
                    Node *Bin, ConnNum Bnum,
                    Node *Cin, ConnNum Cnum,
                    Tri tri,
                    Trans trans,
                    Coef alpha, Coef beta, Type type,
                    Layer layer)
{
  ScaleTrapNode *scal = new ScaleTrapNode(layer, LEFT, tri, beta);
  scal->AddInput(Cin, Cnum);
  
  SplitSingleIter *splitA = new SplitSingleIter(trans==NORMAL ? PARTDOWN : PARTRIGHT, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(trans==NORMAL ? PARTDOWN : PARTRIGHT, POSSTUNIN);
  splitB->AddInput(Bin, Bnum);
  splitB->SetAllStats(FULLUP);
  splitB->SetIndepIters();
  
  SplitSingleIter *splitC = new SplitSingleIter(PARTDIAG, POSSTUNIN, true);
  splitC->AddInput(scal, 0);
  if (tri == LOWER)
    splitC->SetUpStats(FULLUP, FULLUP,
                       PARTUP, NOTUP);
  else
    splitC->SetUpStats(FULLUP, PARTUP,
                       FULLUP, NOTUP);
  
  Node *gemm1;
  if (trans == NORMAL)
    gemm1 = new Gemm(layer, NORMAL, type == REAL ? TRANS : CONJTRANS, alpha, COEFONE, type);
  else
    gemm1 = new Gemm(layer, type == REAL ? TRANS : CONJTRANS, NORMAL, alpha, COEFONE, type);
  
  if (tri == LOWER)
    gemm1->AddInputs(6,
                     splitA, 1,
                     splitB, 0,
                     splitC, 1);
  else
    gemm1->AddInputs(6,
                     splitA, 0,
                     splitB, 1,
                     splitC, 3);
	
  
  Node *gemm2;
  if (trans == NORMAL)
    gemm2 = new Gemm(layer, NORMAL, type == REAL ? TRANS : CONJTRANS, alpha, COEFONE, type);
  else
    gemm2 = new Gemm(layer, type == REAL ? TRANS : CONJTRANS, NORMAL, alpha, COEFONE, type);
  
  if (tri == LOWER)
    gemm2->AddInputs(6,
                     splitB, 2,
                     splitA, 1,
                     splitC, 5);
  else
    gemm2->AddInputs(6,
                     splitB, 1,
                     splitA, 2,
                     splitC, 7);
  
  Node *her2k;
  her2k = new Her2k(layer, tri, trans, alpha, COEFONE, type);
  
  her2k->AddInputs(6, splitA, 1,
                   splitB, 1,
                   splitC, 4);
  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB = splitB->CreateMatchingCombine(0);
  
  CombineSingleIter *comC;
  if (tri == LOWER)
    comC = splitC->CreateMatchingCombine(3,
                                         1, gemm1, 0,
                                         4, her2k, 0,
                                         5, gemm2, 0);
  else
    comC = splitC->CreateMatchingCombine(3,
                                         3, gemm1, 0,
                                         4, her2k, 0,
                                         7, gemm2, 0);
  
  Poss *loopPoss = new Poss(3, comA, comB, comC);
  RealLoop *loop;
#if DOELEM
  if (layer==DMLAYER)
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  else
    throw;
#elif DOBLIS
    loop = new RealLoop(BLISLOOP, loopPoss, BlisKC);
#endif
  
  return loop;
}



RealLoop* Her2kLoopVar4(Node *Ain, ConnNum Anum,
                    Node *Bin, ConnNum Bnum,
                    Node *Cin, ConnNum Cnum,
                    Tri tri,
                    Trans trans,
                    Coef alpha, Coef beta, Type type,
                    Layer layer)
{
  ScaleTrapNode *scal = new ScaleTrapNode(layer, LEFT, tri, beta);
  scal->AddInput(Cin, Cnum);
  
  SplitSingleIter *splitA = new SplitSingleIter(trans==NORMAL ? PARTDOWN : PARTRIGHT, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(trans==NORMAL ? PARTDOWN : PARTRIGHT, POSSTUNIN);
  splitB->AddInput(Bin, Bnum);
  splitB->SetAllStats(FULLUP);
  splitB->SetIndepIters();
  
  SplitSingleIter *splitC = new SplitSingleIter(PARTDIAG, POSSTUNIN, true);
  splitC->AddInput(scal, 0);
  if (tri == LOWER)
    splitC->SetUpStats(FULLUP, FULLUP,
                       FULLUP, NOTUP);
  else
    splitC->SetUpStats(FULLUP, NOTUP,
                       FULLUP, NOTUP);
  splitC->SetIndepIters();
  
  Node *gemm1;
  if (trans == NORMAL)
    gemm1 = new Gemm(layer, NORMAL, type == REAL ? TRANS : CONJTRANS, alpha, COEFONE, type);
  else
    gemm1 = new Gemm(layer, type == REAL ? TRANS : CONJTRANS, NORMAL, alpha, COEFONE, type);
  
  if (tri == LOWER)
    gemm1->AddInputs(6,
                     splitA, 2,
                     splitB, 1,
                     splitC, 5);
  else
    gemm1->AddInputs(6,
                     splitA, 0,
                     splitB, 1,
                     splitC, 3);
	
  
  Node *gemm2;
  if (trans == NORMAL)
    gemm2 = new Gemm(layer, NORMAL, type == REAL ? TRANS : CONJTRANS, alpha, COEFONE, type);
  else
    gemm2 = new Gemm(layer, type == REAL ? TRANS : CONJTRANS, NORMAL, alpha, COEFONE, type);
  
  if (tri == LOWER)
    gemm2->AddInputs(6,
                     splitA, 2,
                     splitB, 1,
                     gemm1, 0);
  else
    gemm2->AddInputs(6,
                     splitB, 0,
                     splitA, 1,
                     gemm1, 0);
  
  Node *her2k;
  her2k = new Her2k(layer, tri, trans, alpha, COEFONE, type);
  
  her2k->AddInputs(6, splitA, 1,
                   splitB, 1,
                   splitC, 4);
  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB = splitB->CreateMatchingCombine(0);
  
  CombineSingleIter *comC;
  if (tri == LOWER)
    comC = splitC->CreateMatchingCombine(2,
                                         4, her2k, 0,
                                         5, gemm2, 0);
  else
    comC = splitC->CreateMatchingCombine(2,
                                         3, gemm2, 0,
                                         4, her2k, 0);
  
  Poss *loopPoss = new Poss(3, comA, comB, comC);
  RealLoop *loop;
#if DOELEM
  if (layer == DMLAYER)
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  else
    throw;
#elif DOBLIS
    loop = new RealLoop(BLISLOOP, loopPoss, BlisKC);
#endif
  
  return loop;
}


RealLoop* Her2kLoopVar9(Node *Ain, ConnNum Anum,
                    Node *Bin, ConnNum Bnum,
                    Node *Cin, ConnNum Cnum,
                    Tri tri,
                    Trans trans,
                    Coef alpha, Coef beta, Type type,
                    Layer layer)
{
  ScaleTrapNode *scal = new ScaleTrapNode(layer, LEFT, tri, beta);
  scal->AddInput(Cin, Cnum);
  
  SplitSingleIter *splitA = new SplitSingleIter(trans==NORMAL ? PARTRIGHT : PARTDOWN, POSSTUNIN, true);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(trans==NORMAL ? PARTRIGHT : PARTDOWN, POSSTUNIN);
  splitB->AddInput(Bin, Bnum);
  splitB->SetAllStats(FULLUP);
  splitB->SetIndepIters();
  
  LoopTunnel *CTun = new LoopTunnel(POSSTUNIN);
  CTun->AddInput(scal, 0);
  CTun->SetAllStats(PARTUP);
  
  Node *her2k;
  her2k = new Her2k(layer, tri, trans, alpha, COEFONE, type);
  
  her2k->AddInputs(6, splitA, 1,
                   splitB, 1,
                   CTun, 0);
  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB = splitB->CreateMatchingCombine(0);
  
  LoopTunnel *CtunOut = new LoopTunnel(POSSTUNOUT);
  CtunOut->AddInput(her2k,0);
  CtunOut->AddInput(CTun,0);
  CtunOut->CopyTunnelInfo(CTun);
  
  Poss *loopPoss = new Poss(3, comA, comB, CtunOut);
  RealLoop *loop;
#if DOELEM
  if (layer == DMLAYER)
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  else
    throw;
#elif DOBLIS
    loop = new RealLoop(BLISLOOP, loopPoss, BlisKC);
#endif
  
  return loop;
}

RealLoop* Tri2kLoopVar10(Node *Ain, ConnNum Anum,
                     Node *Bin, ConnNum Bnum,
                     Node *Cin, ConnNum Cnum,
                     Node *Din, ConnNum Dnum,
                     Node *Ein, ConnNum Enum,
                     Tri tri,
                     Coef alpha, Coef beta, Type type,
                     Layer layer)
#if DOBLIS
{
  LoopTunnel *Atun = new LoopTunnel(POSSTUNIN);
  Atun->AddInput(Ain, Anum);
  Atun->SetAllStats(FULLUP);
  Atun->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(tri==LOWER ? PARTRIGHT : PARTLEFT, POSSTUNIN, false);
  splitB->AddInput(Bin, Bnum);
  splitB->SetAllStats(FULLUP);
  splitB->SetIndepIters();
  
  LoopTunnel *Ctun = new LoopTunnel(POSSTUNIN);
  Ctun->AddInput(Cin, Cnum);
  Ctun->SetAllStats(FULLUP);
  Ctun->SetIndepIters();
  
  SplitSingleIter *splitD = new SplitSingleIter(tri==LOWER ? PARTRIGHT : PARTLEFT, POSSTUNIN, false);
  splitD->AddInput(Din, Dnum);
  splitD->SetAllStats(FULLUP);
  splitD->SetIndepIters();
  
  SplitSingleIter *splitE = new SplitSingleIter(tri==LOWER ? PARTRIGHT : PARTLEFT, POSSTUNIN, true);
  splitE->AddInput(Ein, Enum);
  if (tri==LOWER)
    splitE->SetUpStats(FULLUP, NOTUP,
                       FULLUP, NOTUP);
  else
    splitE->SetUpStats(NOTUP, FULLUP,
                       NOTUP, FULLUP);
  splitE->SetIndepIters();

  GetUpToDiag *diag = new GetUpToDiag(tri, PARTRIGHT);
  if (tri == LOWER)
    diag->AddInputs(8,
                    splitE, 0,
                    splitE, 1,
                    Atun, 0,
                    Ctun, 0);
  else
    diag->AddInputs(8,
                    splitE, 2,
                    splitE, 1,
                    Atun, 0,
                    Ctun, 0);

  
  Node *tri2k = new Tri2k(layer, tri, NORMAL, alpha, beta, type);
  tri2k->AddInputs(10,
                   diag, 1,
                   splitB, 1,
                   diag, 2,
                   splitD, 1,
                   diag, 0);
  
  LoopTunnel *AtunOut = new LoopTunnel(POSSTUNOUT);
  AtunOut->AddInput(Atun,0);
  AtunOut->AddInput(Atun,0);
  AtunOut->CopyTunnelInfo(Atun);
  
  CombineSingleIter *comB = splitB->CreateMatchingCombine(0);
  
  LoopTunnel *CtunOut = new LoopTunnel(POSSTUNOUT);
  CtunOut->AddInput(Ctun,0);
  CtunOut->AddInput(Ctun,0);
  CtunOut->CopyTunnelInfo(Ctun);
  
  CombineSingleIter *comD = splitD->CreateMatchingCombine(0);
  
  CombineSingleIter *comE = splitE->CreateMatchingCombine(1,
                                                1, tri2k, 0);
  
  Poss *loopPoss = new Poss(5, AtunOut, comB, CtunOut, comD, comE);
  RealLoop *loop;
#if DOELEM
  if (layer == DMLAYER)
    throw;
  else
#endif
    loop = new RealLoop(BLISLOOP, loopPoss, BlisNC);

  loop->SetDimName(DIMN);
  
  return loop;
}
#else
{
  throw;
}
#endif

RealLoop* Tri2kLoopVar9(Node *Ain, ConnNum Anum,
                     Node *Bin, ConnNum Bnum,
                     Node *Cin, ConnNum Cnum,
                     Node *Din, ConnNum Dnum,
                     Node *Ein, ConnNum Enum,
                     Tri tri,
                     Coef alpha, Coef beta, Type type,
                     Layer layer)
{
  ScaleTrapNode *scal = new ScaleTrapNode(layer, LEFT, tri, beta);
  scal->AddInput(Ein, Enum);
  
  SplitSingleIter *splitA = new SplitSingleIter(PARTRIGHT, POSSTUNIN, true);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(PARTDOWN, POSSTUNIN);
  splitB->AddInput(Bin, Bnum);
  splitB->SetAllStats(FULLUP);
  splitB->SetIndepIters();

  SplitSingleIter *splitC = new SplitSingleIter(PARTRIGHT, POSSTUNIN);
  splitC->AddInput(Cin, Cnum);
  splitC->SetAllStats(FULLUP);
  splitC->SetIndepIters();
  
  SplitSingleIter *splitD = new SplitSingleIter(PARTDOWN, POSSTUNIN);
  splitD->AddInput(Din, Dnum);
  splitD->SetAllStats(FULLUP);
  splitD->SetIndepIters();
  
  LoopTunnel *ETun = new LoopTunnel(POSSTUNIN);
  ETun->AddInput(scal, 0);
  ETun->SetAllStats(PARTUP);
  
  Node *tri2k = new Tri2k(layer, tri, NORMAL, alpha, COEFONE, type);
  
  tri2k->AddInputs(10, splitA, 1,
                   splitB, 1,
		   splitC, 1,
		   splitD, 1,
                   ETun, 0);

#if DOELEM
    throw;
#elif DOBLIS  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB = splitB->CreateMatchingCombine(0);

  CombineSingleIter *comC = splitC->CreateMatchingCombine(0);
  
  CombineSingleIter *comD = splitD->CreateMatchingCombine(0);
  
  LoopTunnel *EtunOut = new LoopTunnel(POSSTUNOUT);
  EtunOut->AddInput(tri2k,0);
  EtunOut->AddInput(ETun,0);
  EtunOut->CopyTunnelInfo(ETun);
  
  RealLoop *loop;

  Poss *loopPoss = new Poss(5, comA, comB, comC, comD, EtunOut);
  loop = new RealLoop(BLISLOOP, loopPoss, BlisKC);
  loop->SetDimName(DIMK);
  
  return loop;
#endif

}


string Tri2kToTriRK::GetType() const
{
  return "Tri2K to 2 TriRKs " + LayerNumToStr(m_fromLayer);
}

bool Tri2kToTriRK::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Tri2k::GetClass()
      && ((Her2k*)node)->GetLayer() == m_fromLayer)
  {
    return true;
  }
  else
    return false;
}

void Tri2kToTriRK::Apply(Node *node) const
{
  Tri2k *orig = (Tri2k*)node;
  
  if (orig->m_beta != COEFONE)
    throw;
  
  TriRK *one = new TriRK(m_toLayer, orig->m_tri,
                         NORMAL, NORMAL, orig->m_alpha, COEFONE,
                         orig->m_type);
  
  
  TriRK *two = new TriRK(m_toLayer, orig->m_tri,
                         NORMAL, NORMAL, orig->m_alpha, COEFONE,
                         orig->m_type);
  
  one->AddInputs(6,
                 orig->Input(0), orig->InputConnNum(0),
                 orig->Input(1), orig->InputConnNum(1),
                 orig->Input(4), orig->InputConnNum(4));
  
  two->AddInputs(6,
                 orig->Input(2), orig->InputConnNum(2),
                 orig->Input(3), orig->InputConnNum(3),
                 one, 0);
  
  node->m_poss->AddNode(one);
  node->m_poss->AddNode(two);
  
  node->RedirectChildren(two,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

bool Tri2kLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Tri2k::GetClass()) {
    const Tri2k *tri2k = (Tri2k*)node;
    if (tri2k->GetLayer() != m_fromLayer)
      return false;
    if (m_dim == DIMK) {
      return (*(tri2k->GetInputN(0)) <= m_bs);
    }
    else if (m_dim == DIMN)
      return (*(tri2k->GetInputN(4)) <= m_bs);
    else
      throw;
  }
  return false;  
}

void Tri2kLowerLayer::Apply(Node *node) const
{
  Tri2k *tri2k = (Tri2k*)node;
  tri2k->SetLayer(m_toLayer);
}

string Tri2kLowerLayer::GetType() const
{ 
  return "Tri2k lower layer " + LayerNumToStr(m_fromLayer) 
  + " to " + LayerNumToStr(m_toLayer);
}

#if DOBLIS
bool Her2kToTri2K::CanApply(const Node *node) const
{
  return (node->GetNodeClass() == Her2k::GetClass() &&
          ((DLANode*)node)->GetLayer() == m_layer);
}

void Her2kToTri2K::Apply(Node *node) const
{
  Her2k *her2k = (Her2k*)node;
  
  Node *Ain = her2k->Input(0);
  ConnNum Anum = her2k->InputConnNum(0);
  
  Node *Bin = her2k->Input(1);
  ConnNum Bnum = her2k->InputConnNum(1);
  
  Transpose *Atrans = new Transpose(her2k->m_type == REAL ? 
                                    TRANS : CONJTRANS,
                                    true);
  Atrans->AddInput(Ain, Anum);
  
  Transpose *Btrans = new Transpose(her2k->m_type == REAL ? 
                                    TRANS : CONJTRANS,
                                    true);
  Btrans->AddInput(Bin, Bnum);
  
  Tri2k *tri2k = new Tri2k(m_layer,
                           her2k->m_tri,
                           NORMAL,
                           her2k->m_alpha,
                           her2k->m_beta,
                           her2k->m_type);
  
  if (her2k->m_trans == NORMAL) {
    tri2k->AddInputs(8,
                     Ain, Anum,
                     Btrans, 0,
                     Bin, Bnum,
                     Atrans, 0);
  }
  else {
    tri2k->AddInputs(8,
                     Btrans, 0,
                     Ain, Anum,
                     Atrans, 0,
                     Bin, Bnum);
  }
  
  tri2k->AddInput(her2k->Input(2), her2k->InputConnNum(2));
  
  node->m_poss->AddNode(tri2k);
  node->m_poss->AddNode(Atrans);
  node->m_poss->AddNode(Btrans);
  
  her2k->RedirectChildren(tri2k, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}
#endif


string Tri2kLoopExp::GetType() const
{
  switch(m_var) {
      //  case(1):
      //    return "Herk Loop Exp - var 1";
      //  case(2):
      //    return "Herk Loop Exp - var 2";
  case(9):
    return "Tri2k Loop Exp - var 9 (dim K)";
    case(10):
      return "Tri2k Loop Exp - var 10 (dim N)";
  default:
    throw;    
  }
}

bool Tri2kLoopExp::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Tri2k::GetClass()
      && ((Tri2k*)node)->GetLayer() == m_fromLayer)
    return true;
  return false;
}

void Tri2kLoopExp::Apply(Node *node) const
{
  Tri2k *trirk = (Tri2k*)node;
  RealLoop *loop;
  
  NodeConn *connA, *connB, *connC, *connD, *connE;
  connA = trirk->m_inputs[0];
  connB = trirk->m_inputs[1];
  connC = trirk->m_inputs[2];
  connD = trirk->m_inputs[3];
  connE = trirk->m_inputs[4];
  
  switch(m_var) {
    /*
      case(1):
      loop = Tri2kLoopVar1(connA->m_n, connA->m_num,
      connC->m_n, connC->m_num,
      trirk->m_tri,
      trirk->m_transA,
      trirk->m_alpha, trirk->m_beta, trirk->m_type, m_toLayer);
      break;
      case(2):
      loop = Tri2kLoopVar2(connA->m_n, connA->m_num,
      connC->m_n, connC->m_num,
      trirk->m_tri,
      trirk->m_transA,
      trirk->m_alpha, trirk->m_beta, trirk->m_type, m_toLayer);
      break;
    */
  case(9):
    loop = Tri2kLoopVar9(connA->m_n, connA->m_num,
			 connB->m_n, connB->m_num,
			 connC->m_n, connC->m_num,
			 connD->m_n, connD->m_num,
			 connE->m_n, connE->m_num,
			 trirk->m_tri,
			 trirk->m_alpha, trirk->m_beta, trirk->m_type, m_toLayer);
    break;
    case(10):
      loop = Tri2kLoopVar10(connA->m_n, connA->m_num,
                            connB->m_n, connB->m_num,
                            connC->m_n, connC->m_num,
                            connD->m_n, connD->m_num,
                            connE->m_n, connE->m_num,
                            trirk->m_tri,
                            trirk->m_alpha, trirk->m_beta, trirk->m_type, m_toLayer);
      break;
    default:
      throw;
  }
  
  node->m_poss->AddPSet(loop);
  
  node->RedirectChildren(loop->OutTun(2),0);
  node->m_poss->DeleteChildAndCleanUp(node);
}
#endif
