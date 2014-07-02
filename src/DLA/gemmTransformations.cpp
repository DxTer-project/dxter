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

#include "gemm.h"

#if DOBLIS||DOELEM||DOLLDLA
#include "elemRedist.h"
#include "string.h"
#include "helperNodes.h"
#include "loopSupport.h"
#include "blas.h"
#include "pack.h"
#include "gemmTransformations.h"
#include "smmul.h"

using namespace std;

GemmLoopExp::GemmLoopExp(Layer fromLayer, Layer toLayer, int dim)
  : 
  m_fromLayer(fromLayer), 
  m_toLayer(toLayer), 
  m_dim(dim) 
{
#if DOELEM
  m_bsSize = USELEMBS;
#elif DOBLIS
  switch(dim) {
  case (0) :
    m_bsSize = USEBLISMC;
      break;
  case (1) :
    m_bsSize = USEBLISKC;
      break;
  case (2) : 
    m_bsSize = USEBLISNC;
  }
#elif DOLLDLA
  throw;
#endif
}

GemmLoopExp::GemmLoopExp(Layer fromLayer, Layer toLayer, int dim, BSSize bsSize)
  : 
  m_fromLayer(fromLayer), 
  m_toLayer(toLayer), 
  m_dim(dim),
  m_bsSize(bsSize)
{
}

string GemmLoopExp::GetType() const
{
  //If these change, match in LLDLAGemmLoopExp
  string str = "Gemm Loop Exp " 
    + LayerNumToStr(m_fromLayer)
    + " + " 
    + LayerNumToStr(m_toLayer)
    + " bs:" + std::to_string(BSSizeToSize(m_bsSize));
  switch(m_dim) {
    case(0):
      return str + " - m";
    case(1):
      return str + " - k";
    case(-1):
      return str + " - k reversed";
    case(2):
      return str + " - n";
    default:
      throw;
  }
}

bool GemmLoopExp::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Gemm::GetClass()) {
    const Gemm *gemm = (Gemm*)node;
    if (gemm->GetLayer() == m_fromLayer)
      return true;
  }
  return false;
}

void GemmLoopExp::Apply(Node *node) const
{
  Gemm *gemm = (Gemm*)node;
  Loop *loop;
  
  NodeConn *connA, *connB, *connC;
  connA = gemm->m_inputs[0];
  connB = gemm->m_inputs[1];
  connC = gemm->m_inputs[2];
  
  switch(m_dim) {
    case(0):
      loop = GemmVar1Loop(connA->m_n, connA->m_num,
                          connB->m_n, connB->m_num,
                          connC->m_n, connC->m_num,
			  m_bsSize,
                          gemm->m_transA, gemm->m_transB,
                          gemm->m_alpha, gemm->m_beta, m_toLayer, gemm->m_type);
      break;
    case(1):
      loop = GemmVar3Loop(connA->m_n, connA->m_num,
                          connB->m_n, connB->m_num,
                          connC->m_n, connC->m_num,
			  m_bsSize,
                          gemm->m_transA, gemm->m_transB,
                          false,
                          gemm->m_alpha, gemm->m_beta, m_toLayer, gemm->m_type);
      break;
#if DOBLIS||DOELEM
    case(-1):
      loop = GemmVar3Loop(connA->m_n, connA->m_num,
                          connB->m_n, connB->m_num,
                          connC->m_n, connC->m_num,
			  m_bsSize,
                          gemm->m_transA, gemm->m_transB,
                          true,
                          gemm->m_alpha, gemm->m_beta, m_toLayer, gemm->m_type);
      break;
#endif //DOBLIS||DOELEM
    case(2):
      loop = GemmVar2Loop(connA->m_n, connA->m_num,
                          connB->m_n, connB->m_num,
                          connC->m_n, connC->m_num,
			  m_bsSize,
                          gemm->m_transA, gemm->m_transB,
                          gemm->m_alpha, gemm->m_beta, m_toLayer, gemm->m_type);
      break;
    default:
      throw;
  }
  
  if (!(loop->GetBS()))
    throw;
  
  node->m_poss->AddLoop(loop);
  
  node->RedirectChildren(loop->OutTun(2),0);
  node->m_poss->DeleteChildAndCleanUp(node);
}



#if DOELEM
bool DistGemmToLocalGemmStatC::CanApply(const Node *node) const
{
  return IsDMGemm(node);
}


void DistGemmToLocalGemmStatC::Apply(Node *node) const
{
  Gemm *orig = (Gemm*)node;
  RedistNode *node1 = new RedistNode(orig->m_transA==NORMAL ? D_MC_STAR : D_STAR_MC);
  RedistNode *node2 = new RedistNode(orig->m_transB==NORMAL ? D_STAR_MR : D_MR_STAR);
  Gemm *node3 = new Gemm(SMLAYER, orig->m_transA, orig->m_transB, orig->m_alpha, orig->m_beta,orig->m_type);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node2->AddInput(node->Input(1),node->InputConnNum(1));
  node3->AddInput(node1,0);
  node3->AddInput(node2,0);
  node3->AddInput(node->Input(2),node->InputConnNum(2));
  node->m_poss->AddNode(node1);
  node->m_poss->AddNode(node2);
  node->m_poss->AddNode(node3);
  node->RedirectChildren(node3,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

Cost DistGemmToLocalGemmStatC::RHSCostEstimate(const Node *node) const
{
  const Gemm *orig = (Gemm*)node;
  bool normA = orig->m_transA==NORMAL;
  bool normB = orig->m_transB==NORMAL;
  
  const DLANode *input = (DLANode*)(orig->Input(0));
  unsigned int inputNum = orig->InputConnNum(0);
  
  const Sizes *A1 = input->GetM(inputNum);
  const Sizes *A2 = input->GetN(inputNum);
  const Sizes *B1 = orig->GetInputM(1);
  const Sizes *B2 = orig->GetInputN(1);
  const Sizes *localC2 = orig->InputLocalN(2);
  DistType aDist = normA ? D_MC_STAR : D_STAR_MC;
  Cost cost = RedistNode::GetCost(D_MC_MR, aDist, A1, A2);
  cost += RedistNode::GetCost(D_MC_MR, normB ? D_STAR_MR : D_MR_STAR, B1, B2);
  Sizes localA1, localA2;
  GetLocalSizes(aDist, A1, A2, localA1, localA2);
  cost += Gemm::GetCost(SMLAYER,&localA1, &localA2, localC2);
  return cost;
}


bool DistGemmToContribLocalGemmStatANonTrans::CanApply(const Node *node) const
{
  if (IsDMGemm(node) &&
      (((Gemm*)node)->m_transA == NORMAL)) {
    return true;
  }
  return false;
}

void DistGemmToContribLocalGemmStatANonTrans::Apply(Node *node) const
{
  Gemm *orig = (Gemm*)node;
  Trans transB = orig->m_transB;
  DistType redistType = D_LASTDIST;
  if (transB == NORMAL)
    redistType = D_MR_STAR;
  else
    redistType = D_STAR_MR;
  
  RedistNode *redistB = new RedistNode(redistType);
  Gemm *localGemm = new Gemm(SMLAYER, orig->m_transA, transB, orig->m_alpha, COEFZERO, orig->m_type);
  SumScatterNode *node3 = new SumScatterNode(orig->m_beta);
  TempVarNode *tmp = new TempVarNode(D_MC_STAR);
  tmp->AddInput(node->Input(2),node->InputConnNum(2));
  redistB->AddInput(node->Input(1),node->InputConnNum(1));
  localGemm->AddInput(node->Input(0),node->InputConnNum(0));
  localGemm->AddInput(redistB,0);
  localGemm->AddInput(tmp,0);
  
  node3->AddInput(localGemm,0);
  node3->AddInput(node->Input(2),node->InputConnNum(2));
  
  node->m_poss->AddNode(redistB);
  node->m_poss->AddNode(localGemm);
  node->m_poss->AddNode(node3);
  node->m_poss->AddNode(tmp);
  orig->RedirectChildren(node3,0);
  orig->m_poss->DeleteChildAndCleanUp(orig);
}

Cost DistGemmToContribLocalGemmStatANonTrans::RHSCostEstimate(const Node *node) const
{
  const Gemm *orig = (Gemm*)node;
  Trans transB = orig->m_transB;
  DistType redistType = D_LASTDIST;
  if (transB == NORMAL)
    redistType = D_MR_STAR;
  else
    redistType = D_STAR_MR;
  const Sizes *localA1 = orig->InputLocalM(0);
  const Sizes *localA2 = orig->InputLocalN(0);
  const Sizes *B1 = orig->GetInputM(1);
  const Sizes *B2 = orig->GetInputN(1);
  const Sizes *C1 = orig->GetInputM(2);
  const Sizes *C2 = orig->GetInputN(2);
  Cost cost = RedistNode::GetCost(D_MC_MR, redistType, B1, B2);
  Sizes localC1, localC2;
  GetLocalSizes(D_MC_STAR, C1, C2, localC1, localC2);
  cost += Gemm::GetCost(SMLAYER, localA1, localA2, &localC2);
  cost += SumScatterNode::GetCost(&localC1, &localC2, D_MC_MR, D_MC_STAR);
  return cost;
}

//e.g. in GemmNNB
bool DistGemmToContribLocalGemmStatBNonTrans::CanApply(const Node *node) const
{
  if (IsDMGemm(node) &&
      (((Gemm*)node)->m_transB != CONJTRANS) &&
      (((Gemm*)node)->m_transB != TRANS)) {
    return true;
  }
  return false;
}

void DistGemmToContribLocalGemmStatBNonTrans::Apply(Node *node) const
{
  Gemm *orig = (Gemm*)node;
  Trans transA = orig->m_transA;
  Trans transB = orig->m_transB;
  DistType redistType = D_LASTDIST;
  if (transA == NORMAL)
    redistType = D_STAR_MC;
  else
    redistType = D_MC_STAR;
  if (!m_trans) {
    RedistNode *node1 = new RedistNode(redistType);
    Gemm *node2 = new Gemm(SMLAYER, orig->m_transA, orig->m_transB, orig->m_alpha, COEFZERO, orig->m_type);
    SumScatterNode *node3 = new SumScatterNode(orig->m_beta);
    TempVarNode *tmp = new TempVarNode(D_STAR_MR);
    tmp->AddInput(node->Input(2),node->InputConnNum(2));
    node1->AddInput(node->Input(0),node->InputConnNum(0));
    node2->AddInput(node1,0);
    if (transB == NORMAL) {
      node2->AddInput(node->Input(1),node->InputConnNum(1));
    }
    else if (transB == TRANS || transB == CONJTRANS) {
      RedistNode *redist = new RedistNode(D_MR_MC);
      redist->AddInput(node->Input(1),node->InputConnNum(1));
      node2->AddInput(redist,0);
      node->m_poss->AddNode(redist);
    }
    node2->AddInput(tmp,0);
    
    node3->AddInput(node2,0);
    node3->AddInput(node->Input(2),node->InputConnNum(2));
    
    node->m_poss->AddNode(node1);
    node->m_poss->AddNode(node2);
    node->m_poss->AddNode(node3);
    node->m_poss->AddNode(tmp);
    node->RedirectChildren(node3,0);
    orig->m_poss->DeleteChildAndCleanUp(orig);
  }
  else {
    RedistNode *node1 = new RedistNode(redistType);
    Gemm *node2 = new Gemm(SMLAYER,
                           SwitchTrans(orig->m_transB,orig->m_type),
                           SwitchTrans(orig->m_transA,orig->m_type),
                           orig->m_alpha, COEFZERO, orig->m_type);
    SumScatterNode *node3 = new SumScatterNode(orig->m_beta);
    TempVarNode *tmp = new TempVarNode(orig->m_type==COMPLEX ? D_MR_STAR_H : D_MR_STAR_T);
    tmp->AddInput(node->Input(2),node->InputConnNum(2));
    node1->AddInput(node->Input(0),node->InputConnNum(0));
    if (transB == NORMAL) {
      node2->AddInput(node->Input(1),node->InputConnNum(1));
    }
    else if (transB == TRANS || transB == CONJTRANS) {
      RedistNode *redist = new RedistNode(D_MR_MC);
      redist->AddInput(node->Input(1),node->InputConnNum(1));
      node2->AddInput(redist,0);
      node->m_poss->AddNode(redist);
    }
    node2->AddInput(node1,0);
    node2->AddInput(tmp,0);
    
    node3->AddInput(node2,0);
    node3->AddInput(node->Input(2),node->InputConnNum(2));
    
    node->m_poss->AddNode(node1);
    node->m_poss->AddNode(node2);
    node->m_poss->AddNode(node3);
    node->m_poss->AddNode(tmp);
    node->RedirectChildren(node3,0);
    orig->m_poss->DeleteChildAndCleanUp(orig);
  }
}

Cost DistGemmToContribLocalGemmStatBNonTrans::RHSCostEstimate(const Node *node) const
{
  Gemm *orig = (Gemm*)node;
  Trans transA = orig->m_transA;
  Trans transB = orig->m_transB;
  DistType redistType = D_LASTDIST;
  if (transA == NORMAL)
    redistType = D_STAR_MC;
  else
    redistType = D_MC_STAR;
  const Sizes *A1 = orig->GetInputM(0);
  const Sizes *A2 = orig->GetInputN(0);
  const Sizes *B1 = orig->GetInputM(1);
  const Sizes *B2 = orig->GetInputN(1);
  const Sizes *C1 = orig->GetInputM(2);
  const Sizes *C2 = orig->GetInputN(2);
  if (!m_trans) {
    Cost cost = RedistNode::GetCost(D_MC_MR, redistType, A1, A2);
    Sizes localA1, localA2;
    GetLocalSizes(redistType, A1, A2, localA1, localA2);
    Sizes localC1, localC2;
    GetLocalSizes(D_STAR_MR, C1, C2, localC1, localC2);
    if (transB != NORMAL) {
      cost += RedistNode::GetCost(D_MC_MR, D_MR_MC, B1, B2);
    }
    cost += Gemm::GetCost(SMLAYER, &localA1, &localA2, &localC2);
    cost += SumScatterNode::GetCost(&localC1, &localC2, D_MC_MR, D_STAR_MR);
    return cost;
  }
  else {
    Cost cost = RedistNode::GetCost(D_MC_MR, redistType, A1, A2);
    Sizes localA1, localA2;
    GetLocalSizes(redistType, A1, A2, localA1, localA2);
    Sizes localC1, localC2;
    GetLocalSizes(orig->m_type==COMPLEX?D_MR_STAR_H:D_MR_STAR_T, C1, C2, localC1, localC2);
    if (transB != NORMAL) {
      cost += RedistNode::GetCost(D_MC_MR, D_MR_MC, B1, B2);
    }
    cost += Gemm::GetCost(SMLAYER,&localC1, &localA1, &localA2);
    cost += SumScatterNode::GetCost(&localC1, &localC2, D_MC_MR, orig->m_type==COMPLEX?D_MR_STAR_H:D_MR_STAR_T);
    return cost;
  }
}


bool DistGemmToContribLocalGemmStatBTrans::CanApply(const Node *node) const
{
  if (IsDMGemm(node) &&
      (((Gemm*)node)->m_transB != NORMAL)) {
    return true;
  }
  return false;
}

void DistGemmToContribLocalGemmStatBTrans::Apply(Node *node) const
{
  Gemm *orig = (Gemm*)node;
  Trans transA = orig->m_transA;
  Trans transB = orig->m_transB;
  DistType redistType = D_LASTDIST;
  if (transA == NORMAL)
    redistType = D_STAR_MR;
  else
    redistType = D_MR_STAR;
  RedistNode *node1 = new RedistNode(redistType);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  Gemm *node2 = new Gemm(SMLAYER, orig->m_transA, orig->m_transB, orig->m_alpha, COEFZERO, orig->m_type);
  node2->AddInput(node1,0);
  TempVarNode *tmp = new TempVarNode(D_STAR_MC);
  tmp->AddInput(node->Input(2),node->InputConnNum(2));
  if (transB == NORMAL)
    throw;
  node2->AddInput(node->Input(1),node->InputConnNum(1));
  node2->AddInput(tmp,0);
  
  TempVarNode *tmp2 = new TempVarNode(D_MR_MC);
  tmp2->AddInput(node->Input(2),node->InputConnNum(2));
  
  SumScatterFrom *node3 = new SumScatterFrom;
  node3->AddInput(node2,0);
  node3->AddInput(tmp2,0);
  
  RedistNode *redist = new RedistNode(D_MC_MR);
  redist->AddInput(node3,0);
  
  
  Axpy *axpy = new Axpy(SMLAYER, COEFONE);
  axpy->AddInputs(4,
                  redist, 0,
                  node->Input(2), node->InputConnNum(2));
  
  node->m_poss->AddNodes(7,
                 node1,
                 node2,
                 node3,
                 tmp,
                 tmp2,
                 redist,
                 axpy);
  node->RedirectChildren(axpy,0);
  orig->m_poss->DeleteChildAndCleanUp(orig);
}

Cost DistGemmToContribLocalGemmStatBTrans::RHSCostEstimate(const Node *node) const
{
  Gemm *orig = (Gemm*)node;
  Trans transA = orig->m_transA;
  DistType redistType = D_LASTDIST;
  if (transA == NORMAL)
    redistType = D_STAR_MR;
  else
    redistType = D_MR_STAR;
  const Sizes *A1 = orig->GetInputM(0);
  const Sizes *A2 = orig->GetInputN(0);
  const Sizes *C1 = orig->GetInputM(2);
  const Sizes *C2 = orig->GetInputN(2);
  Cost cost = RedistNode::GetCost(D_MC_MR, redistType, A1, A2);
  Sizes localA1, localA2;
  GetLocalSizes(redistType, A1, A2, localA1, localA2);
  Sizes localC1, localC2;
  GetLocalSizes(D_STAR_MC, C1, C2, localC1, localC2);
  cost += Gemm::GetCost(SMLAYER, &localA1, &localA2, &localC2);
  localC1.ClearSizes();
  localC2.ClearSizes();
  GetLocalSizes(D_MR_MC, C1, C2, localC1, localC2);
  cost += SumScatterFrom::GetCost(D_MR_MC, D_STAR_MC, &localC1, &localC2);
  cost += RedistNode::GetCost(D_MR_MC, D_MC_MR, C1, C2);
  cost += Axpy::GetCost(SMLAYER, &localC1, &localC2);
  return cost;
}

bool DistGemmToContribLocalGemmStatATrans::CanApply(const Node *node) const
{
  if (IsDMGemm(node) &&
      (((Gemm*)node)->m_transA != NORMAL)) {
    return true;
  }
  return false;
}

void DistGemmToContribLocalGemmStatATrans::Apply(Node *node) const
{
  Gemm *orig = (Gemm*)node;
  Trans transA = orig->m_transA;
  Trans transB = orig->m_transB;
  DistType redistType = D_LASTDIST;
  if (transB == NORMAL)
    redistType = D_MC_STAR;
  else
    redistType = D_STAR_MC;
  RedistNode *redistB = new RedistNode(redistType);
  redistB->AddInput(node->Input(1),node->InputConnNum(1));
  Gemm *localGemm = new Gemm(SMLAYER, transA, transB, orig->m_alpha, COEFZERO, orig->m_type);
  TempVarNode *tmp = new TempVarNode(D_MR_STAR);
  tmp->AddInput(node->Input(2),node->InputConnNum(2));
  if (transA == NORMAL)
    throw;
  localGemm->AddInput(node->Input(0),node->InputConnNum(0));
  localGemm->AddInput(redistB,0);
  localGemm->AddInput(tmp,0);
  
  TempVarNode *tmp2 = new TempVarNode(D_MR_MC);
  tmp2->AddInput(node->Input(2),node->InputConnNum(2));
  
  SumScatterFrom *node3 = new SumScatterFrom;
  node3->AddInput(localGemm,0);
  node3->AddInput(tmp2,0);
  
  RedistNode *redist = new RedistNode(D_MC_MR);
  redist->AddInput(node3,0);
  
  Axpy *axpy = new Axpy(SMLAYER, COEFONE);
  axpy->AddInputs(4,
                  redist, 0,
                  node->Input(2), node->InputConnNum(2));
  
  node->m_poss->AddNodes(7,
                 redistB,
                 localGemm,
                 node3,
                 tmp,
                 tmp2,
                 redist,
                 axpy);
  node->RedirectChildren(axpy,0);
  orig->m_poss->DeleteChildAndCleanUp(orig);
}

Cost DistGemmToContribLocalGemmStatATrans::RHSCostEstimate(const Node *node) const
{
  Gemm *orig = (Gemm*)node;
  Trans transB = orig->m_transB;
  DistType redistType = D_LASTDIST;
  if (transB == NORMAL)
    redistType = D_MC_STAR;
  else
    redistType = D_STAR_MC;
  const Sizes *B1 = orig->GetInputM(1);
  const Sizes *B2 = orig->GetInputN(1);
  const Sizes *C1 = orig->GetInputM(2);
  const Sizes *C2 = orig->GetInputN(2);
  Cost cost = RedistNode::GetCost(D_MC_MR, redistType, B1, B2);
  Sizes localB1, localB2;
  GetLocalSizes(redistType, B1, B2, localB1, localB2);
  Sizes localC1, localC2;
  GetLocalSizes(D_MR_STAR, C1, C2, localC1, localC2);
  cost += Gemm::GetCost(SMLAYER, &localC1, &localB1, &localB2);
  localC1.ClearSizes();
  localC2.ClearSizes();
  GetLocalSizes(D_MR_MC, C1, C2, localC1, localC2);
  cost += SumScatterFrom::GetCost(D_MR_MC, D_MR_STAR, &localC1, &localC2);
  cost += RedistNode::GetCost(D_MR_MC, D_MC_MR, C1, C2);
  cost += Axpy::GetCost(SMLAYER, &localC1, &localC2);
  return cost;
}


string GemmTrans::GetTransType() const
{
  return "Gemm";
}

bool GemmTrans::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Gemm::GetClass())
    return false;
  Gemm *gemm = (Gemm*)node;
  if (gemm->GetLayer() != SMLAYER)
    return false;
  DLANode *source = (DLANode*)gemm->Input(m_argNum);
  Trans trans = (m_argNum==0 ? gemm->m_transA : gemm->m_transB);
  // m_trans==CONJTRANS might work, but this why expend the additional flops when trans
  //  reduces communication just the same
  if (trans == NORMAL && m_trans==TRANS && gemm->m_type != COMPLEX)
    return source->CanTrans();
  else if (trans == NORMAL && m_trans==CONJTRANS && gemm->m_type == COMPLEX)
    return source->CanTrans();
  else if (trans == CONJTRANS && m_trans == CONJTRANS) {
    return source->CanTrans();
  }
  else if (trans == TRANS && m_trans == TRANS)
    return source->CanTrans();
  else
    return false;
}


bool GemmInputReordering::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Gemm::GetClass())
    return false;
  Gemm *gemm = (Gemm*)node;
  if (gemm->GetLayer() != DMLAYER)
    return false;
  unsigned int Anum, Bnum;
  DLANode *A = gemm->FindNonRedistParent(0, Anum);
  DLANode *B = gemm->FindNonRedistParent(1, Bnum);
  
  if (!A || !B)
    throw;
  
  if (m_type == 0) {
    if (IsDMTrxm(A, true)) {
      DLANode *BUser = B->FindSideEffectingUser(Bnum);
      if (!BUser) {
        return false;
      }
      else if (!IsDMTrxm(BUser, true)) {
        return false;
      }
      B = BUser;
    }
    else  {
      return false;
    }
  }
  else if (m_type == 1) {
    if (IsDMTrxm(B, true)) {
      DLANode *AUser = A->FindSideEffectingUser(Anum);
      if (!AUser || !IsDMTrxm(AUser, true))
        return false;
      A = AUser;
    }
    else
      return false;
  }
  else
    throw;
  
  Trxm *ATrsm = (Trxm*)A;
  Trxm *BTrsm = (Trxm*)B;
  if (ATrsm->FindNonRedistParent(0) != BTrsm->FindNonRedistParent(0))
    return false;
  
  if (ATrsm->m_tri != BTrsm->m_tri) {
    cout << "A : " << TriToStr(ATrsm->m_tri) << endl;
    cout << "B : " << TriToStr(BTrsm->m_tri) << endl;
    throw;
  }
  if (gemm->m_transA == NORMAL) {
    if (gemm->m_transB == NORMAL) {
      //N and N
      return ATrsm->m_side == RIGHT &&
      BTrsm->m_side == LEFT &&
      ATrsm->m_trans == BTrsm->m_trans;
    }
    else {
      //N and (C/T)
      if( ATrsm->m_side == RIGHT &&
         BTrsm->m_side == RIGHT)
      {
        if (ATrsm->m_trans == NORMAL) {
          return BTrsm->m_trans == gemm->m_transB;
        }
        else
          return (ATrsm->m_trans == gemm->m_transB) && (BTrsm->m_trans == NORMAL);
      }
      else {
        cout << "A";
        return false;
      }
    }
  }
  else {
    if (gemm->m_transB == NORMAL) {
      //(C/T) and N
      if( ATrsm->m_side == LEFT &&
         BTrsm->m_side == LEFT)
      {
        if (ATrsm->m_trans == NORMAL) {
          return BTrsm->m_trans == gemm->m_transA;
        }
        else
          return (ATrsm->m_trans == gemm->m_transA) && (BTrsm->m_trans == NORMAL);
      }
      else {
        cout << "B";
        return false;
      }
    }
    else {
      //(C/T) and (C/T)
      if (ATrsm->m_trans == BTrsm->m_trans) {
        if( ATrsm->m_side == LEFT &&
           BTrsm->m_side == RIGHT)
        {
          if (ATrsm->m_trans == NORMAL) {
            return BTrsm->m_trans == NORMAL;
          }
          else
            return ATrsm->m_trans == BTrsm->m_trans;
        }
        else {
          cout << "C";
          return false;
        }
      }
      else {
        cout << "D";
        return false;
      }
    }
  }
}

void GemmInputReordering::Apply(Node *node) const
{
  Gemm *gemm = (Gemm*)node;
  gemm->m_inverseOps.insert(m_inverse);
  unsigned int Anum, Bnum;
  DLANode *A = gemm->FindNonRedistParent(0, Anum);
  DLANode *B = gemm->FindNonRedistParent(1, Bnum);
  
  if (!A || !B)
    throw;
  
  if (m_type == 0) {
    if (IsDMTrxm(A, true)) {
      DLANode *BUser = B->FindSideEffectingUser(Bnum);
      if (!BUser || !IsDMTrxm(BUser, true))
        throw;
      B = BUser;
    }
    else
      throw;
  }
  else if (m_type == 1){
    if (IsDMTrxm(B, true)) {
      DLANode *AUser = A->FindSideEffectingUser(Anum);
      if (!AUser || !IsDMTrxm(AUser, true))
        throw;
      A = AUser;
    }
    else
      throw;
  }
  else
    throw;
  
  Trxm *ATrsm = (Trxm*)A;
  Trxm *BTrsm = (Trxm*)B;
  
  if (m_type == 0) {
    gemm->m_alpha = gemm->m_alpha * ATrsm->m_coeff / BTrsm->m_coeff;
    
    Node *oldIn = gemm->Input(0);
    oldIn->RemoveChild(gemm,gemm->InputConnNum(0));
    gemm->ChangeInput1Way(gemm->Input(0),gemm->InputConnNum(0),ATrsm->Input(1),ATrsm->InputConnNum(1));
    if (!oldIn->m_children.size())
      gemm->m_poss->DeleteChildAndCleanUp(oldIn);
    
    oldIn = gemm->Input(1);
    oldIn->RemoveChild(gemm,gemm->InputConnNum(1));
    gemm->ChangeInput1Way(gemm->Input(1),gemm->InputConnNum(1),BTrsm,0);
    if (!oldIn->m_children.size())
      gemm->m_poss->DeleteChildAndCleanUp(oldIn);
  }
  else if (m_type == 1) {
    gemm->m_alpha = gemm->m_alpha * BTrsm->m_coeff / ATrsm->m_coeff;
    
    Node *oldIn = gemm->Input(1);
    oldIn->RemoveChild(gemm,gemm->InputConnNum(1));
    gemm->ChangeInput1Way(gemm->Input(1),gemm->InputConnNum(1),BTrsm->Input(1),BTrsm->InputConnNum(1));
    if (!oldIn->m_children.size())
      gemm->m_poss->DeleteChildAndCleanUp(oldIn);
    
    oldIn = gemm->Input(0);
    oldIn->RemoveChild(gemm,gemm->InputConnNum(0));
    gemm->ChangeInput1Way(gemm->Input(0),gemm->InputConnNum(0),ATrsm,0);
    if (!oldIn->m_children.size())
      gemm->m_poss->DeleteChildAndCleanUp(oldIn);
  }
  else
    throw;
}
#endif

Loop* GemmVar1Loop(Node *Ain, unsigned int Anum,
                   Node *Bin, unsigned int Bnum,
                   Node *Cin, unsigned int Cnum,
		   BSSize bs,
                   Trans transA, Trans transB,
                   Coef alpha, Coef beta,
                   Layer layer, Type type)
{
  SplitSingleIter *splitA = new SplitSingleIter(transA==NORMAL ? PARTDOWN : PARTRIGHT, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  LoopTunnel *Btun = new LoopTunnel(POSSTUNIN);
  Btun->AddInput(Bin, Bnum);
  Btun->SetAllStats(FULLUP);
  Btun->SetIndepIters();
  
  SplitSingleIter *splitC = new SplitSingleIter(PARTDOWN, POSSTUNIN, true);
  splitC->AddInput(Cin, Cnum);
  splitC->SetUpStats(FULLUP, FULLUP,
                     NOTUP, NOTUP);
  splitC->SetIndepIters();

  
  Node *gepp;
  
  gepp = new Gemm(layer, transA, transB, alpha, beta, type);
  gepp->AddInput(splitA, 1);
  gepp->AddInput(Btun, 0);
  gepp->AddInput(splitC, 1);
  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  LoopTunnel *BtunOut = new LoopTunnel(POSSTUNOUT);
  BtunOut->AddInput(Btun, 0);
  BtunOut->AddInput(Btun, 0);
  BtunOut->CopyTunnelInfo(Btun);
  
  CombineSingleIter *comC = splitC->CreateMatchingCombine(1,
                                                1, gepp, 0);
  
  Poss *loopPoss = new Poss(3, comA, BtunOut, comC);
  Loop *loop = NULL;
#if DOELEM
  if (layer == DMLAYER)
    loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  else
    throw;
#elif DOBLIS
  loop = new Loop(BLISLOOP, loopPoss, USEBLISMC);
#elif DOLLDLA
    loop = new Loop(LLDLALOOP, loopPoss, USELLDLAMU);
#endif

  loop->SetDimName(DIMM);
  
  return loop;
}


Loop* GemmVar3Loop(Node *Ain, unsigned int Anum,
                   Node *Bin, unsigned int Bnum,
                   Node *Cin, unsigned int Cnum,
		   BSSize bs,
                   Trans transA, Trans transB,
                   bool reverse,
                   Coef alpha, Coef beta,
                   Layer layer, Type type)
{
  PartDir aDir;
  if (transA == NORMAL) {
    if (reverse) {
      aDir = PARTLEFT;
    }
    else {
      aDir = PARTRIGHT;
    }
  }
  else {
    if (reverse) {
      aDir = PARTUPWARD;
    }
    else {
      aDir = PARTDOWN;
    }
  }
  SplitSingleIter *splitA = new SplitSingleIter(aDir, POSSTUNIN, true);
  splitA->AddInput(Ain, Anum);
  splitA->SetUpStats(FULLUP, FULLUP,
                     FULLUP, FULLUP);
  splitA->SetIndepIters();
  
  PartDir bDir;
  if (transB == NORMAL) {
    if (reverse) {
      bDir = PARTUPWARD;
    }
    else {
      bDir = PARTDOWN;
    }
  }
  else {
    if (reverse) {
      bDir = PARTLEFT;
    }
    else {
      bDir = PARTRIGHT;
    }
  }
  SplitSingleIter *splitB = new SplitSingleIter(bDir, POSSTUNIN);
  splitB->AddInput(Bin, Bnum);
  splitB->SetUpStats(FULLUP, FULLUP,
                     FULLUP, FULLUP);
  splitB->SetIndepIters();

#if DOLLDLA
  SMMul *scale = NULL;
  if (beta != COEFONE) {
    ConstVal *constVal = new ConstVal(beta.LLDLAStr(),beta);
    constVal->AddInput(Cin, Cnum);
    
    scale = new SMMul(type, layer);
    scale->SetLayer(layer);
    scale->AddInputs(4, 
		     constVal, 0,
		     Cin, Cnum);
  }
#else  
  ScaleNode *scale = NULL;
  if (beta != COEFONE) {
    scale = new ScaleNode(layer, beta);
    scale->AddInput(Cin, Cnum);
  }
#endif
  
  LoopTunnel *Ctun = new LoopTunnel(POSSTUNIN);
  if (beta != COEFONE) {
    Ctun->AddInput(scale, 0);
  }
  else {
    Ctun->AddInput(Cin, Cnum);
  }
  Ctun->SetUpStats(PARTUP, PARTUP,
                   PARTUP, PARTUP);
  
  Node *gepp;
  
  gepp = new Gemm(layer, transA, transB, alpha, COEFONE, type);
  gepp->AddInput(splitA, 1);
  gepp->AddInput(splitB, 1);
  gepp->AddInput(Ctun, 0);
  
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  
  CombineSingleIter *comB = splitB->CreateMatchingCombine(0);
  
  LoopTunnel *CtunOut = new LoopTunnel(POSSTUNOUT);
  CtunOut->AddInput(gepp, 0);
  CtunOut->AddInput(Ctun, 1);
  CtunOut->CopyTunnelInfo(Ctun);
  
  Poss *loopPoss = new Poss(3, comA, comB, CtunOut);
  Loop *loop;
#if DOELEM
  if (layer == DMLAYER)
    loop = new Loop(ELEMLOOP, loopPoss, bs);
  else
    throw;
#elif DOBLIS
    loop = new Loop(BLISLOOP, loopPoss, bs);
#elif DOLLDLA
    loop = new Loop(LLDLALOOP, loopPoss, bs);
#endif

  loop->SetDimName(DIMK);
  
  return loop;
}


Loop* GemmVar2Loop(Node *Ain, unsigned int Anum,
                   Node *Bin, unsigned int Bnum,
                   Node *Cin, unsigned int Cnum,
		   BSSize bs,
                   Trans transA, Trans transB,
                   Coef alpha, Coef beta,
                   Layer layer, Type type)
{
  LoopTunnel *Atun = new LoopTunnel(POSSTUNIN);
  Atun->AddInput(Ain, Anum);
  Atun->SetAllStats(FULLUP);
  Atun->SetIndepIters();
  
  SplitSingleIter *splitB = new SplitSingleIter(transB==NORMAL ? PARTRIGHT : PARTDOWN, POSSTUNIN);
  splitB->AddInput(Bin, Bnum);
  splitB->SetAllStats(FULLUP);
  splitB->SetIndepIters();
  
  SplitSingleIter *splitC = new SplitSingleIter(PARTRIGHT, POSSTUNIN, true);
  splitC->AddInput(Cin, Cnum);
  splitC->SetUpStats(FULLUP, NOTUP,
                     FULLUP, NOTUP);
  splitC->SetIndepIters();
  
  Node *gepp;
  
  gepp = new Gemm(layer, transA, transB, alpha, beta, type);
  gepp->AddInput(Atun, 0);
  gepp->AddInput(splitB, 1);
  gepp->AddInput(splitC, 1);
  
  
  LoopTunnel *AtunOut = new LoopTunnel(POSSTUNOUT);
  AtunOut->AddInput(Atun, 0);
  AtunOut->AddInput(Atun, 1);
  AtunOut->CopyTunnelInfo(Atun);
  
  CombineSingleIter *comB = splitB->CreateMatchingCombine(0);
  
  CombineSingleIter *comC = splitC->CreateMatchingCombine(1,
                                                1, gepp, 0);
  
  Poss *loopPoss = new Poss(3, AtunOut, comB, comC);
  Loop *loop = NULL;
#if DOELEM
  if (layer == DMLAYER)
    loop = new Loop(ELEMLOOP, loopPoss, bs);
  else
    throw;
#elif DOBLIS
    loop = new Loop(BLISLOOP, loopPoss, bs);
#else
    loop = new Loop(LLDLALOOP, loopPoss, bs);
#endif

  loop->SetDimName(DIMN);
  
  return loop;
}

#if DOBLIS
string BLISGemmLoopExp::GetType() const
{
  return "BLISGemmLoopExp";
}

bool BLISGemmLoopExp::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Gemm::GetClass()) {
    const Gemm *gemm = (Gemm*)node;
    if (gemm->GetLayer() == m_fromLayer)
      return true;
  }
  return false;
}

void BLISGemmLoopExp::Apply(Node *node) const
{
  Gemm *gemm = (Gemm*)node;
  
  NodeConn *connA, *connB, *connC;
  connA = gemm->m_inputs[0];
  connB = gemm->m_inputs[1];
  connC = gemm->m_inputs[2];
  
  Node *Bin = connB->m_n;
  unsigned int Bnum = connB->m_num;
  Node *Cin = connC->m_n;
  unsigned int Cnum = connC->m_num;
  
  SplitSingleIter *splitA = new SplitSingleIter(gemm->m_transA == NORMAL ? PARTDOWN : PARTRIGHT, POSSTUNIN);
  splitA->AddInput(connA->m_n, connA->m_num);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  if (gemm->m_transB != NORMAL) {
    Bin = AddTranspose(gemm->m_transB, true, Bin, Bnum, true);
    Bnum = 0;
  }
  PackBuff *bBuff = new PackBuff(Bin->GetName(Bnum).m_name,
                                 PACKCOLPANS, PACKBPANEL, NOTTRI, NOTTRIDIAG, GEN,
                                 false, false, false, false,
                                 USEKRSIZE, USENRSIZE );
  Pack *bPack = new Pack(PACKCOLPANS, 2, false, false, false, false, false);
  bBuff->AddInput(Bin,Bnum);
  bPack->AddInput(Bin,Bnum);
  bPack->AddInput(bBuff, 0);
  
  node->m_poss->AddNode(bBuff);
  node->m_poss->AddNode(bPack);
  
  LoopTunnel *Btun = new LoopTunnel(POSSTUNIN);
  Btun->AddInput(bPack, 0);
  Btun->SetAllStats(FULLUP);
  Btun->SetIndepIters();
  
  SplitSingleIter *splitC = new SplitSingleIter(PARTDOWN, POSSTUNIN, true);
  splitC->AddInput(Cin, Cnum);
  splitC->SetUpStats(FULLUP, FULLUP,
                     PARTUP, PARTUP);
  splitC->SetIndepIters();
  
  Node *Ain = splitA;
  unsigned int Anum = 1;
  string name;
  
  if (gemm->m_transA != NORMAL) {
    Ain = AddTranspose(gemm->m_transA, true, Ain, Anum, false);
    Anum = 0;
    name = Ain->GetName(0).str();
  }
  else {
    name = splitA->GetName(1, BLISLOOP).str();
  }
  
  PackBuff *aBuff = new PackBuff(name,
                                 PACKROWPANS, PACKABLOCK, NOTTRI, NOTTRIDIAG,
                                 GEN,
                                 false, false, false, false,
                                 USEMRSIZE, USEKRSIZE);
  Pack *aPack = new Pack(PACKROWPANS, 2, false, false, false, false, false);
  aBuff->AddInput(Ain, Anum);
  aPack->AddInput(Ain, Anum);
  aPack->AddInput(aBuff, 0);
  
  Gemm *gebp = new Gemm(m_toLayer, NORMAL, NORMAL,
                        gemm->m_alpha, gemm->m_beta, gemm->m_type);
  gebp->AddInput(aPack, 0);
  gebp->AddInput(Btun, 0);
  gebp->AddInput(splitC, 1);
  
  Combine *comA = splitA->CreateMatchingCombine(0);
  
  LoopTunnel *BtunOut = new LoopTunnel(POSSTUNOUT);
  BtunOut->AddInput(Btun, 0);
  BtunOut->AddInput(Btun, 0);
  BtunOut->CopyTunnelInfo(Btun);
  
  Combine *comC = splitC->CreateMatchingCombine(1,
                                                1, gebp, 0);
  
  Poss *loopPoss = new Poss(3, comA, BtunOut, comC);
  Loop *loop = new Loop(BLISLOOP, loopPoss, bs);

  loop->SetDimName(DIMM);
  
  
  
  node->m_poss->AddLoop(loop);
  
  node->RedirectChildren(loop->OutTun(2), 0);
  
  node->m_poss->DeleteChildAndCleanUp(node);
}
#endif

bool GemmLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Gemm::GetClass()) {
    const Gemm *gemm = (Gemm*)node;
    if (gemm->GetLayer() != m_fromLayer)
      return false;
    if (m_dim == DIMK) {
#if !DOLLDLA
      if (gemm->m_transA == NORMAL) {
        return (*(gemm->GetInputN(0)) <= m_bs);
      }
      else {
        return (*(gemm->GetInputM(0)) <= m_bs);
      }
#else
      return (*(gemm->GetInputN(0)) <= m_bs);
#endif
    }
    else if (m_dim == DIMN) {
#if !DOLLDLA
      if (gemm->m_transB == NORMAL) {
        return (*(gemm->GetInputN(1)) <= m_bs);
      }
      else {
        return (*(gemm->GetInputM(1)) <= m_bs);
      }
#else
      return (*(gemm->GetInputN(1)) <= m_bs);
#endif
    }
    else
      throw;
  }
  return false;
  
}

void GemmLowerLayer::Apply(Node *node) const
{
  Gemm *gemm = (Gemm*)node;
  gemm->SetLayer(m_toLayer);
}

string GemmLowerLayer::GetType() const
{ 
  return "Gemm lower layer " + LayerNumToStr(m_fromLayer) 
  + " to " + LayerNumToStr(m_toLayer);
}

#if DOBLIS||DOELEM
string SplitSingleIterGemm::GetType() const
{ 
  return "SplitSingleIter Gemm " + LayerNumToStr(m_layer);
}

bool SplitSingleIterGemm::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Gemm::GetClass()) {
    const Gemm *gemm = (Gemm*)node;
    if (gemm->GetLayer() != m_layer)
      return false;
    return gemm->m_beta == COEFONE;
  }
  return false;
}

void SplitSingleIterGemm::Apply(Node *node) const
{
  Gemm *gemm = (Gemm*)node;
  Node *Cin = gemm->Input(2);
  unsigned int Cnum = gemm->InputConnNum(2);

#if DOELEM
  TempVarNode *CTmp = new TempVarNode(D_MC_MR, "CTemp");
#else
  TempVarNode *CTmp = new TempVarNode("CTemp");
#endif
  CTmp->AddInput(Cin, Cnum);

  Axpy *axpy = new Axpy(m_layer, COEFONE);
  gemm->RedirectChildren(axpy, 0);

  gemm->ChangeInput2Way(gemm->Input(2), gemm->InputConnNum(2),
			CTmp, 0);
  gemm->m_beta = COEFZERO;

  axpy->AddInputs(4,
		  gemm, 0,
		  Cin, Cnum);

  node->m_poss->AddNode(CTmp);
  node->m_poss->AddNode(axpy);
}
#endif //DOBLIS||DOELEM
#endif
