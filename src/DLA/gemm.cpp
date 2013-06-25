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

#include "gemm.h"
#include "distributions.h"
#include "string.h"
#include "helperNodes.h"
#include "loopSupport.h"
#include "blas.h"
#include "pack.h"

using namespace std;

bool IsDMGemm(const Node *node)
{
  if (node->GetNodeClass() != Gemm::GetClass())
    return false;
  const Gemm *gemm = (Gemm*)node;
  if (gemm->m_layer != DMLAYER)
    return false;
  else
    return true;
}

Gemm::Gemm(Layer layer, Trans transA, Trans transB, Coef alpha, Coef beta, Type type)
: m_transA(transA),
m_transB(transB),
m_alpha(alpha),
m_beta(beta),
  m_type(type),
  m_comm(CORECOMM)
{
  SetLayer(layer);
}

bool Gemm::IsBLISParallelizable() const
{
  return GetLayer() == S3LAYER;
}

void Gemm::Parallelize(Comm comm)
{
  if (GetLayer() == S3LAYER)
    m_comm = comm;
  else
    throw;
}


Node* Gemm::BlankInst()
{
  return new Gemm(ABSLAYER, NORMAL, NORMAL, COEFONE, COEFONE, REAL);
}

NodeType Gemm::GetType() const
{
  return "Gemm "
  + TransToStr(m_transA)
  + TransToStr(m_transB) + " "
  + LayerNumToStr(GetLayer())
    + " " + CommToStr(m_comm);
}

void Gemm::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  const Gemm *gemm = (Gemm*)orig;
  m_transA = gemm->m_transA;
  m_transB = gemm->m_transB;
  m_alpha = gemm->m_alpha;
  m_beta = gemm->m_beta;
  m_type = gemm->m_type;
  m_comm = gemm->m_comm;
}

void Gemm::FlattenCore(ofstream &out) const
{
  DLAOp<3,1>::FlattenCore(out);
  WRITE(m_transA);
  WRITE(m_transB);
  WRITE(m_alpha);
  WRITE(m_beta);
  WRITE(m_type);
  WRITE(m_comm);
}

void Gemm::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<3,1>::UnflattenCore(in, info);
  READ(m_transA);
  READ(m_transB);
  READ(m_alpha);
  READ(m_beta);
  READ(m_type);
  READ(m_comm);
}

DistType Gemm::GetDistType(unsigned int num) const
{
#if DODPPHASE
  switch(GetLayer()) {
    case (ABSLAYER):
    case (DMLAYER):
      return D_MC_MR;
    case (SMLAYER):
      return InputDistType(2);
    default:
      throw;
  }
#else
  return InputDistType(2);
#endif
}

Phase Gemm::MaxPhase() const
{  switch(GetLayer()) {
  case (ABSLAYER):
#if DODPPHASE
  case (DMLAYER):
    return DPPHASE;
  case (SMLAYER):
    return NUMPHASES;
#else
    return SR1PHASE;
  case (S1LAYER):
    return SR2PHASE;
  case (S2LAYER):
    return SR3PHASE;
  case (S3LAYER):
    return NUMPHASES;
#endif
  default:
    throw;
}
}

bool Gemm::DoNotCullDP() const
{
  return GetLayer() == DMLAYER;
}

bool Gemm::ShouldCullSR() const
{
  if (GetLayer() == SMLAYER)
    return m_hasRefined;
  else
    return false;
}

bool Gemm::CanTransposeInputs() const
{
  return GetLayer() == SMLAYER;
}

Cost Gemm::GetCost(Layer layer, const Sizes *localDim1, const Sizes *localDim2, const Sizes *localDim3)
{
  if (layer == SMLAYER || layer == S1LAYER || layer == S2LAYER || layer == S3LAYER)
    return TWO * GAMMA * localDim1->SumProds111(*localDim2, *localDim3);
  else
    throw;
}

void Gemm::SanityCheck()
{
  if (GetLayer() == ABSLAYER || GetLayer() == DMLAYER) {
    if (InputDistType(2) != D_MC_MR) {
      cout << "input not D_MC_MR 7";
      throw;
    }
  }
}

void Gemm::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3,1>::Prop();
    switch(GetLayer()) {
      case(ABSLAYER):
      case (DMLAYER):
        m_cost = ZERO;
        break;
      case (SMLAYER):
      {
        DLANode *in0 = (DLANode*)Input(0);
        unsigned int num0 = InputConnNum(0);
        DLANode *in2 = (DLANode*)Input(2);
        unsigned int num2 = InputConnNum(2);
        const Sizes *size1 = in0->LocalM(num0);
        const Sizes *size2 = in0->LocalN(num0);
        const Sizes *size3 = in2->LocalN(num2);
        m_cost = GetCost(SMLAYER, size1, size2, size3);
        break;
      }
      case (S1LAYER):
      case (S2LAYER):
      case (S3LAYER):
        DLANode *in0 = (DLANode*)Input(0);
        unsigned int num0 = InputConnNum(0);
        DLANode *in2 = (DLANode*)Input(2);
        unsigned int num2 = InputConnNum(2);
        const Sizes *size1 = in0->LocalM(num0);
        const Sizes *size2 = in0->LocalN(num0);
        const Sizes *size3 = in2->LocalN(num2);
        if (size1->NumSizes() != size2->NumSizes())
          throw;
        if (size2->NumSizes() != size3->NumSizes())
          throw;
        m_cost = GetCost(S3LAYER, size1, size2, size3);
    }
  }
}


void LocalGemmTransUpdate(DistType t0, DistType t1, Trans &transA, Trans &transB)
{
  if (transA == NORMAL) {
    if (t0 == D_MC_MR || t0 == D_MC_STAR || t0 == D_STAR_MC || t0 == D_STAR_MR)
      transA = NORMAL;
    else if (t0 == D_STAR_MC_H || t0 == D_STAR_MR_H || t0 == D_MR_STAR_H)
      transA = CONJTRANS;
    else if (t0 == D_STAR_MC_T || t0 == D_STAR_MR_T || t0 == D_MR_STAR_T)
      transA = TRANS;
    else if (t0 == D_STAR_STAR)
      transA = NORMAL;
    else {
      cout << "BAD 1!!!!!!\n";
      cout << TransToStr(transA) << endl;
      cout << DistTypeToStr(t0) << endl;
      throw;
    }
  }
  else if (transA == TRANS) {
    if (t0 == D_STAR_MC || t0 == D_STAR_MR || t0 == D_MC_MR || t0 == D_MC_STAR || t0 == D_MR_STAR)
      transA = TRANS;
    else if (t0 == D_MC_STAR_T || t0 == D_MR_STAR_T || t0 == D_STAR_MR_T || t0 == D_STAR_MC_T)
      transA = NORMAL;
    else if (t0 == D_STAR_STAR)
      transA = TRANS;
    else {
      cout << "BAD2!!!!!!\n";
      cout << DistTypeToStr(t0) << endl;
      throw;
    }
  }
  else if (transA == CONJTRANS) {
    if (t0 == D_STAR_MC || t0 == D_STAR_MR || t0 == D_MC_STAR || t0 == D_MC_MR)
      transA = CONJTRANS;
    else if (t0 == D_MC_STAR_H || t0 == D_MR_STAR_H || t0 == D_STAR_MC_H)
      transA = NORMAL;
    else if (t0 == D_STAR_STAR)
      transA = CONJTRANS;
    else {
      cout << "BAD 3!!!!!!\n";
      throw;
    }
  }
  else {
    cout << "BAD 4!!!!!!\n";
    throw;
  }
  if (t0 == D_STAR_STAR) {
    transB = transB;
  }
  else if (t0 == D_MC_MR) {
    if (transB == NORMAL) {
      if (t1 == D_MR_STAR || t1 == D_MC_STAR)
        transB = (NORMAL);
      else if (t1 == D_STAR_MR_T || t1 == D_STAR_MC_T)
        transB = (TRANS);
      else if (t1 == D_STAR_MR_H || t1 == D_STAR_MC_H)
        transB = (CONJTRANS);
      else {
        cout << "BAD 5!!!!!!\n";
        cout << DistTypeToStr(t1) << endl;
        throw;
      }
    }
    else if (transB == TRANS) {
      if (t1 == D_STAR_MR || t1 == D_STAR_MC)
        transB = (TRANS);
      else if (t1 == D_MR_STAR_T || t1 == D_MC_STAR_T)
        transB = (NORMAL);
      else {
        cout << "BAD 6!!!!!!\n";
        throw;
      }
    }
    else if (transB == CONJTRANS) {
      if (t1 == D_STAR_MR || t1 == D_STAR_MC)
        transB = (CONJTRANS);
      else if (t1 == D_MR_STAR_H || t1 == D_MC_STAR_H)
        transB = (NORMAL);
      else {
        cout << "BAD 7!!!!!!\n";
        throw;
      }
    }
    else  {
      cout << "BAD 8!!!!!!\n";
      throw;
    }
  }
  else if (transB == NORMAL) {
    if (t1 == D_STAR_MR || t1 == D_MC_MR)
      transB = (NORMAL);
    else if (t1 == D_MR_STAR_T)
      transB = (TRANS);
    else if (t1 == D_MR_STAR_H)
      transB = (CONJTRANS);
    else {
      cout << "BAD 9!!!!!!\n";
      throw;
    }
  }
  else if (transB == TRANS) {
    if (t1 == D_STAR_MR_T)
      transB = (NORMAL);
    else if (t1 == D_MR_STAR)
      transB = (TRANS);
    else if (t1 == D_MR_MC || t1 == D_MC_MR)
      transB = (TRANS);
    else {
      cout << "BAD 10!!!!!!\n";
      throw;
    }
  }
  else if (transB == CONJTRANS) {
    if (t1 == D_STAR_MR_H)
      transB = (NORMAL);
    else if (t1 == D_MR_STAR)
      transB = (CONJTRANS);
    else if (t1 == D_MR_MC || t1 == D_MC_MR)
      transB = (CONJTRANS);
    else {
      cout << "BAD 11!!!!!!" << DistTypeToStr(t0) << ", " << DistTypeToStr(t1) << "\n";
      throw;
    }
  }
  else
    cout << "BAD 12!!!!!!\n";
  
}

void Gemm::PrintCode(IndStream &out)
{
  out.Indent();
  if (GetLayer() == ABSLAYER || GetLayer() == DMLAYER) {
    if (GetLayer() == ABSLAYER)
      *out << "AbsGemm( ";
    else
      *out << "DistGemm( ";
    *out << TransToStr(m_transA) << ", " << TransToStr(m_transB)
    << ", \n\t";
    out << m_alpha;
    *out << ", "
    << GetInputName(0).str() << ", " << GetInputName(1).str()
    << ", ";
    out << m_beta;
    *out << ", " << GetInputName(2).str() << " );\n";
  }
  else if (GetLayer() == SMLAYER) {
    string transAStr, transBStr;
    DistType t0 = InputDistType(0);
    DistType t1 = InputDistType(1);
    Trans transA = m_transA;
    Trans transB = m_transB;
    LocalGemmTransUpdate(t0, t1, transA, transB);
    transAStr = TransToStr(transA);
    transBStr = TransToStr(transB);
    
    *out << "internal::LocalGemm( " << transAStr << ", " << transBStr << ", \n" << out.Tabs(1);
    out << m_alpha;
    *out << ","
    << GetInputName(0).str() << ", " << GetInputName(1).str() << ", \n" << out.Tabs(1);
    out << m_beta;
    *out << ", " << GetInputName(2).str() << " );\n";
  }
  else if (GetLayer() == S1LAYER ||
           GetLayer() == S2LAYER ||
           GetLayer() == S3LAYER) {
    string transAStr = TransToStr(m_transA);
    string transBStr = TransToStr(m_transB);
    
    if (GetLayer() == S1LAYER) {
      *out << "BLISGemmLimitedN( ";
      *out << transAStr << ", " << transBStr << ", \n" << out.Tabs(1);
      out << m_alpha;
      *out << ","
      << GetInputName(0).str() << ", " << GetInputName(1).str() << ", \n" << out.Tabs(1);
      out << m_beta;
      *out << ", " << GetInputName(2).str() << " );\n";
    }
    else if (GetLayer() == S2LAYER) {
      *out << "GemmRankKUpdate( ";
      *out << transAStr << ", " << transBStr << ", \n" << out.Tabs(1);
      out << m_alpha;
      *out << ","
      << GetInputName(0).str() << ", " << GetInputName(1).str() << ", \n" << out.Tabs(1);
      out << m_beta;
      *out << ", " << GetInputName(2).str() << " );\n";
    }
    else if (GetLayer() == S3LAYER) {
      if (m_comm == CORECOMM) 
	*out << "bli_gemm_ker_var2( ";
      else
	*out << "bli_gemm_ker_var2_par( " << CommToStr(m_comm) << ", ";
      out << m_alpha;
      *out<< ", &"
      << GetInputName(0).str() << ", &" << GetInputName(1).str() << ", \n"
      << out.Tabs(2);
      out << m_beta;
      *out << ", &" << GetInputName(2).str() << ", (gemm_t*)NULL );\n";
    }
    else
      throw;
    
  }
  else {
    throw;
  }
}

void Gemm::UpdateInnerPackingMultiple(PackSize size)
{
  Node *input = Input(0);
  if (input->GetNodeClass() != Pack::GetClass())
    throw;
  Node *packInput = input->Input(1);
  if (packInput->GetNodeClass() != PackBuff::GetClass())
    throw;
  PackBuff *buff = (PackBuff*)packInput;
  if (buff->m_children.size() > 1)
    throw;
  buff->m_n = size;
}


string GemmLoopExp::GetType() const
{
  switch(m_dim) {
    case(0):
      return "Gemm Loop Exp - m";
    case(1):
      return "Gemm Loop Exp - k";
    case(-1):
      return "Gemm Loop Exp - k reversed";
    case(2):
      return "Gemm Loop Exp - n";
    default:
      throw;
  }
}

bool GemmLoopExp::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == Gemm::GetClass()) {
    const Gemm *gemm = (Gemm*)node;
    if (gemm->GetLayer() == m_fromLayer)
      return true;
  }
  return false;
}

void GemmLoopExp::Apply(Poss *poss, Node *node) const
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
                          gemm->m_transA, gemm->m_transB,
                          gemm->m_alpha, gemm->m_beta, m_toLayer, gemm->m_type);
      break;
    case(1):
      loop = GemmVar3Loop(connA->m_n, connA->m_num,
                          connB->m_n, connB->m_num,
                          connC->m_n, connC->m_num,
                          gemm->m_transA, gemm->m_transB,
                          false,
                          gemm->m_alpha, gemm->m_beta, m_toLayer, gemm->m_type);
      break;
    case(-1):
      loop = GemmVar3Loop(connA->m_n, connA->m_num,
                          connB->m_n, connB->m_num,
                          connC->m_n, connC->m_num,
                          gemm->m_transA, gemm->m_transB,
                          true,
                          gemm->m_alpha, gemm->m_beta, m_toLayer, gemm->m_type);
      break;
    case(2):
      loop = GemmVar2Loop(connA->m_n, connA->m_num,
                          connB->m_n, connB->m_num,
                          connC->m_n, connC->m_num,
                          gemm->m_transA, gemm->m_transB,
                          gemm->m_alpha, gemm->m_beta, m_toLayer, gemm->m_type);
      break;
    default:
      throw;
  }
  
  poss->AddLoop(loop);
  
  node->RedirectChildren(loop->OutTun(2),0);
  node->m_poss->DeleteChildAndCleanUp(node);
}




bool DistGemmToLocalGemmStatC::CanApply(const Poss *poss, const Node *node) const
{
  return IsDMGemm(node);
}

void DistGemmToLocalGemmStatC::Apply(Poss *poss, Node *node) const
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
  poss->AddNode(node1);
  poss->AddNode(node2);
  poss->AddNode(node3);
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


bool DistGemmToContribLocalGemmStatANonTrans::CanApply(const Poss *poss, const Node *node) const
{
  if (IsDMGemm(node) &&
      (((Gemm*)node)->m_transA == NORMAL)) {
    return true;
  }
  return false;
}

void DistGemmToContribLocalGemmStatANonTrans::Apply(Poss *poss, Node *node) const
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
  
  poss->AddNode(redistB);
  poss->AddNode(localGemm);
  poss->AddNode(node3);
  poss->AddNode(tmp);
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
bool DistGemmToContribLocalGemmStatBNonTrans::CanApply(const Poss *poss, const Node *node) const
{
  if (IsDMGemm(node) &&
      (((Gemm*)node)->m_transB != CONJTRANS) &&
      (((Gemm*)node)->m_transB != TRANS)) {
    return true;
  }
  return false;
}

void DistGemmToContribLocalGemmStatBNonTrans::Apply(Poss *poss, Node *node) const
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
      poss->AddNode(redist);
    }
    node2->AddInput(tmp,0);
    
    node3->AddInput(node2,0);
    node3->AddInput(node->Input(2),node->InputConnNum(2));
    
    poss->AddNode(node1);
    poss->AddNode(node2);
    poss->AddNode(node3);
    poss->AddNode(tmp);
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
      poss->AddNode(redist);
    }
    node2->AddInput(node1,0);
    node2->AddInput(tmp,0);
    
    node3->AddInput(node2,0);
    node3->AddInput(node->Input(2),node->InputConnNum(2));
    
    poss->AddNode(node1);
    poss->AddNode(node2);
    poss->AddNode(node3);
    poss->AddNode(tmp);
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


bool DistGemmToContribLocalGemmStatBTrans::CanApply(const Poss *poss, const Node *node) const
{
  if (IsDMGemm(node) &&
      (((Gemm*)node)->m_transB != NORMAL)) {
    return true;
  }
  return false;
}

void DistGemmToContribLocalGemmStatBTrans::Apply(Poss *poss, Node *node) const
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
  
  poss->AddNodes(7,
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

bool DistGemmToContribLocalGemmStatATrans::CanApply(const Poss *poss, const Node *node) const
{
  if (IsDMGemm(node) &&
      (((Gemm*)node)->m_transA != NORMAL)) {
    return true;
  }
  return false;
}

void DistGemmToContribLocalGemmStatATrans::Apply(Poss *poss, Node *node) const
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
  
  poss->AddNodes(7,
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

bool GemmTrans::CanApply(const Poss *poss, const Node *node) const
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

bool GemmInputReordering::CanApply(const Poss *poss, const Node *node) const
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

void GemmInputReordering::Apply(Poss *poss, Node *node) const
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

Loop* GemmVar1Loop(Node *Ain, unsigned int Anum,
                   Node *Bin, unsigned int Bnum,
                   Node *Cin, unsigned int Cnum,
                   Trans transA, Trans transB,
                   Coef alpha, Coef beta,
                   Layer layer, Type type)
{
  Split *splitA = new Split(transA==NORMAL ? PARTDOWN : PARTRIGHT, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  LoopTunnel *Btun = new LoopTunnel(POSSTUNIN);
  Btun->AddInput(Bin, Bnum);
  Btun->SetAllStats(FULLUP);
  Btun->SetIndepIters();
  
  Split *splitC = new Split(PARTDOWN, POSSTUNIN, true);
  splitC->AddInput(Cin, Cnum);
  splitC->SetUpStats(FULLUP, FULLUP,
                     PARTUP, PARTUP);
  splitC->SetIndepIters();
  
  Node *gepp;
  
  gepp = new Gemm(layer, transA, transB, alpha, beta, type);
  gepp->AddInput(splitA, 1);
  gepp->AddInput(Btun, 0);
  gepp->AddInput(splitC, 1);
  
  Combine *comA = splitA->CreateMatchingCombine(0);
  
  LoopTunnel *BtunOut = new LoopTunnel(POSSTUNOUT);
  BtunOut->AddInput(Btun, 0);
  BtunOut->AddInput(Btun, 0);
  BtunOut->CopyTunnelInfo(Btun);
  
  Combine *comC = splitC->CreateMatchingCombine(1,
                                                1, gepp, 0);
  
  Poss *loopPoss = new Poss(3, comA, BtunOut, comC);
  Loop *loop;
  if (layer == DMLAYER)
    loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  else
    loop = new Loop(BLISLOOP, loopPoss, USEBLISMC);
  
  return loop;
}


Loop* GemmVar3Loop(Node *Ain, unsigned int Anum,
                   Node *Bin, unsigned int Bnum,
                   Node *Cin, unsigned int Cnum,
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
  Split *splitA = new Split(aDir, POSSTUNIN, true);
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
  Split *splitB = new Split(bDir, POSSTUNIN);
  splitB->AddInput(Bin, Bnum);
  splitB->SetUpStats(FULLUP, FULLUP,
                     FULLUP, FULLUP);
  splitB->SetIndepIters();
  
  ScaleNode *scale = new ScaleNode(layer, beta);
  scale->AddInput(Cin, Cnum);
  
  LoopTunnel *Ctun = new LoopTunnel(POSSTUNIN);
  Ctun->AddInput(scale, 0);
  Ctun->SetUpStats(PARTUP, PARTUP,
                   PARTUP, PARTUP);
  
  Node *gepp;
  
  gepp = new Gemm(layer, transA, transB, alpha, beta, type);
  gepp->AddInput(splitA, 1);
  gepp->AddInput(splitB, 1);
  gepp->AddInput(Ctun, 0);
  
  Combine *comA = splitA->CreateMatchingCombine(0);
  
  Combine *comB = splitB->CreateMatchingCombine(0);
  
  LoopTunnel *CtunOut = new LoopTunnel(POSSTUNOUT);
  CtunOut->AddInput(gepp, 0);
  CtunOut->AddInput(Ctun, 0);
  CtunOut->CopyTunnelInfo(Ctun);
  
  Poss *loopPoss = new Poss(3, comA, comB, CtunOut);
  Loop *loop;
  if (layer == DMLAYER)
    loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  else
    loop = new Loop(BLISLOOP, loopPoss, USEBLISKC);
  
  return loop;
}


Loop* GemmVar2Loop(Node *Ain, unsigned int Anum,
                   Node *Bin, unsigned int Bnum,
                   Node *Cin, unsigned int Cnum,
                   Trans transA, Trans transB,
                   Coef alpha, Coef beta,
                   Layer layer, Type type)
{
  LoopTunnel *Atun = new LoopTunnel(POSSTUNIN);
  Atun->AddInput(Ain, Anum);
  Atun->SetAllStats(FULLUP);
  Atun->SetIndepIters();
  
  Split *splitB = new Split(transB==NORMAL ? PARTRIGHT : PARTDOWN, POSSTUNIN);
  splitB->AddInput(Bin, Bnum);
  splitB->SetAllStats(FULLUP);
  splitB->SetIndepIters();
  
  Split *splitC = new Split(PARTRIGHT, POSSTUNIN, true);
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
  AtunOut->AddInput(Atun, 0);
  AtunOut->CopyTunnelInfo(Atun);
  
  Combine *comB = splitB->CreateMatchingCombine(0);
  
  Combine *comC = splitC->CreateMatchingCombine(1,
                                                1, gepp, 0);
  
  Poss *loopPoss = new Poss(3, AtunOut, comB, comC);
  Loop *loop;
  if (layer == DMLAYER)
    loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  else
    loop = new Loop(BLISLOOP, loopPoss, USEBLISNC);
  
  return loop;
}

string BLISGemmLoopExp::GetType() const
{
  return "BLISGemmLoopExp";
}

bool BLISGemmLoopExp::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == Gemm::GetClass()) {
    const Gemm *gemm = (Gemm*)node;
    if (gemm->GetLayer() == m_fromLayer)
      return true;
  }
  return false;
}

void BLISGemmLoopExp::Apply(Poss *poss, Node *node) const
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
  
  Split *splitA = new Split(gemm->m_transA == NORMAL ? PARTDOWN : PARTRIGHT, POSSTUNIN);
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
  
  poss->AddNode(bBuff);
  poss->AddNode(bPack);
  
  LoopTunnel *Btun = new LoopTunnel(POSSTUNIN);
  Btun->AddInput(bPack, 0);
  Btun->SetAllStats(FULLUP);
  Btun->SetIndepIters();
  
  Split *splitC = new Split(PARTDOWN, POSSTUNIN, true);
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
  Loop *loop = new Loop(BLISLOOP, loopPoss, USEBLISMC);
  
  
  
  poss->AddLoop(loop);
  
  node->RedirectChildren(loop->OutTun(2), 0);
  
  node->m_poss->DeleteChildAndCleanUp(node);
}


bool GemmLowerLayer::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == Gemm::GetClass()) {
    const Gemm *gemm = (Gemm*)node;
    if (gemm->GetLayer() != m_fromLayer)
      return false;
    if (m_dim == DIMK) {
      if (gemm->m_transA == NORMAL) {
        return (*(gemm->InputLocalN(0)) <= m_bs);
      }
      else {
        return (*(gemm->InputLocalM(0)) <= m_bs);
      }
    }
    else if (m_dim == DIMN) {
      if (gemm->m_transB == NORMAL) {
        return (*(gemm->InputLocalN(1)) <= m_bs);
      }
      else {
        return (*(gemm->InputLocalM(1)) <= m_bs);
      }
    }
    else
      throw;
  }
  return false;
  
}

void GemmLowerLayer::Apply(Poss *poss, Node *node) const
{
  Gemm *gemm = (Gemm*)node;
  gemm->SetLayer(m_toLayer);
}

string GemmLowerLayer::GetType() const
{ 
  return "Gemm lower layer " + LayerNumToStr(m_fromLayer) 
  + " to " + LayerNumToStr(m_toLayer);
}
