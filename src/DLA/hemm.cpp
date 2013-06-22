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



#include "hemm.h"
#include "blas.h"
#include "distributions.h"
#include "string.h"
#include "helperNodes.h"
#include "loopSupport.h"
#include "pack.h"

using namespace std;

Hemm::Hemm(Layer layer, Side side, Tri tri, Coef alpha, Coef beta, Type type) 
  : m_side(side), 
    m_tri(tri), 
    m_alpha(alpha), 
    m_beta(beta), 
    m_type(type) 
{
  SetLayer(layer);
}


void Hemm::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  const Hemm *hemm = (Hemm*)orig;
  m_side = hemm->m_side;
  m_tri = hemm->m_tri;
  m_alpha = hemm->m_alpha;
  m_beta = hemm->m_beta;
  m_type = hemm->m_type;
}

void Hemm::FlattenCore(ofstream &out) const
{
  DLAOp<3,1>::FlattenCore(out);
  WRITE(m_side);
  WRITE(m_tri);
  WRITE(m_alpha);
  WRITE(m_beta);
  WRITE(m_type);
}

void Hemm::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<3,1>::UnflattenCore(in, info);
  READ(m_side);
  READ(m_tri);
  READ(m_alpha);
  READ(m_beta);
  READ(m_type);
}

NodeType Hemm::GetType() const
{
  string str = "Hemm " + LayerNumToStr(GetLayer())
    + SideToStr(m_side) + " " + TriToStr(m_tri);
  return str;
}

DistType Hemm::GetDistType(unsigned int num) const
{ 
  switch(GetLayer()) {
  case (ABSLAYER):
  case (DMLAYER):
    return D_MC_MR;
  case (SMLAYER):
  case (S1LAYER):
  case (S2LAYER):
    return InputDistType(2); 
      
  default:
    throw;
  } 
}

Phase Hemm::MaxPhase() const 
{  
  switch(GetLayer()) {
#if DODPPHASE
  case (ABSLAYER):
  case (DMLAYER):
    return DPPHASE;
    case (SMLAYER):
      return NUMPHASES;
  default:
    throw;
#elif DOSR1PHASE
    case (ABSLAYER):
      return SR1PHASE;
    case (S1LAYER):
      return SR2PHASE;
    default:
      throw;
#endif
  }
}

void Hemm::SanityCheck()
{
  DLAOp<3,1>::SanityCheck();
  if (GetLayer() == ABSLAYER || GetLayer() == DMLAYER) {
    if (InputDistType(2) != D_MC_MR) {
      cout << "input not D_MC_MR 7";
      throw;
    }
  }
  else if (GetLayer() == SMLAYER) {
    DistType t0 = InputDistType(0);
    DistType t1 = InputDistType(1);
    if (t0 != D_STAR_STAR)
      m_poss->MarkInsane();
    if ((m_side == RIGHT && t1 != D_VC_STAR && t1 != D_VR_STAR && t1 != D_STAR_STAR && t1 != D_MC_STAR && t1 != D_MR_STAR)
	|| 
	(m_side == LEFT && t1 != D_STAR_VC && t1 != D_STAR_VR && t1 != D_STAR_STAR && t1 != D_STAR_MC && t1 != D_STAR_MR))
      m_poss->MarkInsane();
    if (t1 != InputDistType(2))
      m_poss->MarkInsane();
  }
  else if (GetLayer() == S1LAYER || GetLayer() == S2LAYER) {
  }
  else
    throw;
}

bool Hemm::ShouldCullDP() const 
{
  switch (GetLayer()) {
  case (ABSLAYER):
  case (DMLAYER):
    return true;
  case (SMLAYER):
    return false;
  default:
    throw;
  }
}

void Hemm::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3,1>::Prop();
    if (GetLayer() == ABSLAYER || GetLayer() == DMLAYER)
      m_cost = ZERO;
    else if (GetLayer() == SMLAYER)
      m_cost = GetCost(SMLAYER, m_side,LocalM(0),LocalN(0));
    else if (GetLayer() == S1LAYER || GetLayer() == S2LAYER)
      m_cost = 0;
    else
      throw;
  }
}

Cost Hemm::GetCost(Layer layer, Side side, const Sizes *localMs, const Sizes *localNs)
{
  if (layer != SMLAYER)
    throw;
  if (side == LEFT)
    return TWO * GAMMA * localMs->SumProds21(*localNs);
  else
    return TWO * GAMMA * localNs->SumProds21(*localMs);
}

void Hemm::PrintCode(IndStream &out)
{
  out.Indent();
  if (GetLayer() == ABSLAYER)
    *out << "AbsHemm( ";
  else if (GetLayer() == DMLAYER)
    *out << "DistHemm( ";
  else if (GetLayer() == SMLAYER)
    *out << "blas::Hemm( ";
  else
    throw;
  *out << SideToStr(m_side) << ", " << TriToStr(m_tri) 
       << ", \n" << out.Tabs(1);
  out << m_alpha;
  *out << ","
       << GetInputName(0).str() << ", " 
       << GetInputName(1).str() << ", " << endl
       << out.Tabs(1);
  out << m_beta;
  *out << ", " << GetInputName(2).str() << " );\n";
}

string HemmLoopExp::GetType() const
{
  switch(m_var) {
  case(-8):
      return "Hemm Loop Exp var 8 (altered)";
    case(4):
      return "Hemm Loop Exp var 4";
    case(8):
      return "Hemm Loop Exp var 8";
    default:
      throw;    
  }
}

bool HemmLoopExp::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == Hemm::GetClass()
      && ((Hemm*)node)->GetLayer() == m_fromLayer)
    {
      return true;
    }
  else
    return false;
}

void HemmLoopExp::Apply(Poss *poss, Node *node) const
{
  Hemm *hemm = (Hemm*)node;
  Loop *loop;
  
  NodeConn *connA, *connB, *connC;
  connA = hemm->m_inputs[0];
  connB = hemm->m_inputs[1];
  connC = hemm->m_inputs[2];
  
  switch(m_var) {
    case(-8):
      loop = HemmLoopVar8Altered(connA->m_n, connA->m_num,
				 connB->m_n, connB->m_num,
				 connC->m_n, connC->m_num,
				 hemm->m_side, hemm->m_tri,
				 hemm->m_alpha, hemm->m_beta,
				 hemm->m_type,
				 m_toLayer);
      break;
    case(4):
      loop = HemmLoopVar4(connA->m_n, connA->m_num,
			  connB->m_n, connB->m_num,
			  connC->m_n, connC->m_num,
			  hemm->m_side, hemm->m_tri,
			  hemm->m_alpha, hemm->m_beta,
			  hemm->m_type,
			  m_toLayer);
      break;
    case(8):
      loop = HemmLoopVar8(connA->m_n, connA->m_num,
			  connB->m_n, connB->m_num,
			  connC->m_n, connC->m_num,
			  hemm->m_side, hemm->m_tri,
			  hemm->m_alpha, hemm->m_beta,
			  hemm->m_type,
			  m_toLayer);
      break;
    default:
      throw;
  }
  
  poss->AddLoop(loop);
  
  node->RedirectChildren(loop->OutTun(2),0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

bool DistHemmToLocalHemm::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == Hemm::GetClass()
      && ((Hemm*)node)->GetLayer() == DMLAYER) 
    {
      return true;
    }
  else 
    return false;
}

void DistHemmToLocalHemm::Apply(Poss *poss, Node *node) const
{
  Hemm *orig = (Hemm*)node;
  bool left = orig->m_side == LEFT;
  RedistNode *redist1 = new RedistNode(D_STAR_STAR);
  redist1->AddInput(node->Input(0),node->InputConnNum(0));
  RedistNode *redist2 = new RedistNode(left ? m_leftType : m_rightType);
  redist2->AddInput(node->Input(1),node->InputConnNum(1));
  RedistNode *redist3 = new RedistNode(left ? m_leftType : m_rightType);
  redist3->AddInput(node->Input(2),node->InputConnNum(2));
  Hemm *hemm = new Hemm(SMLAYER, orig->m_side, orig->m_tri, orig->m_alpha, orig->m_beta,orig->m_type);
  hemm->AddInput(redist1,0);
  hemm->AddInput(redist2,0);
  hemm->AddInput(redist3,0);
  RedistNode *redist4 = new RedistNode(D_MC_MR);
  redist4->AddInput(hemm,0);
  poss->AddNodes(5, redist1, redist2, redist3, hemm, redist4);
  node->RedirectChildren(redist4,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}


Cost DistHemmToLocalHemm::RHSCostEstimate(const Node *node) const
{
  Hemm *hemm = (Hemm*)node;
  bool left = hemm->m_side == LEFT;
  const Sizes *A1 = hemm->GetInputM(0);
  const Sizes *A2 = hemm->GetInputN(0);
  const Sizes *B1 = hemm->GetInputM(1);
  const Sizes *B2 = hemm->GetInputN(1);
  const Sizes *C1 = hemm->GetInputM(2);
  const Sizes *C2 = hemm->GetInputN(2);
  Cost cost = RedistNode::GetCost(D_MC_MR, D_STAR_STAR, A1, A2);
  DistType type = left ? m_leftType : m_rightType;
  cost += RedistNode::GetCost(D_MC_MR, type, B1, B2);
  cost += RedistNode::GetCost(D_MC_MR, type, C1, C2);
  Sizes localC1, localC2;
  GetLocalSizes(type, C1, C2, localC1, localC2);
  cost += Hemm::GetCost(SMLAYER, hemm->m_side, &localC1, &localC2);
  cost += RedistNode::GetCost(type, D_MC_MR, C1, C2);
  return cost;
}


Loop* HemmLoopVar4(Node *Ain, unsigned int Anum,
		   Node *Bin, unsigned int Bnum,
		   Node *Cin, unsigned int Cnum,
		   Side side, Tri tri,
		   Coef alpha, Coef beta,
		   Type type,
		   Layer layer)
{
  
  LoopTunnel *Atun = new LoopTunnel(POSSTUNIN);
  Atun->AddInput(Ain, Anum);
  Atun->SetAllStats(FULLUP);
  Atun->SetIndepIters();
  
  Split *splitB = new Split(side==LEFT ? PARTRIGHT : PARTDOWN, POSSTUNIN);
  splitB->AddInput(Bin, Bnum);
  splitB->SetAllStats(FULLUP);
  splitB->SetIndepIters();
  
  Split *splitC = new Split(side==LEFT ? PARTRIGHT : PARTDOWN, POSSTUNIN, true);
  splitC->AddInput(Cin, Cnum);
  if (side == LEFT) {
    splitC->SetUpStats(FULLUP, NOTUP,
		       FULLUP, NOTUP);
  }
  else {
    splitC->SetUpStats(FULLUP, FULLUP,
		       NOTUP, NOTUP);
  }
  splitC->SetIndepIters();

  Node *hemm;
  
  hemm = new Hemm(layer, side, tri, alpha, beta, type);
  hemm->AddInput(Atun, 0);
  hemm->AddInput(splitB, 1);
  hemm->AddInput(splitC, 1);
  
  
  LoopTunnel *AtunOut = new LoopTunnel(POSSTUNOUT);
  AtunOut->AddInput(Atun, 0);
  AtunOut->AddInput(Atun, 0);
  AtunOut->CopyTunnelInfo(Atun);
  
  Combine *comB = splitB->CreateMatchingCombine(0);
  
  Combine *comC = splitC->CreateMatchingCombine(1,
                                                1, hemm, 0);
  
  Poss *loopPoss = new Poss(3, AtunOut, comB, comC);
  Loop *loop;
  if (layer == DMLAYER)
    loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  else
    loop = new Loop(BLISLOOP, loopPoss, USEBLISKC);
  
  return loop;
}


Loop* HemmLoopVar8(Node *Ain, unsigned int Anum,
		   Node *Bin, unsigned int Bnum,
		   Node *Cin, unsigned int Cnum,
		   Side side, Tri tri,
		   Coef alpha, Coef beta,
		   Type type,
		   Layer layer)
{
  Split *splitA = new Split(PARTDIAG, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  Split *splitB = new Split(side==LEFT ? PARTDOWN : PARTRIGHT, POSSTUNIN);
  splitB->AddInput(Bin, Bnum);
  splitB->SetAllStats(FULLUP);
  splitB->SetIndepIters();

  ScaleNode *scale = new ScaleNode(layer, beta);
  scale->AddInput(Cin, Cnum);
  
  Split *splitC = new Split(side==LEFT ? PARTDOWN : PARTRIGHT, POSSTUNIN, true);
  splitC->AddInput(scale, 0);
  splitC->SetAllStats(PARTUP);

  Node *gemm1;
  if ((side == LEFT && tri == UPPER) || (side == RIGHT && tri == LOWER)) {
    gemm1 = new Gemm(layer, NORMAL, NORMAL, alpha, COEFONE, type);
  }
  else if (side == LEFT && tri == LOWER) {
    gemm1 = new Gemm(layer, type == COMPLEX ? CONJTRANS : TRANS, NORMAL, alpha, COEFONE, type);
  }
  else if (side == RIGHT && tri == UPPER) {
    gemm1 = new Gemm(layer, NORMAL, type == COMPLEX ? CONJTRANS : TRANS, alpha, COEFONE, type);
  }
  else
    throw;

  if (side == LEFT) {
    if (tri == UPPER) {
      gemm1->AddInput(splitA, 3);
    }
    else {
      gemm1->AddInput(splitA, 1);
    }
    gemm1->AddInput(splitB, 1);
  }
  else {
    gemm1->AddInput(splitB, 1);
    if (tri == UPPER) {
      gemm1->AddInput(splitA, 3);
    }
    else {
      gemm1->AddInput(splitA, 1);
    }
  }
  gemm1->AddInput(splitC, 0);

  Node *hemm;
  
  hemm = new Hemm(layer, side, tri, alpha, COEFONE, type);
  hemm->AddInput(splitA, 4);
  hemm->AddInput(splitB, 1);
  hemm->AddInput(splitC, 1);

  Node *gemm2;

  if ((side == RIGHT && tri == UPPER) || (side == LEFT && tri == LOWER)) {
    gemm2 = new Gemm(layer, NORMAL, NORMAL, alpha, COEFONE, type);
  }
  else if (side == LEFT && tri == UPPER) {
    gemm2 = new Gemm(layer, type == COMPLEX ? CONJTRANS : TRANS, NORMAL, alpha, COEFONE, type);
  }
  else if (side == RIGHT && tri == LOWER) {
    gemm2 = new Gemm(layer, NORMAL, type == COMPLEX ? CONJTRANS : TRANS, alpha, COEFONE, type);
  }
  else
    throw;


  if (side == LEFT) {
    if (tri == UPPER) {
      gemm2->AddInput(splitA, 7);
    }
    else {
      gemm2->AddInput(splitA, 5);
    }
    gemm2->AddInput(splitB, 1);
  }
  else {
    gemm2->AddInput(splitB, 1);
    if (tri == UPPER) {
      gemm2->AddInput(splitA, 7);
    }
    else {
      gemm2->AddInput(splitA, 5);
    }
  }
  gemm2->AddInput(splitC, 2);

  Combine *comA = splitA->CreateMatchingCombine(0);
  
  Combine *comB = splitB->CreateMatchingCombine(0);
  
  Combine *comC = splitC->CreateMatchingCombine(3,
                                                0, gemm1, 0,
						1, hemm, 0,
						2, gemm2, 0);
  
  Poss *loopPoss = new Poss(3, comA, comB, comC);
  Loop *loop;
  if (layer == DMLAYER)
    loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  else
    loop = new Loop(BLISLOOP, loopPoss, USEBLISKC);
  
  return loop;
}



Loop* HemmLoopVar8Altered(Node *Ain, unsigned int Anum,
		   Node *Bin, unsigned int Bnum,
		   Node *Cin, unsigned int Cnum,
		   Side side, Tri tri,
		   Coef alpha, Coef beta,
		   Type type,
			  Layer layer)
{  
  Split *splitA = new Split(PARTDIAG, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  Split *splitB = new Split(side==LEFT ? PARTDOWN : PARTRIGHT, POSSTUNIN);
  splitB->AddInput(Bin, Bnum);
  splitB->SetAllStats(FULLUP);
  splitB->SetIndepIters();

  ScaleNode *scale = new ScaleNode(layer, beta);
  scale->AddInput(Cin, Cnum);
  
  Split *splitC = new Split(side==LEFT ? PARTDOWN : PARTRIGHT, POSSTUNIN, true);
  splitC->AddInput(scale, 0);
  splitC->SetAllStats(PARTUP);

  ViewPan *ARow, *ACol;
  ARow = new ViewPan(false, "ARowPan");
  ACol = new ViewPan(true, "AColPan");

  if (tri == LOWER) {
    ARow->AddInputs(4,
		    splitA, 1,
		    splitA, 4);
    ACol->AddInputs(4,
		    splitA, 4,
		    splitA, 5);
  }
  else {
    ARow->AddInputs(4,
		    splitA, 4,
		    splitA, 7);
    ACol->AddInputs(4,
		    splitA, 3,
		    splitA, 4);
  }
  
  ViewAroundDiag *view;
  view = new ViewAroundDiag(side == LEFT ? true: false, "C");

  view->AddInputs(6,
		  splitC, 0,
		  splitC, 1,
		  splitC, 2);

  ViewAroundDiagCombine *com = new ViewAroundDiagCombine(side == LEFT ? true: false);

  MakeTrapNode *trapACol, *trapARow;
  if ((side == LEFT && tri == LOWER) || (side == RIGHT && tri == UPPER)) {
    trapACol = new MakeTrapNode(LEFT, LOWER, 0);
    trapARow = new MakeTrapNode(RIGHT, LOWER, -1);
  }
  else {
    trapACol = new MakeTrapNode(RIGHT, UPPER, 0);
    trapARow = new MakeTrapNode(LEFT, UPPER, 1);
  }

  trapACol->AddInput(ACol, 0);
  trapARow->AddInput(ARow, 0);

  Node *gemm1, *gemm2;
  gemm1 =  new Gemm(layer, NORMAL, NORMAL, alpha, COEFONE, type);
  if (side == LEFT)
    gemm2 = new Gemm(layer, type == COMPLEX? CONJTRANS : TRANS, NORMAL, alpha, COEFONE, type);
  else
    gemm2 = new Gemm(layer, NORMAL, type == COMPLEX? CONJTRANS : TRANS, alpha, COEFONE, type);

  if (side == LEFT) {
    if (tri == LOWER) {
      gemm1->AddInputs(6,
		      trapACol, 0,
		      splitB, 1, 
		       view, 1);
      gemm2->AddInputs(6,
		       trapARow, 0,
		       splitB, 1,
		       view, 0);
      com->AddInputs(4,
		      gemm2, 0,
		      gemm1, 0);
    }
    else {
      gemm1->AddInputs(6,
		       trapACol, 0,
		       splitB, 1,
		       view, 0);
      gemm2->AddInputs(6,
		       trapARow, 0,
		       splitB, 1,
		       view, 1);
      com->AddInputs(4,
		      gemm1, 0,
		      gemm2, 0);
    }
  }
  else {
    if (tri == LOWER) {
      gemm1->AddInputs(6,
		       splitB, 1,
		       trapARow, 0,
		       view, 0);
      gemm2->AddInputs(6,
		       splitB, 1,
		       trapACol, 0,
		       view, 1);
      com->AddInputs(4,
		      gemm1, 0,
		      gemm2, 0);

    }
    else {
      gemm1->AddInputs(6,
		       splitB, 1,
		       trapARow, 0,
		       view, 1);
      gemm2->AddInputs(6,
		       splitB, 1,
		       trapACol, 0,
		       view, 0);
      com->AddInputs(4,
		      gemm2, 0,
		      gemm1, 0);
    }
  }

  com->AddInputs(6,
		 splitC, 0,
		 splitC, 1,
		 splitC, 2);


  Combine *comA = splitA->CreateMatchingCombine(0);
  
  Combine *comB = splitB->CreateMatchingCombine(0);
  
  Combine *comC = splitC->CreateMatchingCombine(3,
                                                0, com, 0,
						1, com, 1,
						2, com, 2);
  
  Poss *loopPoss = new Poss(3, comA, comB, comC);
  Loop *loop;
  if (layer == DMLAYER)
    loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  else
    loop = new Loop(BLISLOOP, loopPoss, USEBLISKC);
  
  return loop;
}


LocalSymmAcc::LocalSymmAcc(Side side, Tri tri, Type type, Coef alpha)
  :m_side(side), m_tri(tri), m_type(type), m_alpha(alpha)
{
  
}

void LocalSymmAcc::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<5,2>::Duplicate(orig,shallow,possMerging);
  const LocalSymmAcc *acc = (LocalSymmAcc*)orig;
  m_side = acc->m_side;
  m_tri = acc->m_tri;
  m_type = acc->m_type;
}

void LocalSymmAcc::SanityCheck()
{
  DLAOp<5,2>::SanityCheck();
  DistType t0 = InputDistType(0);
  DistType t1 = InputDistType(1);
  DistType t2 = InputDistType(2);
  DistType t3 = InputDistType(3);
  DistType t4 = InputDistType(4);
  if (m_tri == LOWER) {
    if (t0 != D_MC_MR)
      throw;
    if (t1 != D_MC_STAR)
      throw;
    if ((t2 != D_STAR_MR_T && m_type == REAL) || (t2 != D_STAR_MR_H && m_type == COMPLEX))
      throw;
    if (t3 != D_MC_STAR)
      throw;
    if (t4 != D_MR_STAR)
      throw;
  }
  else {
    if (t0 != D_MC_MR)
      throw;
    if (t1 != D_STAR_MC)
      throw;
    if ((t2 != D_MR_STAR_T && m_type == REAL) || (t2 != D_MR_STAR_H && m_type == COMPLEX))
      throw;
    if ((t3 != D_MC_STAR_T && m_type == REAL) || (t3 != D_MC_STAR_H && m_type == COMPLEX))
      throw;
    if ((t4 != D_MR_STAR_T && m_type == REAL) || (t4 != D_MR_STAR_H && m_type == COMPLEX))
      throw;
  }
}

void LocalSymmAcc::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<5,2>::Prop();
    m_cost = GetCost(m_side, InputLocalM(1), InputLocalN(1));
  }
}

Cost LocalSymmAcc::GetCost(Side side, const Sizes *localMs, const Sizes *localNs)
{
  if (side == LEFT)
    return TWO * GAMMA * localMs->SumProds21(*localNs);
  else
    return TWO * GAMMA * localNs->SumProds21(*localMs);
}


void LocalSymmAcc::PrintCode(IndStream &out)
{
  out.Indent();
  *out << "LocalSymmetricAccumulate";

  if (m_side == LEFT)
    *out << "L";
  else
    *out << "R";

  if (m_tri == LOWER)
    *out << "L";
  else
    *out << "U";

  if (m_type == REAL)
    *out << "( TRANSPOSE, ";
  else
    *out << "( ADJOINT, ";
  
  *out << endl << out.Tabs(1);
 
  out << m_alpha;

  for (unsigned int i = 0; i < 5; ++i)
    *out << ", " << GetInputNameStr(i);

  *out << " );\n";
}

void LocalSymmAcc::FlattenCore(ofstream &out) const
{
  DLAOp<5,2>::FlattenCore(out);
  WRITE(m_side);
  WRITE(m_tri);
  WRITE(m_type);
  WRITE(m_alpha);
}


void LocalSymmAcc::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<5,2>::UnflattenCore(in, info);
  READ(m_side);
  READ(m_tri);
  READ(m_type);
  READ(m_alpha);
}


bool DistHemmToLocalHemmStatA::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == Hemm::GetClass()
      && ((Hemm*)node)->GetLayer() == DMLAYER) 
    {
      return true;
    }
  else
    return false;
}

void DistHemmToLocalHemmStatA::Apply(Poss *poss, Node *node) const
{
  Hemm *orig = (Hemm*)node;
  Type type = orig->m_type;
  bool complex = type == COMPLEX;
  
  NodeConn *Aconn = orig->m_inputs[0];
  NodeConn *Bconn = orig->m_inputs[1];
  NodeConn *Cconn = orig->m_inputs[2];

  if (orig->m_side == LEFT) {
    RedistNode *redist1 = new RedistNode(D_MC_STAR);
    redist1->AddInput(Bconn->m_n, Bconn->m_num);

    RedistNode *redist2 = new RedistNode(D_VR_STAR);
    redist2->AddInput(Bconn->m_n, Bconn->m_num);

    RedistNode *redist3 = new RedistNode(complex ? D_STAR_MR_H : D_STAR_MR_T);
    redist3->AddInput(redist2, 0);    


    TempVarNode *tmp1 = new TempVarNode(D_MC_STAR);
    tmp1->AddInput(Cconn->m_n, Cconn->m_num);
    TempVarNode *tmp2 = new TempVarNode(D_MR_STAR);
    tmp2->AddInput(Cconn->m_n, Cconn->m_num);

    LocalSymmAcc *acc = new LocalSymmAcc(LEFT, orig->m_tri, type, orig->m_alpha);

    acc->AddInputs(10, 
		   Aconn->m_n, Aconn->m_num,
		   redist1, 0,
		   redist3, 0,
		   tmp1, 0,
		   tmp2, 0);

    SumScatterFrom *from = new SumScatterFrom;
    TempVarNode *tmp3 = new TempVarNode(D_MR_MC);
    tmp3->AddInput(Cconn->m_n, Cconn->m_num);
    from->AddInputs(4,
		    acc, 1,
		    tmp3, 0);
    RedistNode *redist4 = new RedistNode(D_MC_MR);
    redist4->AddInput(from, 0);
    SumScatterNode *up = new SumScatterNode(COEFONE);
    up->AddInputs(4,
		  acc, 0,
		  redist4, 0);

    Axpy *axpy = new Axpy(SMLAYER, COEFONE);
    axpy->AddInputs(4,
		    up, 0,
		    Cconn->m_n, Cconn->m_num);

    poss->AddNodes(11, redist1, redist2, redist3, tmp1, tmp2, acc,
		   from, tmp3, redist4, up, axpy);
    node->RedirectChildren(axpy,0);
    node->m_poss->DeleteChildAndCleanUp(node);
  }
  else {
    RedistNode *redist1 = new RedistNode(complex ? D_MR_STAR_H : D_MR_STAR_T);
    redist1->AddInput(Bconn->m_n, Bconn->m_num);

    RedistNode *redist2 = new RedistNode(D_STAR_MC);
    redist2->AddInput(Bconn->m_n, Bconn->m_num);
    
    TempVarNode *tmp1 = new TempVarNode(complex ? D_MC_STAR_H : D_MC_STAR_T);
    tmp1->AddInput(Cconn->m_n, Cconn->m_num);

    TempVarNode *tmp2 = new TempVarNode(complex ? D_MR_STAR_H : D_MR_STAR_T);
    tmp2->AddInput(Cconn->m_n, Cconn->m_num);

    LocalSymmAcc *acc = new LocalSymmAcc(RIGHT, orig->m_tri, type, orig->m_alpha);

    acc->AddInputs(10, 
		   Aconn->m_n, Aconn->m_num,
		   redist2, 0,
		   redist1, 0,
		   tmp1, 0,
		   tmp2, 0);

    TempVarNode *tmp3 = new TempVarNode(complex ? D_MC_MR_H : D_MC_MR_T);
    tmp3->AddInput(Cconn->m_n, Cconn->m_num);

    SumScatterFrom *from = new SumScatterFrom;
    from->AddInputs(4,
		    acc, 0,
		    tmp3, 0);

    RedistNode *redist3 = new RedistNode(complex ? D_MR_MC_H : D_MR_MC_T);
    redist3->AddInput(from, 0);

    SumScatterNode *up = new SumScatterNode(COEFONE);
    up->AddInputs(4,
		  acc, 1,
		  redist3, 0);

    RedistNode *redist4 = new RedistNode(D_MC_MR);
    redist4->AddInput(up, 0);

    Axpy *axpy = new Axpy(SMLAYER, COEFONE);
    axpy->AddInputs(4,
		    redist4, 0,
		    Cconn->m_n, Cconn->m_num);

    poss->AddNodes(11, redist1, redist2, tmp1, tmp2, acc,
		   from, tmp3, redist3, up, redist4, axpy);
    node->RedirectChildren(axpy,0);
    node->m_poss->DeleteChildAndCleanUp(node);
  }
}

Cost DistHemmToLocalHemmStatA::RHSCostEstimate(const Node *node) const
{
  Hemm *orig = (Hemm*)node;
  Type type = orig->m_type;
  bool complex = type == COMPLEX;
  const Sizes *B1 = orig->GetInputM(1);
  const Sizes *B2 = orig->GetInputN(1);
  const Sizes *C1 = orig->GetInputM(2);
  const Sizes *C2 = orig->GetInputN(2);
  if (orig->m_side == LEFT) {
    Cost cost = RedistNode::GetCost(D_MC_MR, D_MC_STAR, B1, B2);
    cost += RedistNode::GetCost(D_MC_MR, D_VR_STAR, B1, B2);
    cost += RedistNode::GetCost(D_VR_STAR, complex ? D_STAR_MR_H : D_STAR_MR_T, B1, B2);
    Sizes localB1, localB2;
    GetLocalSizes(D_MC_STAR, B1, B2, localB1, localB2);
    cost += LocalSymmAcc::GetCost(LEFT, &localB1, &localB2);

    Sizes localC1, localC2;
    GetLocalSizes(D_MR_MC, C1, C2, localC1, localC2);
    cost+=SumScatterFrom::GetCost(D_MR_MC, D_MR_STAR, &localC1, &localC2);

    cost+= RedistNode::GetCost(D_MR_MC, D_MC_MR, C1, C2);
    cost += SumScatterNode::GetCost(orig->InputLocalM(2), orig->InputLocalN(2), D_MC_MR, D_MC_STAR);
    cost += Axpy::GetCost(SMLAYER, orig->InputLocalM(2), orig->InputLocalN(2));
    return cost;
  }
  else {
    DistType type = complex ? D_MR_STAR_H : D_MR_STAR_T;
    Cost cost = RedistNode::GetCost(D_MC_MR, type, B1, B2);
    cost += RedistNode::GetCost(D_MC_MR, D_STAR_MC, B1, B2);

    Sizes localB1, localB2;
    GetLocalSizes(D_STAR_MC, B1, B2, localB1, localB2);
    cost += LocalSymmAcc::GetCost(RIGHT, &localB1, &localB2);

    type = complex ? D_MC_MR_H : D_MC_MR_T;
    Sizes localC1, localC2;
    GetLocalSizes(type, C1, C2, localC1, localC2);
    cost+=SumScatterFrom::GetCost(type, complex ? D_MC_STAR_H : D_MC_STAR_T, &localC1, &localC2);
    cost+= RedistNode::GetCost(type, complex ? D_MR_MC_H : D_MR_MC_T, C1, C2);

    cost += SumScatterNode::GetCost(orig->InputLocalN(2), orig->InputLocalM(2),
				    complex ? D_MR_MC_H : D_MR_MC_T,
				    complex ? D_MR_STAR_H : D_MR_STAR_T);
    cost+= RedistNode::GetCost(complex ? D_MR_MC_H : D_MR_MC_T, D_MC_MR, C1, C2);
    cost += Axpy::GetCost(SMLAYER, orig->InputLocalM(2), orig->InputLocalN(2));
    return cost;
  }
}


string BLISHemmLoopExp::GetType() const
{
  return "BLIS Hemm Loop Exp";
}

bool BLISHemmLoopExp::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == Hemm::GetClass()) {
    const Hemm *hemm = (Hemm*)node;
    if (hemm->GetLayer() == m_fromLayer)
      return true;
  }
  return false;
}

void BLISHemmLoopExp::Apply(Poss *poss, Node *node) const
{
  Hemm *hemm = (Hemm*)node;

  NodeConn *connA, *connB, *connC;
  connA = hemm->m_inputs[0];
  connB = hemm->m_inputs[1];
  connC = hemm->m_inputs[2];
  
  Node *Ain = connA->m_n;
  unsigned int Anum = connA->m_num;
  Node *Bin = connB->m_n;
  unsigned int Bnum = connB->m_num;
  Node *Cin = connC->m_n;
  unsigned int Cnum = connC->m_num;

  if (hemm->m_side == RIGHT) {
    throw;
    Ain = AddTranspose(TRANS, false, Ain, Anum, true);
    Anum = 0;
    Bin = AddTranspose(TRANS, false, Bin, Bnum, true);
    Bnum = 0;
    Cin = AddTranspose(TRANS, true, Cin, Cnum, true);
    Cnum = 0;
  }

  Split *splitA = new Split(PARTDOWN, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();

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
                     NOTUP, NOTUP);
  splitC->SetIndepIters();
  
  PackBuff *aBuff = new PackBuff(splitA->GetName(1, BLISLOOP).m_name,
				 PACKROWPANS, PACKABLOCK, hemm->m_tri, NOTTRIDIAG,
				 hemm->m_type == REAL ? SYMM : HERM,
				 true, false, false, false,
				 USEMRSIZE, USEKRSIZE);
  Pack *aPack = new Pack(PACKROWPANS, 2, false, true, false, false, false);
  aBuff->AddInput(splitA, 1);
  aPack->AddInput(splitA, 1);
  aPack->AddInput(aBuff, 0);  

  Gemm *gebp = new Gemm(m_toLayer, NORMAL, NORMAL,
			hemm->m_alpha, hemm->m_beta, hemm->m_type);
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

bool HemmRightToLeft::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == Hemm::GetClass()) {
    const Hemm *hemm = (Hemm*)node;
    if (hemm->GetLayer() != m_layer)
      return false;
    return hemm->m_side == RIGHT;
  }
  return false;
}

void HemmRightToLeft::Apply(Poss *poss, Node *node) const
{
  Hemm *hemm = (Hemm*)node;
  hemm->m_side = LEFT;

  InsertTranspose(CONJ, true, hemm, 0, true);
  InsertTranspose(TRANS, true, hemm, 1, true);
  InsertTranspose(TRANS, false, hemm, 2, true);

  Transpose *newTrans = new Transpose(TRANS, false);
  poss->AddNode(newTrans);
  hemm->RedirectAllChildren(newTrans);
  newTrans->AddInput(hemm, 0);
}


bool HemmLowerLayer::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == Hemm::GetClass()) {
    const Hemm *hemm = (Hemm*)node;
    if (hemm->GetLayer() != m_fromLayer)
      return false;
    if (m_dim == DIMK)
      return (*(hemm->InputLocalM(0)) <= m_bs
	      && *(hemm->InputLocalN(0)) <= m_bs);
    else if (m_dim == DIMN)
      return (*(hemm->InputLocalN(1)) <= m_bs);
    else
      throw;
  }
  return false;
}

void HemmLowerLayer::Apply(Poss *poss, Node *node) const
{
  Hemm *hemm = (Hemm*)node;
  hemm->SetLayer(m_toLayer);
}
