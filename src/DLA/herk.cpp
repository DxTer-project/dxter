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

#include "layers.h"

#if DOBLIS||DOELEM

#include "herk.h"
#include "blas.h"
#include "elemRedist.h"
#include "string.h"
#include "helperNodes.h"
#include "loopSupport.h"
#include "blis.h"
#include "pack.h"

using namespace std;

HerkProps::HerkProps(Tri tri, Trans transA, Trans transB, Coef alpha, Coef beta, Type type)
:m_alpha(alpha), m_beta(beta)
{
  m_tri = tri;
  m_transA = transA;
  m_transB = transB;
  m_type = type;
}

void HerkProps::Duplicate(const HerkProps *orig)
{
  m_tri = orig->m_tri;
  m_transA = orig->m_transA;
  m_transB = orig->m_transB;
  m_alpha = orig->m_alpha;
  m_beta = orig->m_beta;
  m_type = orig->m_type;
}

void HerkProps::FlattenCore(ofstream &out) const
{
  WRITE(m_transA);
  WRITE(m_transB);
  WRITE(m_tri);
  WRITE(m_alpha);
  WRITE(m_beta);
  WRITE(m_type);
}

void HerkProps::UnflattenCore(ifstream &in, SaveInfo &info)
{
  READ(m_transA);
  READ(m_transB);
  READ(m_tri);
  READ(m_alpha);
  READ(m_beta);
  READ(m_type);
}


Herk::Herk(Layer layer, Tri tri, Trans trans, Coef alpha, Coef beta, Type type)
: HerkProps(tri, trans, trans, alpha, beta, type)
{
  SetLayer(layer);
  if (trans == NORMAL) {
    m_transA = NORMAL;
    m_transB = (type == REAL ? TRANS : CONJTRANS);
  }
  else {
    m_transA = trans;
    m_transB = NORMAL;
  }
}

#if DOELEM
const DistType& Herk::GetDistType(unsigned int num) const
{
  switch (GetLayer()) {
    case (ABSLAYER):
    case (DMLAYER):
      return MC_MR;
    default:
      throw;
  }  
}
#endif

void Herk::SanityCheck()
{
  DLAOp<2,1>::SanityCheck();
  if (m_inputs.size() != 2) {
    throw;
  }
  if (m_type == REAL) {
    if (m_transA == CONJTRANS || m_transB == CONJTRANS)
      throw;
  }
  
#if DOELEM
  if (GetLayer() != ABSLAYER && GetLayer() != DMLAYER)
    throw;
  if (InputDistType(0) != D_MC_MR)
    throw;
  if (InputDistType(1) != D_MC_MR)
    throw;
#else
  
#endif
}



void Herk::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<2,1>::Duplicate(orig, shallow, possMerging);
  const Herk *herk = (Herk*)orig;
  HerkProps::Duplicate(herk);
}


NodeType Herk::GetType() const
{
  string str = "Herk " + TransToStr(m_transA) + TransToStr(m_transB) + TriToStr(m_tri) + LayerNumToStr(GetLayer());
  return str;
}

void Herk::FlattenCore(ofstream &out) const
{
  DLAOp<2,1>::FlattenCore(out);
  HerkProps::FlattenCore(out);
}

void Herk::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<2,1>::UnflattenCore(in, info);
  HerkProps::UnflattenCore(in,info);
}

void Herk::PrintCode(IndStream &out)
{
  out.Indent();
  if (GetLayer() == ABSLAYER) {
    if (m_type == REAL)
      *out << "AbsSyrk( " ;
    else
      *out << "AbsHerk( " ;
  }
  else if(GetLayer() == DMLAYER) {
    if (m_type == REAL)
      *out << "DistSyrk( " ;
    else
      *out << "DistHerk( " ;
  }
  *out << TriToStr(m_tri) << ", "
  << TransToStr(m_transA) << ", " << TransToStr(m_transB);
  out << m_alpha;
  *out << ", "
  << GetInputName(0).str() << ","
  << "\n" << out.Tabs(1);
  out << m_beta;
  *out << ", " << GetInputName(1).str() << ","
  << " );\n";
}

Phase Herk::MaxPhase() const
{
#if DODPPHASE
  switch (GetLayer()) {
    case (ABSLAYER):
    case (DMLAYER):
      return DPPHASE;
    case (SMLAYER):
      return NUMPHASES;
    default:
      throw;
  }
#elif DOBLIS
  switch (GetLayer()) {
    case (ABSLAYER):
      return SR1PHASE;
    case (S1LAYER):
      return SR2PHASE;
    case (S2LAYER):
      return SR3PHASE;
    case (S3LAYER):
      return NUMPHASES;
    default:
      throw;
  }
  throw;
#else 
  throw;
#endif
}

bool Herk::ShouldCullDP() const
{
  switch (GetLayer()) {
    case (ABSLAYER):
      return false;
    case (DMLAYER):
      return true;
    default:
      throw;
  }
}


string HerkLoopExp::GetType() const
{
  switch(m_var) {
    case(1):
      return "Herk Loop Exp - var 1";
    case(2):
      return "Herk Loop Exp - var 2";
    case(5):
      return "Herk Loop Exp - var 5";
    default:
      throw;
  }
}

bool HerkLoopExp::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == Herk::GetClass()
      && ((Herk*)node)->GetLayer() == m_fromLayer)
    return true;
  return false;
}

void HerkLoopExp::Apply(Poss *poss, Node *node) const
{
  Herk *herk = (Herk*)node;
  Loop *loop;
  
  NodeConn *connA, *connC;
  connA = herk->m_inputs[0];
  connC = herk->m_inputs[1];
  
  switch(m_var) {
    case(1):
      loop = HerkLoopVar1(connA->m_n, connA->m_num,
                          connC->m_n, connC->m_num,
                          herk->m_tri,
                          herk->m_transA,
                          herk->m_alpha, herk->m_beta, herk->m_type, m_toLayer);
      break;
    case(2):
      loop = HerkLoopVar2(connA->m_n, connA->m_num,
                          connC->m_n, connC->m_num,
                          herk->m_tri,
                          herk->m_transA,
                          herk->m_alpha, herk->m_beta, herk->m_type, m_toLayer);
      break;
    case(5):
      loop = HerkLoopVar5(connA->m_n, connA->m_num,
                          connC->m_n, connC->m_num,
                          herk->m_tri,
                          herk->m_transA,
                          herk->m_alpha, herk->m_beta, herk->m_type, m_toLayer);
      break;
    default:
      throw;
  }
  
  poss->AddLoop(loop);
  
  
  node->RedirectChildren(loop->OutTun(1),0);
  node->m_poss->DeleteChildAndCleanUp(node);
}


string TriRKLoopExp::GetType() const
{
  switch(m_var) {
      //  case(1):
      //    return "Herk Loop Exp - var 1";
      //  case(2):
      //    return "Herk Loop Exp - var 2";
    case(5):
      return "TriRK Loop Exp - var 5 (dim K)";
    case(7):
      return "TriRK Loop Exp - var 7 (dim N)";
    default:
      throw;
  }
}

bool TriRKLoopExp::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == TriRK::GetClass()
      && ((TriRK*)node)->GetLayer() == m_fromLayer)
    return true;
  return false;
}

void TriRKLoopExp::Apply(Poss *poss, Node *node) const
{
  TriRK *trirk = (TriRK*)node;
  Loop *loop;
  
  NodeConn *connA, *connB, *connC;
  connA = trirk->m_inputs[0];
  connB = trirk->m_inputs[1];
  connC = trirk->m_inputs[2];
  
  switch(m_var) {
      /*
       case(1):
       loop = TriRKLoopVar1(connA->m_n, connA->m_num,
       connC->m_n, connC->m_num,
       trirk->m_tri,
       trirk->m_transA,
       trirk->m_alpha, trirk->m_beta, trirk->m_type, m_toLayer);
       break;
       case(2):
       loop = TriRKLoopVar2(connA->m_n, connA->m_num,
       connC->m_n, connC->m_num,
       trirk->m_tri,
       trirk->m_transA,
       trirk->m_alpha, trirk->m_beta, trirk->m_type, m_toLayer);
       break;
       */
    case(5):
      loop = TriRKLoopVar5(connA->m_n, connA->m_num,
                           connB->m_n, connB->m_num,
                           connC->m_n, connC->m_num,
                           trirk->m_tri,
                           trirk->m_alpha, trirk->m_beta, trirk->m_type, m_toLayer);
      break;
    case(7):
      loop = TriRKLoopVar7(connA->m_n, connA->m_num,
                           connB->m_n, connB->m_num,
                           connC->m_n, connC->m_num,
                           trirk->m_tri,
                           trirk->m_alpha, trirk->m_beta, trirk->m_type, m_toLayer);
      break;
    default:
      throw;
  }
  
  poss->AddLoop(loop);
  
  node->RedirectChildren(loop->OutTun(2),0);
  node->m_poss->DeleteChildAndCleanUp(node);
}



TriRK::TriRK(Layer layer, Tri tri, Trans transA, Trans transB, Coef alpha, Coef beta, Type type)
: HerkProps(tri, transA, transB, alpha, beta, type),
m_comm(CORECOMM)
{
  SetLayer(layer);
}


Phase TriRK::MaxPhase() const
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
    case (S3LAYER):
      return NUMPHASES;
    default:
      throw;
  }
#else
  throw;
#endif
}

bool TriRK::IsBLISParallelizable() const
{
  return GetLayer() == S3LAYER;
}

bool TriRK::IsParallel() const
{
  return m_comm != CORECOMM;
}


bool TriRK::RemoveParallelization()
{
  m_comm = CORECOMM;
  return false;
}



void TriRK::Parallelize(Comm comm)
{
  if (GetLayer() == S3LAYER)
    m_comm = comm;
  else
    throw;
}

#if DOELEM
const DistType& TriRK::GetDistType(unsigned int num) const
{
  if (GetLayer() == SMLAYER)
    return InputDistType(2);
  else if (GetLayer() == S1LAYER || GetLayer() == S2LAYER || GetLayer() == S3LAYER)
    return InputDistType(2);
  else
    throw;
}
#endif

NodeType TriRK::GetType() const
{
  string str = "TriRK " 
    + TriToStr(m_tri) 
    + LayerNumToStr(GetLayer())
    + " " + CommToStr(m_comm);
  return str;
}

void TriRK::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  const TriRK *triRK = (TriRK*)orig;
  HerkProps::Duplicate(triRK);
  DLAOp<3,1>::Duplicate(orig, shallow, possMerging);
  m_comm = triRK->m_comm;
}


void TriRK::SanityCheck()
{
  DLAOp<3,1>::SanityCheck();
  if (m_inputs.size() != 3) {
    cout << "m_inputs.size() != 3 2\n";
    throw;
  }
#if DOELEM
  if (GetLayer() == SMLAYER) {
    if (InputDistType(2) != D_MC_MR)
      throw;
    if (m_type == REAL) {
      if (m_transA == CONJTRANS || m_transB == CONJTRANS)
        throw;
    }
  }
#elif DOBLIS
  if (GetLayer() == S1LAYER || GetLayer() == S2LAYER || GetLayer() == S3LAYER) {
    if (m_transA != NORMAL || m_transB != NORMAL)
      throw;
    
    if (*InputLocalM(0) != *InputLocalM(2))
      throw;
    
    if (*InputLocalN(0) != *InputLocalM(1))
      throw;
    
    if (*InputLocalN(1) != *InputLocalN(2))
      throw;
  }
#endif
}

void TriRK::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3,1>::Prop();
    
#if DOELEM
    if (GetLayer() == SMLAYER) {
      DistType t0 = InputDistType(0);
      DistType t1 = InputDistType(1);
      
      if (UpdateTrans(m_transB,t1) == NORMAL)
        m_cost = GAMMA * InputLocalM(0)->SumProds111(*InputLocalN(0),*InputLocalN(1));
      else
        m_cost = GAMMA * InputLocalM(0)->SumProds111(*InputLocalN(0),*InputLocalM(1));
      
      if (t0 == D_STAR_STAR && t1 == D_STAR_STAR) {
        
      }
      else if (m_transA == NORMAL && m_transB == CONJTRANS) {
        if (t0 != D_MC_STAR && t0 != D_STAR_MC_T && t0 != D_STAR_MC_H) {
          m_poss->MarkInsane();
        }
        if (t1 != D_STAR_MR_H && t1 != D_MR_STAR) {
          m_poss->MarkInsane();
        }
      }
      else if (m_transA == CONJTRANS && m_transB == NORMAL) {
        if (t0 != D_STAR_MC && t0 != D_MC_STAR_H) {
          //    cout << "Marked insane - first arg : " << DistTypeToStr(t0);
          m_poss->MarkInsane();
        }
        if (t1 != D_MR_STAR_H && t1 != D_STAR_MR && t1 != D_MR_STAR_T) {
          //    cout << "Marked insane - second arg : " << DistTypeToStr(t1);
          m_poss->MarkInsane();
        }
      }
      else if (m_transA == NORMAL && m_transB == TRANS) {
        if (t0 != D_MC_STAR && t0 != D_STAR_MC_T && t0 != D_STAR_MC_H) {
          m_poss->MarkInsane();
        }
        if (t1 != D_STAR_MR_T && t1 != D_MR_STAR) {
          m_poss->MarkInsane();
        }
      }
      else if (m_transA == TRANS && m_transB == NORMAL) {
        if (t0 != D_STAR_MC && t0 != D_MC_STAR_T) {
          m_poss->MarkInsane();
        }
        if (t1 != D_MR_STAR_H && t1 != D_STAR_MR && t1 != D_MR_STAR_T) {
          m_poss->MarkInsane();
        }
      }
      else
        throw;
    }
#elif DOBLIS
     if (GetLayer() == S1LAYER || GetLayer() == S2LAYER || GetLayer() == S3LAYER) {
      const Sizes *sizes1 = InputLocalM(0);
      const Sizes *sizes2 = InputLocalN(0);
      const Sizes *sizes3 = InputLocalN(1);
      unsigned int parFactor = NumCoresInComm(m_comm);
      m_cost = GAMMA * sizes1->SumProds111(*sizes2,*sizes3) / parFactor;
      if (GetLayer() == S3LAYER) {
	//BAM, really, only a portion of A might be brough into the L2 since the computation
	// kernel is parallelized in some way.  It's not just by the number of cores in the
	// comm since, they're split over the m and n dimensions
	m_cost += AdditionalCostForBringingIntoL2(this, 0, sizes1->SumProds11(*sizes2), m_comm);
      }
    }
#endif
  }
}

void TriRK::PrintCode(IndStream &out)
{
  string transAStr, transBStr;

#if DOELEM
  DistType t0 = InputDistType(0);
  DistType t1 = InputDistType(1);
  
  if (GetLayer() == SMLAYER) {
    if (m_transA == NORMAL) {
      if (t0 == D_STAR_MC_T)
        transAStr = TransToStr(TRANS) + ", ";
      else if (t0 == D_STAR_MC_H)
        transAStr = TransToStr(CONJTRANS) + ", ";
      else if (t0 != D_MC_STAR)
        throw;
    }
    else {
      if (t0 == D_STAR_MC)
        transAStr = TransToStr(m_transA) + ", ";
      else if ((m_transA == CONJTRANS && t0 != D_MC_STAR_H)
               || (m_transA == TRANS && t0 != D_MC_STAR_T))
        throw;
    }
    if (m_transB == NORMAL) {
      if (t1 == D_MR_STAR_T)
        transBStr = TransToStr(TRANS) + ", ";
      else if (t1 == D_MR_STAR_H)
        transBStr = TransToStr(CONJTRANS) + ", ";
      else if (t1 != D_STAR_MR)
        throw;
    }
    else {
      if (t1 == D_MR_STAR)
        transBStr = TransToStr(m_transB) + ", ";
      else if ((m_transB == CONJTRANS && t1 != D_STAR_MR_H)
               || (m_transB == TRANS && t1 != D_STAR_MR_T))
        throw;
    }
    
    out.Indent();
    *out << "internal::LocalTrrk( " << TriToStr(m_tri) << ", "
    << transAStr << transBStr;
    out << m_alpha;
    *out << ", "
    << "\n" << out.Tabs(1)
    << GetInputName(0).str() << ","
    << "\n" << out.Tabs(1) << GetInputName(1).str() << ","
    << "\n" << out.Tabs(1);
    out << m_beta;
    *out << ", " << GetInputName(2).str()
    << " );\n";
  }
#elif DOBLIS  
  if (GetLayer() == S1LAYER || GetLayer() == S2LAYER) {
    out.Indent();
    *out << "BlisTrrk" << LayerNumToStr(GetLayer()) << "( ";
    out << m_alpha;
    *out << ", "
    << "\n" << out.Tabs(1)
    << GetInputName(0).str() << ","
    << "\n" << out.Tabs(1) << GetInputName(1).str() << ","
    << "\n" << out.Tabs(1);
    out << m_beta;
    *out << ", " << GetInputName(2).str()
    << " );\n";
    
  }
  else if (GetLayer() == S3LAYER) {
    out.Indent();
    *out << "bli_herk_";
    if (m_tri == UPPER)
      *out << "u";
    else
      *out << "l";
    if (m_comm == CORECOMM)
      *out << "_ker_var2( ";
    else
      *out << "_ker_var2_par( ";
    out << m_alpha;
    *out << ", &"
    << GetInputName(0).str() << ", &" << GetInputName(1).str() << ", \n"
    << out.Tabs(2)
	 << "&BLIS_ONE, &" << GetInputName(2).str() << ", (herk_t*)NULL ";
    if (m_comm != CORECOMM)
      *out << ", " << CommToStr(GetSubComm(m_comm));
    *out << ");\n";

  }
#endif
}

void TriRK::FlattenCore(ofstream &out) const
{
  DLAOp<3,1>::FlattenCore(out);
  HerkProps::FlattenCore(out);
  WRITE(m_comm);
}

void TriRK::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<3,1>::UnflattenCore(in, info);
  TriRK::UnflattenCore(in,info);
  READ(m_comm);
}

bool TriRK::CanTransposeInputs() const
{
  if (GetLayer() == SMLAYER)
    return true;
  else
    return false;
}

#if DOELEM
string TriRKTrans::GetTransType() const
{
  return "TriRK";
}

bool TriRKTrans::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != TriRK::GetClass())
    return false;
  TriRK *trirk = (TriRK*)node;
  if (trirk->GetLayer() != SMLAYER)
    return false;
  DLANode *source = (DLANode*)(trirk->Input(m_argNum));
  Trans trans = (m_argNum==0 ? trirk->m_transA : trirk->m_transB);
  // m_trans==CONJTRANS might work, but this why expend the additional flops when trans
  //  reduces communication just the same
  if (trans == NORMAL && m_trans == TRANS)
    return source->CanTrans();
  else if (trans == NORMAL && m_trans == CONJTRANS)
    return source->CanTrans();
  else if (trans == CONJTRANS && m_trans == CONJTRANS)
    return source->CanTrans();
  else if (trans == TRANS && m_trans == TRANS)
    return source->CanTrans();
  else
    return false;
}
#endif

#if DOELEM
bool DistHerkToLocalTriRK::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == Herk::GetClass()
      && ((Herk*)node)->GetLayer() == DMLAYER)
    return true;
  return false;
}


void DistHerkToLocalTriRK::Apply(Poss *poss, Node *node) const
{
  Herk *orig = (Herk*)node;
  RedistNode *node1;
  RedistNode *node2;
  if (orig->m_transA == NORMAL)
    node1 = new RedistNode(D_MC_STAR);
  else
    node1 = new RedistNode(D_STAR_MC);
  if (orig->m_transB == NORMAL)
    node2 = new RedistNode(D_STAR_MR);
  else
    node2 = new RedistNode(D_MR_STAR);
  TriRK *node3 = new TriRK(SMLAYER, orig->m_tri, orig->m_transA,
                           orig->m_transB, orig->m_alpha, orig->m_beta,
                           orig->m_type);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node2->AddInput(node->Input(0),node->InputConnNum(0));
  node3->AddInput(node1,0);
  node3->AddInput(node2,0);
  node3->AddInput(node->Input(1),node->InputConnNum(1));
  poss->AddNodes(3, node1, node2, node3);
  node->RedirectChildren(node3,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}
#endif

Loop* HerkLoopVar1(Node *Ain, unsigned int Anum,
                   Node *Cin, unsigned int Cnum,
                   Tri tri,
                   Trans trans,
                   Coef alpha, Coef beta, Type type,
                   Layer layer)
{
  Split *splitA = new Split(trans==NORMAL ? PARTDOWN : PARTRIGHT, POSSTUNIN, true);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  Split *splitC = new Split(PARTDIAG, POSSTUNIN);
  splitC->AddInput(Cin, Cnum);
  if (tri == LOWER)
    splitC->SetUpStats(FULLUP, FULLUP,
                       NOTUP, NOTUP);
  else
    splitC->SetUpStats(FULLUP, FULLUP,
                       FULLUP, NOTUP);
  splitC->SetIndepIters();
  
  Node *gemm;
  if (trans == NORMAL)
    gemm = new Gemm(layer, NORMAL, type == REAL ? TRANS : CONJTRANS, alpha, beta, type);
  else
    gemm = new Gemm(layer, type == REAL ? TRANS : CONJTRANS, NORMAL, alpha, beta, type);
  
  if (tri == LOWER)
    gemm->AddInputs(6,
                    splitA, 1,
                    splitA, 0,
                    splitC, 1);
  else
    gemm->AddInputs(6,
                    splitA, 1,
                    splitA, 2,
                    splitC, 7);
  
  Node *herk;
  herk = new Herk(layer, tri, trans, alpha, beta, type);
  
  herk->AddInputs(4, splitA, 1,
                  splitC, 4);
  
  Combine *comA = splitA->CreateMatchingCombine(0);
  
  Combine *comC;
  if (tri == LOWER)
    comC = splitC->CreateMatchingCombine(2,
                                         1, gemm, 0,
                                         4, herk, 0);
  else
    comC = splitC->CreateMatchingCombine(2,
                                         4, herk, 0,
                                         7, gemm, 0);
  
  Poss *loopPoss = new Poss(2, comA, comC);
  Loop *loop;
#if DOELEM
  if (layer == DMLAYER)
    loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  else
    throw;
#elif DOBLIS
    loop = new Loop(BLISLOOP, loopPoss, USEBLISMC);
#endif
  
  return loop;
}

Loop* HerkLoopVar2(Node *Ain, unsigned int Anum,
                   Node *Cin, unsigned int Cnum,
                   Tri tri,
                   Trans trans,
                   Coef alpha, Coef beta, Type type,
                   Layer layer)
{
  Split *splitA = new Split(trans==NORMAL ? PARTDOWN : PARTRIGHT, POSSTUNIN, true);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  Split *splitC = new Split(PARTDIAG, POSSTUNIN);
  splitC->AddInput(Cin, Cnum);
  if (tri == LOWER)
    splitC->SetUpStats(FULLUP, FULLUP,
                       FULLUP, NOTUP);
  else
    splitC->SetUpStats(FULLUP, NOTUP,
                       FULLUP, NOTUP);
  splitC->SetIndepIters();
  
  Node *gemm;
  if (trans == NORMAL)
    gemm = new Gemm(layer, NORMAL, type == REAL ? TRANS : CONJTRANS, alpha, beta, type);
  else
    gemm = new Gemm(layer, type == REAL ? TRANS : CONJTRANS, NORMAL, alpha, beta, type);
  
  if (tri == LOWER)
    gemm->AddInputs(6,
                    splitA, 2,
                    splitA, 1,
                    splitC, 5);
  else
    gemm->AddInputs(6,
                    splitA, 0,
                    splitA, 1,
                    splitC, 3);
  
  Node *herk = new Herk(layer, tri, trans, alpha, beta, type);
  
  herk->AddInputs(4, splitA, 1,
                  splitC, 4);
  
  Combine *comA = splitA->CreateMatchingCombine(0);
  
  Combine *comC;
  if (tri == LOWER)
    comC = splitC->CreateMatchingCombine(2,
                                         4, herk, 0,
                                         5, gemm, 0);
  else
    comC = splitC->CreateMatchingCombine(2,
                                         3, gemm, 0,
                                         4, herk, 0);
  
  Poss *loopPoss = new Poss(2, comA, comC);
  Loop *loop;
#if DOELEM
  if (layer == DMLAYER)
    loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  else
    throw;
#elif DOBLIS
    loop = new Loop(BLISLOOP, loopPoss, USEBLISKC);
#endif
  return loop;
}

Loop* HerkLoopVar5(Node *Ain, unsigned int Anum,
                   Node *Cin, unsigned int Cnum,
                   Tri tri,
                   Trans trans,
                   Coef alpha, Coef beta, Type type,
                   Layer layer)
{
  ScaleTrapNode *scal = new ScaleTrapNode(layer, LEFT, tri, beta);
  scal->AddInput(Cin, Cnum);
  
  Split *splitA = new Split(trans==NORMAL ? PARTRIGHT : PARTDOWN, POSSTUNIN, true);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  LoopTunnel *Ctun = new LoopTunnel(POSSTUNIN);
  Ctun->AddInput(scal, 0);
  Ctun->SetAllStats(PARTUP);
  
  Node *herk;
  herk = new Herk(layer, tri, trans, alpha, beta, type);
  herk->AddInputs(4, splitA, 1,
                  Ctun, 0);
  
  Combine *comA = splitA->CreateMatchingCombine(0);
  
  LoopTunnel *CtunOut = new LoopTunnel(POSSTUNOUT);
  CtunOut->AddInput(herk,0);
  CtunOut->AddInput(Ctun,0);
  CtunOut->CopyTunnelInfo(Ctun);
  
  Poss *loopPoss = new Poss(2, comA, CtunOut);
  Loop *loop;
#if DOELEM
    loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  #else
    loop = new Loop(BLISLOOP, loopPoss, USEBLISKC);
#endif
  
  return loop;
}

Loop* TriRKLoopVar5(Node *Ain, unsigned int Anum,
                    Node *Bin, unsigned int Bnum,
                    Node *Cin, unsigned int Cnum,
                    Tri tri,
                    Coef alpha, Coef beta, Type type,
                    Layer layer)
{
  ScaleTrapNode *scal = new ScaleTrapNode(layer, LEFT, tri, beta);
  scal->AddInput(Cin, Cnum);
  
  Split *splitA = new Split(PARTRIGHT, POSSTUNIN, true);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  Split *splitB = new Split(PARTDOWN, POSSTUNIN, false);
  splitB->AddInput(Bin, Bnum);
  splitB->SetAllStats(FULLUP);
  splitB->SetIndepIters();
  
  LoopTunnel *Ctun = new LoopTunnel(POSSTUNIN);
  Ctun->AddInput(scal, 0);
  Ctun->SetAllStats(PARTUP);
  
  Node *herk;
  herk = new TriRK(layer, tri, NORMAL, NORMAL, alpha, COEFONE, type);
  herk->AddInputs(6,
                  splitA, 1,
                  splitB, 1,
                  Ctun, 0);

#if DOELEM
    throw;
#else  
  Combine *comA = splitA->CreateMatchingCombine(0);
  
  Combine *comB = splitB->CreateMatchingCombine(0);
  
  LoopTunnel *CtunOut = new LoopTunnel(POSSTUNOUT);
  CtunOut->AddInput(herk,0);
  CtunOut->AddInput(Ctun,0);
  CtunOut->CopyTunnelInfo(Ctun);
  
  Loop *loop;

    Poss *loopPoss = new Poss(3, comA, comB, CtunOut);
    loop = new Loop(BLISLOOP, loopPoss, USEBLISKC);
  
  loop->SetDim(DIMK);
  
  return loop;
#endif
}



Loop* TriRKLoopVar7(Node *Ain, unsigned int Anum,
                    Node *Bin, unsigned int Bnum,
                    Node *Cin, unsigned int Cnum,
                    Tri tri,
                    Coef alpha, Coef beta, Type type,
                    Layer layer)
#if DOBLIS
{
  LoopTunnel *Atun = new LoopTunnel(POSSTUNIN);
  Atun->AddInput(Ain, Anum);
  Atun->SetAllStats(FULLUP);
  Atun->SetIndepIters();
  
  Split *splitB = new Split(tri==LOWER ? PARTRIGHT : PARTLEFT, POSSTUNIN, false);
  splitB->AddInput(Bin, Bnum);
  splitB->SetAllStats(FULLUP);
  splitB->SetIndepIters();
  
  Split *splitC = new Split(tri==LOWER ? PARTRIGHT : PARTLEFT, POSSTUNIN, true);
  splitC->AddInput(Cin, Cnum);
  if (tri==LOWER)
    splitC->SetUpStats(FULLUP, NOTUP,
                       FULLUP, NOTUP);
  else
    splitC->SetUpStats(NOTUP, FULLUP,
                       NOTUP, FULLUP);
  splitC->SetIndepIters();
  
  GetUpToDiag *diag = new GetUpToDiag(tri, PARTRIGHT);
  if (tri == LOWER)
    diag->AddInputs(6,
                    splitC, 0,
                    splitC, 1,
                    Atun, 0);
  else
    diag->AddInputs(6,
                    splitC, 2,
                    splitC, 1,
                    Atun, 0);
  
  Node *trirk;
  trirk = new TriRK(layer, tri, NORMAL, NORMAL, alpha, beta, type);
  trirk->AddInputs(6,
                   diag, 1,
                   splitB, 1,
                   diag, 0);
  
  LoopTunnel *AtunOut = new LoopTunnel(POSSTUNOUT);
  AtunOut->AddInput(Atun,0);
  AtunOut->AddInput(Atun,0);
  AtunOut->CopyTunnelInfo(Atun);
  
  Combine *comB = splitB->CreateMatchingCombine(0);
  
  Combine *comC = splitC->CreateMatchingCombine(1,
                                                1, trirk, 0);
  
  Poss *loopPoss = new Poss(3, AtunOut, comB, comC);
  Loop *loop;
  if (layer == DMLAYER)
    throw;
  else
    loop = new Loop(BLISLOOP, loopPoss, USEBLISNC);

  loop->SetDim(DIMN);
  
  return loop;
}
#else
{
  throw;
}
#endif


#if DOBLIS
string BLISTriRKLoopExp::GetType() const
{
  return "BLISTriRKLoopExp " + LayerNumToStr(m_fromLayer)
  + " " + LayerNumToStr(m_toLayer);
}


bool BLISTriRKLoopExp::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() != TriRK::GetClass())
    return false;
  const TriRK *trirk = (TriRK*)node;
  if (trirk->m_layer != m_fromLayer)
    return false;
  else
    return true;
  
}

void BLISTriRKLoopExp::Apply(Poss *poss, Node *node) const
{
  TriRK *orig = (TriRK*)node;
  
  
  Node *aSrc = node->Input(0);
  unsigned int aSrcNum = node->InputConnNum(0);
  if (orig->m_transA != NORMAL) {
    throw;
    //    aSrc = AddTranspose(orig->m_transA, true, aSrc, aSrcNum, true);
    //    aSrcNum = 0;
  }
  
  Node *bSrc = node->Input(1);
  unsigned int bSrcNum = node->InputConnNum(1);
  if (orig->m_transB != NORMAL) {
    throw;
    //    bSrc = AddTranspose(orig->m_transB, true, bSrc, bSrcNum, true);
    //    bSrcNum = 0;
  }
  
  Node *cSrc = node->Input(2);
  unsigned int cSrcNum = node->InputConnNum(2);
  
  Loop *loop = BLISHerkLoop(aSrc, aSrcNum,
                            bSrc, bSrcNum,
                            cSrc, cSrcNum,
                            orig->m_tri,
                            orig->m_alpha,
                            orig->m_type,
                            m_toLayer);
  
  poss->AddLoop(loop);
  
  node->RedirectChildren(loop->OutTun(2),0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

Loop* BLISHerkLoop(Node *Ain, unsigned int Anum,
                   Node *Bin, unsigned int Bnum,
                   Node *Cin, unsigned int Cnum,
                   Tri tri,
                   Coef alpha, Type type,
                   Layer layer)
{
  Split *splitA = new Split(tri==LOWER ? PARTUPWARD : PARTDOWN, POSSTUNIN);
  splitA->AddInput(Ain, Anum);
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  PackBuff *bBuff = new PackBuff(Bin->GetName(Bnum).m_name,
                                 PACKCOLPANS, PACKBPANEL, NOTTRI, NOTTRIDIAG, GEN,
                                 false, false, false, false,
                                 USEKRSIZE, USENRSIZE );
  bBuff->AddInput(Bin, Bnum);
  
  Pack *bPack = new Pack(PACKCOLPANS, 2, false, false, false, false, false);
  bPack->AddInput(Bin, Bnum);
  bPack->AddInput(bBuff, 0);
  
  
  LoopTunnel *Btun = new LoopTunnel(POSSTUNIN);
  Btun->AddInput(bPack, 0);
  Btun->SetAllStats(FULLUP);
  Btun->SetIndepIters();
  
  Split *splitC = new Split(tri==LOWER ? PARTUPWARD : PARTDOWN, POSSTUNIN, true);
  splitC->AddInput(Cin, Cnum);
  if (tri==LOWER)
    splitC->SetUpStats(NOTUP, NOTUP,
                       FULLUP, FULLUP);
  else
    splitC->SetUpStats(FULLUP, FULLUP,
                       NOTUP, NOTUP);
  splitC->SetIndepIters();
  
  SetObjProps *props = new SetObjProps(tri, NOTTRIDIAG, type == REAL ? SYMM : HERM);
  props->AddInput(splitC, 1);
  
  
  GetUpToDiag *diag = new GetUpToDiag(tri, PARTDOWN);
  if (tri == LOWER)
    diag->AddInputs(6,
                    splitC, 0,
                    props, 0,
                    Btun, 0);
  else
    diag->AddInputs(6,
                    splitC, 2,
                    props, 0,
                    Btun, 0);
  
  
  PackBuff *aBuff = new PackBuff(splitA->GetName(1,BLISLOOP).m_name,
                                 PACKROWPANS, PACKABLOCK, NOTTRI, NOTTRIDIAG, GEN,
                                 false, false, false, false,
                                 USEMRSIZE, USEKRSIZE );
  aBuff->AddInput(splitA, 1);
  
  Pack *aPack = new Pack(PACKROWPANS, 2, false, false, false, false, false);
  aPack->AddInput(splitA, 1);
  aPack->AddInput(aBuff, 0);
  
  TriRK *herkbp = new TriRK(layer, tri, NORMAL, NORMAL, alpha, COEFONE, type);
  herkbp->AddInputs(6,
                    aPack, 0,
                    diag, 1,
                    diag, 0);
  /*
   CombineDiag *triCombine = new CombineDiag;
   triCombine->AddInputs(4,
   herkbp, 0,
   splitC, 1);
   */
  
  Combine *comA = splitA->CreateMatchingCombine(0);
  
  LoopTunnel *BtunOut = new LoopTunnel(POSSTUNOUT);
  BtunOut->AddInput(Btun, 0);
  BtunOut->AddInput(Btun, 0);
  BtunOut->CopyTunnelInfo(Btun);
  
  Combine *comC = splitC->CreateMatchingCombine(1,
                                                1, herkbp, 0);
  
  Poss *loopPoss = new Poss(3, comA, BtunOut, comC);
  Loop *loop;
  if (layer == DMLAYER)
    throw;
  else
    loop = new Loop(BLISLOOP, loopPoss, USEBLISMC);

  loop->SetDim(DIMM);
  
  return loop;
}
#endif

bool TriRKLowerLayer::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == TriRK::GetClass()) {
    const TriRK *trirk = (TriRK*)node;
    if (trirk->GetLayer() != m_fromLayer)
      return false;
    if (trirk->m_transA != NORMAL || trirk->m_transB != NORMAL)
      throw;
    if (m_dim == DIMK) {
      return (*(trirk->InputLocalN(0)) <= m_bs);
    }
    else if (m_dim == DIMN) {
      return (*(trirk->InputLocalN(2)) <= m_bs);
    }
    else
      throw;
  }
  return false;
  
}

void TriRKLowerLayer::Apply(Poss *poss, Node *node) const
{
  TriRK *trirk = (TriRK*)node;
  trirk->SetLayer(m_toLayer);
}


string TriRKLowerLayer::GetType() const
{ 
  return "TriRK lower layer " + LayerNumToStr(m_fromLayer) 
  + " to " + LayerNumToStr(m_toLayer);
}

#if DOBLIS
bool HerkToTriRK::CanApply(const Poss *poss, const Node *node) const
{
  return (node->GetNodeClass() == Herk::GetClass() &&
          ((DLANode*)node)->GetLayer() == m_layer);
}

void HerkToTriRK::Apply(Poss *poss, Node *node) const
{
  Herk *herk = (Herk*)node;
  Transpose *Atrans = new Transpose(herk->m_type == REAL ? 
                                    TRANS : CONJTRANS,
                                    true);
  Node *Ain = herk->Input(0);
  unsigned int Anum = herk->InputConnNum(0);
  Atrans->AddInput(Ain, Anum);
  
  TriRK *trirk = new TriRK(m_layer,
                           herk->m_tri,
                           NORMAL,
                           NORMAL,
                           herk->m_alpha,
                           herk->m_beta,
                           herk->m_type);
  
  if (herk->m_transA != NORMAL) {
    trirk->AddInputs(4,
                     Atrans, 0,
                     Ain, Anum);
  }
  else if (herk->m_transB != NORMAL) {
    trirk->AddInputs(4,
                     Ain, Anum,
                     Atrans, 0);
  }
  else
    throw;
  
  trirk->AddInput(herk->Input(1), herk->InputConnNum(1));
  
  poss->AddNode(trirk);
  poss->AddNode(Atrans);
  
  herk->RedirectChildren(trirk, 0);
  poss->DeleteChildAndCleanUp(node);
}
#endif
#endif
