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
#if DOBLIS||DOELEM


#include "elemRedist.h"
#include "loopSupport.h"
#include "blas.h"
#include <stdio.h>
#include <string>
#include <iomanip>
#include <sstream>
#include "TriInv.h"

void TriInv::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  m_tri = ((TriInv*)orig)->m_tri;
}

void TriInv::FlattenCore(ofstream &out) const
{
  DLAOp<1,1>::FlattenCore(out);
  WRITE(m_tri);
}

void TriInv::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<1,1>::UnflattenCore(in, info);
  READ(m_tri);
}

Phase TriInv::MaxPhase() const
{
#if DODPPHASE
  if (m_layer == ABSLAYER || m_layer == DMLAYER)
    return DPPHASE;
  else
    return NUMPHASES;
#else
  throw;
#endif
}

void TriInv::PrintCode(IndStream &out)
{
  out.Indent();
  if (m_layer == ABSLAYER)
    *out << "AbsTriInv( " << TriToStr(m_tri) << ", NON_UNIT" << GetInputName(0).str() << " );\n";
  else if (m_layer == DMLAYER)
    *out << "DistTriInv( " << TriToStr(m_tri) << ", NON_UNIT, " << GetInputName(0).str() << " );\n";
  else if (m_layer == SMLAYER)
    *out << "internal::LocalTriangularInverse( " 
	 << TriToStr(m_tri) << ", NON_UNIT, " 
	 << GetInputName(0).str() << " );\n";
  else
    throw;  
}

void TriInv::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<1,1>::Prop();

#if DOELEM
  if (m_layer == ABSLAYER || m_layer == DMLAYER) {
    if (InputDataType(0).m_dist != D_MC_MR) {
      cout << "input not D_MC_MR";
      throw;
    }
  }
  else if (m_layer == SMLAYER) {
    if (InputDataType(0).m_dist != D_STAR_STAR) {
      cout << "input not D_STAR_STAR";
      throw;
    }
  }
#else
  throw;
#endif


    if (m_layer == ABSLAYER || m_layer == DMLAYER)
      m_cost = ZERO;
    else if (m_layer == SMLAYER)
      m_cost = GAMMA * 1.0/3 * GetM(0)->SumCubes();
    else 
      throw;
  }
}

bool TriInvLoopExp::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != TriInv::GetClass())
    return false;
  const TriInv *triInv = (TriInv*)node;
  return triInv->GetLayer() == ABSLAYER;
}

string TriInvLoopExp::GetType() const
{
  switch(m_var) {
    case(1):
      return "TriInv Loop Exp var 1";
    case(2):
      return "TriInv Loop Exp var 2";
    case(3):
      return "TriInv Loop Exp var 3";
    case(8):
      return "TriInv Loop Exp var 8";
    default:
      throw;    
  }
}

void TriInvLoopExp::Apply(Node *node) const
{
  TriInv *triInv = (TriInv*)node;
  Tri tri = triInv->m_tri;
  Loop *loop = NULL;
  switch(m_var) {
  case(1):
    if (tri == LOWER)
      loop = TriInvAlgVar1Lower(triInv->Input(0), triInv->InputConnNum(0));
    else
      loop = TriInvAlgVar1Upper(triInv->Input(0), triInv->InputConnNum(0));
    break;
  case(2):
    if (tri == UPPER)
      loop = TriInvAlgVar2Upper(triInv->Input(0), triInv->InputConnNum(0));
    else
      loop = TriInvAlgVar2Lower(triInv->Input(0), triInv->InputConnNum(0));
    break;
  case(3):
    if (tri == UPPER)
      loop = TriInvAlgVar3Upper(triInv->Input(0), triInv->InputConnNum(0));
    else
      loop = TriInvAlgVar3Lower(triInv->Input(0), triInv->InputConnNum(0));
    break;
  case(8):
    if (tri != UPPER)
      loop = TriInvAlgVar8Lower(triInv->Input(0), triInv->InputConnNum(0));
    break;
  default:
    throw;
  }

  node->m_poss->AddLoop(loop);
  
  node->RedirectChildren(loop->OutTun(0),0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

#if DOELEM
bool DistTriInvToLocalTriInv::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != TriInv::GetClass())
    return false;
  const TriInv *triInv = (TriInv*)node;
  return triInv->GetLayer() == DMLAYER;
}

void DistTriInvToLocalTriInv::Apply(Node *node) const
{
  TriInv *triInv = (TriInv*)node;
  RedistNode *node1 = new RedistNode(D_STAR_STAR);
  TriInv *node2 = new TriInv(SMLAYER, triInv->m_tri);
  RedistNode *node3 = new RedistNode(D_MC_MR);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node2->AddInput(node1,0);
  node3->AddInput(node2,0);
  node->m_poss->AddNodes(3, node1, node2, node3);
  node->RedirectChildren(node3,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}
#endif

Loop* TriInvAlgVar1Lower(Node *in, unsigned int num)
{

  Split *split = new Split(PARTDIAG, POSSTUNIN, true);
  split->AddInput(in, num);
  split->SetUpStats(FULLUP, FULLUP,
		    NOTUP, NOTUP);

  Trxm* trmm1 = new Trxm(false, DMLAYER, RIGHT, LOWER, NONUNIT, NORMAL, COEFONE, REAL);
  trmm1->AddInputs(4, split, 0, split, 1);
  Poss *poss5 = new Poss(trmm1,false);
  PSet *set5 = new PSet(poss5);

  Trxm *trsm1 = new Trxm(true, DMLAYER, LEFT, LOWER, NONUNIT, NORMAL, COEFNEGONE, REAL);
  trsm1->AddInputs(4, split, 4, set5->OutTun(0), 0);
  Poss *poss7 = new Poss(trsm1,false);
  PSet *set7 = new PSet(poss7);

  TriInv *triInv = new TriInv(DMLAYER, LOWER);
  triInv->AddInput(split, 4);
  Poss *poss4 = new Poss(triInv,false);
  PSet *set4 = new PSet(poss4);

  Combine *comA2 = new Combine(PARTDIAG, POSSTUNOUT);
  comA2->AddInput(split,0);
  comA2->AddInput(set7->OutTun(0),0);
  comA2->AddInput(split,2);
  comA2->AddInput(split,3);
  comA2->AddInput(set4->OutTun(0),0);
  comA2->AddInput(split,5);
  comA2->AddInput(split,6);
  comA2->AddInput(split,7);
  comA2->AddInput(split,8);
  comA2->AddInput(split,9);

  comA2->CopyTunnelInfo(split);

  Poss *loopPoss = new Poss(1, comA2);
#if DOELEM
  Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);

  return loop;
#else
  throw;
#endif
}

Loop* TriInvAlgVar1Upper(Node *in, unsigned int num)
{
  Split *split = new Split(PARTDIAG, POSSTUNIN, true);
  split->AddInput(in,num);
  split->SetUpStats(FULLUP, NOTUP,
		    FULLUP, NOTUP);

  Trxm* trmm = new Trxm(false, DMLAYER, LEFT, UPPER, NONUNIT, NORMAL, COEFNEGONE, REAL);
  trmm->AddInputs(4, split, 0, split, 3);
  Poss *poss5 = new Poss(trmm,false);
  PSet *set5 = new PSet(poss5);

  Trxm *trsm = new Trxm(true, DMLAYER, RIGHT, UPPER, NONUNIT, NORMAL, COEFONE, REAL);
  trsm->AddInputs(4, split, 4, set5->OutTun(0), 0);
  Poss *poss7 = new Poss(trsm,false);
  PSet *set7 = new PSet(poss7);

  TriInv *triInv = new TriInv(DMLAYER, UPPER);
  triInv->AddInput(split, 4);
  Poss *poss4 = new Poss(triInv,false);
  PSet *set4 = new PSet(poss4);


  Combine *comA2 = new Combine(PARTDIAG, POSSTUNOUT);
  comA2->AddInput(split,0);
  comA2->AddInput(split,1);
  comA2->AddInput(split,2);
  comA2->AddInput(set7->OutTun(0),0);
  comA2->AddInput(set4->OutTun(0),0);
  comA2->AddInput(split,5);
  comA2->AddInput(split,6);
  comA2->AddInput(split,7);
  comA2->AddInput(split,8);
  comA2->AddInput(split,9);

  comA2->CopyTunnelInfo(split);

  Poss *loopPoss = new Poss(1, comA2);
#if DOELEM
  Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);

  return loop;
#else
  throw;
#endif
}

Loop* TriInvAlgVar2Lower(Node *in, unsigned int num)
{

  Split *split = new Split(PARTDIAG, POSSTUNIN, true);
  split->AddInput(in, num);
  split->SetUpStats(FULLUP, FULLUP,
		    FULLUP, NOTUP);

  Trxm *trsm1 = new Trxm(true, DMLAYER, LEFT, LOWER, NONUNIT, NORMAL, COEFONE, REAL);
  trsm1->AddInputs(4, split, 8, split, 5);
  Poss *poss5 = new Poss(trsm1,false);
  PSet *set5 = new PSet(poss5);

  Trxm *trsm2 = new Trxm(true, DMLAYER, RIGHT, LOWER, NONUNIT, NORMAL, COEFNEGONE, REAL);
  trsm2->AddInputs(4, split, 4, set5->OutTun(0), 0);
  Poss *poss7 = new Poss(trsm2,false);
  PSet *set7 = new PSet(poss7);

  TriInv *triInv = new TriInv(DMLAYER, LOWER);
  triInv->AddInput(split, 4);
  Poss *poss4 = new Poss(triInv,false);
  PSet *set4 = new PSet(poss4);

  Combine *comA2 = new Combine(PARTDIAG, POSSTUNOUT);
  comA2->AddInput(split,0);
  comA2->AddInput(split,1);
  comA2->AddInput(split,2);
  comA2->AddInput(split,3);
  comA2->AddInput(set4->OutTun(0),0);
  comA2->AddInput(set7->OutTun(0),0);
  comA2->AddInput(split,6);
  comA2->AddInput(split,7);
  comA2->AddInput(split,8);
  comA2->AddInput(split,9);

  comA2->CopyTunnelInfo(split);

  Poss *loopPoss = new Poss(1, comA2);

#if DOELEM
  Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);

  return loop;
#else
  throw;
#endif
}


Loop* TriInvAlgVar2Upper(Node *in, unsigned int num)
{
  Split *split = new Split(PARTDIAG, POSSTUNIN, true);
  split->AddInput(in,num);
  split->SetUpStats(FULLUP, FULLUP,
		    FULLUP, NOTUP);

  Trxm *Trsm1 = new Trxm(true, DMLAYER, RIGHT, UPPER, NONUNIT, NORMAL, COEFONE, REAL);
  Trsm1->AddInputs(4, split, 8, split, 7);
  Poss *poss5 = new Poss(Trsm1,false);
  PSet *set5 = new PSet(poss5);

  Trxm *Trsm2 = new Trxm(true, DMLAYER, LEFT, UPPER, NONUNIT, NORMAL, COEFNEGONE, REAL);
  Trsm2->AddInputs(4, split, 4, set5->OutTun(0), 0);
  Poss *poss7 = new Poss(Trsm2,false);
  PSet *set7 = new PSet(poss7);

  TriInv *triInv = new TriInv(DMLAYER, UPPER);
  triInv->AddInput(split, 4);
  Poss *poss4 = new Poss(triInv,false);
  PSet *set4 = new PSet(poss4);


  Combine *comA2 = new Combine(PARTDIAG, POSSTUNOUT);
  comA2->AddInput(split,0);
  comA2->AddInput(split,1);
  comA2->AddInput(split,2);
  comA2->AddInput(split,3);
  comA2->AddInput(set4->OutTun(0),0);
  comA2->AddInput(split,5);
  comA2->AddInput(split,6);
  comA2->AddInput(set7->OutTun(0),0);
  comA2->AddInput(split,8);
  comA2->AddInput(split,9);

  comA2->CopyTunnelInfo(split);

  Poss *loopPoss = new Poss(1, comA2);

#if DOELEM
  Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);

  return loop;
#else
  throw;
#endif
}

Loop* TriInvAlgVar8Lower(Node *in, unsigned int num)
{
  Split *split = new Split(PARTDIAGBACK, POSSTUNIN, true);
  split->AddInput(in, num);
  split->SetUpStats(NOTUP, FULLUP,
                    NOTUP, FULLUP);

						
  //L_{21} = L_{21} L_{22}
  Trxm *trmm = new Trxm(false, DMLAYER, RIGHT, LOWER, NONUNIT, NORMAL, COEFONE, REAL);
  trmm->AddInputs(4, split, 8, split, 5);
  Poss *poss = new Poss(trmm,false);
  PSet *set = new PSet(poss);

  //L_{21} = - L_{11}^{-1} L_{21}
  Trxm *Trsm1 = new Trxm(true, DMLAYER, LEFT, LOWER, NONUNIT, NORMAL, COEFNEGONE, REAL);
  Trsm1->AddInputs(4, split, 4, set->OutTun(0), 0);
  Poss *poss5 = new Poss(Trsm1,false);
  PSet *set5 = new PSet(poss5);

  //L_{11} = L_{11}^{-1}
  TriInv *triInv = new TriInv(DMLAYER, LOWER);
  triInv->AddInput(split, 4);
  Poss *poss4 = new Poss(triInv,false);
  PSet *set4 = new PSet(poss4);

  Combine *comA2 = new Combine(PARTDIAGBACK, POSSTUNOUT);
  comA2->AddInput(split,0);
  comA2->AddInput(split,1);
  comA2->AddInput(split,2);
  comA2->AddInput(split,3);
  comA2->AddInput(set4->OutTun(0),0);
  comA2->AddInput(set5->OutTun(0),0);
  comA2->AddInput(split,6);
  comA2->AddInput(split,7);
  comA2->AddInput(split,8);
  comA2->AddInput(split,9);

  comA2->CopyTunnelInfo(split);

  Poss *loopPoss = new Poss(1, comA2);
#if DOELEM
  Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);

  return loop;
#else
  throw;
#endif
}



Loop* TriInvAlgVar3Lower(Node *in, unsigned int num)
{

  Split *split = new Split(PARTDIAG, POSSTUNIN, true);
  split->AddInput(in, num);
  split->SetUpStats(FULLUP, FULLUP,
		    PARTUP, NOTUP);

  Trxm *trsm1 = new Trxm(true, DMLAYER, RIGHT, LOWER, NONUNIT, NORMAL, COEFNEGONE, REAL);
  trsm1->AddInputs(4, split, 4, split, 5);
  Poss *poss5 = new Poss(trsm1,false);
  PSet *set5 = new PSet(poss5);

  Gemm *gemm = new Gemm(DMLAYER, NORMAL, NORMAL, COEFONE, COEFONE,REAL);
  gemm->AddInputs(6, set5->OutTun(0), 0, split, 1, split, 2);
  Poss *poss6 = new Poss(gemm,false);
  PSet *set6 = new PSet(poss6);

  Trxm *trsm2 = new Trxm(true, DMLAYER, LEFT, LOWER, NONUNIT, NORMAL, COEFONE, REAL);
  trsm2->AddInputs(4, split, 4, split, 1);
  Poss *poss7 = new Poss(trsm2,false);
  PSet *set7 = new PSet(poss7);

  TriInv *triInv = new TriInv(DMLAYER, LOWER);
  triInv->AddInput(split, 4);
  Poss *poss4 = new Poss(triInv,false);
  PSet *set4 = new PSet(poss4);

  Combine *comA2 = new Combine(PARTDIAG, POSSTUNOUT);
  comA2->AddInput(split,0);
  comA2->AddInput(set7->OutTun(0),0);
  comA2->AddInput(set6->OutTun(0),0);
  comA2->AddInput(split,3);
  comA2->AddInput(set4->OutTun(0),0);
  comA2->AddInput(set5->OutTun(0),0);
  comA2->AddInput(split,6);
  comA2->AddInput(split,7);
  comA2->AddInput(split,8);
  comA2->AddInput(split,9);

  comA2->CopyTunnelInfo(split);

  Poss *loopPoss = new Poss(1, comA2);
#if DOELEM
  Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);

  return loop;
#else
  throw;
#endif
}


Loop* TriInvAlgVar3Upper(Node *in, unsigned int num)
{
  Split *split = new Split(PARTDIAG, POSSTUNIN, true);
  split->AddInput(in,num);
  split->SetUpStats(FULLUP, PARTUP,
		    FULLUP, NOTUP);

  Trxm *Trsm1 = new Trxm(true, DMLAYER, LEFT, UPPER, NONUNIT, NORMAL, COEFNEGONE, REAL);
  Trsm1->AddInputs(4, split, 4, split, 7);
  Poss *poss5 = new Poss(Trsm1,false);
  PSet *set5 = new PSet(poss5);

  Gemm *gemm = new Gemm(DMLAYER, NORMAL,NORMAL, COEFONE, COEFONE,REAL);
  gemm->AddInputs(6, split, 3, set5->OutTun(0), 0, split, 6);
  Poss *poss6 = new Poss(gemm,false);
  PSet *set6 = new PSet(poss6);

  Trxm *Trsm2 = new Trxm(true, DMLAYER, RIGHT, UPPER, NONUNIT, NORMAL, COEFONE, REAL);
  Trsm2->AddInputs(4, split, 4, split, 3);
  Poss *poss7 = new Poss(Trsm2,false);
  PSet *set7 = new PSet(poss7);

  TriInv *triInv = new TriInv(DMLAYER, UPPER);
  triInv->AddInput(split, 4);
  Poss *poss4 = new Poss(triInv,false);
  PSet *set4 = new PSet(poss4);

  Combine *comA2 = new Combine(PARTDIAG, POSSTUNOUT);
  comA2->AddInput(split,0);
  comA2->AddInput(split,1);
  comA2->AddInput(split,2);
  comA2->AddInput(set7->OutTun(0),0);
  comA2->AddInput(set4->OutTun(0),0);
  comA2->AddInput(split,5);
  comA2->AddInput(set6->OutTun(0),0);
  comA2->AddInput(set5->OutTun(0),0);
  comA2->AddInput(split,8);
  comA2->AddInput(split,9);

  comA2->CopyTunnelInfo(split);

  Poss *loopPoss = new Poss(1, comA2);
#if DOELEM
  Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);

  return loop;
#else
  throw;
#endif
}
#endif
