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
#if DOELEM||DOBLIS

#include "blas.h"
#include "hetrmm.h"
#include "elemRedist.h"
#include "string.h"
#include "helperNodes.h"
#include "loopSupport.h"

using namespace std;

NodeType Hetrmm::GetType() const 
{
  return "HetrmmLoop" + LayerNumToStr(GetLayer());
}

void Hetrmm::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<1,1>::Duplicate(orig, shallow, possMerging);
  m_tri = ((Hetrmm*)orig)->m_tri;
}

void Hetrmm::FlattenCore(ofstream &out) const
{
  DLAOp<1,1>::FlattenCore(out);
  WRITE(m_tri);
}

void Hetrmm::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<1,1>::UnflattenCore(in, info);
  READ(m_tri);
}


void Hetrmm::Prop()
{
#if DOELEM
  if (!IsValidCost(m_cost)) {
    DLAOp<1,1>::Prop();

    if (m_layer == ABSLAYER || m_layer == DMLAYER) {
      if (InputDistType(0) != D_MC_MR)
	throw;
    }
    else if (m_layer == SMLAYER) {
      if (InputDistType(0) != D_STAR_STAR)
	throw;
    }
    else
      throw;

    if (m_layer == SMLAYER)
      m_cost = GAMMA * 1.0/3 * GetM(0)->SumCubes();
    else
      m_cost = 0;
  }
#else
  throw;
#endif
}


Phase Hetrmm::MaxPhase() const 
{
#if DODPPHASE
  if (m_layer == DMLAYER || m_layer == ABSLAYER)
    return DPPHASE;
  else if (m_layer == SMLAYER)
    return NUMPHASES;
  else
    throw;
#else
  throw;
#endif
}

#if DOELEM
const DistType& Hetrmm::GetDistType(unsigned int num) const 
{ 
  if (m_layer == DMLAYER || m_layer == ABSLAYER)
    return MC_MR; 
  else if (m_layer == SMLAYER)
    return STAR_STAR;
  else
    throw;
}
#endif

bool Hetrmm::ShouldCullDP() const 
{
#if DODPPHASE
  return m_layer == ABSLAYER || m_layer == DMLAYER;
#else
  throw;
#endif
}


void Hetrmm::PrintCode(IndStream &out)
{
  out.Indent();
#if DODPPHASE
  if (m_layer == ABSLAYER)
    *out << "HetrmmLoop( " << TransToStr(TRANS) << ", " << TriToStr(m_tri) << ", " << GetInputName(0).str() << " );\n";
  else if (m_layer == DMLAYER)
    *out << "Hetrmm( " << TransToStr(TRANS) << ", " << TriToStr(m_tri) << ", " << GetInputName(0).str() << " );\n";
  else if (m_layer == SMLAYER)
  *out << "internal::LocalHetrmm( " << TransToStr(TRANS) << ", " << TriToStr(m_tri) << ", " << GetInputName(0).str() << " );\n";
  else
    throw;
#else
  throw;
#endif
}

string HetrmmLoopExp::GetType() const
{
  switch(m_var) {
    case(1):
      return "Hetrmm Loop Exp var 1";
    case(2):
      return "Hetrmm Loop Exp var 2";
    case(3):
      return "Hetrmm Loop Exp var 3";
    default:
      throw;    
  }
}


bool HetrmmLoopExp::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == Hetrmm::GetClass())
    return ((Hetrmm*)node)->GetLayer() == ABSLAYER;
    return true;
  return false;
}

void HetrmmLoopExp::Apply(Poss *poss, Node *node) const
{
  Hetrmm *hetrmm = (Hetrmm*)node;
  Loop *loop=NULL;
  switch(m_var) {
  case(1):
    if (hetrmm->m_tri == LOWER) 
      loop = HetrmmAlgVar1Lower(hetrmm->Input(0), hetrmm->InputConnNum(0));
    else 
      loop = HetrmmAlgVar1Upper(hetrmm->Input(0), hetrmm->InputConnNum(0));
    break;
  case(2):
    if (hetrmm->m_tri == LOWER) 
      loop = HetrmmAlgVar2Lower(hetrmm->Input(0), hetrmm->InputConnNum(0));
    else
      loop = HetrmmAlgVar2Upper(hetrmm->Input(0), hetrmm->InputConnNum(0)); 
    break;
  case(3):
    if (hetrmm->m_tri == LOWER) 
      loop = HetrmmAlgVar3Lower(hetrmm->Input(0), hetrmm->InputConnNum(0));
    else 
      loop = HetrmmAlgVar3Upper(hetrmm->Input(0), hetrmm->InputConnNum(0));
    break;
  default:
    throw;
  }

  poss->AddLoop(loop);
  
  node->RedirectChildren(loop->OutTun(0),0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

#if DOELEM
bool DistHetrmmToLocalHetrmm::CanApply(const Poss *poss, const Node *node) const
{
  if (node->GetNodeClass() == Hetrmm::GetClass()) {
    return ((Hetrmm*)node)->GetLayer() == DMLAYER;
  }
  return false;
}

void DistHetrmmToLocalHetrmm::Apply(Poss *poss, Node *node) const
{
  Hetrmm *hetrmm = (Hetrmm*)node;
  RedistNode *node1 = new RedistNode(D_STAR_STAR);
  Hetrmm *node2 = new Hetrmm(SMLAYER, hetrmm->m_tri);
  RedistNode *node3 = new RedistNode(D_MC_MR);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node2->AddInput(node1,0);
  node3->AddInput(node2,0);
  poss->AddNodes(3, node1, node2, node3);
  node->RedirectChildren(node3,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}
#endif

Loop* HetrmmAlgVar1Lower(Node *in, unsigned int num)
{
  Split *split = new Split(PARTDIAG, POSSTUNIN, true);
  split->AddInput(in,num);
  split->SetUpStats(PARTUP, FULLUP,
		    NOTUP, NOTUP);

  // L_{10} = L_{11}^T L_{10}
  Trxm *trmm = new Trxm(false, DMLAYER, LEFT, LOWER, NONUNIT, CONJTRANS, COEFONE, COMPLEX);
  trmm->AddInputs(4, split, 4, split, 1);
  Poss *poss9 = new Poss(trmm,false);
  PSet *set9 = new PSet(poss9);

  // L_{11} = L_{11}^T L_{11}
  Hetrmm *hetrmm = new Hetrmm(DMLAYER, LOWER);
  hetrmm->AddInput(split, 4);
  Poss *poss10 = new Poss(hetrmm,false);
  PSet *set10 = new PSet(poss10);

  // L_{00} = L_{10}^T L_{10} + L_{00}
  Herk *tri2 = new Herk(DMLAYER, LOWER, CONJTRANS, COEFONE, COEFONE, COMPLEX);
  tri2->AddInputs(4, 
		 split, 1,
		 split, 0);
  Poss *poss8 = new Poss(tri2,false);
  PSet *set8 = new PSet(poss8);


  Combine *comA3 = new Combine(PARTDIAG, POSSTUNOUT);
  comA3->AddInput(set8->OutTun(0),0);
  comA3->AddInput(set9->OutTun(0),0);
  comA3->AddInput(split,2);
  comA3->AddInput(split,3);
  comA3->AddInput(set10->OutTun(0),0);
  comA3->AddInput(split,5);
  comA3->AddInput(split,6);
  comA3->AddInput(split,7);
  comA3->AddInput(split,8);
  comA3->AddInput(split,9);


  comA3->CopyTunnelInfo(split);

  Poss *loopPoss = new Poss(1, comA3);
#if DOELEM
  Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  return loop;
#else
  throw;
#endif
}



Loop* HetrmmAlgVar1Upper(Node *in, unsigned int num)
{
  Split *split = new Split(PARTDIAG, POSSTUNIN, true);
  split->AddInput(in, num);
  split->SetUpStats(PARTUP, NOTUP,
		    FULLUP, NOTUP);

  Herk *tri2 = new Herk(DMLAYER, UPPER, NORMAL, COEFONE, COEFONE, COMPLEX);
  tri2->AddInputs(4,
		 split, 3,
		 split, 0);
  Poss *poss8 = new Poss(tri2,false);
  PSet *set8 = new PSet(poss8);

  Trxm *trmm = new Trxm(false, DMLAYER, RIGHT, UPPER, NONUNIT, CONJTRANS, COEFONE, COMPLEX);
  trmm->AddInputs(4, split, 4, split, 3);
  Poss *poss9 = new Poss(trmm,false);
  PSet *set9 = new PSet(poss9);

  Hetrmm *hetrmm = new Hetrmm(DMLAYER, UPPER);
  hetrmm->AddInput(split, 4);
  Poss *poss10 = new Poss(hetrmm,false);
  PSet *set10 = new PSet(poss10);

  Combine *comA3 = new Combine(PARTDIAG, POSSTUNOUT);
  comA3->AddInput(set8->OutTun(0),0);
  comA3->AddInput(split,1);
  comA3->AddInput(split,2);
  comA3->AddInput(set9->OutTun(0),0);
  comA3->AddInput(set10->OutTun(0),0);
  comA3->AddInput(split,5);
  comA3->AddInput(split,6);
  comA3->AddInput(split,7);
  comA3->AddInput(split,8);
  comA3->AddInput(split,9);

  comA3->CopyTunnelInfo(split);

  Poss *loopPoss = new Poss(1, comA3);
#if DOELEM
  Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  return loop;
#else
  throw;
#endif
}

Loop* HetrmmAlgVar2Lower(Node *in, unsigned int num)
{
  Split *split = new Split(PARTDIAG, POSSTUNIN, true);
  split->AddInput(in,num);
  split->SetUpStats(FULLUP, FULLUP,
		    NOTUP, NOTUP);

  // L_{10} = L_{11}^T L_{10}
  Trxm *trmm = new Trxm(false, DMLAYER, LEFT, LOWER, NONUNIT, CONJTRANS, COEFONE, COMPLEX);
  trmm->AddInputs(4, split, 4, split, 1);
  Poss *poss9 = new Poss(trmm,false);
  PSet *set9 = new PSet(poss9);

  // L_{10} = L_{21}^T L_{20} + L_{10}
  Gemm *gemm = new Gemm(DMLAYER, TRANS, NORMAL, COEFONE, COEFONE, COMPLEX);
  gemm->AddInputs(6, split, 5, split, 2, set9->OutTun(0), 0);
  Poss *poss7 = new Poss(gemm,false);
  PSet *set7 = new PSet(poss7);

  // L_{11} = L_{11}^T L_{11}
  Hetrmm *hetrmm = new Hetrmm(DMLAYER, LOWER);
  hetrmm->AddInput(split, 4);
  Poss *poss10 = new Poss(hetrmm,false);
  PSet *set10 = new PSet(poss10);

  // L_{11} = L_{21}^L_{21} + L_{11}
  Herk *tri2 = new Herk(DMLAYER, LOWER, CONJTRANS, COEFONE, COEFONE, COMPLEX);
  tri2->AddInputs(4,
		  split, 5,
		  set10->OutTun(0), 0);
  Poss *poss8 = new Poss(tri2,false);
  PSet *set8 = new PSet(poss8);

  Combine *comA3 = new Combine(PARTDIAG, POSSTUNOUT);
  comA3->AddInput(split,0);
  comA3->AddInput(set7->OutTun(0),0);
  comA3->AddInput(split,2);
  comA3->AddInput(split,3);
  comA3->AddInput(set8->OutTun(0),0);
  comA3->AddInput(split,5);
  comA3->AddInput(split,6);
  comA3->AddInput(split,7);
  comA3->AddInput(split,8);
  comA3->AddInput(split,9);

  comA3->CopyTunnelInfo(split);

  Poss *loopPoss = new Poss(1, comA3);
#if DOELEM
  Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  return loop;
#else
  throw;
#endif
}

Loop* HetrmmAlgVar2Upper(Node *in, unsigned int num)
{
  Split *split = new Split(PARTDIAG, POSSTUNIN, true);
  split->AddInput(in, num);
  split->SetUpStats(FULLUP, NOTUP,
		    FULLUP, NOTUP);

  Trxm *trmm = new Trxm(false, DMLAYER, RIGHT, UPPER, NONUNIT, CONJTRANS, COEFONE, COMPLEX);
  trmm->AddInputs(4, split, 4, split, 3);
  Poss *poss9 = new Poss(trmm,false);
  PSet *set9 = new PSet(poss9);

  Gemm *gemm = new Gemm(DMLAYER, NORMAL, CONJTRANS, COEFONE, COEFONE, COMPLEX);
  gemm->AddInputs(6, split, 6, split, 7, set9->OutTun(0), 0);
  Poss *poss1 = new Poss(gemm,false);
  PSet *set1 = new PSet(poss1);

  Hetrmm *hetrmm = new Hetrmm(DMLAYER, UPPER);
  hetrmm->AddInput(split, 4);
  Poss *poss10 = new Poss(hetrmm,false);
  PSet *set10 = new PSet(poss10);

  Herk *tri2 = new Herk(DMLAYER, UPPER, NORMAL, COEFONE, COEFONE, COMPLEX);
  tri2->AddInputs(4,
		 split, 7,
		  set10->OutTun(0), 0);
  Poss *poss8 = new Poss(tri2,false);
  PSet *set8 = new PSet(poss8);


  Combine *comA3 = new Combine(PARTDIAG, POSSTUNOUT);
  comA3->AddInput(split,0);
  comA3->AddInput(split,1);
  comA3->AddInput(split,2);
  comA3->AddInput(set1->OutTun(0),0);
  comA3->AddInput(set8->OutTun(0),0);
  comA3->AddInput(split,5);
  comA3->AddInput(split,6);
  comA3->AddInput(split,7);
  comA3->AddInput(split,8);
  comA3->AddInput(split,9);

  comA3->CopyTunnelInfo(split);

  Poss *loopPoss = new Poss(1, comA3);
#if DOELEM
  Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  return loop;
#else
  throw;
#endif
}

Loop* HetrmmAlgVar3Lower(Node *in, unsigned int num)
{
  Split *split = new Split(PARTDIAG, POSSTUNIN, true);
  split->AddInput(in,num);
  split->SetUpStats(FULLUP, FULLUP,
		      FULLUP, NOTUP);

  // L_{21} = L_{22}^T L_{21}
  Trxm *trmm = new Trxm(false, DMLAYER, LEFT, LOWER, NONUNIT, CONJTRANS, COEFONE, COMPLEX);
  trmm->AddInputs(4, split, 8, split, 5);
  Poss *poss9 = new Poss(trmm,false);
  PSet *set9 = new PSet(poss9);

  // L_{11} = L_{11}^T L_{11}
  Hetrmm *hetrmm = new Hetrmm(DMLAYER, LOWER);
  hetrmm->AddInput(split, 4);
  Poss *poss10 = new Poss(hetrmm,false);
  PSet *set10 = new PSet(poss10);

  // L_{11} = L_{21}^TL_{21} + L_{11}
  Herk *tri2 = new Herk(DMLAYER, LOWER, CONJTRANS, COEFONE, COEFONE, COMPLEX);
  tri2->AddInputs(4, 
		  split, 5,
		  set10->OutTun(0), 0);
  Poss *poss8 = new Poss(tri2,false);
  PSet *set8 = new PSet(poss8);



  Combine *comA3 = new Combine(PARTDIAG, POSSTUNOUT);
  comA3->AddInput(split,0);
  comA3->AddInput(split,1);
  comA3->AddInput(split,2);
  comA3->AddInput(split,3);
  comA3->AddInput(set8->OutTun(0),0);
  comA3->AddInput(set9->OutTun(0),0);
  comA3->AddInput(split,6);
  comA3->AddInput(split,7);
  comA3->AddInput(split,8);
  comA3->AddInput(split,9);


  comA3->CopyTunnelInfo(split);

  Poss *loopPoss = new Poss(1, comA3);
#if DOELEM
  Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  return loop;
#else
  throw;
#endif
}


Loop* HetrmmAlgVar3Upper(Node *in, unsigned int num)
{
  Split *split = new Split(PARTDIAG, POSSTUNIN, true);
  split->AddInput(in, num);
  split->SetUpStats(FULLUP, FULLUP,
		    FULLUP, NOTUP);


  Hetrmm *hetrmm = new Hetrmm(DMLAYER, UPPER);
  hetrmm->AddInput(split, 4);
  Poss *poss10 = new Poss(hetrmm,false);
  PSet *set10 = new PSet(poss10);

  Herk *tri2 = new Herk(DMLAYER, UPPER, NORMAL, COEFONE, COEFONE, COMPLEX);
  tri2->AddInputs(4,
		 split, 7,
		  set10->OutTun(0), 0);
  Poss *poss8 = new Poss(tri2,false);
  PSet *set8 = new PSet(poss8);

  Trxm *trmm = new Trxm(false, DMLAYER, RIGHT, UPPER, NONUNIT, CONJTRANS, COEFONE, COMPLEX);
  trmm->AddInputs(4, split, 8, split, 7);
  Poss *poss9 = new Poss(trmm,false);
  PSet *set9 = new PSet(poss9);

  Combine *comA3 = new Combine(PARTDIAG, POSSTUNOUT);
  comA3->AddInput(split,0);
  comA3->AddInput(split,1);
  comA3->AddInput(split,2);
  comA3->AddInput(split,3);
  comA3->AddInput(set8->OutTun(0),0);
  comA3->AddInput(split,5);
  comA3->AddInput(split,6);
  comA3->AddInput(set9->OutTun(0),0);
  comA3->AddInput(split,8);
  comA3->AddInput(split,9);

  comA3->CopyTunnelInfo(split);

  Poss *loopPoss = new Poss(1, comA3);
#if DOELEM
  Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);
  return loop;
#else
  throw;
#endif
}
#endif
