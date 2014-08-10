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


#include "twoSidedTrxm.h"
#ifndef SKIPTWOSIDED
#include "elemRedist.h"
#include "realLoop.h"
#include "blas.h"
#include <stdio.h>
#include <string>
#include <iomanip>
#include <sstream>
#include "loopSupport.h"
#include "helperNodes.h"
#include "elemRedist.h"


void TwoSidedTrxm::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<2,1>::Duplicate(orig, shallow, possMerging);
  m_tri = ((TwoSidedTrxm*)orig)->m_tri;
  m_invert = ((TwoSidedTrxm*)orig)->m_invert;
}

NodeType TwoSidedTrxm::GetType() const 
{
  if (m_invert)
      return "TwoSidedTrsm" + TriToStr(m_tri) + " " + LayerNumToStr(GetLayer());
  else 
    return "TwoSidedTrmm" + TriToStr(m_tri) + " " + LayerNumToStr(GetLayer());
}

Phase TwoSidedTrxm::MaxPhase() const 
{  switch(GetLayer()) {
  case (ABSLAYER):
#if DODPPHASE
  case (DMLAYER):
    return DPPHASE;
  case (SMLAYER):
    return NUMPHASES;
#else
  case (S1LAYER):
    return SR1PHASE;
  case (S2LAYER):
    return NUMPHASES;
#endif
  default:
    throw;
  }
}

bool TwoSidedTrxm::ShouldCullDP() const 
{
#if DOELEM
  return m_layer == DMLAYER;
#else
  throw;
#endif
}

void TwoSidedTrxm::FlattenCore(ofstream &out) const
{
  DLAOp<2,1>::FlattenCore(out);
  WRITE(m_tri);
  WRITE(m_invert);
}

void TwoSidedTrxm::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<2,1>::UnflattenCore(in, info);
  READ(m_tri);
  READ(m_invert);
}

void TwoSidedTrxm::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();
#if DOELEM
    if (m_layer == DMLAYER) {
      if (InputDataType(0).m_dist != D_MC_MR)
	throw;
      else if (InputDataType(1).m_dist != D_MC_MR)
	throw;
    }
    else if (m_layer == SMLAYER) {
      if (InputDataType(0).m_dist != D_STAR_STAR)
	throw;
      else if (InputDataType(1).m_dist != D_STAR_STAR)
	throw;
    }
#endif

    m_cost = ZERO;
  }
}

void TwoSidedTrxm::PrintCode(IndStream &out)
{
  out.Indent();
#if DOBLIS
  if (m_layer == S2LAYER) {
    *out << "libflame_Hegst( &"
	 << GetInputName(1).str()
      << ", &" << GetInputName(0).str() << " );\n";
    return;
  }
#elif DOELEM
  if (m_layer == DMLAYER) {
    if (m_invert)
      *out << "TwoSidedTrsm( ";
    else
      *out << "TwoSidedTrmm( ";
  }
  else if (m_layer == SMLAYER) {
    if (m_invert)
      *out << "internal::TwoSidedTrsm( ";
    else
      *out << "internal::TwoSidedTrmm( ";
  }
#endif
  else
    throw;

    *out << TriToStr(m_tri) 
	 << ",\n" << out.Tabs(1) << GetInputName(1).str()
	 << ",\n" << out.Tabs(1) << GetInputName(0).str() << " );\n";
}

#if DODPPHASE
bool DistTwoSidedTrxmToLocalTwoSidedTrxm::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == TwoSidedTrxm::GetClass()) {
    const TwoSidedTrxm *hegst = (TwoSidedTrxm*)node;
    return hegst->GetLayer() == DMLAYER;
  }
  return false;
}

void DistTwoSidedTrxmToLocalTwoSidedTrxm::Apply(Node *node) const
{
  TwoSidedTrxm *hegst = (TwoSidedTrxm*)node;
  RedistNode *redist1 = new RedistNode(D_STAR_STAR);
  RedistNode *redist2 = new RedistNode(D_STAR_STAR);
  TwoSidedTrxm *node2 = new TwoSidedTrxm(SMLAYER, hegst->m_invert, hegst->m_tri);
  RedistNode *node3 = new RedistNode(D_MC_MR);
  redist1->AddInput(node->Input(0),node->InputConnNum(0));
  redist2->AddInput(node->Input(1),node->InputConnNum(1));
  node2->AddInput(redist1,0);
  node2->AddInput(redist2,0);
  node3->AddInput(node2,0);
  node->m_poss->AddNodes(4, redist1, redist2, node2, node3);
  node->RedirectChildren(node3,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}
#endif


string TwoSidedTrxmLoopExp::GetType() const 
{
  std::stringstream str;
  str << "TwoSidedTrxmLoop expansion variant "
      << m_varNum;
  return str.str();

}

bool TwoSidedTrxmLoopExp::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != TwoSidedTrxm::GetClass()) 
    return false;
  const TwoSidedTrxm *hegst = (TwoSidedTrxm*)node;
  if (hegst->GetLayer() != m_fromLayer)
    return false;
  return true;
}

void TwoSidedTrxmLoopExp::Apply(Node *node) const
{
  TwoSidedTrxm *hegst = (TwoSidedTrxm*)node;
  RealLoop *loop;
  if (hegst->m_tri != LOWER)
    throw;

  if (!hegst->m_invert) {
    if (m_varNum == 1) {
      loop = TwoSidedTrmmLowerVar1Alg(hegst->Input(0), hegst->InputConnNum(0), 
				   hegst->Input(1), hegst->InputConnNum(1), 
				   m_toLayerBLAS, m_toLayerTwoSidedTrxm);
    }
    else if (m_varNum == 2) {
      loop = TwoSidedTrmmLowerVar2Alg(hegst->Input(0), hegst->InputConnNum(0), 
				   hegst->Input(1), hegst->InputConnNum(1), 
				   m_toLayerBLAS, m_toLayerTwoSidedTrxm);
    }
    else if (m_varNum == 4) {
      loop = TwoSidedTrmmLowerVar4Alg(hegst->Input(0), hegst->InputConnNum(0), 
				   hegst->Input(1), hegst->InputConnNum(1), 
				   m_toLayerBLAS, m_toLayerTwoSidedTrxm);
    }
    else
      throw;
  }
  else {
    if (m_varNum == 1) {
      loop = TwoSidedTrsmLowerVar1Alg(hegst->Input(0), hegst->InputConnNum(0), 
				   hegst->Input(1), hegst->InputConnNum(1), 
				   m_toLayerBLAS, m_toLayerTwoSidedTrxm);
    }
    else if (m_varNum == 2) {
      loop = TwoSidedTrsmLowerVar2Alg(hegst->Input(0), hegst->InputConnNum(0), 
				   hegst->Input(1), hegst->InputConnNum(1), 
				   m_toLayerBLAS, m_toLayerTwoSidedTrxm);
    }
    else if (m_varNum == 4) {
      loop = TwoSidedTrsmLowerVar4Alg(hegst->Input(0), hegst->InputConnNum(0), 
				   hegst->Input(1), hegst->InputConnNum(1), 
				   m_toLayerBLAS, m_toLayerTwoSidedTrxm);
    }
    else
      throw;
  }

  node->m_poss->AddPSet(loop);
  
  node->RedirectChildren(loop->OutTun(0),0);
  node->m_poss->DeleteChildAndCleanUp(node);
} 



RealLoop* TwoSidedTrsmLowerVar1Alg(
			     Node *Lin, ConnNum Lnum,
			     Node *Ain, ConnNum Anum,
			     Layer layerBLAS, Layer layerTwoSidedTrxm)
{
  SplitSingleIter *splitA = new SplitSingleIter(PARTDIAG, POSSTUNIN, true);
  splitA->AddInput(Ain, Anum);
  splitA->SetUpStats(FULLUP, FULLUP,
		     NOTUP, NOTUP);

  SplitSingleIter *splitL = new SplitSingleIter(PARTDIAG, POSSTUNIN);
  splitL->AddInput(Lin, Lnum);
  splitL->SetAllStats(FULLUP);

#if DOELEM
  TempVarNode *Yin = new TempVarNode(D_MC_MR, "Y10");
#else
  TempVarNode *Yin = new TempVarNode("Y10");
#endif
  Yin->SetLayer(layerBLAS);
  Yin->AddInput(splitA, 5);

  // Y10 = B10 * A00;
  Hemm *hemm = new Hemm(layerBLAS, RIGHT, LOWER, COEFONE, COEFZERO, COMPLEX);
  hemm->AddInput(splitA,0);
  hemm->AddInput(splitL,1);
  hemm->AddInput(Yin,0);

  // A10 = A10 * inv( tril( B00 )' );   
  Trxm *trsm = new Trxm(true, layerBLAS, RIGHT, LOWER, NONUNIT, CONJTRANS, COEFONE, COMPLEX);
  trsm->AddInput(splitL,0);
  trsm->AddInput(splitA,1);

  // A10 = A10 - 1/2 * Y10;
  Axpy *axpy1 = new Axpy(layerBLAS, COEFNEGONEHALF);
  axpy1->AddInput(hemm,0);
  axpy1->AddInput(trsm,0);

  // A11 = A11 - A10 * B10' - B10 * A10';
  Her2k *her2k = new Her2k(layerBLAS, LOWER, NORMAL, COEFNEGONE, COEFONE, COMPLEX);
  her2k->AddInput(axpy1,0);
  her2k->AddInput(splitL,1);
  her2k->AddInput(splitA,4);

  // A11 = inv( tril( B11 ) ) * A11 * inv( tril( B11 )' );
  TwoSidedTrxm *hegst = new TwoSidedTrxm(layerTwoSidedTrxm, true, LOWER);
  hegst->AddInput(splitL,4);
  hegst->AddInput(her2k,0);

  // A10 = A10 - 1/2 * Y10;
  Axpy *axpy2 = new Axpy(layerBLAS, COEFNEGONEHALF);
  axpy2->AddInput(hemm,0);
  axpy2->AddInput(axpy1,0);

  // A10 = inv( tril( B11 ) ) * A10;
  Trxm *trsm2 = new Trxm(true, layerBLAS, LEFT, LOWER, NONUNIT, NORMAL, COEFONE, COMPLEX);
  trsm2->AddInput(splitL,4);
  trsm2->AddInput(axpy2,0);

  CombineSingleIter *comA = new CombineSingleIter(PARTDIAG, POSSTUNOUT);
  comA->AddInput(splitA,0);
  comA->AddInput(trsm2,0);
  comA->AddInput(splitA,2);
  comA->AddInput(splitA,3);
  comA->AddInput(hegst,0);
  comA->AddInput(splitA,5);
  comA->AddInput(splitA,6);
  comA->AddInput(splitA,7);
  comA->AddInput(splitA,8);
  comA->AddInput(splitA,9);
  
  comA->CopyTunnelInfo(splitA);

  CombineSingleIter *comL = splitL->CreateMatchingCombine(0);

  Poss *loopPoss = new Poss(2,
			    comA,
			    comL);
  RealLoop *loop;
#if DOELEM
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  #else
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISOUTERBS);
#endif

  return loop;
}

RealLoop* TwoSidedTrsmLowerVar2Alg(
			     Node *Lin, ConnNum Lnum,
			     Node *Ain, ConnNum Anum,
			     Layer layerBLAS, Layer layerTwoSidedTrxm)
{
  SplitSingleIter *splitA = new SplitSingleIter(PARTDIAG, POSSTUNIN, true);
  splitA->AddInput(Ain, Anum);
  splitA->SetUpStats(FULLUP, FULLUP,
		     PARTUP, NOTUP);

  SplitSingleIter *splitL = new SplitSingleIter(PARTDIAG, POSSTUNIN);
  splitL->AddInput(Lin, Lnum);
  splitL->SetAllStats(FULLUP);

#if DOELEM
  TempVarNode *Yin = new TempVarNode(D_MC_MR, "Y10");
#else
  TempVarNode *Yin = new TempVarNode("Y10");
#endif
  Yin->SetLayer(layerBLAS);
  Yin->AddInput(splitA, 5);

  // Y10 = 1/2 * L10 * A00
  Hemm *hemm = new Hemm(layerBLAS, RIGHT, LOWER, COEFONEHALF, COEFZERO, COMPLEX);
  hemm->AddInput(splitA,0);
  hemm->AddInput(splitL,1);
  hemm->AddInput(Yin,0);

  // A10 = A10 - Y10; 
  Axpy *axpy1 = new Axpy(layerBLAS, COEFNEGONE);
  axpy1->AddInput(hemm,0);
  axpy1->AddInput(splitA,1);

  // A11 = A11 - A10 * L10' - B10 * L10';
  Her2k *her2k = new Her2k(layerBLAS, LOWER, NORMAL, COEFNEGONE, COEFONE, COMPLEX);
  her2k->AddInput(axpy1,0);
  her2k->AddInput(splitL,1);
  her2k->AddInput(splitA,4);

  // A11 = inv( tril( B11 ) ) * A11 * inv( tril( B11 )' );
  TwoSidedTrxm *hegst = new TwoSidedTrxm(layerTwoSidedTrxm, true, LOWER);
  hegst->AddInput(splitL,4);
  hegst->AddInput(her2k,0);

  // A21 = A21 - A20 * B10'; 
  Gemm *gemm = new Gemm(layerBLAS, NORMAL, TRANS, COEFNEGONE, COEFONE, COMPLEX);
  gemm->AddInput(splitA,2);
  gemm->AddInput(splitL,1);
  gemm->AddInput(splitA,5);

  // A21 = A21 * inv( tril( B11 )' );
  Trxm *trsm = new Trxm(true, layerBLAS, RIGHT, LOWER, NONUNIT, CONJTRANS, COEFONE, COMPLEX);
  trsm->AddInput(splitL,4);
  trsm->AddInput(gemm,0);

  // A10 = A10 - Y10;
  Axpy *axpy2 = new Axpy(layerBLAS, COEFNEGONE);
  axpy2->AddInput(hemm,0);
  axpy2->AddInput(axpy1,0);

  // A10 = inv( tril( B11 ) ) * A10;
  Trxm *trsm2 = new Trxm(true, layerBLAS, LEFT, LOWER, NONUNIT, NORMAL, COEFONE,COMPLEX);
  trsm2->AddInput(splitL,4);
  trsm2->AddInput(axpy2,0);

  CombineSingleIter *comA = new CombineSingleIter(PARTDIAG, POSSTUNOUT);
  comA->AddInput(splitA,0);
  comA->AddInput(trsm2,0);
  comA->AddInput(splitA,2);
  comA->AddInput(splitA,3);
  comA->AddInput(hegst,0);
  comA->AddInput(trsm,0);
  comA->AddInput(splitA,6);
  comA->AddInput(splitA,7);
  comA->AddInput(splitA,8);
  comA->AddInput(splitA,9);
  
  comA->CopyTunnelInfo(splitA);

  CombineSingleIter *comL = splitL->CreateMatchingCombine(0);

  Poss *loopPoss = new Poss(2,
			    comA,
			    comL);
  RealLoop *loop;
#if DOELEM
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  #else
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISOUTERBS);
#endif

  return loop;
}


RealLoop* TwoSidedTrsmLowerVar4Alg(
			     Node *Lin, ConnNum Lnum,
			     Node *Ain, ConnNum Anum,
			     Layer layerBLAS, Layer layerTwoSidedTrxm)
{
  SplitSingleIter *splitA = new SplitSingleIter(PARTDIAG, POSSTUNIN, true);
  splitA->AddInput(Ain, Anum);
  splitA->SetUpStats(PARTUP, FULLUP,
		     PARTUP, PARTUP);

  SplitSingleIter *splitL = new SplitSingleIter(PARTDIAG, POSSTUNIN);
  splitL->AddInput(Lin, Lnum);
  splitL->SetAllStats(FULLUP);

#if DOELEM
  TempVarNode *Yin = new TempVarNode(D_MC_MR, "Y21");
#else
  TempVarNode *Yin = new TempVarNode("Y21");
#endif
  Yin->SetLayer(layerBLAS);
  Yin->AddInput(splitA, 5);

  Trxm *trsm = new Trxm(true, layerBLAS, LEFT, LOWER, NONUNIT, NORMAL, COEFONE,COMPLEX);
  trsm->AddInput(splitL,4);
  trsm->AddInput(splitA,1);

  TwoSidedTrxm *hegst = new TwoSidedTrxm(layerTwoSidedTrxm, true, LOWER);
  hegst->AddInput(splitL,4);
  hegst->AddInput(splitA,4);

  Gemm *gemm = new Gemm(layerBLAS, NORMAL, NORMAL, COEFNEGONE, COEFONE, COMPLEX);
  gemm->AddInput(splitL,5);
  gemm->AddInput(trsm,0);
  gemm->AddInput(splitA,2);

  Hemm *hemm = new Hemm(layerBLAS, RIGHT, LOWER, COEFONE, COEFZERO, COMPLEX);
  hemm->AddInput(hegst,0);
  hemm->AddInput(splitL,5);
  hemm->AddInput(Yin,0);

  Trxm *trsm2 = new Trxm(true, layerBLAS, RIGHT, LOWER, NONUNIT, CONJTRANS, COEFONE,COMPLEX);
  trsm2->AddInput(splitL,4);
  trsm2->AddInput(splitA,5);

  Axpy *axpy1 = new Axpy(layerBLAS, COEFONEHALF);
  axpy1->AddInput(hemm,0);
  axpy1->AddInput(trsm2,0);
  
  Her2k *her2k = new Her2k(layerBLAS, LOWER, NORMAL, COEFNEGONE, COEFONE, COMPLEX);
  her2k->AddInput(splitL,5);
  her2k->AddInput(axpy1,0);
  her2k->AddInput(splitA,8);

  Axpy *axpy2 = new Axpy(layerBLAS, COEFONEHALF);
  axpy2->AddInput(hemm,0);
  axpy2->AddInput(axpy1,0);

  CombineSingleIter *comA = new CombineSingleIter(PARTDIAG, POSSTUNOUT);
  comA->AddInput(splitA,0);
  comA->AddInput(trsm,0);
  comA->AddInput(gemm,0);
  comA->AddInput(splitA,3);
  comA->AddInput(hegst,0);
  comA->AddInput(axpy2,0);
  comA->AddInput(splitA,6);
  comA->AddInput(splitA,7);
  comA->AddInput(her2k,0);
  comA->AddInput(splitA,9);
  
  comA->CopyTunnelInfo(splitA);

  CombineSingleIter *comL = new CombineSingleIter(PARTDIAG, POSSTUNOUT);
  comL->AddInput(splitL,0);
  comL->AddInput(splitL,1);
  comL->AddInput(splitL,2);
  comL->AddInput(splitL,3);
  comL->AddInput(splitL,4);
  comL->AddInput(splitL,5);
  comL->AddInput(splitL,6);
  comL->AddInput(splitL,7);
  comL->AddInput(splitL,8);
  comL->AddInput(splitL,9);

  comL->CopyTunnelInfo(splitL);

  Poss *loopPoss = new Poss(2,
			    comA,
			    comL);
  RealLoop *loop;
  #if DOELEM
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  #else
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISOUTERBS);
#endif

  return loop;
}

RealLoop* TwoSidedTrmmLowerVar1Alg(
			     Node *Lin, ConnNum Lnum,
			     Node *Ain, ConnNum Anum,
			     Layer layerBLAS, Layer layerTwoSidedTrxm)
{
  SplitSingleIter *splitA = new SplitSingleIter(PARTDIAG, POSSTUNIN, true);
  splitA->AddInput(Ain, Anum);
  splitA->SetUpStats(FULLUP, FULLUP,
		     PARTUP, NOTUP);

  SplitSingleIter *splitL = new SplitSingleIter(PARTDIAG, POSSTUNIN);
  splitL->AddInput(Lin, Lnum);
  splitL->SetAllStats(FULLUP);

#if DOELEM
  TempVarNode *Yin = new TempVarNode(D_MC_MR, "Y21");
#else
  TempVarNode *Yin = new TempVarNode("Y21");
#endif
  Yin->SetLayer(layerBLAS);
  Yin->AddInput(splitA, 7);

  // Y21 = A22 * B21;
  Hemm *hemm = new Hemm(layerBLAS, RIGHT, LOWER, COEFONE, COEFZERO, COMPLEX);
  hemm->AddInput(splitA,8);
  hemm->AddInput(splitL,5);
  hemm->AddInput(Yin,0);

  // A21 = A21 * tril( B11 );
  Trxm *trmm = new Trxm(false, layerBLAS, RIGHT, LOWER, NONUNIT, NORMAL, COEFONE, COMPLEX);
  trmm->AddInput(splitL,4);
  trmm->AddInput(splitA,5);

  // A21 = A21 + 1/2 * Y21;
  Axpy *axpy1 = new Axpy(layerBLAS, COEFONEHALF);
  axpy1->AddInput(hemm,0);
  axpy1->AddInput(trmm,0);

  // A11 = tril( B11 )' * A11 * tril( B11 );  
  TwoSidedTrxm *hegst = new TwoSidedTrxm(layerTwoSidedTrxm, true, LOWER);
  hegst->AddInput(splitL,4);
  hegst->AddInput(splitA,4);

  // A11 = A11 + A21' * B21 + B21' * A21;  
  Her2k *her2k = new Her2k(layerBLAS, LOWER, NORMAL, COEFNEGONE, COEFONE, COMPLEX);
  her2k->AddInput(axpy1,0);
  her2k->AddInput(splitL,5);
  her2k->AddInput(hegst,0);

  // A21 = A21 + 1/2 * Y21;
  Axpy *axpy2 = new Axpy(layerBLAS, COEFONEHALF);
  axpy2->AddInput(hemm,0);
  axpy2->AddInput(axpy1,0);

  // A21 = tril( B22 )' * A21; 
  Trxm *trmm2 = new Trxm(false, layerBLAS, LEFT, LOWER, NONUNIT, CONJTRANS, COEFONE, COMPLEX);
  trmm2->AddInput(splitL,8);
  trmm2->AddInput(axpy2,0);

  CombineSingleIter *comA = new CombineSingleIter(PARTDIAG, POSSTUNOUT);
  comA->AddInput(splitA,0);
  comA->AddInput(splitA,1);
  comA->AddInput(splitA,2);
  comA->AddInput(splitA,3);
  comA->AddInput(her2k,0);
  comA->AddInput(trmm2,0);
  comA->AddInput(splitA,6);
  comA->AddInput(splitA,7);
  comA->AddInput(splitA,8);
  comA->AddInput(splitA,9);
  
  comA->CopyTunnelInfo(splitA);

  CombineSingleIter *comL = splitL->CreateMatchingCombine(0);

  Poss *loopPoss = new Poss(2,
			    comA,
			    comL);
  RealLoop *loop;
#if DOELEM
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  #else
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISOUTERBS);
#endif

  return loop;
}

RealLoop* TwoSidedTrmmLowerVar2Alg(
			     Node *Lin, ConnNum Lnum,
			     Node *Ain, ConnNum Anum,
			    Layer layerBLAS, Layer layerTwoSidedTrxm)
{
  SplitSingleIter *splitA = new SplitSingleIter(PARTDIAG, POSSTUNIN, true);
  splitA->AddInput(Ain, Anum);
  splitA->SetUpStats(FULLUP, NOTUP,
		     PARTUP, NOTUP);

  SplitSingleIter *splitL = new SplitSingleIter(PARTDIAG, POSSTUNIN);
  splitL->AddInput(Lin, Lnum);
  splitL->SetAllStats(FULLUP);

#if DOELEM
  TempVarNode *Yin = new TempVarNode(D_MC_MR, "Y21");
#else
  TempVarNode *Yin = new TempVarNode("Y21");
#endif
  Yin->SetLayer(layerBLAS);
  Yin->AddInput(splitA, 5);

  Trxm *trmm = new Trxm(false, layerBLAS, LEFT, LOWER, NONUNIT, CONJTRANS, COEFONE, COMPLEX);
  trmm->AddInputs(4, splitL, 4, splitA, 1);

  Gemm *gemm = new Gemm(layerBLAS, CONJTRANS, NORMAL, COEFONE, COEFONE, COMPLEX);
  gemm->AddInput(splitL,5);
  gemm->AddInput(splitA,2);
  gemm->AddInput(trmm,0);

  Hemm *hemm = new Hemm(layerBLAS, LEFT, LOWER, COEFONE, COEFZERO, COMPLEX);
  hemm->AddInput(splitA,8);
  hemm->AddInput(splitL,5);
  hemm->AddInput(Yin,0);

  Trxm *trmm2 = new Trxm(false, layerBLAS, RIGHT, LOWER, NONUNIT, NORMAL, COEFONE, COMPLEX);
  trmm2->AddInput(splitL,4);
  trmm2->AddInput(splitA,5);

  Axpy *axpy1 = new Axpy(layerBLAS, COEFONEHALF);
  axpy1->AddInput(hemm,0);
  axpy1->AddInput(trmm2,0);


  TwoSidedTrxm *hegst = new TwoSidedTrxm(layerTwoSidedTrxm, false, LOWER);
  hegst->AddInput(splitL,4);
  hegst->AddInput(splitA,4);

  Her2k *her2k = new Her2k(layerBLAS, LOWER, CONJTRANS, COEFONE, COEFONE, COMPLEX);
  her2k->AddInput(axpy1,0);
  her2k->AddInput(splitL,5);
  her2k->AddInput(hegst,0);

  Axpy *axpy2 = new Axpy(layerBLAS, COEFONEHALF);
  axpy2->AddInput(hemm,0);
  axpy2->AddInput(axpy1,0);


  CombineSingleIter *comA = new CombineSingleIter(PARTDIAG, POSSTUNOUT);
  comA = splitA->CreateMatchingCombine(3, 
				       1, gemm, 0,
				       4, her2k, 0,
				       5, axpy2, 0);

  CombineSingleIter *comL = new CombineSingleIter(PARTDIAG, POSSTUNOUT);
  comL = splitL->CreateMatchingCombine(0);

  Poss *loopPoss = new Poss(2,
			    comA,
			    comL);
  RealLoop *loop;
  #if DOELEM
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  #else
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISOUTERBS);
#endif

  return loop;
}

RealLoop* TwoSidedTrmmLowerVar4Alg(
			    Node *Lin, ConnNum Lnum,
			    Node *Ain, ConnNum Anum,
			    Layer layerBLAS, Layer layerTwoSidedTrxm)
{
  SplitSingleIter *splitA = new SplitSingleIter(PARTDIAG, POSSTUNIN, true);
  splitA->AddInput(Ain, Anum);
  splitA->SetUpStats(PARTUP, FULLUP,
		     PARTUP, NOTUP);


  SplitSingleIter *splitL = new SplitSingleIter(PARTDIAG, POSSTUNIN);
  splitL->AddInput(Lin, Lnum);
  splitL->SetAllStats(FULLUP);

#if DOELEM
  TempVarNode *Yin = new TempVarNode(D_MC_MR, "Y10");
  Yin->SetLayer(ABSLAYER);
#else
  TempVarNode *Yin = new TempVarNode("Y10");
  Yin->SetLayer(S3LAYER);
#endif
  Yin->AddInput(splitA, 1);

  //  InputNode *Yin = new InputNode("Y10 input", bigSize, BLIS_KC_BSVAL, "Y10");

  Hemm *hemm = new Hemm(layerBLAS, LEFT, LOWER, COEFONE, COEFZERO, COMPLEX);
  hemm->AddInput(splitA,4);
  hemm->AddInput(splitL,1);
  hemm->AddInput(Yin,0);

  Axpy *axpy1 = new Axpy(layerBLAS, COEFONEHALF);
  axpy1->AddInput(hemm,0);
  axpy1->AddInput(splitA,1);

  Her2k *her2k = new Her2k(layerBLAS, LOWER, CONJTRANS, COEFONE, COEFONE, COMPLEX);
  her2k->AddInput(axpy1,0);
  her2k->AddInput(splitL,1);
  her2k->AddInput(splitA,0);

  Axpy *axpy2 = new Axpy(layerBLAS, COEFONEHALF);
  axpy2->AddInput(hemm,0);
  axpy2->AddInput(axpy1,0);

  Trxm *trmm = new Trxm(false, layerBLAS, LEFT, LOWER, NONUNIT, CONJTRANS, COEFONE, COMPLEX);
  trmm->AddInputs(4, splitL, 4, axpy2, 0);
  
  TwoSidedTrxm *hegst = new TwoSidedTrxm(layerTwoSidedTrxm, false, LOWER);
  hegst->AddInput(splitL,4);
  hegst->AddInput(splitA,4);

  Gemm *gemm = new Gemm(layerBLAS, NORMAL, NORMAL, COEFONE, COEFONE, COMPLEX);
  gemm->AddInput(splitA, 5);
  gemm->AddInput(splitL, 1);
  gemm->AddInput(splitA, 2);

  Trxm *trmm2 = new Trxm(false, layerBLAS, RIGHT, LOWER, NONUNIT, NORMAL, COEFONE, COMPLEX);
  trmm2->AddInput(splitL,4);
  trmm2->AddInput(splitA,5);

  CombineSingleIter *comA;
  comA = splitA->CreateMatchingCombine(5,
				       0,her2k,0,
				       1,trmm,0,
				       2,gemm,0,
				       4,hegst,0,
				       5,trmm2,0);

  CombineSingleIter *comL = splitL->CreateMatchingCombine(0);

  Poss *loopPoss = new Poss(2,
			    comA,
			    comL);
  RealLoop *loop;
  #if DOELEM
    loop = new RealLoop(ELEMLOOP, loopPoss, ElemBS);
  #else
    loop = new RealLoop(BLISLOOP, loopPoss, USEBLISOUTERBS);
#endif
  return loop;
}



bool TwoSidedTrxmLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == TwoSidedTrxm::GetClass()) {
    const TwoSidedTrxm *hegst = (TwoSidedTrxm*)node;
    if (hegst->GetLayer() != m_fromLayer)
      return false;
    return true;
  }
  return false;
  
}

void TwoSidedTrxmLowerLayer::Apply(Node *node) const
{
  TwoSidedTrxm *hegst = (TwoSidedTrxm*)node;
  hegst->SetLayer(m_toLayer);
}

#endif //!SKIPTWOSIDED
