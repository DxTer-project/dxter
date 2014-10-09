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
#if DOTENSORS
#include "base.h"
#include "yaxppx.h"
#include "string.h"
#include "helperNodes.h"
#include "loopSupport.h"

using namespace std;

YAxpPx::YAxpPx(Layer layer, Coef alpha, Coef beta, string start, string end)
  : m_alpha(alpha), m_beta(beta), m_permutation(start, end)
{ 
  SetLayer(layer);
}

YAxpPx::YAxpPx(Layer layer, Coef alpha, Coef beta)
  : m_alpha(alpha), m_beta(beta)
{ 
  SetLayer(layer);
}

YAxpPx::YAxpPx(Layer layer, Coef alpha, Coef beta, const Permutation &perm)
  : m_alpha(alpha), m_beta(beta)
{ 
  SetLayer(layer);
  m_permutation = perm;
}

NodeType YAxpPx::GetType() const
{
  return "YAxpPx " + LayerNumToStr(GetLayer());
}

void YAxpPx::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<3,1>::Duplicate(orig, shallow, possMerging);
  const YAxpPx *axpy = (YAxpPx*)orig;
  m_alpha = axpy->m_alpha;
  m_beta = axpy->m_beta;
  m_permutation = axpy->m_permutation;
}

void YAxpPx::FlattenCore(ofstream &out) const
{
  DLAOp<3,1>::FlattenCore(out);
  WRITE(m_alpha);
  WRITE(m_alpha);
  throw;
  //m_permutation
}

void YAxpPx::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLAOp<3,1>::UnflattenCore(in, info);
  throw;
  //m_permutation
}

Phase YAxpPx::MaxPhase() const 
{
  if (m_layer == DMLAYER || m_layer == ABSLAYER)
    return DPTENSORPHASE;
  else if (m_layer == SMLAYER)
    return NUMPHASES;
  else
    throw;
}
/*
bool YAxpPx::ShouldCullDP() const 
{
#if DODPTENSORPHASE
  return m_layer == DMLAYER;
#else
  throw;
#endif
}
*/
bool YAxpPx::DoNotCullDP() const 
{
#if DODPTENSORPHASE
  return m_layer == DMLAYER;
#else
  throw;
#endif
}

void YAxpPx::PrintCode(IndStream &out)
{
  out.Indent();
  *out << "YAxpPx( ";
  out << m_alpha;
  *out << ", " << GetInputName(0).str() << ", ";
  out << m_beta;
  *out  << ", " << GetInputName(1).str() << ", "
	<< PermutationVarName(m_permutation.m_permutation) << ", " 
	<< GetInputName(2).str() << " );\n";
}



void YAxpPx::AddVariables(VarSet &set) const
{
  DLAOp<3,1>::AddVariables(set);
  Var var(PermutationVarType,m_permutation.m_permutation);
  set.insert(var);
}

void YAxpPx::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3,1>::Prop();

    if (m_layer == ABSLAYER || m_layer == DMLAYER) {
      m_cost = 3 * TotalNumberOfElements(0);
      Dim numDims = InputNumDims(0);
      if (InputNumDims(1) != numDims || InputNumDims(2) != numDims)
	throw;
      for (Dim dim = 0; dim < numDims; ++dim) {
	if (*InputLen(0,dim) != *InputLen(2,dim))
	  throw;
	if (*InputLen(0,dim) != *InputLen(1,m_permutation.MapFinishToStart(dim)))
	  throw;
      }
    }
    else if (m_layer == SMLAYER) {
      m_cost = 3 * TotalNumberOfLocalElements(0);
      Dim numDims = InputNumDims(0);
      if (InputNumDims(1) != numDims || InputNumDims(2) != numDims)
	throw;
      if (m_permutation.Size() != numDims)
	throw;
      const DistType in0Type = InputDataType(0).GetEffectiveDist();
      const DistType in1Type = InputDataType(1).GetEffectiveDist();
      const DistType in2Type = InputDataType(2).GetEffectiveDist();
      for (Dim dim = 0; dim < numDims; ++dim) {
	if (in0Type.m_dists[dim] != in2Type.m_dists[dim])
	  throw;
	if (*InputLen(0,dim) != *InputLen(2,dim))
	  throw;

	Dim mapping = m_permutation.MapFinishToStart(dim);
	if (in0Type.m_dists[dim] != in1Type.m_dists[mapping]) {
	  cout << dim << endl;
	  cout << in0Type.PrettyStr() << endl;
	  cout << mapping << endl;
	  cout << in1Type.PrettyStr() << endl;
	  throw;
	}
	if (*InputLen(0,dim) != *InputLen(1,mapping)) {
	  cout << "Input 0:\n";
	  InputLen(0,dim)->Print();
	  cout << "Input 1:\n";
	  InputLen(1,mapping)->Print();
	  throw;
	}
      }
    }
    else {
      cout << LayerNumToStr(m_layer) << endl;
      throw;
    }
  }
}



bool DistYAxpPxToDefaultLocalYAxpPx::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != YAxpPx::GetClass())
    return false;
  if (((YAxpPx*)node)->GetLayer() != DMLAYER)
    return false;
  return true;
}

void DistYAxpPxToDefaultLocalYAxpPx::Apply(Node *node) const
{
  YAxpPx *orig = (YAxpPx*)node;
  YAxpPx *newYAxpPx = new YAxpPx(SMLAYER, orig->m_alpha, orig->m_beta, orig->m_permutation);

  newYAxpPx->AddInput(node->Input(0),node->InputConnNum(0));

  if (orig->m_permutation.HasPerm()) {
    const DataTypeInfo &inputType = orig->InputDataType(0);
    if (inputType.HasPerm())
      throw;
    DistType newType = inputType.GetDist();
    for(Dim dim = 0; dim < newType.m_numDims; ++dim) {
      newType.m_dists[dim] = inputType.GetDist().m_dists[orig->m_permutation.MapFinishToStart(dim)];
    }
    RedistNode *redist = new RedistNode(newType);
    redist->AddInput(node->Input(1),node->InputConnNum(1));

    Poss *poss = new Poss(redist, false);
    RealPSet *set = new RealPSet(poss);
    node->m_poss->AddPSet(set,true,true);
    if (set->m_inTuns.empty())
      throw;
    if (set->m_outTuns.empty())
      throw;

    newYAxpPx->AddInput(set->OutTun(0), 0);
    //    node->m_poss->AddNode(redist);
  }
  else
    newYAxpPx->AddInput(node->Input(1),node->InputConnNum(1));
  
  newYAxpPx->AddInput(node->Input(2),node->InputConnNum(2));

  node->m_poss->AddNode(newYAxpPx);
  
  node->RedirectChildren(newYAxpPx,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

YAxpPxLoopExp::YAxpPxLoopExp(Layer fromLayer, Layer toLayer)
  : 
  m_fromLayer(fromLayer), 
  m_toLayer(toLayer)
{
}

string YAxpPxLoopExp::GetType() const
{
  string str = "YAxpPx Loop Exp " 
    + LayerNumToStr(m_fromLayer)
    + " + " 
    + LayerNumToStr(m_toLayer);
  return str;
}

bool YAxpPxLoopExp::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == YAxpPx::GetClass()) {
    const YAxpPx *axpy = (YAxpPx*)node;
    if (axpy->GetLayer() == m_fromLayer) {
      Dim numDims = axpy->InputNumDims(0);
      for(Dim dim = 0; dim < numDims; ++dim)
	if (!(*(axpy->InputLen(0,dim)) <= TensorBS.GetSize()))
	  return true;
    }
  }
  return false;
}

void YAxpPxLoopExp::Apply(Node *node) const
{
  YAxpPx *axpy = (YAxpPx*)node;

  Dim m_dim = 0;
  Size size = (*(axpy->InputLen(0,0)))[0];

  Dim numDims = axpy->InputNumDims(0);
  for(Dim dim = 1; dim < numDims; ++dim) {
    Size newSize = (*(axpy->InputLen(0,dim)))[0];
    if (newSize > size) {
      m_dim = dim;
      size = newSize;
    }
  }

  if (size <= TensorBS.GetSize())
    throw;
  
  
  NodeConn *connA, *connB, *connC;
  connA = axpy->m_inputs[0];
  connB = axpy->m_inputs[1];
  connC = axpy->m_inputs[2];

  SplitSingleIter *xTun = new SplitSingleIter(m_dim, POSSTUNIN, false);
  xTun->AddInput(connA->m_n, connA->m_num);
  xTun->SetAllStats(FULLUP);
  xTun->SetIndepIters();


  Dim mappedDim = axpy->m_permutation.MapFinishToStart(m_dim);

  SplitSingleIter *pxTun = NULL;
  if (mappedDim != m_dim) {
    pxTun = new SplitSingleIter(mappedDim, POSSTUNIN, false);
    pxTun->AddInput(connB->m_n, connB->m_num);
    pxTun->SetAllStats(FULLUP);
    pxTun->SetIndepIters();
  }


  SplitSingleIter *yTun = new SplitSingleIter(m_dim, POSSTUNIN, true);
  yTun->AddInput(connC->m_n, connC->m_num);
  yTun->SetUpStats(FULLUP,FULLUP,
		     NOTUP,NOTUP);
  yTun->SetIndepIters();
  
  YAxpPx *newAxpy = new YAxpPx(m_toLayer, 
			       axpy->m_alpha, 
			       axpy->m_beta,
			       axpy->m_permutation);
  newAxpy->AddInput(xTun, 1);
  if (mappedDim != m_dim) {
    newAxpy->AddInput(pxTun, 1);
  }
  else {
    newAxpy->AddInput(xTun, 1);
  }
  newAxpy->AddInput(yTun, 1);
  
  CombineSingleIter *xOut = xTun->CreateMatchingCombine(0);
  CombineSingleIter *pxOut = NULL;
  if (mappedDim != m_dim) {
    pxOut = pxTun->CreateMatchingCombine(0);
  }
  CombineSingleIter *yOut = yTun->CreateMatchingCombine(1,
							 1, newAxpy, 0);
					
  if (mappedDim != m_dim) {
    Poss *loopPoss = new Poss(3, xOut, pxOut, yOut);
    RealLoop *loop = new RealLoop(TENSORLOOP, loopPoss, TensorBS);
    node->m_poss->AddPSet(loop);
    node->RedirectChildren(loop->OutTun(2),0);
  }
  else {
    Poss *loopPoss = new Poss(2, xOut, yOut);
    RealLoop *loop = new RealLoop(TENSORLOOP, loopPoss, TensorBS);
    node->m_poss->AddPSet(loop);
    node->RedirectChildren(loop->OutTun(1),0);
  }

  node->m_poss->DeleteChildAndCleanUp(node, false, false, true);
}

#endif // DOTENSORS
