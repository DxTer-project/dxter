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

#include "DLAOp.h"
#include "mvmul.h"
#include "loopSupport.h"
#include "LLDLA.h"

#if DOLLDLA

MVMul::MVMul(Type type, Layer layer)
{
  m_type = type;
  m_layer = layer;
}

void MVMul::PrintCode(IndStream &out)
{
  if (m_layer != LLDLAPRIMITIVELAYER) {
    cout << "ERROR: Attempt to generate code from non-primitive scalar vector multiply\n";
    throw;
  }
  const DataTypeInfo &inInfo = InputDataType(1);
  const Stride rowStride = inInfo.m_rowStride;
  const Stride colStride = inInfo.m_colStride;
  
  out.Indent();
  if (rowStride == NONUNITSTRIDE && colStride == NONUNITSTRIDE) {
    PrintGeneralStride(out);
  } else if (rowStride == UNITSTRIDE && colStride == NONUNITSTRIDE) {
    PrintColStride(out);
  } else if (rowStride == NONUNITSTRIDE && colStride == UNITSTRIDE) {
    PrintRowStride(out);
  } else {
    *out << "ERROR: BAD STRIDE\n";
  }
}


void MVMul::PrintRowStride(IndStream &out)
{
  *out << "row_stride_mmul_2x2_2x1( " <<
    GetInputName(0).str() << ", " <<
    GetInputName(1).str() << ", " <<
    GetInputName(2).str() << ");\n";
}

void MVMul::PrintColStride(IndStream &out)
{
  *out << "col_stride_mmul_2x2_2x1( " <<
    GetInputName(0).str() << ", " <<
    GetInputName(1).str() << ", " <<
    GetInputName(2).str() << ");\n";
}

void MVMul::PrintGeneralStride(IndStream &out)
{
  *out << "gen_stride_mmul_2x2_2x1( " <<
    GetInputName(0).str() << ", " <<
    GetInputName(1).str() << ", " <<
    GetInputName(2).str() << ");\n";
}

Node* MVMul::BlankInst()
{
  return new MVMul(REAL, LLDLAPRIMITIVELAYER);
}

Phase MVMul::MaxPhase() const 
{ 
  switch (m_layer)
    { 
    case(ABSLAYER):
      return LLDLALOOPPHASE;
    case(LLDLAMIDLAYER):
      return LLDLAPRIMPHASE;
    case (LLDLAPRIMITIVELAYER):
      return NUMPHASES; 
    default:
      throw;
    }
}

void MVMul::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3, 1>::Prop();

    if (*GetInputN(0) != *GetInputM(1)) {
      cout << "ERROR: Input dimensions don't match\n";
      throw;
    } else if (*GetInputM(0) != *GetInputM(2)) {
      cout << "ERROR: Input dimensions don't match\n";
      throw;
    } else if (*GetInputN(2) != 1 || *GetInputN(1) != 1) {
      cout << "Error: Input vectors have more than 1 column\n";
      throw;
    }
    
    if (m_layer == LLDLAPRIMITIVELAYER) {
      if (*GetInputM(0) != LLDLA_MU || *GetInputN(0) != LLDLA_MU) {
	cout << "ERROR: Primitive matrix must be 2 x 2\n";
	throw;
      }
      if (*GetInputM(1) != LLDLA_MU) {
	cout << "ERROR: Primitive vector must be 2 x 1\n";
	throw;
      }
      if (*GetInputM(2) != LLDLA_MU) {
	cout << "ERROR: Primitive vector must be 2 x 1\n";
	throw;
      }
    }
    m_cost = ZERO;
  }
}

NodeType MVMul::GetType() const
{
  return "MVMul" +  LayerNumToStr(GetLayer());
}

void MVMul::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<3,1>::Duplicate(orig, shallow, possMerging);
  const MVMul *rhs = (MVMul*)orig;
  m_type = rhs->m_type;
  m_layer = rhs->m_layer;
  return;
}

bool MVMulLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == MVMul::GetClass()) {
    const MVMul *mvmul = (MVMul*) node;
    if (mvmul->GetLayer() != m_fromLayer) {
      return false;
    }
    if (*(mvmul->GetInputM(0)) <= m_bs &&
	*(mvmul->GetInputN(0)) <= m_bs) {
      return true;
    } else {
      return false;
    }
  }
  else {
    throw;
  }
}

void MVMulLowerLayer::Apply(Node *node) const
{
  MVMul *svmul = (MVMul*) node;
  svmul->SetLayer(m_toLayer);
}

string MVMulLowerLayer::GetType() const
{
  return "MVMul lower layer " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer);
}

string MVMulLoopRef::GetType() const
{
  switch (m_dim)
    {
    case (DIMM):
      return "MVMulLoopRef - dim m";
    case (DIMN):
      return "MVMulLoopRef - dim n";
    default:
      throw;
    }  
}

bool MVMulLoopRef::CanApply(const Node *node) const
{
  const MVMul *mul = (MVMul*) node;
  if (mul->GetLayer() != m_fromLayer) {
    return false;
  }

  if (m_dim == DIMM) {
    return !(*(mul->GetInputM(0)) <= BSSizeToSize(m_bs));
  } else {
    return !(*(mul->GetInputN(0)) <= BSSizeToSize(m_bs));
  }
  return false;
}

void MVMulLoopRef::Apply(Node *node) const
{
  if (m_dim == DIMM) {
    ApplyRowSplit(node);
  } else if (m_dim == DIMN) {
    ApplyColSplit(node);
  } else {
    cout << "ERROR: Not refining mvmul loop around m or n dimension\n";
    throw;
  }
}

void MVMulLoopRef::ApplyRowSplit(Node *node) const
{
  MVMul *mul = (MVMul*) node;
  
  // Split A on m dimension
  SplitSingleIter *splitA = new SplitSingleIter(PARTDOWN, POSSTUNIN, true);
  splitA->AddInput(mul->Input(0), mul->InputConnNum(0));
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();

  // Split y (the result vector) on m dimension
  SplitSingleIter *splitY = new SplitSingleIter(PARTDOWN, POSSTUNIN, false);
  splitY->AddInput(mul->Input(2), mul->InputConnNum(2));
  splitY->SetUpStats(FULLUP, FULLUP,
		     NOTUP, NOTUP);

  // Create tunnel for x
  LoopTunnel *xTun = new LoopTunnel(POSSTUNIN);
  xTun->AddInput(mul->Input(1), mul->InputConnNum(1));
  xTun->SetAllStats(FULLUP);

  // Create new mvmul for loop body
  MVMul *newMul = new MVMul(mul->m_type, m_toLayer);
  newMul->SetLayer(m_toLayer);

  // Attach the new multiply to arguments
  newMul->AddInput(splitA, 1);
  newMul->AddInput(xTun, 0);
  newMul->AddInput(splitY, 1);

  // Create outputs
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  CombineSingleIter *comY = splitY->CreateMatchingCombine(1, 1, newMul, 0);
  
  LoopTunnel *xOut = new LoopTunnel(POSSTUNOUT);
  xOut->AddInput(xTun, 0);
  xOut->AddInput(xTun, 1);
  xOut->CopyTunnelInfo(xTun);

  // Create poss
  Poss *loopPoss = new Poss(3, comA, xTun, comY);
  RealLoop *loop = new RealLoop(LLDLALOOP, loopPoss, m_bs);
  loop->SetDimName(m_dim);

  // Set up graph to correctly incorporate this loop
  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(2), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

void MVMulLoopRef::ApplyColSplit(Node *node) const
{
  MVMul *mul = (MVMul*) node;

  // Split a on n dimension
  SplitSingleIter *splitA = new SplitSingleIter(PARTRIGHT, POSSTUNIN, true);
  splitA->AddInput(mul->Input(0), mul->InputConnNum(0));
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();

  // Split x on m dimension
  SplitSingleIter *splitX = new SplitSingleIter(PARTDOWN, POSSTUNIN, false);
  splitX->AddInput(mul->Input(1), mul->InputConnNum(1));
  splitX->SetAllStats(FULLUP);
  splitX->SetIndepIters();

  // Create tunnel for y (the result vector)
  LoopTunnel *tunY = new LoopTunnel(POSSTUNIN);
  tunY->AddInput(mul->Input(2), mul->InputConnNum(2));
  tunY->SetAllStats(PARTUP);

  // Create new mvmul for loop body
  MVMul *newMul = new MVMul(mul->m_type, m_toLayer);
  newMul->SetLayer(m_toLayer);

  // Attach to arguments
  newMul->AddInput(splitA, 1);
  newMul->AddInput(splitX, 1);
  newMul->AddInput(tunY, 0);

  // Create outputs
  CombineSingleIter *comA = splitA->CreateMatchingCombine(0);
  CombineSingleIter *comX = splitX->CreateMatchingCombine(0);
  
  LoopTunnel *outY = new LoopTunnel(POSSTUNOUT);
  outY->AddInput(newMul, 0);
  outY->AddInput(tunY, 1);
  outY->CopyTunnelInfo(tunY);

  // Create the poss
  Poss *loopPoss = new Poss(3, comA, comX, outY);
  RealLoop *loop = new RealLoop(LLDLALOOP, loopPoss, m_bs);
  loop->SetDimName(m_dim);

  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(2), 0);
  node->m_poss->DeleteChildAndCleanUp(node);

}

#endif // DOLLDLA
