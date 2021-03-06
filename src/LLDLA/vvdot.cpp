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

#include "LLDLA.h"
#include "vvdot.h"
#include "DLAReg.h"
#include "regArith.h"
#include "regLoadStore.h"

#if DOLLDLA

VVDot::VVDot(Layer layer)
{
  m_layer = layer;
}

void VVDot::PrintCode(IndStream &out)
{
  const DataTypeInfo &inInfo = InputDataType(1);
  const Stride rowStride = inInfo.m_rowStride;
  const Stride colStride = inInfo.m_colStride;

  out.Indent();

  if (m_layer == ABSLAYER) {
    if (GetDataType() == REAL_DOUBLE) {
      *out << "simple_mmul( ";
    } else if (GetDataType() == REAL_SINGLE) {
      *out << "simple_mmul_float( ";
    }
    *out << "1, " <<
      "1, " <<
      InputDataType(0).m_numColsVar << ", " <<
      GetInputName(0).str() << ", " <<
      InputDataType(0).m_rowStrideVar << ", " <<
      InputDataType(0).m_colStrideVar << ", " <<
      GetInputName(1).str() << ", " <<
      InputDataType(1).m_rowStrideVar << ", " <<
      InputDataType(1).m_colStrideVar << ", " <<
      GetInputName(2).str() << ", " <<
      InputDataType(2).m_rowStrideVar << ", " <<
      InputDataType(2).m_colStrideVar << " );\n";
    return;
  }

  if (m_layer != LLDLAPRIMITIVELAYER) {
    cout << "ERROR: Attempt to generate code from non-primitive dot product\n";
    LOG_FAIL("replacement for throw call");
  }

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

void VVDot::PrintRowStride(IndStream &out)
{
  *out << "row_stride_mmul_1x2_2x1( " <<
    GetInputName(0).str() << ", " <<
    InputDataType(0).m_rowStrideVar << ", " <<
    GetInputName(1).str() << ", " <<
    InputDataType(1).m_rowStrideVar << ", " <<
    GetInputName(2).str() << ", " <<
    InputDataType(2).m_rowStrideVar << " );\n";
}

void VVDot::PrintColStride(IndStream &out)
{
  *out << "col_stride_mmul_1x2_2x1( " <<
    GetInputName(0).str() << ", " <<
    InputDataType(0).m_colStrideVar << ", " <<
    GetInputName(1).str() << ", " <<
    InputDataType(0).m_colStrideVar << ", " <<
    GetInputName(2).str() << ", " <<
    InputDataType(0).m_colStrideVar << " );\n";
}

void VVDot::PrintGeneralStride(IndStream &out)
{
  *out << "gen_stride_mmul_1x2_2x1( " <<
    GetInputName(0).str() << ", " <<
    GetInputName(1).str() << ", " <<
    GetInputName(2).str() << " );\n";
}

void VVDot::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3, 1>::Prop();
    if (*GetInputM(2) != 1 || *GetInputN(2) != 1) {
      cout << "ERROR: Result of dot product must be a scalar\n";
    }
    
    if (*GetInputN(0) != *GetInputM(1)) {
      cout << "ERROR: Dot product argument dimensions do not match\n";
    }

    if (m_layer == LLDLAPRIMITIVELAYER) {
      //      if (*GetInputN(0) != m_regWidth) {
      //cout << "ERROR: First argument to primitive dot prod must be 1x2\n";
      //LOG_FAIL("replacement for throw call");
      //      } else if (*GetInputN(1) != 1) {
      //cout << "ERROR: Second argument to primitive dot prod must be 2x1\n";
      //LOG_FAIL("replacement for throw call");
      //      }
    } else {
      if (*GetInputM(0) != 1) {
	cout << "ERROR: First argument to dot prod must be a row vector\n";
	LOG_FAIL("replacement for throw call");
      } else if (*GetInputN(1) != 1) {
	cout << "ERROR: Second argument to primitive dot prod must be a column vector\n";
	LOG_FAIL("replacement for throw call");
      }
    }
    if (m_layer == ABSLAYER) {
      m_cost = TWO * GetInputN(0)->Sum() - 1;
    } else {
      m_cost = ZERO;
    }
  }
  return;
}

Node* VVDot::BlankInst()
{
  return new VVDot(LLDLAPRIMITIVELAYER);
}

Phase VVDot::MaxPhase() const 
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
      LOG_FAIL("replacement for throw call");
      throw;
    }
}

NodeType VVDot::GetType() const
{
  return "VVDot" +  LayerNumToStr(GetLayer());
}

VVDotLoopRef::VVDotLoopRef(Layer fromLayer, Layer toLayer, BSSize bs)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_bs = bs;
}

string VVDotLoopRef::GetType() const
{
  return "VVDotLoopRef";
}

bool VVDotLoopRef::CanApply(const Node *node) const
{
  const VVDot *dot = static_cast<const VVDot*>(node);
  if (dot->GetLayer() != m_fromLayer) {
    return false;
  }
  if (*(dot->GetInputN(0)) <= m_bs.GetSize()) {
    return false;
  } else if (*(dot->GetInputM(1)) <= m_bs.GetSize()) {
    return false;
  } else {
    return true;
  }
}

void VVDotLoopRef::Apply(Node *node) const
{
  VVDot *dot = static_cast<VVDot*>(node);
  Type dataType = dot->GetDataType();

  // Split for row vector
  SplitSingleIter *splitRow = new SplitSingleIter(PARTRIGHT, POSSTUNIN, true);
  splitRow->AddInput(dot->Input(0), dot->InputConnNum(0));
  splitRow->SetAllStats(FULLUP);
  splitRow->SetIndepIters();

  // Split for col vector
  SplitSingleIter *splitCol = new SplitSingleIter(PARTDOWN, POSSTUNIN, false);
  splitCol->AddInput(dot->Input(1), dot->InputConnNum(1));
  splitCol->SetAllStats(FULLUP);
  splitCol->SetIndepIters();

  LoopTunnel *scalarTun = new LoopTunnel(POSSTUNIN);
  scalarTun->AddInput(dot->Input(2), dot->InputConnNum(2));
  scalarTun->SetAllStats(PARTUP);

  // Create new dot product for interior of loop
  VVDot *newDot = new VVDot(m_toLayer);
  newDot->SetLayer(m_toLayer);

  // Attach inputs to dot product
  newDot->AddInput(splitRow, 1);
  newDot->AddInput(splitCol, 1);
  newDot->AddInput(scalarTun, 0);

  // Create outputs
  CombineSingleIter *rowCom = splitRow->CreateMatchingCombine(0);
  CombineSingleIter *colCom = splitCol->CreateMatchingCombine(0);

  LoopTunnel *scalarTunOut = new LoopTunnel(POSSTUNOUT);
  scalarTunOut->AddInput(newDot, 0);
  scalarTunOut->AddInput(scalarTun, 1);
  scalarTunOut->CopyTunnelInfo(scalarTun);

  // Create the poss
  Poss *loopPoss = new Poss(3, rowCom, colCom, scalarTunOut);
  RealLoop* loop;
  if (dataType == REAL_SINGLE) {
    loop = new RealLoop(LLDLALOOP, loopPoss, LLDLAMuSingle);
  } else if (dataType == REAL_DOUBLE) {
    loop = new RealLoop(LLDLALOOP, loopPoss, LLDLAMuDouble);
  } else {
    cout << "Error: Bad dataType in vadd apply\n";
    LOG_FAIL("replacement for throw call");
    throw;
  }
  loop->SetDimName(DIMK);
  
  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(2), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;
}

VVDotLowerLayer::VVDotLowerLayer(Layer fromLayer, Layer toLayer, Size bs)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_bs = bs;
}

bool VVDotLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == VVDot::GetClass()) {
    const VVDot *vvdot = static_cast<const VVDot*>(node);
    if (vvdot->GetLayer() != m_fromLayer) {
      return false;
    }

    if (*(vvdot->GetInputN(0)) <= m_bs &&
	*(vvdot->GetInputM(1)) <= m_bs) {
      return true;
    } else {
      return false;
    }
  }
  else {
    LOG_FAIL("replacement for throw call");
    throw;
  }
}

void VVDotLowerLayer::Apply(Node *node) const
{
 VVDot *svmul = static_cast<VVDot*>(node);
  svmul->SetLayer(m_toLayer);
  return;
}

string VVDotLowerLayer::GetType() const
{
  return "VVDot lower layer " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer);
}

VVDotToRegArith::VVDotToRegArith(Layer fromLayer, Layer toLayer)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
}

bool VVDotToRegArith::CanApply(const Node* node) const
{
  if (node->GetNodeClass() == VVDot::GetClass()) {
    const VVDot* vvdot = static_cast<const VVDot*>(node);
    //    int regWidth = arch->VecRegWidth(vvdot->GetDataType());
    if (vvdot->GetLayer() != m_fromLayer) {
      return false;
    }
    if (!(vvdot->GetInputN(0)->EvenlyDivisibleBy(vvdot->GetVecRegWidth()))) {
      return false;
    }
    if (*(vvdot->GetInputN(0)) > 0 &&
	*(vvdot->GetInputM(1)) > 0) {
      return true;
    } else {
      return false;
    }
  }
  return false;
}

void VVDotToRegArith::Apply(Node *node) const
{
  VVDot* vvdot = static_cast<VVDot*>(node);
  Type dataType = vvdot->GetDataType();

  // Split A on N dimension
  SplitSingleIter* splitA = new SplitSingleIter(PARTRIGHT, POSSTUNIN, true);
  splitA->AddInput(vvdot->Input(0), vvdot->InputConnNum(0));
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();
  
  // Split B on M dimension
  SplitSingleIter* splitB = new SplitSingleIter(PARTDOWN, POSSTUNIN, false);
  splitB->AddInput(vvdot->Input(1), vvdot->InputConnNum(1));
  splitB->SetAllStats(FULLUP);
  splitB->SetIndepIters();
  
  // Create new node to accumulate data in the loop
  TempVecReg* accum = new TempVecReg();
  accum->AddInput(vvdot->Input(2), vvdot->InputConnNum(2));

  node->m_poss->AddNode(accum);

  // Create tunnel for accumulator
  LoopTunnel* accTun = new LoopTunnel(POSSTUNIN);
  accTun->AddInput(accum, 0);
  accTun->SetAllStats(PARTUP);

  // Create loads for A and B
  LoadToRegs* loadA = new LoadToRegs();
  loadA->AddInput(splitA, 1);

  LoadToRegs* loadB = new LoadToRegs();
  loadB->AddInput(splitB, 1);

  // Create FMA operation that updates accum
  FMAdd* fmadd = new FMAdd();
  fmadd->AddInput(loadA, 0);
  fmadd->AddInput(loadB, 0);
  fmadd->AddInput(accTun, 0);

  // Create output tunnel for accumulate result
  LoopTunnel* accOut = new LoopTunnel(POSSTUNOUT);
  accOut->AddInput(fmadd, 0);
  accOut->AddInput(accTun, 1);
  accOut->CopyTunnelInfo(accTun);

  // Create combine outputs for A and B
  CombineSingleIter* comA = splitA->CreateMatchingCombine(0);
  CombineSingleIter* comB = splitB->CreateMatchingCombine(0);

  // Create the poss
  Poss* loopPoss = new Poss(3, comA, comB, accOut);
  RealLoop* loop;
  if (dataType == REAL_SINGLE) {
    loop = new RealLoop(LLDLALOOP, loopPoss, LLDLAMuSingle);
  } else if (dataType == REAL_DOUBLE) {
    loop = new RealLoop(LLDLALOOP, loopPoss, LLDLAMuDouble);
  } else {
    cout << "Error: Bad dataType in vadd apply\n";
    LOG_FAIL("replacement for throw call");
    throw;
  }

  // Add needed components to poss
  node->m_poss->AddPSet(loop);

  // Accumulate result of operation in C
  AccumReg* accumInC = new AccumReg();
  accumInC->AddInput(loop->OutTun(2), 0);
  accumInC->AddInput(vvdot->Input(2), vvdot->InputConnNum(2));

  node->m_poss->AddNode(accumInC);
  node->RedirectChildren(accumInC, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;
}

string VVDotToRegArith::GetType() const
{
  return "VVDot register arith " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer);
}

#endif // DOLLDLA
