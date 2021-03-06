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

#include "DLAOp.h"
#include "mvmul.h"
#include "svmul.h"
#include "loopSupport.h"
#include "regLoadStore.h"
#include "regArith.h"
#include "LLDLA.h"

#if DOLLDLA

MVMul::MVMul(Layer layer)
{
  m_layer = layer;
  return;
}

void MVMul::PrintCode(IndStream &out)
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
    *out << InputDataType(0).m_numRowsVar << ", " <<
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
    cout << "ERROR: Attempt to generate code from non-primitive matrix vector multiply\n";
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


void MVMul::PrintRowStride(IndStream &out)
{
  *out << "row_stride_mmul_2x2_2x1( " <<
    GetInputName(0).str() << ", " <<
    InputDataType(0).m_rowStrideVar << ", " <<
    GetInputName(1).str() << ", " <<
    InputDataType(1).m_rowStrideVar << ", " <<
    GetInputName(2).str() << ", " <<
    InputDataType(2).m_rowStrideVar << ");\n";
  return;
}

void MVMul::PrintColStride(IndStream &out)
{
  *out << "col_stride_mmul_2x2_2x1( " <<
    GetInputName(0).str() << ", " <<
    InputDataType(0).m_colStrideVar << ", " <<
    GetInputName(1).str() << ", " <<
    InputDataType(1).m_colStrideVar << ", " <<
    GetInputName(2).str() << ", " <<
    InputDataType(2).m_colStrideVar << ");\n";
}

void MVMul::PrintGeneralStride(IndStream &out)
{
  *out << "gen_stride_mmul_2x2_2x1( " <<
    GetInputName(0).str() << ", " <<
    GetInputName(1).str() << ", " <<
    GetInputName(2).str() << " );\n";
}

Node* MVMul::BlankInst()
{
  return new MVMul(ABSLAYER);
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
      LOG_FAIL("replacement for throw call");
      throw;
    }
}

void MVMul::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3, 1>::Prop();

    if (!InputDimsAreOneRepeatedSizeEach(0)
	|| !InputDimsAreOneRepeatedSizeEach(1)
	|| !InputDimsAreOneRepeatedSizeEach(2)) {
      cout << "ERROR: MVMul input dimensions are not single, repeated sizes" << endl;
      LOG_FAIL("replacement for throw call");
    }

    if (*GetInputN(0) != *GetInputM(1)) {
      cout << "ERROR: Input dimensions don't match\n";
      LOG_FAIL("replacement for throw call");
    } else if (*GetInputM(0) != *GetInputM(2)) {
      cout << "ERROR: Input dimensions don't match\n";
      LOG_FAIL("replacement for throw call");
    } else if (!IsInputColVector(1) || !IsInputColVector(2)) {
      cout << "Error: MVMul inputs that should be column vectors are not\n";
      LOG_FAIL("replacement for throw call");
    }
    
    if (m_layer == LLDLAPRIMITIVELAYER) {
      if (*GetInputM(0) != GetVecRegWidth() || *GetInputN(0) != GetVecRegWidth()) {
	cout << "ERROR: Primitive matrix must be 2 x 2 for MVMul\n";
	LOG_FAIL("replacement for throw call");
      }
      if (*GetInputM(1) != GetVecRegWidth()) {
	cout << "ERROR: Primitive vector must be 2 x 1 MVMul\n";
	LOG_FAIL("replacement for throw call");
      }
      if (*GetInputM(2) != GetVecRegWidth()) {
	cout << "ERROR: Primitive vector must be 2 x 1 MVMul\n";
	LOG_FAIL("replacement for throw call");
      }
    }

    if (m_layer == ABSLAYER) {
      m_cost = TWO * GetInputM(0)->SumProds11(*GetInputN(0));
    } else {
      m_cost = ZERO;
    }
  }
}

NodeType MVMul::GetType() const
{
  return "MVMul" +  LayerNumToStr(GetLayer());
}

void MVMul::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<3,1>::Duplicate(orig, shallow, possMerging);
  const MVMul *rhs = static_cast<const MVMul*>(orig);
  m_layer = rhs->m_layer;
  return;
}

MVMulLowerLayer::MVMulLowerLayer(Layer fromLayer, Layer toLayer, Size bs)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_bs = bs;
}

bool MVMulLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == MVMul::GetClass()) {
    const MVMul *mvmul = static_cast<const MVMul*>(node);
    if (mvmul->GetLayer() != m_fromLayer) {
      return false;
    }

    int regWidth = mvmul->GetVecRegWidth();
    if (m_toLayer == LLDLAPRIMITIVELAYER) {
      if (*mvmul->GetInputM(0) != regWidth ||
	  *mvmul->GetInputN(0) != regWidth ||
	  *mvmul->GetInputM(1) != regWidth ||
	  *mvmul->GetInputN(1) != 1        ||
	  *mvmul->GetInputM(2) != regWidth ||
	  *mvmul->GetInputN(2) != 1) {
	return false;
      } else {
	return true;
      }
    } else if (*mvmul->GetInputM(0) <= m_bs &&
	       *mvmul->GetInputN(0) <= m_bs) {
      return true;
    } else {
      return false;
    }
  }
  cout << "Error: Applying MVMulLowerLayer to non-MVMul node\n";
  LOG_FAIL("replacement for throw call");
  throw;
}

void MVMulLowerLayer::Apply(Node *node) const
{
  MVMul *svmul = static_cast<MVMul*>(node);
  svmul->SetLayer(m_toLayer);
}

string MVMulLowerLayer::GetType() const
{
  return "lower layer " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer);
}

MVMulLoopRef::MVMulLoopRef(Layer fromLayer, Layer toLayer, DimName dim, BSSize bs)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_dim = dim;
  m_bs = bs;
}

string MVMulLoopRef::GetType() const
{
  switch (m_dim)
    {
    case (DIMM):
      return "MVMulLoopRef - dim m " + std::to_string((long long int) m_bs.GetSize());
    case (DIMN):
      return "MVMulLoopRef - dim n " + std::to_string((long long int) m_bs.GetSize());
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }  
}

bool MVMulLoopRef::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == MVMul::GetClass()) {
    const MVMul *mvmul = static_cast<const MVMul*>(node);

    if ((mvmul->Input(0)->GetVecRegWidth() > ((int) m_bs.GetSize()))) {
      return false;
    }
    if (mvmul->GetLayer() != m_fromLayer) {
      return false;
    }

    if (m_dim == DIMM) {
      return mvmul->InputMIsMultipleOf(0, m_bs.GetSize())
	//	&& *(mvmul->GetInputN(0)) > m_bs.GetSize();
      	&& *(mvmul->GetInputM(0)) > m_bs.GetSize();
    } else {
      return mvmul->InputNIsMultipleOf(0, m_bs.GetSize())
      	&& *(mvmul->GetInputN(0)) > m_bs.GetSize();
	//	&& mvmul->InputMIsMultipleOf(0, m_bs.GetSize());
    }
  }
  LOG_FAIL("replacement for throw call");
  throw;
}

void MVMulLoopRef::Apply(Node *node) const {
  if (m_dim == DIMM) {
    ApplyRowSplit(node);
  } else if (m_dim == DIMN) {
    ApplyColSplit(node);
  } else {
    cout << "ERROR: Not refining mvmul loop around m or n dimension\n";
    LOG_FAIL("replacement for throw call");
  }
}

void MVMulLoopRef::ApplyRowSplit(Node *node) const {
  MVMul *mul = static_cast<MVMul*>(node);
  
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
  splitY->SetIndepIters();

  // Create tunnel for x
  LoopTunnel *xTun = new LoopTunnel(POSSTUNIN);
  xTun->AddInput(mul->Input(1), mul->InputConnNum(1));
  xTun->SetAllStats(FULLUP);

  // Create new mvmul for loop body
  MVMul *newMul = new MVMul(m_toLayer);
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
  Poss *loopPoss = new Poss(3, comA, xOut, comY);
  RealLoop *loop = new RealLoop(LLDLALOOP, loopPoss, m_bs);
  loop->SetDimName(m_dim);

  // Set up graph to correctly incorporate this loop
  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(2), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

void MVMulLoopRef::ApplyColSplit(Node *node) const
{
  MVMul *mul = static_cast<MVMul*>(node);

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

  // Create new svmul for loop body
  MVMul *newMul = new MVMul(m_toLayer);
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

MVMulToRegArith::MVMulToRegArith(Layer fromLayer, Layer toLayer)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
}

string MVMulToRegArith::GetType() const
{
  return "MVMulToRegArith from " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer);
}

bool MVMulToRegArith::CanApply(const Node* node) const
{
  if (node->GetNodeClass() == MVMul::GetClass()) {
    const MVMul* mvmul = static_cast<const MVMul*>(node);
    return *mvmul->GetInputM(0) == mvmul->GetVecRegWidth();
  }
  return false;
}

void MVMulToRegArith::Apply(Node* node) const
{
  MVMul* mvmul = static_cast<MVMul*>(node);
  
  // Split matrix A into mu x 1 column vectors
  SplitSingleIter* splitA = new SplitSingleIter(PARTRIGHT, POSSTUNIN, true);
  splitA->AddInput(mvmul->Input(0), mvmul->InputConnNum(0));
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();

  // Split x into individual elements
  SplitSingleIter* splitX = new SplitSingleIter(PARTDOWN, POSSTUNIN, false);
  splitX->AddInput(mvmul->Input(1), mvmul->InputConnNum(1));
  splitX->SetAllStats(FULLUP);
  splitX->SetIndepIters();

  // Create vector register with elements of y
  LoadToRegs* loadY = new LoadToRegs();
  loadY->AddInput(mvmul->Input(2), mvmul->InputConnNum(2));
  
  node->m_poss->AddNode(loadY);
  
  // Create tunnel for y register
  LoopTunnel* yTun = new LoopTunnel(POSSTUNIN);
  yTun->AddInput(loadY, 0);
  yTun->SetAllStats(PARTUP);
  
  // Create load for elements of A in loop body
  LoadToRegs* loadA = new LoadToRegs();
  loadA->AddInput(splitA, 1);

  // Create duplicate load for x
  DuplicateRegLoad* loadX = new DuplicateRegLoad();
  loadX->AddInput(splitX, 1);

  // Create new FMA instruction for loop body
  FMAdd* fmadd = new FMAdd();
  fmadd->AddInput(loadA, 0);
  fmadd->AddInput(loadX, 0);
  fmadd->AddInput(yTun, 0);

  // Create output tunnel for FMA result
  LoopTunnel* fmaOut = new LoopTunnel(POSSTUNOUT);
  fmaOut->AddInput(fmadd, 0);
  fmaOut->AddInput(yTun, 1);
  fmaOut->CopyTunnelInfo(yTun);

  // Create combines for A and x
  CombineSingleIter* comA = splitA->CreateMatchingCombine(0);
  CombineSingleIter* comX = splitX->CreateMatchingCombine(0);

  // Create the poss
  Poss* loopPoss = new Poss(3, comA, comX, fmaOut);
  RealLoop* loop = new RealLoop(LLDLALOOP, loopPoss, UnitBS);

  node->m_poss->AddPSet(loop);

  // Store result of computation back to y
  StoreFromRegs* storeToY = new StoreFromRegs();
  storeToY->AddInput(loop->OutTun(2), 0);
  storeToY->AddInput(mvmul->Input(2), mvmul->InputConnNum(2));

  node->m_poss->AddNode(storeToY);
  
  node->RedirectChildren(storeToY);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;
}

#endif // DOLLDLA
