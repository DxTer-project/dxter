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
#include "vmmul.h"
#include "LLDLA.h"
#include "loopSupport.h"
#include "regArith.h"
#include "regLoadStore.h"

#if DOLLDLA

VMMul::VMMul(Layer layer, Type type)
{
  m_type = type;
  m_layer = layer;
  return;
}

void VMMul::PrintCode(IndStream &out) {

  const DataTypeInfo &inInfo = InputDataType(1);
  const Stride rowStride = inInfo.m_rowStride;
  const Stride colStride = inInfo.m_colStride;

  out.Indent();

  if (m_layer == ABSLAYER) {
    *out << "simple_mmul( " <<
      "1, " <<
      InputDataType(1).m_numRowsVar << ", " <<
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

  if (rowStride == NONUNITSTRIDE && colStride == NONUNITSTRIDE) {
    PrintGeneralStride(out);
  } else if (rowStride == UNITSTRIDE && colStride == NONUNITSTRIDE) {
    PrintColStride(out);
  } else if (rowStride == NONUNITSTRIDE && colStride == UNITSTRIDE) {
    PrintRowStride(out);
  } else {
    *out << "ERROR: BAD STRIDE\n";
  }
  return;
}

void VMMul::PrintRowStride(IndStream &out)
{
  *out << "row_stride_mmul_1x2_2x2( " <<
    GetInputName(0).str() << ", " <<
    InputDataType(0).m_rowStrideVar << ", " <<
    GetInputName(1).str() << ", " <<
    InputDataType(1).m_rowStrideVar << ", " <<
    GetInputName(2).str() << ", " <<
    InputDataType(2).m_rowStrideVar << ");\n";
  return;
}

void VMMul::PrintColStride(IndStream &out)
{
  *out << "col_stride_mmul_1x2_2x2( " <<
    GetInputName(0).str() << ", " <<
    InputDataType(0).m_colStrideVar << ", " <<
    GetInputName(1).str() << ", " <<
    InputDataType(1).m_colStrideVar << ", " <<
    GetInputName(2).str() << ", " <<
    InputDataType(2).m_colStrideVar << ");\n";
}

void VMMul::PrintGeneralStride(IndStream &out)
{
  *out << "gen_stride_mmul_1x2_2x2( " <<
    GetInputName(0).str() << ", " <<
    GetInputName(1).str() << ", " <<
    GetInputName(2).str() << ");\n";
}

Node* VMMul::BlankInst()
{
  return new VMMul(LLDLAPRIMITIVELAYER, REAL);
}

Phase VMMul::MaxPhase() const
{
  switch(m_layer)
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

void VMMul::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3, 1>::Prop();

    if (*GetInputN(0) != *GetInputM(1)) {
      cout << "ERROR: Input dimensions for vmmul don't match\n";
      throw;
    } else if (*GetInputM(0) != *GetInputM(2)) {
      cout << "ERROR: Input dimensions don't match\n";
      throw;
    } else if (*GetInputM(2) != 1) {
      cout << "Error: Input vectors have more than 1 row\n";
      throw;
    }
    
    if (m_layer == LLDLAPRIMITIVELAYER) {
      if (*GetInputM(0) != 1 || *GetInputN(0) != LLDLA_MU) {
	cout << "ERROR: Primitive vector for vmmul must be 1 x 2\n";
	throw;
      }
      if (*GetInputN(1) != LLDLA_MU) {
	cout << "ERROR: Primitive matrix must be 2 x 2\n";
	throw;
      }
      if (*GetInputM(2) != 1) {
	cout << "ERROR: Primitive vector must be 1 x 2\n";
	throw;
      }
    }
    m_cost = ZERO;

  }
}

NodeType VMMul::GetType() const
{
  return "VMMul " + LayerNumToStr(GetLayer());
}

void VMMul::Duplicate(const Node* orig, bool shallow, bool possMerging)
{
  DLAOp<3, 1>::Duplicate(orig, shallow, possMerging);
  const VMMul* rhs = (VMMul*) orig;
  m_type = rhs->m_type;
  m_layer = rhs->m_layer;
  return;
}

string VMMulLoopRef::GetType() const
{
  switch (m_dim) {
  case(DIMN):
    return "VMMulLoopRef N dim";
  case(DIMM):
    cout << "Error: DIMN is not valid dimension for  VMMulLoopRef\n";
    throw;
  case(DIMK):
    return "VMMulLoopRef K dim";
  case(BADDIM):
    cout << "Error: VMMulLoopRef has BADDIM\n";
    throw;
  }
  cout << "No matching dimension case for VMMulLoopRef\n";
  throw;
}

bool VMMulLoopRef::CanApply(const Node* node) const
{
  if (node->GetNodeClass() == VMMul::GetClass()) {
    const VMMul* vmmul = (VMMul*) node;
    if (vmmul->GetLayer() != m_fromLayer) {
      return false;
    }
    if (m_dim == DIMN) {
    return !(*(vmmul->GetInputN(1)) <= m_bs.GetSize());
    } else if (m_dim == DIMK) {
      return !((*vmmul->GetInputM(1)) <= m_bs.GetSize());
    } else {
      cout << "Error: No matching VMMulLoopRef for given dimension\n";
      throw;
    }
  }
  return false;
}

void VMMulLoopRef::Apply(Node* node) const
{
  if (m_dim == DIMN) {
    ApplyDimN(node);
  } else {
    ApplyDimK(node);
  }
  return;
}

void VMMulLoopRef::ApplyDimK(Node* node) const
{
  VMMul* vmmul = (VMMul*) node;

  SplitSingleIter* splitA = new SplitSingleIter(PARTDOWN, POSSTUNIN, true);
  splitA->AddInput(vmmul->Input(1), vmmul->InputConnNum(1));
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();

  SplitSingleIter* splitX = new SplitSingleIter(PARTRIGHT, POSSTUNIN, false);
  splitX->AddInput(vmmul->Input(0), vmmul->InputConnNum(0));
  splitX->SetAllStats(FULLUP);
  splitX->SetIndepIters();

  LoopTunnel* inY = new LoopTunnel(POSSTUNIN);
  inY->AddInput(vmmul->Input(2), vmmul->InputConnNum(2));
  inY->SetAllStats(PARTUP);
  inY->SetIndepIters();

  VMMul* newVMMul = new VMMul(m_toLayer, vmmul->m_type);
  newVMMul->AddInput(splitX, 1);
  newVMMul->AddInput(splitA, 1);
  newVMMul->AddInput(inY, 0);

  CombineSingleIter* comA = splitA->CreateMatchingCombine(0);
  CombineSingleIter* comX = splitX->CreateMatchingCombine(0);

  LoopTunnel* outY = new LoopTunnel(POSSTUNOUT);
  outY->AddInput(newVMMul, 0);
  outY->AddInput(inY, 1);
  outY->CopyTunnelInfo(inY);

  Poss* loopPoss = new Poss(3, comX, comA, outY);
  RealLoop* loop = new RealLoop(LLDLALOOP, loopPoss, m_bs);

  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(2), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;  
}

void VMMulLoopRef::ApplyDimN(Node* node) const
{
  VMMul* vmmul = (VMMul*) node;

  SplitSingleIter* splitA = new SplitSingleIter(PARTRIGHT, POSSTUNIN, true);
  splitA->AddInput(vmmul->Input(1), vmmul->InputConnNum(1));
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();

  LoopTunnel* tunX = new LoopTunnel(POSSTUNIN);
  tunX->AddInput(vmmul->Input(0), vmmul->InputConnNum(0));
  tunX->SetAllStats(FULLUP);
  tunX->SetIndepIters();

  SplitSingleIter* splitY = new SplitSingleIter(PARTRIGHT, POSSTUNIN, false);
  splitY->AddInput(vmmul->Input(2), vmmul->InputConnNum(2));
  splitY->SetUpStats(FULLUP, NOTUP,
		     FULLUP, NOTUP);
  splitY->SetIndepIters();

  VMMul* newMul = new VMMul(m_toLayer, vmmul->m_type);

  // Attach to arguments
  newMul->AddInput(tunX, 0);
  newMul->AddInput(splitA, 1);
  newMul->AddInput(splitY, 1);

  CombineSingleIter* comA = splitA->CreateMatchingCombine(0);
  CombineSingleIter* comY = splitY->CreateMatchingCombine(1, 1, newMul, 0);

  LoopTunnel* outX = new LoopTunnel(POSSTUNOUT);
  outX->AddInput(tunX, 0);
  outX->AddInput(tunX, 1);
  outX->CopyTunnelInfo(tunX);

  Poss* loopPoss = new Poss(3, comA, outX, comY);
  RealLoop* loop = new RealLoop(LLDLALOOP, loopPoss, m_bs);

  loop->SetDimName(DIMN);

  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(2), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;
}

bool VMMulLowerLayer::CanApply(const Node* node) const
{
  if (node->GetNodeClass() == VMMul::GetClass()) {
    const VMMul* vmmul = (VMMul*) node;
    if (vmmul->GetLayer() != m_fromLayer) {
      return false;
    }
    if (m_toLayer == LLDLAPRIMITIVELAYER) {
      if (*(vmmul->GetInputM(1)) == m_bs &&
	  *(vmmul->GetInputN(1)) == m_bs) {
	return true;
      } else {
	return false;
      }
    } else {
      if (*(vmmul->GetInputM(1)) <= m_bs &&
	  *(vmmul->GetInputN(1)) <= m_bs) {
	return true;
      } else {
	return false;
      }
    }
  } else {
    throw;
  }
}

void VMMulLowerLayer::Apply(Node* node) const
{
  VMMul* vmmul = (VMMul*) node;
  vmmul->SetLayer(m_toLayer);
  return;
}

string VMMulLowerLayer::GetType() const
{
  return "VMMul lower layer from " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer);
}

string VMMulToRegArith::GetType() const
{
  return "VMMulToRegArith from " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer);
}

bool VMMulToRegArith::CanApply(const Node* node) const
{
  if (node->GetNodeClass() == VMMul::GetClass()) {
    VMMul* vmmul = (VMMul*) node;
    return *vmmul->GetInputN(1) == LLDLA_MU;
  }
  return false;
}

void VMMulToRegArith::Apply(Node* node) const
{
  VMMul* vmmul = (VMMul*) node;

  SplitSingleIter* splitA = new SplitSingleIter(PARTDOWN, POSSTUNIN, false);
  splitA->AddInput(vmmul->Input(1), vmmul->InputConnNum(1));
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();

  SplitSingleIter* splitX = new SplitSingleIter(PARTRIGHT, POSSTUNIN, true);
  splitX->AddInput(vmmul->Input(0), vmmul->InputConnNum(0));
  splitX->SetAllStats(FULLUP);
  splitX->SetIndepIters();

  LoadToRegs* loadY = new LoadToRegs();
  loadY->AddInput(vmmul->Input(2), vmmul->InputConnNum(2));

  node->m_poss->AddNode(loadY);

  LoopTunnel* yTun = new LoopTunnel(POSSTUNIN);
  yTun->AddInput(loadY, 0);
  yTun->SetAllStats(PARTUP);

  LoadToRegs* loadA = new LoadToRegs();
  loadA->AddInput(splitA, 1);

  DuplicateRegLoad* loadX = new DuplicateRegLoad();
  loadX->AddInput(splitX, 1);

  FMAdd* fmadd = new FMAdd();
  fmadd->AddInput(loadX, 0);
  fmadd->AddInput(loadA, 0);
  fmadd->AddInput(yTun, 0);

  LoopTunnel* fmaOut = new LoopTunnel(POSSTUNOUT);
  fmaOut->AddInput(fmadd, 0);
  fmaOut->AddInput(yTun, 1);
  fmaOut->CopyTunnelInfo(yTun);

  CombineSingleIter* comA = splitA->CreateMatchingCombine(0);
  CombineSingleIter* comX = splitX->CreateMatchingCombine(0);

  Poss* loopPoss = new Poss(3, comA, comX, fmaOut);
  RealLoop* loop = new RealLoop(LLDLALOOP, loopPoss, UnitBS);

  loop->SetDimName(DIMM);

  node->m_poss->AddPSet(loop);

  StoreFromRegs* storeToY = new StoreFromRegs();
  storeToY->AddInput(loop->OutTun(2), 0);
  storeToY->AddInput(vmmul->Input(2), vmmul->InputConnNum(2));

  node->m_poss->AddNode(storeToY);
  
  node->RedirectChildren(storeToY);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;
}

#endif // DOLLDLA
