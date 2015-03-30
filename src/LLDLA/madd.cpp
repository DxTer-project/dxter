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

#include "madd.h"
#include "vadd.h"
#include "regArith.h"
#include "regLoadStore.h"

#if DOLLDLA

#include "logging.h"

MAdd::MAdd(Layer layer)
{
  m_layer = layer;
}

void MAdd::PrintCode(IndStream &out)
{
  if (m_layer == ABSLAYER) {
    if (GetDataType() == REAL_DOUBLE) {
      out.Indent();
      *out << "simple_add( ";
    } else if (GetDataType() == REAL_SINGLE) {
      out.Indent();
      *out << "simple_add_float( ";
    }
    *out << InputDataType(0).m_numRowsVar << ", " <<
      InputDataType(0).m_numColsVar << ", " <<
      GetInputName(0).str() << ", " <<
      InputDataType(0).m_rowStrideVar << ", " <<
      InputDataType(0).m_colStrideVar << ", " <<
      GetInputName(1).str() << ", " <<
      InputDataType(1).m_rowStrideVar << ", " <<
      InputDataType(1).m_colStrideVar << " );\n";
    
    return;
  }
  if (m_layer != LLDLAPRIMITIVELAYER) {
    cout << "ERROR: Attempt to generate code from non-primitive matrix add\n";
    LOG_FAIL("Replacement for call to throw");
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


void MAdd::PrintRowStride(IndStream &out)
{
  *out << "row_stride_add_2x2( " <<
    GetInputName(0).str() << ", " <<
    InputDataType(0).m_rowStrideVar << ", " <<
    GetInputName(1).str() << ", " <<
    InputDataType(1).m_rowStrideVar << " );\n";
}

void MAdd::PrintColStride(IndStream &out)
{
  *out << "col_stride_add_2x2( " <<
    GetInputName(0).str() << ", " <<
    InputDataType(0).m_colStrideVar << ", " <<
    GetInputName(1).str() << ", " <<
    InputDataType(1).m_colStrideVar << " );\n";
}

void MAdd::PrintGeneralStride(IndStream &out)
{
  *out << "gen_stride_add_2x2( " <<
    GetInputName(0).str() << ", " <<
    InputDataType(0).m_rowStrideVar << ", " <<
    InputDataType(0).m_colStrideVar << ", " <<
    GetInputName(1).str() << ", " <<
    InputDataType(0).m_rowStrideVar << ", " <<
    InputDataType(0).m_colStrideVar << " );\n";
}

void MAdd::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2, 1>::Prop();

    if ((*GetInputM(0) != *GetInputM(1)) || (*GetInputN(0) != *GetInputN(1))) {
      cout << "ERROR: Cannot MAdd two matrices of different dimension\n";
      LOG_FAIL("Replacement for call to throw;");
    }

    if (m_layer == LLDLAPRIMITIVELAYER) {
      if ((*GetInputM(0) != GetVecRegWidth()) || (*GetInputN(0) != GetVecRegWidth())
	  || (*GetInputM(1) != GetVecRegWidth()) || (*GetInputN(1) != GetVecRegWidth())) {
	cout << "ERROR: MAdd of matrices that do not have GetVecRegWidth() dimensions in LLDLAPRIMITIVELAYER\n";
      }
    }

    if (m_layer == ABSLAYER) {
      m_cost = GetInputM(0)->SumProds11(*GetInputN(0));
    } else {
      m_cost = ZERO;
    }
  }
}

Node* MAdd::BlankInst()
{
  return new MAdd(LLDLAPRIMITIVELAYER);
}

NodeType MAdd::GetType() const
{
  return "MAdd" + LayerNumToStr(GetLayer());
}

Phase MAdd::MaxPhase() const
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
      LOG_FAIL("Replacement for call to throw;");
    }
}

MAddLoopRef::MAddLoopRef(Layer fromLayer, Layer toLayer, DimName dim, BSSize bs)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_dim = dim;
  m_bs = bs;
}

string MAddLoopRef::GetType() const
{
  switch (m_dim) {
  case (DIMM):
    return "MAddM" + std::to_string((long long int) m_bs.GetSize());
  case(DIMN):
    return "MAddN" + std::to_string((long long int) m_bs.GetSize());
  default:
    return "ERROR: Bad dimension in MAddLoopRef transform";
  }
  return "MAdd";
}

bool MAddLoopRef::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == MAdd::GetClass()) {
    const MAdd *madd = static_cast<const MAdd*>(node);
    if (madd->GetLayer() != m_fromLayer) {
      return false;
    }
    if (m_dim == DIMM) {
      return (*(madd->GetInputM(0)) > m_bs.GetSize()) &&
	madd->GetInputM(0)->EvenlyDivisibleBy(m_bs.GetSize());
    } else if (m_dim == DIMN) {
      return (*(madd->GetInputN(0)) > m_bs.GetSize()) &&
	madd->GetInputN(0)->EvenlyDivisibleBy(m_bs.GetSize());
    } else {
      return false;
    }
  }
  return false;
}

void MAddLoopRef::Apply(Node *node) const
{
  MAdd *madd = static_cast<MAdd*>(node);
  
  SplitSingleIter *split0 = new SplitSingleIter(m_dim == DIMM ? PARTDOWN : PARTRIGHT, POSSTUNIN, true);
  split0->AddInput(madd->Input(0), madd->InputConnNum(0));

  SplitSingleIter *split1 = new SplitSingleIter(m_dim == DIMM ? PARTDOWN : PARTRIGHT, POSSTUNIN, false);
  split1->AddInput(madd->Input(1), madd->InputConnNum(1));

  split0->SetAllStats(FULLUP);
  if (m_dim == DIMM) {
    split1->SetUpStats(FULLUP, FULLUP,
		       NOTUP, NOTUP);
  } else {
    split1->SetUpStats(FULLUP, NOTUP,
		       FULLUP, NOTUP);
  }

  split0->SetIndepIters();
  split1->SetIndepIters();

  MAdd *newMAdd = new MAdd(m_toLayer);
  newMAdd->AddInput(split0, 1);
  newMAdd->AddInput(split1, 1);

  CombineSingleIter *com0 = split0->CreateMatchingCombine(0);
  CombineSingleIter *com1 = split1->CreateMatchingCombine(1, 1, newMAdd, 0);

  Poss *loopPoss = new Poss(2, com0, com1);
  RealLoop *loop = new RealLoop(LLDLALOOP, loopPoss, m_bs);
  loop->SetDimName(m_dim);

  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(1), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;
}

MAddToVAddLoopRef::MAddToVAddLoopRef(Layer fromLayer, Layer toLayer, VecType vecType)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_vecType = vecType;
}

string MAddToVAddLoopRef::GetType() const
{
  switch (m_vecType) {
  case(COLVECTOR):
    return "MAddToVAddLoopRef from " + LayerNumToStr(m_fromLayer) + " to " +
      LayerNumToStr(m_toLayer) + " of type COLVECTOR";
  case(ROWVECTOR):
    return "MAddToVAddLoopRef from " + LayerNumToStr(m_fromLayer) + " to " +
      LayerNumToStr(m_toLayer) + " of type ROWVECTOR";    
  default:
    LOG_FAIL("Replacement for call to throw;");
  }
}

bool MAddToVAddLoopRef::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == MAdd::GetClass()) {
    const MAdd *madd = static_cast<const MAdd*>(node);
    if (madd->GetLayer() != m_fromLayer) {
      return false;
    }
    if (m_vecType == ROWVECTOR) {
      return *(madd->GetInputM(0)) > 1
	&& *(madd->GetInputN(0)) > madd->GetVecRegWidth();
    } else {
      return *(madd->GetInputM(0)) > madd->GetVecRegWidth()
	&& *(madd->GetInputN(0)) > 1;
    }
  }
  cout << "ERROR: Applying MAddToVAddLoopRef to non MAdd node\n";
  cout << "Node has class " << node->GetNodeClass() << endl;
  LOG_FAIL("Replacement for call to throw;");
}

void MAddToVAddLoopRef::Apply(Node *node) const
{
  MAdd *madd = static_cast<MAdd*>(node);

  SplitSingleIter *split0 = new SplitSingleIter(m_vecType == ROWVECTOR ? PARTDOWN : PARTRIGHT, POSSTUNIN, true);
  split0->AddInput(madd->Input(0), madd->InputConnNum(0));

  SplitSingleIter *split1 = new SplitSingleIter(m_vecType == ROWVECTOR ? PARTDOWN : PARTRIGHT, POSSTUNIN, false);
  split1->AddInput(madd->Input(1), madd->InputConnNum(1));

  split0->SetAllStats(FULLUP);
  if (m_vecType == ROWVECTOR) {
    split1->SetUpStats(FULLUP, FULLUP,
		       NOTUP, NOTUP);
  } else {
    split1->SetUpStats(FULLUP, NOTUP,
		       FULLUP, NOTUP);
  }

  split0->SetIndepIters();
  split1->SetIndepIters();

  VAdd *newVAdd = new VAdd(m_toLayer, m_vecType);
  newVAdd->AddInput(split0, 1);
  newVAdd->AddInput(split1, 1);

  CombineSingleIter *com0 = split0->CreateMatchingCombine(0);
  CombineSingleIter *com1 = split1->CreateMatchingCombine(1, 1, newVAdd, 0);
  
  Poss *loopPoss = new Poss(2, com0, com1);
  RealLoop *loop = new RealLoop(LLDLALOOP, loopPoss, UnitBS);

  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(1), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
  return;
}

MAddLowerLayer::MAddLowerLayer(Layer fromLayer, Layer toLayer, Size bs)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_bs = bs;
}

bool MAddLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == MAdd::GetClass()) {
    const MAdd *madd = static_cast<const MAdd*>(node);
    if (madd->GetLayer() != m_fromLayer) {
      return false;
    }
    if (*(madd->GetInputM(1)) <= m_bs &&
	*(madd->GetInputN(1)) <= m_bs &&
	*(madd->GetInputM(0)) <= m_bs &&
	*(madd->GetInputN(0)) <= m_bs) {
      return true;
    } else {
      return false;
    }
  }
  else {
    LOG_FAIL("Replacement for call to throw;");
  }
}

void MAddLowerLayer::Apply(Node *node) const
{
  MAdd *madd = static_cast<MAdd*>(node);
  madd->SetLayer(m_toLayer);
}

string MAddLowerLayer::GetType() const
{
  return "MAdd lower layer " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer);
}

string MAddToRegArith::GetType() const
{
  return "MAddToRegArith " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer);
}

MAddToRegArith::MAddToRegArith(Layer fromLayer, Layer toLayer)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
}

bool MAddToRegArith::CanApply(const Node* node) const
{
  if (node->GetNodeClass() == MAdd::GetClass()) {
    const MAdd* madd = static_cast<const MAdd*>(node);
    if (madd->GetLayer() != m_fromLayer) {
      return false;
    }
    if ((*(madd->GetInputM(0)) == madd->GetVecRegWidth()) ||
	(*(madd->GetInputN(0)) == madd->GetVecRegWidth())) {
      return true;
    } else {
      return false;
    }
  }
  return false;
}

void MAddToRegArith::Apply(Node* node) const
{
  MAdd* madd = static_cast<MAdd*>(node);

  // Set direction of split
  bool splitIntoRows = *(madd->GetInputM(0)) > madd->GetVecRegWidth();

  // Split matrices A and B
  SplitSingleIter* splitA;
  SplitSingleIter* splitB;
  if (splitIntoRows) {
    splitA = new SplitSingleIter(PARTDOWN, POSSTUNIN, true);
    splitB = new SplitSingleIter(PARTDOWN, POSSTUNIN, false);
  } else {
    splitA = new SplitSingleIter(PARTRIGHT, POSSTUNIN, true);
    splitB = new SplitSingleIter(PARTRIGHT, POSSTUNIN, false);
  }

  splitA->AddInput(madd->Input(0), madd->InputConnNum(0));
  splitA->SetAllStats(FULLUP);
  splitA->SetIndepIters();

  splitB->AddInput(madd->Input(1), madd->InputConnNum(1));
  if (splitIntoRows) {
    splitB->SetUpStats(FULLUP, FULLUP,
		       NOTUP, NOTUP);
  } else {
    splitB->SetUpStats(FULLUP, NOTUP,
		       FULLUP, NOTUP);
  }
  splitB->SetIndepIters();

  // Create loads for A and B
  LoadToRegs* loadA = new LoadToRegs();
  loadA->AddInput(splitA, 1);

  LoadToRegs* loadB = new LoadToRegs();
  loadB->AddInput(splitB, 1);

  // Create new add node
  Add* add = new Add();
  add->AddInput(loadA, 0);
  add->AddInput(loadB, 0);

  // Create store to write back to B
  StoreFromRegs* storeToB = new StoreFromRegs();
  storeToB->AddInput(add, 0);
  storeToB->AddInput(splitB, 1);

  // Create combines for A and B
  CombineSingleIter* comA = splitA->CreateMatchingCombine(0);
  CombineSingleIter* comB = splitB->CreateMatchingCombine(1, 1, storeToB, 0);

  Poss* loopPoss = new Poss(2, comA, comB);
  RealLoop* loop = new RealLoop(LLDLALOOP, loopPoss, UnitBS);
  loop->SetDimName(splitIntoRows ? DIMM : DIMN);

  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(1), 0);
  node->m_poss->DeleteChildAndCleanUp(node);

  return;
}

#endif // DOLLDLA
