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

#if DOLLDLA

#include "LLDLA.h"
#include "regArith.h"
#include "regLoadStore.h"
#include "vadd.h"

VAdd::VAdd(Layer layer, VecType vecType)
{
  m_layer = layer;
  m_vecType = vecType;
}

void VAdd::PrintCode(IndStream &out)
{
  out.Indent();
  if (m_layer == ABSLAYER) {
    if (GetDataType() == REAL_DOUBLE) {
      *out << "simple_add( ";
    } else if (GetDataType() == REAL_SINGLE) {
      *out << "simple_add_float( ";
    }
    if (GetVecType() == COLVECTOR) {
      *out << InputDataType(1).m_numRowsVar << ", " <<
	" 1, " <<
	GetInputName(0).str() << ", " <<
	InputDataType(0).m_rowStrideVar << ", " <<
	InputDataType(0).m_colStrideVar << ", " <<
	GetInputName(1).str() << ", " <<
	InputDataType(1).m_rowStrideVar << ", " <<
	InputDataType(1).m_colStrideVar << ");\n";
    } else {
      *out << " 1, " <<
	InputDataType(1).m_numColsVar << ", " <<
	GetInputName(0).str() << ", " <<
	InputDataType(0).m_rowStrideVar << ", " <<
	InputDataType(0).m_colStrideVar << ", " <<
	GetInputName(1).str() << ", " <<
	InputDataType(1).m_rowStrideVar << ", " <<
	InputDataType(1).m_colStrideVar << ");\n";
    }
    return;
  }

  if (m_layer != LLDLAPRIMITIVELAYER) {
    cout << "ERROR: Attempt to generate code from non-primitive vector add\n";
    LOG_FAIL("replacement for throw call");
  }
  const DataTypeInfo &inInfo = InputDataType(1);
  const Stride rowStride = inInfo.m_rowStride;
  const Stride colStride = inInfo.m_colStride;

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

void VAdd::PrintRowStride(IndStream &out)
{
  if (GetVecType() == COLVECTOR) {
    *out << "row_stride_add_2x1( " <<
      GetInputName(0).str() << ", " <<
      InputDataType(0).m_rowStrideVar << ", " <<
      GetInputName(1).str() << ", " <<
      InputDataType(0).m_rowStrideVar << " );\n";
  } else {
    *out << "row_stride_add_1x2( " <<
      GetInputName(0).str() << ", " <<
      InputDataType(0).m_rowStrideVar << ", " <<
      GetInputName(1).str() << ", " <<
      InputDataType(0).m_rowStrideVar << " );\n";
  }
}

void VAdd::PrintColStride(IndStream &out)
{
  if (GetVecType() == COLVECTOR) {
    *out << "col_stride_add_2x1( " <<
      GetInputName(0).str() << ", " <<
      InputDataType(0).m_colStrideVar << ", " <<
      GetInputName(1).str() << ", " <<
      InputDataType(0).m_colStrideVar << " );\n";
  } else {
    *out << "col_stride_add_1x2( " <<
      GetInputName(0).str() << ", " <<
      InputDataType(0).m_colStrideVar << ", " <<
      GetInputName(1).str() << ", " <<
      InputDataType(0).m_colStrideVar << " );\n";
  }

}

void VAdd::PrintGeneralStride(IndStream &out)
{
  if (GetVecType() == COLVECTOR) {
    *out << "gen_stride_add_2x1( " <<
      GetInputName(0).str() << ", " <<
      GetInputName(1).str() << ");\n";
  } else {
    *out << "gen_stride_add_1x2( " <<
      GetInputName(0).str() << ", " <<
      GetInputName(1).str() << ");\n";
  }
}

void VAdd::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2, 1>::Prop();
    
    VectorOpInputDimensionCheck(0);
    VectorOpInputDimensionCheck(1);

    if ((*GetInputM(0) != *GetInputM(1)) || (*GetInputN(0) != *GetInputN(1))) {
      cout << "ERROR: Cannot VAdd two vectors of different dimension\n";
      LOG_FAIL("replacement for throw call");
    }

    if (m_layer == ABSLAYER) {
      if (GetVecType() == ROWVECTOR) {
	m_cost = GetInputN(1)->Sum();
      } else {
	m_cost = GetInputM(1)->Sum();
      }
    } else {
      m_cost = ZERO;
    }

  }
}

void VAdd::VectorOpInputDimensionCheck(ConnNum inputNum)
{
  if (GetVecType() == ROWVECTOR && *GetInputM(inputNum) != 1) {
    cout << "ERROR: " << GetType() << " input # " << inputNum << " has more than 1 row\n";
    LOG_FAIL("replacement for throw call");
    throw;
  } else if (GetVecType() == COLVECTOR && *GetInputN(inputNum) != 1) {
    cout << "Input num is " << std::to_string((long long int) inputNum) << endl;
    GetInputN(inputNum)->Print();
    this->m_poss->PrintTransVecUp();
    cout << "ERROR: " << GetType() << " input # " << inputNum  << " has more than 1 column\n";
    LOG_FAIL("replacement for throw call");
    throw;
  }
  int regWidth = GetVecRegWidth();
  if (m_layer == LLDLAPRIMITIVELAYER) {
    if (GetVecType() == ROWVECTOR && *GetInputN(inputNum) != regWidth) {
      cout << "ERROR: " << GetType() << " input # " << inputNum << " does not have m_regWidth columns\n";
      LOG_FAIL("replacement for throw call");
      throw;
    } else if(GetVecType() == COLVECTOR && *GetInputM(inputNum) != regWidth) {
      cout << "ERROR: " << GetType() << " input # " << inputNum << " does not have m_regWidth rows\n";
      LOG_FAIL("replacement for throw call");
      throw;
    }
  }
}

Node* VAdd::BlankInst()
{
  return new VAdd(LLDLAPRIMITIVELAYER, COLVECTOR);
}

Node* VAdd::GetNewInst()
{
  return new VAdd(m_layer, GetVecType());
}

void VAdd::Duplicate(const Node *orig, bool shallow, bool possMerging) {
  DLAOp<2, 1>::Duplicate(orig, shallow, possMerging);
  const VAdd* vadd = static_cast<const VAdd*>(orig);
  m_layer = vadd->m_layer;
  m_vecType = vadd->GetVecType();
}

string VecTypeToString(VecType vType) {
  if (vType == COLVECTOR) {
    return "COLVECTOR";
  } else {
    return "ROWVECTOR";
  }
}

NodeType VAdd::GetType() const
{
  return "VADD" + LayerNumToStr(GetLayer()) + VecTypeToString(GetVecType());
}

VecType VAdd::GetVecType() const {
  /*  if (IsInputRowVector(0) && IsInputRowVector(1)) {
    return ROWVECTOR;
  }
  if (IsInputColVector(0) && IsInputColVector(1)) {
    return COLVECTOR;
    }*/
  ///  LOG_FAIL("replacement for throw call");
  return m_vecType;
}

Phase VAdd::MaxPhase() const
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

VAddLoopRef::VAddLoopRef(Layer fromLayer, Layer toLayer, VecType vtype)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_vtype = vtype;
}

string VAddLoopRef::GetType() const
{
  switch(m_vtype)
    {
    case(ROWVECTOR):
      return "VAddLoopRef - row vector";
    case(COLVECTOR):
      return "VAddLoopRef - column vector";
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }
}

bool VAddLoopRef::CheckRowVectorDimension(const VAdd* vadd) const {
  return *(vadd->GetInputN(0)) > vadd->GetVecRegWidth() &&
    vadd->GetInputN(0)->EvenlyDivisibleBy(vadd->GetVecRegWidth());
}

bool VAddLoopRef::CheckColVectorDimension(const VAdd* vadd) const {
  return *(vadd->GetInputM(0)) > vadd->GetVecRegWidth() &&
    vadd->GetInputM(0)->EvenlyDivisibleBy(vadd->GetVecRegWidth());
}

bool VAddLoopRef::CanApply(const Node *node) const {
  if (node->GetNodeClass() == VAdd::GetClass()) {
    const VAdd *vadd = static_cast<const VAdd*>(node);
    if (vadd->GetLayer() != m_fromLayer) {
      return false;
    }
    if (vadd->GetVecType() != m_vtype) {
      return false;
    }
    if (m_vtype == ROWVECTOR) {
      return CheckRowVectorDimension(vadd);
    } else if (m_vtype == COLVECTOR) {
      return CheckColVectorDimension(vadd);
    } else {
      LOG_FAIL("bad vector type in VAddLoopRef::CanApply");
      throw;
    }
    return false;
  }
  LOG_FAIL("bad node in VAddLoopRef::CanApply");
  throw;
}

void VAddLoopRef::Apply(Node *node) const
{
  VAdd *vadd = static_cast<VAdd*>(node);

  SplitSingleIter *split1 = new SplitSingleIter(m_vtype == COLVECTOR ? PARTDOWN : PARTRIGHT, POSSTUNIN, true);
  split1->AddInput(vadd->Input(1), vadd->InputConnNum(1));

  SplitSingleIter *split0 = new SplitSingleIter(m_vtype == COLVECTOR ? PARTDOWN : PARTRIGHT, POSSTUNIN, false);
  split0->AddInput(vadd->Input(0), vadd->InputConnNum(0));

  split0->SetAllStats(FULLUP);
  if (m_vtype == COLVECTOR) {
    split1->SetUpStats(FULLUP, FULLUP,
		      NOTUP, NOTUP);
  } else {
    split1->SetUpStats(FULLUP, NOTUP,
		       FULLUP, NOTUP);
  }

  split0->SetIndepIters();
  split1->SetIndepIters();

  VAdd *newVAdd = new VAdd(vadd->m_layer, m_vtype);
  newVAdd->SetLayer(m_toLayer);

  newVAdd->AddInput(split0, 1);
  newVAdd->AddInput(split1, 1);

  CombineSingleIter *com0 = split0->CreateMatchingCombine(0);
  CombineSingleIter *com1 = split1->CreateMatchingCombine(1, 1, newVAdd, 0);

  Poss *loopPoss = new Poss(2, com0, com1);
  RealLoop* loop;
  if (vadd->GetDataType() == REAL_SINGLE) {
    loop = new RealLoop(LLDLALOOP, loopPoss, LLDLAMuSingle);
  } else if (vadd->GetDataType() == REAL_DOUBLE) {
    loop = new RealLoop(LLDLALOOP, loopPoss, LLDLAMuDouble);
  } else {
    LOG_FAIL("Error: Bad GetDataType in vadd apply\n");
    throw;
  }
  loop->SetDimName(m_vtype == COLVECTOR ? DIMM : DIMN);
  
  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(1), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

VAddLowerLayer::VAddLowerLayer(Layer fromLayer, Layer toLayer, Size bs)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_bs = bs;
}

bool VAddLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == VAdd::GetClass()) {
    const VAdd *vadd = static_cast<const VAdd*>(node);
    if (vadd->GetLayer() != m_fromLayer) {
      return false;
    }

    if (vadd->GetVecType() == ROWVECTOR) {
      if (*vadd->GetInputN(0) <= m_bs &&
	  *vadd->GetInputN(1) <= m_bs) {
	return true;
      } else {
	return false;
      }
    } else if (vadd->GetVecType() == COLVECTOR) {
      if (*vadd->GetInputM(0) <= m_bs &&
	  *vadd->GetInputM(1) <= m_bs) {
	return true;
      } else {
	return false;
      }
    }
  }
  cout << "Error: Applying VAddLowerLayer to non-VAdd node\n";
  LOG_FAIL("replacement for throw call");
  throw;
}

void VAddLowerLayer::Apply(Node *node) const
{
  VAdd *vadd = static_cast<VAdd*>(node);
  vadd->SetLayer(m_toLayer);
  return;
}

string VAddLowerLayer::GetType() const
{
  return "VAdd lower layer " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer);
}

VAddToRegArith::VAddToRegArith(Layer fromLayer, Layer toLayer)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
}

bool VAddToRegArith::CanApply(const Node* node) const
{
  if (node->GetNodeClass() == VAdd::GetClass()) {
    const VAdd* vadd = static_cast<const VAdd*>(node);
    if (*(vadd->GetInputM(0)) > 1 &&
	vadd->GetInputM(0)->EvenlyDivisibleBy(vadd->GetVecRegWidth())) {
      return true;
    } else if (*(vadd->GetInputN(0)) > 1 &&
	       vadd->GetInputN(0)->EvenlyDivisibleBy(vadd->GetVecRegWidth())) {
      return true;
    }
    return false;
  }
  cout << "ERROR: Applying VAddToRegArith to non-VAdd node\n";
  LOG_FAIL("replacement for throw call");
  throw;
}

void VAddToRegArith::Apply(Node* node) const
{
  VAdd* vadd = static_cast<VAdd*>(node);

  bool splitDown = *vadd->GetInputN(0) == ONE;

  SplitSingleIter* splitX;
  SplitSingleIter* splitY;

  if (splitDown) {
    splitX = new SplitSingleIter(PARTDOWN, POSSTUNIN, true);
    splitY = new SplitSingleIter(PARTDOWN, POSSTUNIN, false);
  } else {
    splitX = new SplitSingleIter(PARTRIGHT, POSSTUNIN, true);
    splitY = new SplitSingleIter(PARTRIGHT, POSSTUNIN, false);
  }

  splitX->AddInput(vadd->Input(0), vadd->InputConnNum(0));
  splitX->SetAllStats(FULLUP);
  splitX->SetIndepIters();

  splitY->AddInput(vadd->Input(1), vadd->InputConnNum(1));
  if (splitDown) {
    splitY->SetUpStats(FULLUP, FULLUP,
		       NOTUP, NOTUP);
  } else {
    splitY->SetUpStats(FULLUP, NOTUP,
		       FULLUP, NOTUP);
  }
  splitY->SetIndepIters();

  LoadToRegs* loadX = new LoadToRegs();
  loadX->AddInput(splitX, 1);

  LoadToRegs* loadY = new LoadToRegs();
  loadY->AddInput(splitY, 1);

  Add* add = new Add();
  add->AddInput(loadX, 0);
  add->AddInput(loadY, 0);

  StoreFromRegs* storeToY = new StoreFromRegs();
  storeToY->AddInput(add, 0);
  storeToY->AddInput(splitY, 1);

  CombineSingleIter* comX = splitX->CreateMatchingCombine(0);
  CombineSingleIter* comY = splitY->CreateMatchingCombine(1, 1, storeToY, 0);
  
  Poss* loopPoss = new Poss(2, comX, comY);
  RealLoop* loop;
  if (vadd->GetDataType() == REAL_SINGLE) {
    loop = new RealLoop(LLDLALOOP, loopPoss, LLDLAMuSingle);
  } else if (vadd->GetDataType() == REAL_DOUBLE) {
    loop = new RealLoop(LLDLALOOP, loopPoss, LLDLAMuDouble);
  } else {
    cout << "Error: Bad GetDataType() in vadd apply\n";
    LOG_FAIL("replacement for throw call");
    throw;
  }
  loop->SetDimName(splitDown ? DIMM : DIMN);

  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(1), 0);
  node->m_poss->DeleteChildAndCleanUp(node);

  return;
}

string VAddToRegArith::GetType() const
{
  return "VAddToRegArith " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_fromLayer);
}

#endif // DOLLDLA
