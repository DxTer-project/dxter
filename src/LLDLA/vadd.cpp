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

#include "LLDLA.h"

#if DOLLDLA

#include "LLDLA.h"
#include "regArith.h"
#include "regLoadStore.h"
#include "vadd.h"

VAdd::VAdd(VecType vecType, Layer layer, Type type)
{
  m_vecType = vecType;
  m_layer = layer;
  m_type = type;
  m_regWidth = arch->VecRegWidth(m_type);
}

void VAdd::PrintCode(IndStream &out)
{
  out.Indent();
  if (m_layer == ABSLAYER) {
    if (m_type == REAL_DOUBLE) {
      *out << "simple_add( ";
    } else if (m_type == REAL_SINGLE) {
      *out << "simple_add_float( ";
    }
    if (m_vecType == COLVECTOR) {
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
    throw;
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
  if (m_vecType == COLVECTOR) {
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
  if (m_vecType == COLVECTOR) {
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
  if (m_vecType == COLVECTOR) {
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
      throw;
    }

    m_cost = ZERO;
  }
}

void VAdd::VectorOpInputDimensionCheck(ConnNum inputNum)
{
  if (m_vecType == ROWVECTOR && *GetInputM(inputNum) != 1) {
    cout << "ERROR: " << GetType() << " input # " << inputNum << " has more than 1 row\n";
    throw;
  } else if (m_vecType == COLVECTOR && *GetInputN(inputNum) != 1) {
    cout << "ERROR: " << GetType() << " input # " << inputNum  << " has more than 1 column\n";
    throw;
  }
  
  if (m_layer == LLDLAPRIMITIVELAYER) {
    if (m_vecType == ROWVECTOR && *GetInputN(inputNum) != m_regWidth) {
      cout << "ERROR: " << GetType() << " input # " << inputNum << " does not have m_regWidth columns\n";
      throw;
    } else if(m_vecType == COLVECTOR && *GetInputM(inputNum) != m_regWidth) {
      cout << "ERROR: " << GetType() << " input # " << inputNum << " does not have m_regWidth rows\n";
      throw;
    }
  }
}

Node* VAdd::BlankInst()
{
  return new VAdd(COLVECTOR, LLDLAPRIMITIVELAYER, REAL_SINGLE);
}

Node* VAdd::GetNewInst()
{
  return new VAdd(m_vecType, m_layer, m_type);
}

NodeType VAdd::GetType() const
{
  return "VADD" + LayerNumToStr(GetLayer()) + (char)m_type + (char)m_vecType;
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
      throw;
    }
}

VAddLoopRef::VAddLoopRef(Layer fromLayer, Layer toLayer, VecType vtype, BSSize bs, Type type)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_vtype = vtype;
  m_bs = bs;
  m_type = type;
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
      throw;
    }
}

bool VAddLoopRef::CanApply(const Node *node) const
{
  const VAdd *vadd = (VAdd*) node;
  if (vadd->GetLayer() != m_fromLayer || m_type != vadd->m_type) {
    return false;
  }
  if (m_vtype == ROWVECTOR) {
    if (!(*(vadd->GetInputN(0)) <= m_bs.GetSize())
	&& !(*(vadd->GetInputN(1)) <= m_bs.GetSize())) {
      return true;
    } else {
      return false;
    }
  } else if (m_vtype == COLVECTOR) {
    if (!(*(vadd->GetInputM(0)) <= m_bs.GetSize())
	&& !(*(vadd->GetInputM(1)) <= m_bs.GetSize())) {
      return true;
    } else {
      return false;
    }
  } else {
    throw;
  }
  return false;
}

void VAddLoopRef::Apply(Node *node) const
{
  VAdd *vadd = (VAdd*) node;

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

  VAdd *newVAdd = new VAdd(vadd->m_vecType, vadd->m_layer, vadd->m_type);
  newVAdd->SetLayer(m_toLayer);

  newVAdd->AddInput(split0, 1);
  newVAdd->AddInput(split1, 1);

  CombineSingleIter *com0 = split0->CreateMatchingCombine(0);
  CombineSingleIter *com1 = split1->CreateMatchingCombine(1, 1, newVAdd, 0);

  Poss *loopPoss = new Poss(2, com0, com1);
  RealLoop* loop;
  if (m_type == REAL_SINGLE) {
    loop = new RealLoop(LLDLALOOP, loopPoss, LLDLAMuSingle);
  } else if (m_type == REAL_DOUBLE) {
    loop = new RealLoop(LLDLALOOP, loopPoss, LLDLAMuDouble);
  } else {
    cout << "Error: Bad m_type in vadd apply\n";
    throw;
  }
  loop->SetDimName(m_vtype == COLVECTOR ? DIMM : DIMN);
  
  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(1), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

VAddLowerLayer::VAddLowerLayer(Layer fromLayer, Layer toLayer, Size bs, Type type)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_bs = bs;
  m_type = type;
  m_regWidth = arch->VecRegWidth(m_type);
}

bool VAddLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == VAdd::GetClass()) {
    const VAdd *vadd = (VAdd*) node;
    if (vadd->GetLayer() != m_fromLayer || m_type != vadd->m_type) {
      return false;
    }

    if (vadd->m_vecType == ROWVECTOR) {
      if (*vadd->GetInputN(0) <= m_bs &&
	  *vadd->GetInputN(1) <= m_bs) {
	return true;
      } else {
	return false;
      }
    } else if (vadd->m_vecType == COLVECTOR) {
      if (*vadd->GetInputM(0) <= m_bs &&
	  *vadd->GetInputM(1) <= m_bs) {
	return true;
      } else {
	return false;
      }
    }
  }
  cout << "Error: Applying VAddLowerLayer to non-VAdd node\n";
  throw;
}

void VAddLowerLayer::Apply(Node *node) const
{
  VAdd *vadd = (VAdd*) node;
  vadd->SetLayer(m_toLayer);
  return;
}

string VAddLowerLayer::GetType() const
{
  return "VAdd lower layer " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer) + " type " + std::to_string((long long int) m_type);
}

VAddToRegArith::VAddToRegArith(Layer fromLayer, Layer toLayer, Type type)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_type = type;
  m_regWidth = arch->VecRegWidth(m_type);
}


bool VAddToRegArith::CanApply(const Node* node) const
{
  if (node->GetNodeClass() == VAdd::GetClass()) {
    VAdd* vadd = (VAdd*) node;
    return m_type == vadd->m_type;
  }
  return false;
}

void VAddToRegArith::Apply(Node* node) const
{
  VAdd* vadd = (VAdd*) node;

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

  LoadToRegs* loadX = new LoadToRegs(m_type);
  loadX->AddInput(splitX, 1);

  LoadToRegs* loadY = new LoadToRegs(m_type);
  loadY->AddInput(splitY, 1);

  Add* add = new Add(m_type);
  add->AddInput(loadX, 0);
  add->AddInput(loadY, 0);

  StoreFromRegs* storeToY = new StoreFromRegs(m_type);
  storeToY->AddInput(add, 0);
  storeToY->AddInput(splitY, 1);

  CombineSingleIter* comX = splitX->CreateMatchingCombine(0);
  CombineSingleIter* comY = splitY->CreateMatchingCombine(1, 1, storeToY, 0);
  
  Poss* loopPoss = new Poss(2, comX, comY);
  RealLoop* loop;
  if (m_type == REAL_SINGLE) {
    loop = new RealLoop(LLDLALOOP, loopPoss, LLDLAMuSingle);
  } else if (m_type == REAL_DOUBLE) {
    loop = new RealLoop(LLDLALOOP, loopPoss, LLDLAMuDouble);
  } else {
    cout << "Error: Bad m_type in vadd apply\n";
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
    + " to " + LayerNumToStr(m_fromLayer) + " type " + std::to_string((long long int) m_type);
}

#endif // DOLLDLA
