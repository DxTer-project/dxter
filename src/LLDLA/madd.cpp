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

#include "madd.h"
#include "vadd.h"

#if DOLLDLA

MAdd::MAdd(Type type, Layer layer)
{
  m_type = type;
  m_layer = layer;
}

void MAdd::PrintCode(IndStream &out)
{
  if (m_layer == ABSLAYER) {
    *out << "simple_add( " <<
      InputDataType(0).m_numRowsVar << ", " <<
      InputDataType(0).m_numColsVar << ", " <<
      GetInputName(0).str() << ", " <<
      InputDataType(0).m_rowStrideVar << ", " <<
      InputDataType(0).m_colStrideVar << ", " <<
      GetInputName(1).str() << ", " <<
      InputDataType(1).m_rowStrideVar << ", " <<
      InputDataType(1).m_colStrideVar << ");\n";
    
    return;
  }
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


void MAdd::PrintRowStride(IndStream &out)
{
  *out << "row_stride_add_2x2( " <<
    GetInputName(0).str() << ", " <<
    GetInputName(1).str() << ");\n";
}

void MAdd::PrintColStride(IndStream &out)
{
  *out << "col_stride_add_2x2( " <<
    GetInputName(0).str() << ", " <<
    GetInputName(1).str() << ");\n";
}

void MAdd::PrintGeneralStride(IndStream &out)
{
  *out << "gen_stride_add_2x2( " <<
    GetInputName(0).str() << ", " <<
    GetInputName(1).str() << ");\n";
}

void MAdd::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp::Prop();
    
    if ((*GetInputM(0) != *GetInputM(1)) || (*GetInputN(0) != *GetInputN(1))) {
      cout << "ERROR: Cannot MAdd two matrices of different dimension\n";
      throw;
    }

    if (m_layer == LLDLAPRIMITIVELAYER) {
      if ((*GetInputM(0) != LLDLA_MU) || (*GetInputN(0) != LLDLA_MU)
	  || (*GetInputM(1) != LLDLA_MU) || (*GetInputN(1) != LLDLA_MU)) {
	cout << "ERROR: MAdd of matrices that do not have LLDLA_MU dimensions in LLDLAPRIMITIVELAYER\n";
      }
    }

    m_cost = ZERO;
  }
}

Node* MAdd::BlankInst()
{
  return new MAdd(REAL, LLDLAPRIMITIVELAYER);
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
      throw;
    }
}

string MAddLoopRef::GetType() const
{
  return "MAdd";
}

bool MAddLoopRef::CanApply(const Node *node) const
{
  const MAdd *madd = (MAdd*) node;
  if (madd->GetLayer() != m_fromLayer) {
    return false;
  }
  if (m_dim == DIMM) {
    return !(*(madd->GetInputM(0)) <= m_bs) && !(*(madd->GetInputM(1)) <= m_bs);
  } else if (m_dim == DIMN) {
    return !(*(madd->GetInputN(0)) <= m_bs) && !(*(madd->GetInputN(1)) <= m_bs);
  } else {
    return false;
  }
}

void MAddLoopRef::Apply(Node *node) const
{
  MAdd *madd = (MAdd*) node;
  
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

  MAdd *newMAdd = new MAdd(madd->m_type, m_toLayer);
  newMAdd->AddInput(split0, 1);
  newMAdd->AddInput(split1, 1);

  CombineSingleIter *com0 = split0->CreateMatchingCombine(0);
  CombineSingleIter *com1 = split1->CreateMatchingCombine(1, 1, newMAdd, 0);

  Poss *loopPoss = new Poss(2, com0, com1);

  RealLoop *loop = new RealLoop(LLDLALOOP, loopPoss, USELLDLAMU);
  loop->SetDimName(m_dim == DIMM ? DIMM : DIMN);

  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(1), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

string MAddToVAddLoopRef::GetType() const
{
  return "MAddToVAddLoopRef";
}

bool MAddToVAddLoopRef::CanApply(const Node *node) const
{
  const MAdd *madd = (MAdd*) node;
  if (madd->GetLayer() != m_fromLayer) {
    return false;
  }
  if (m_dim == DIMM) {
    return (*(madd->GetInputN(0)) == LLDLA_MU) && (*(madd->GetInputN(1)) == LLDLA_MU) &&
      !(*(madd->GetInputM(0)) <= LLDLA_MU) && !(*(madd->GetInputM(1)) <= LLDLA_MU);
  } else if (m_dim == DIMN) {
    return (*(madd->GetInputM(0)) == LLDLA_MU) && (*(madd->GetInputM(0)) == LLDLA_MU) &&
      !(*(madd->GetInputN(0)) <= LLDLA_MU) && !(*(madd->GetInputN(1)) <= LLDLA_MU);
  } else {
    return false;
  }  
}

void MAddToVAddLoopRef::Apply(Node *node) const
{
  MAdd *madd = (MAdd*) node;

  SplitSingleIter *split0 = new SplitSingleIter(m_dim == DIMM ? PARTDOWN : PARTRIGHT, POSSTUNIN, true);
  split0->AddInput(madd->Input(0), madd->InputConnNum(0));

  SplitSingleIter *split1 = new SplitSingleIter(m_dim == DIMN ? PARTDOWN : PARTRIGHT, POSSTUNIN, false);
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

  VAdd *newVAdd = new VAdd(m_dim == DIMM ? ROWVECTOR : COLVECTOR, m_toLayer, madd->m_type);
  newVAdd->AddInput(split0, 1);
  newVAdd->AddInput(split1, 1);

  CombineSingleIter *com0 = split0->CreateMatchingCombine(0);
  CombineSingleIter *com1 = split1->CreateMatchingCombine(1, 1, newVAdd, 0);

  
  Poss *loopPoss = new Poss(2, com0, com1);

  RealLoop *loop = new RealLoop(LLDLALOOP, loopPoss, USEUNITBS);
  loop->SetDimName(m_dim == DIMM ? DIMM : DIMN);

  node->m_poss->AddPSet(loop);
  node->RedirectChildren(loop->OutTun(1), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

bool MAddLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == MAdd::GetClass()) {
    const MAdd *madd = (MAdd*) node;
    if (madd->GetLayer() != m_fromLayer) {
      return false;
    }
    if (*(madd->GetInputM(1)) <= m_bs &&
	*(madd->GetInputN(1)) <= m_bs) {
      return true;
    } else {
      return false;
    }
  }
  else {
    throw;
  }
}

void MAddLowerLayer::Apply(Node *node) const
{
  MAdd *madd = (MAdd*) node;
  madd->SetLayer(m_toLayer);
}

string MAddLowerLayer::GetType() const
{
  return "VAdd lower layer " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer);
}

#endif // DOLLDLA
