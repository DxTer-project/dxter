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
#include "svmul.h"

#if DOLLDLA

SVMul::SVMul(VecType vecType, Layer layer, Type type)
{
  m_vecType = vecType;
  m_layer = layer;
  m_type = type;
}

void SVMul::PrintCode(IndStream &out)
{
  /*  if (m_layer != LLDLAPRIMITIVELAYER) {
    cout << "ERROR: Attempt to generate code from non-primitive scalar vector multiply\n";
    throw;
    }*/
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


void SVMul::PrintRowStride(IndStream &out)
{
  if (m_vecType == COLVECTOR) {
    *out << "row_stride_smul_2x1( " <<
      GetInputName(0).str() << ", " <<
      GetInputName(1).str() << ");\n";
  } else {
    *out << "row_stride_smul_1x2( " <<
      GetInputName(0).str() << ", " <<
      GetInputName(1).str() << ");\n";
  }
}

void SVMul::PrintColStride(IndStream &out)
{
  if (m_vecType == COLVECTOR) {
    *out << "col_stride_smul_2x1( " <<
      GetInputName(0).str() << ", " <<
      GetInputName(1).str() << ");\n";
  } else {
    *out << "col_stride_smul_1x2( " <<
      GetInputName(0).str() << ", " <<
      GetInputName(1).str() << ");\n";
  }

}

void SVMul::PrintGeneralStride(IndStream &out)
{
  if (m_vecType == COLVECTOR) {
    *out << "gen_stride_smul_2x1( " <<
      GetInputName(0).str() << ", " <<
      GetInputName(1).str() << ");\n";
  } else {
    *out << "gen_stride_smul_1x2( " <<
      GetInputName(0).str() << ", " <<
      GetInputName(1).str() << ");\n";
  }
}

void SVMul::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();

    if (!((DLANode*) Input(0))->IsScalar(InputConnNum(0))) {
      cout << "ERROR: SVMul input 0 is not a scalar\n";
      throw;
    }

    VectorOpInputDimensionCheck(1);

    m_cost = ZERO;
  }
}

void SVMul::VectorOpInputDimensionCheck(unsigned int inputNum)
{
  if (m_vecType == ROWVECTOR && *GetInputM(inputNum) != 1) {
    cout << "ERROR: " << GetType() << " input # " << inputNum << " has more than 1 row\n";
    throw;
  } else if (m_vecType == COLVECTOR && *GetInputN(inputNum) != 1) {
    cout << "ERROR: " << GetType() << " input # " << inputNum  << " has more than 1 column\n";
  }
  
  if (m_layer == LLDLAPRIMITIVELAYER) {
    if (m_vecType == ROWVECTOR && *GetInputN(inputNum) != LLDLA_MU) {
      cout << "ERROR: " << GetType() << " input # " << inputNum << " does not have LLDLA_MU columns\n";
      throw;
    } else if(m_vecType == COLVECTOR && *GetInputM(inputNum) != LLDLA_MU) {
      cout << "ERROR: " << GetType() << " input # " << inputNum << " does not have LLDLA_MU rows\n";
      throw;
    }
  }
}

Node* SVMul::BlankInst()
{
  return new SVMul(ROWVECTOR, LLDLAPRIMITIVELAYER, REAL);
}

NodeType SVMul::GetType() const
{
  return "SVMul" +  LayerNumToStr(GetLayer());
}

void SVMul::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLAOp<2,1>::Duplicate(orig, shallow, possMerging);
  const SVMul *rhs = (SVMul*)orig;
  m_type = rhs->m_type;
  m_vecType = rhs->m_vecType;
}

Phase SVMul::MaxPhase() const 
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

string SVMulLoopRef::GetType() const
{
  switch(m_vtype)
    {
    case(ROWVECTOR):
      return "SVMulLoopRef - row vector";
    case(COLVECTOR):
      return "SVMulLoopRef - column vector";
    default:
      throw;
    }
}


bool SVMulLoopRef::CanApply(const Node *node) const
{
  const SVMul *svmul = (SVMul*) node;
  if (svmul->GetLayer() != m_fromLayer) {
    return false;
  }
  if (m_vtype == ROWVECTOR) {
    return !(*(svmul->GetInputN(1)) <= BSSizeToSize(m_bs));
  } 
  else if (m_vtype == COLVECTOR) {
    return !(*(svmul->GetInputM(1)) <= BSSizeToSize(m_bs));
  } 
  else {
    throw;
  }
  return false;
}

void SVMulLoopRef::Apply(Node *node) const
{
  SVMul *svmul = (SVMul*) node;

  SplitSingleIter *split = new SplitSingleIter(m_vtype == COLVECTOR ? PARTDOWN : PARTRIGHT, POSSTUNIN, true);
  split->AddInput(svmul->Input(1), svmul->InputConnNum(1));

  if (m_vtype == COLVECTOR) {
    split->SetUpStats(FULLUP, FULLUP,
		      NOTUP, NOTUP);		     
  } else {
    split->SetUpStats(FULLUP, NOTUP,
		       FULLUP, NOTUP);
  }

  split->SetIndepIters();

  LoopTunnel *scalarTun = new LoopTunnel(POSSTUNIN);
  scalarTun->AddInput(svmul->Input(0),svmul->InputConnNum(0));
  scalarTun->SetAllStats(FULLUP);
  scalarTun->SetIndepIters();

  SVMul *newMul = new SVMul(svmul->m_vecType, svmul->m_layer, svmul->m_type);
  newMul->SetLayer(m_toLayer);

  newMul->AddInput(scalarTun, 0);
  newMul->AddInput(split, 1);

  LoopTunnel *scalarTunOut = new LoopTunnel(POSSTUNOUT);
  scalarTunOut->AddInput(scalarTun, 0);
  scalarTunOut->AddInput(scalarTun, 0);
  scalarTunOut->CopyTunnelInfo(scalarTun);

  CombineSingleIter *com = split->CreateMatchingCombine(1, 
					      1, newMul, 0);
  
  Poss *loopPoss = new Poss(2, scalarTunOut, com);
  
  Loop *loop = new Loop(LLDLALOOP, loopPoss, USELLDLAMU);
  // Row vectors are partitioned in the N dimension, column vectors in the M dimension
  loop->SetDimName(m_vtype == COLVECTOR ? DIMM : DIMN);
  
  node->m_poss->AddLoop(loop);
  node->RedirectChildren(loop->OutTun(1), 0);
  node->m_poss->DeleteChildAndCleanUp(node);

}

bool SVMulLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == SVMul::GetClass()) {
    const SVMul *smul = (SVMul*) node;
    if (smul->GetLayer() != m_fromLayer) {
      return false;
    }
    if (*(smul->GetInputM(1)) <= m_bs &&
	*(smul->GetInputN(1)) <= m_bs) {
      return true;
    } else {
      return false;
    }
  }
  else {
    throw;
  }
}

void SVMulLowerLayer::Apply(Node *node) const
{
  SVMul *svmul = (SVMul*) node;
  svmul->SetLayer(m_toLayer);
}

string SVMulLowerLayer::GetType() const
{
  return "SVMul lower layer " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer);
}

#endif // DOLLDLA
