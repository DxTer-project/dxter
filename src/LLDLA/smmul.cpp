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
#include "smmul.h"
#include "svmul.h"

#if DOLLDLA

SMMul::SMMul(Layer layer)
{
  SetLayer(layer);
}

void SMMul::PrintCode(IndStream &out)
{

  if (GetLayer() == ABSLAYER) {
    if (GetDataType() == REAL_DOUBLE) {
      out.Indent();
      *out << "simple_smul( ";
    } else if (GetDataType() == REAL_SINGLE) {
      out.Indent();
      *out << "simple_smul_float( ";
    }
    *out << InputDataType(1).m_numRowsVar << ", " <<
      InputDataType(1).m_numColsVar << ", " <<
      GetInputName(0).str() << ", " <<
      GetInputName(1).str() << ", " <<
      InputDataType(1).m_rowStrideVar << ", " <<
      InputDataType(1).m_colStrideVar << ");\n";
    return;
  }

  if (m_layer != LLDLAPRIMITIVELAYER) {
    cout << "ERROR: Attempt to generate code from non-primitive scalar matrix multiply\n";
    LOG_FAIL("replacement for throw call");
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


void SMMul::PrintRowStride(IndStream &out)
{
  *out << "row_stride_smul_2x2( " <<
    GetInputName(0).str() << ", " <<
    GetInputName(1).str() << ", " <<
    InputDataType(1).m_rowStrideVar << ");\n";
}

void SMMul::PrintColStride(IndStream &out)
{
  *out << "col_stride_smul_2x2( " <<
    GetInputName(0).str() << ", " <<
    GetInputName(1).str() << ", " <<
    InputDataType(1).m_rowStrideVar << ");\n";

}

void SMMul::PrintGeneralStride(IndStream &out)
{
  *out << "gen_stride_mmul_2x2_2x2( " <<
    GetInputName(0).str() << ", " <<
    GetInputName(1).str() << ", " <<
    InputDataType(1).m_rowStrideVar << ", " <<
    InputDataType(1).m_colStrideVar << ");\n";
}


void SMMul::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2, 1>::Prop();

    if (!((DLANode*) Input(0))->IsScalar(InputConnNum(0))) {
      cout << "ERROR: SMMul input 0 is not a scalar\n";
    }
    
    if (GetLayer() == LLDLAPRIMITIVELAYER || GetLayer() == LLDLAMIDLAYER) {
      if (*GetInputM(1) != GetVecRegWidth() || *GetInputN(1) != GetVecRegWidth()) {
	GetInputM(1)->Print();
	cout << endl;
	GetInputN(1)->Print();
	cout << "ERROR: SMMul input 1 must be an LLDLAMU by LLDLAMU matrix\n";
      }
    }      

    if (m_layer == ABSLAYER) {
      m_cost = GetInputM(1)->SumProds11(*GetInputN(1));
    } else {
      m_cost = ZERO;
    }
  }
}

Phase SMMul::MaxPhase() const
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
    }
}

Node* SMMul::BlankInst()
{
  return new SMMul(ABSLAYER);
}

NodeType SMMul::GetType() const
{
  return "SMMul" + LayerNumToStr(GetLayer());
}


SMulToSVMul::SMulToSVMul(Layer fromLayer, Layer toLayer, VecType vecType)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_vecType = vecType;
}

string SMulToSVMul::GetType() const
{
  switch (m_vecType)
    {
    case (ROWVECTOR):
      return "SMulToSVMul - dim m";
    case (COLVECTOR):
      return "SMulToSVMul - dim n";
    default:
      LOG_FAIL("replacement for throw call");
    }  
}

bool SMulToSVMul::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == SMMul::GetClass()) {
    const SMMul *smmul = static_cast<const SMMul*>(node);
    if (smmul->GetLayer() != m_fromLayer) {
      return false;
    }
    if (m_vecType == ROWVECTOR) {
      return *(smmul->GetInputM(1)) > 1;
    } else {
      return *(smmul->GetInputN(1)) > 1;
    }
  }
  cout << "ERROR: Applying SMMulToSVMul to non SMMul node\n";
  cout << "Node has class " << node->GetNodeClass() << endl;
  LOG_FAIL("replacement for throw call");
}

void SMulToSVMul::Apply(Node *node) const
{
  SMMul* mul = (SMMul*)node;


  //Create a split for the input matrix
  // If we're splitting on the m dimension, then the split moves down
  //    (i.e., it's horizontal)
  // If we're splitting on the n dimension, then the split moves right
  //    (i.e., it's vertical)
  SplitSingleIter* split = new SplitSingleIter(m_vecType == ROWVECTOR ? PARTDOWN : PARTRIGHT, POSSTUNIN, true);
  // Add input, which is the matrix input to mul
  split->AddInput(mul->Input(1), mul->InputConnNum(1));
  //Set the update statuses
  //If we're moving down, then everything above the thick line is fully updated
  if (m_vecType == ROWVECTOR) {
    split->SetUpStats(FULLUP, FULLUP,
		      NOTUP, NOTUP);
  }
  //If we're moving right, then everything left of the thick line is fully updated
  else {
    split->SetUpStats(FULLUP, NOTUP,
		       FULLUP, NOTUP);
  }

  //Independent iterations
  split->SetIndepIters();

  //Create a loop tunnel - the scalar is input to each iteration of the loop
  LoopTunnel* scalarTun = new LoopTunnel(POSSTUNIN);
  //Wire up the scalar input
  scalarTun->AddInput(mul->Input(0), mul->InputConnNum(0));
  //Always updated since it doesn't change
  scalarTun->SetAllStats(FULLUP);
  //Independent iterations
  scalarTun->SetIndepIters();
  

  //Create a new SMul or the same type and in my m_toLayer layer
  SVMul* newMul = new SVMul(m_toLayer);
  newMul->SetLayer(m_toLayer);
  newMul->AddInput(scalarTun, 0);
  newMul->AddInput(split, 1);

  //Create an output tunnel for the scalar just for consistency
  LoopTunnel* scalarTunOut = new LoopTunnel(POSSTUNOUT);
  //Wire two inputs, the first is for the data that is passed out
  // and the second is to the input tunnel to which this one
  // should match
  scalarTunOut->AddInput(scalarTun, 0);
  scalarTunOut->AddInput(scalarTun, 1);
  scalarTunOut->CopyTunnelInfo(scalarTun);

  //Create an output tunnel for the matrix and overwrite
  // the 1st (0-based) partition of the output matrix
  CombineSingleIter* com = split->CreateMatchingCombine(1,
							1, newMul, 0);
  
  //Put all of this into single poss (this constructor
  // will recursively move up the flow of data, adding
  // all nodes it finds
  Poss* loopPoss = new Poss(2, scalarTunOut, com);
  //Put that poss into a loop - it's LLDLALOOP type and
  // uses the m_regWidth blocksize
  RealLoop* loop = new RealLoop(LLDLALOOP, loopPoss, UnitBS);

  //Set the dimension over which this loop iterates
  if (m_vecType == ROWVECTOR) {
  loop->SetDimName(DIMM);
  } else {
  loop->SetDimName(DIMN);
  }
  
  //Add the loop to the node's owning Poss
  node->m_poss->AddPSet(loop);
  //Redirect children of the node to use the loop's output
  // (it's the 1st (0-base) output since the 0th is the scalar)
  node->RedirectChildren(loop->OutTun(1), 0);
  //Delete the node and recursively move up to delete any stragling
  // nodes (in this case, there shouldn't be any)
  node->m_poss->DeleteChildAndCleanUp(node);
}

SMulToScalarArith::SMulToScalarArith(Layer fromLayer, Layer toLayer, DimName dim)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_dim = dim;
}

string SMulToScalarArith::GetType() const
{
  switch (m_dim) {
  case (DIMM):
    return "SMulToScalarArithMDIM";
  case (DIMN):
    return "SMulToScalarArithNDIM";
  default:
    cout << "Invalid dimension in SMulToScalarArith::GetType()" << endl;
    LOG_FAIL("replacement for throw call");
  }
}

bool SMulToScalarArith::CanApply(const Node* node) const
{
  if (node->GetNodeClass() == SMMul::GetClass()) {
    const SMMul *mul = (SMMul*)node;
    if (mul->GetLayer() != m_fromLayer) {
      return false;
    }
    return true;
  }
  cout << "ERROR: Called SMulToScalarArith::CanApply on non SMul node" << endl;
  LOG_FAIL("replacement for throw call");
}

void SMulToScalarArith::Apply(Node* node) const
{
  
}

SMulLowerLayer::SMulLowerLayer(Layer fromLayer, Layer toLayer, Size bs)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_bs = bs;
}

bool SMulLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == SMMul::GetClass()) {
    const SMMul *smul = (SMMul*)node;
    if (smul->GetLayer() != m_fromLayer)
      return false;
    if (*(smul->GetInputM(1)) <= m_bs &&
	*(smul->GetInputN(1)) <= m_bs)
      return true;
    else
      return false;
  }
  else
    LOG_FAIL("replacement for throw call");  
}

void SMulLowerLayer::Apply(Node *node) const
{
  SMMul *smul = (SMMul*)node;
  smul->SetLayer(m_toLayer);
}

string SMulLowerLayer::GetType() const
{ 
  return "SMul lower layer " + LayerNumToStr(m_fromLayer) 
  + " to " + LayerNumToStr(m_toLayer);
}


#endif //DOLLDLA
