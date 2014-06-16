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
#include "primitiveSMul.h"

#if DOLLDLA

PrimitiveSMul::PrimitiveSMul(Type type)
{
  m_layer = LLDLAPRIMITIVELAYER;
  m_type = type;
}

void PrimitiveSMul::PrintCode(IndStream &out)
{
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


void PrimitiveSMul::PrintRowStride(IndStream &out)
{
  *out << "row_stride_smul_2x2( " <<
    GetInputName(0).str() << ", " <<
    GetInputName(1).str() << ");\n";
}

void PrimitiveSMul::PrintColStride(IndStream &out)
{
  *out << "COL STRIDE not yet implemented\n";
}

void PrimitiveSMul::PrintGeneralStride(IndStream &out)
{
  *out << "gen_stride_smul_2x2( " <<
    GetInputName(0).str() << ", " <<
    GetInputName(1).str() << ");\n";
}


void PrimitiveSMul::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp::Prop();

    if (m_layer != LLDLAPRIMITIVELAYER) {
      cout << "ERROR: PrimitiveSMul appears in layer " <<  LayerNumToStr(m_layer) << "\n" ;
      throw;
    }

    if (*GetInputM(1) != LLDLA_MU || *GetInputN(1) != LLDLA_MU) {
      GetInputM(1)->Print();
      cout << endl;
      GetInputN(1)->Print();
      cout << "ERROR: PrimitiveSMul input 1 must be an LLDLAMU by LLDLAMU matrix\n";
    }

    if (!((DLANode*) Input(0))->IsScalar(InputConnNum(0))) {
      cout << "ERROR: PrimitiveSMul input 0 is not a scalar\n";
    }
    
    m_cost = ZERO;
  }
}

Node* PrimitiveSMul::BlankInst()
{
  return new PrimitiveSMul(REAL);
}

NodeType PrimitiveSMul::GetType() const
{
  return "PrimitiveSMul" + LayerNumToStr(GetLayer());
}

string SMulLoopRef::GetType() const
{
  switch (m_dim)
    {
    case (DIMM):
      return "SMulLoopRef - dim m";
    case (DIMN):
      return "SMulLoopRef - dim n";
    default:
      throw;
    }
  
}

bool SMulLoopRef::CanApply(const Node *node) const
{
  const PrimitiveSMul *mul = (PrimitiveSMul*)node;
  if (mul->GetLayer() != m_fromLayer)
    return false;
  
  if (m_dim == DIMM) {
    if (!(*(mul->GetInputM(1)) <= BSSizeToSize(m_bs)))
      return true;
  }
  else if (m_dim == DIMN) {
    if (!(*(mul->GetInputN(1)) <= BSSizeToSize(m_bs)))
      return true;
  }
  else
    throw;
  return false;
}

void SMulLoopRef::Apply(Node *node) const
{
  PrimitiveSMul *mul = (PrimitiveSMul*)node;


  //Create a split for the input matrix
  // If we're splitting on the m dimension, then the split moves down 
  //    (i.e., it's horizontal)
  // If we're splitting on the n dimension, then the split moves right
  //    (i.e., it's vertical)
  Split *split = new Split(m_dim==DIMM ? PARTDOWN : PARTRIGHT, POSSTUNIN);
  // Add input, which is the matrix input to mul
  split->AddInput(mul->Input(1), mul->InputConnNum(1));
  //Set the update statuses
  //If we're moving down, then everything above the thick line is fully updated
  if (m_dim == DIMM) {
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
  LoopTunnel *scalarTun = new LoopTunnel(POSSTUNIN);
  //Wire up the scalar input
  scalarTun->AddInput(mul->Input(0), mul->InputConnNum(0));
  //Always updated since it doesn't change
  scalarTun->SetAllStats(FULLUP);
  //Independent iterations
  scalarTun->SetIndepIters();
  

  //Create a new SMul or the same type and in my m_toLayer layer
  PrimitiveSMul *newMul = new PrimitiveSMul(mul->m_type);
  newMul->SetLayer(m_toLayer);

  //Wire inputs - the scalar loop tunnel and 
  //   the 1st (0-based) partition of the matrix split tunnel
  newMul->AddInput(scalarTun, 0);
  newMul->AddInput(split, 1);

  //Create an output tunnel for the scalar just for consistency
  LoopTunnel *scalarTunOut = new LoopTunnel(POSSTUNOUT);
  //Wire two inputs, the first is for the data that is passed out
  // and the second is to the input tunnel to which this one
  // should match
  scalarTunOut->AddInput(scalarTun, 0);
  scalarTunOut->AddInput(scalarTun, 0);
  scalarTunOut->CopyTunnelInfo(scalarTun);

  
  //Create an output tunnel for the matrix and overwrite
  // the 1st (0-based) partition of the output matrix
  Combine *com = split->CreateMatchingCombine(1,
					     1, newMul, 0);
  
  //Put all of this into single poss (this constructor
  // will recursively move up the flow of data, adding
  // all nodes it finds
  Poss *loopPoss = new Poss(2, scalarTunOut, com);
  //Put that poss into a loop - it's LLDLALOOP type and
  // uses the LLDLA_MU blocksize
  Loop *loop = new Loop(LLDLALOOP, loopPoss, USELLDLAMU);

  //Set the dimension over which this loop iterates
  loop->SetDimName(m_dim);
  
  //Add the loop to the node's owning Poss
  node->m_poss->AddLoop(loop);
  //Redirect children of the node to use the loop's output
  // (it's the 1st (0-base) output since the 0th is the scalar)
  node->RedirectChildren(loop->OutTun(1), 0);
  //Delete the node and recursively move up to delete any stragling
  // nodes (in this case, there shouldn't be any)
  node->m_poss->DeleteChildAndCleanUp(node);
}

#endif //DOLLDLA
