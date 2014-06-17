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
#include "vadd.h"

#if DOLLDLA

VAdd::VAdd(VecType vecType, Layer layer, Type type)
{
  m_vecType = vecType;
  m_layer = layer;
  m_type = type;
}

void VAdd::PrintCode(IndStream &out)
{
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


void VAdd::PrintRowStride(IndStream &out)
{
  if (m_vecType == COLVECTOR) {
    *out << "row_stride_add_2x1( " <<
      GetInputName(0).str() << ", " <<
      GetInputName(1).str() << ");\n";
  } else {
    *out << "row_stride_add_1x2( " <<
      GetInputName(0).str() << ", " <<
      GetInputName(1).str() << ");\n";
  }
}

void VAdd::PrintColStride(IndStream &out)
{
  *out << "COL STRIDE not yet implemented\n";
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
    DLAOp::Prop();
    
    VectorOpInputDimensionCheck(0);
    VectorOpInputDimensionCheck(1);

    if ((*GetInputM(0) != *GetInputM(1)) || (*GetInputN(0) != *GetInputN(1))) {
      cout << "ERROR: Cannot VAdd two vectors of different dimension\n";
      throw;
    }

    m_cost = ZERO;
  }

}

void VAdd::VectorOpInputDimensionCheck(unsigned int inputNum)
{
  if (m_vecType == ROWVECTOR && *GetInputM(inputNum) != 1) {
    cout << "ERROR: " << GetType() << " input # " << inputNum << " has more than 1 row\n";
    throw;
  } else if (m_vecType == COLVECTOR && *GetInputN(inputNum) != 1) {
    cout << "ERROR: " << GetType() << " input # " << inputNum  << " has more than 1 column\n";
    throw;
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

Node* VAdd::BlankInst()
{
  return new VAdd(ROWVECTOR, LLDLAPRIMITIVELAYER, REAL);
}

NodeType VAdd::GetType() const
{
  return "VADD" + LayerNumToStr(GetLayer());
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
  if (vadd->GetLayer() != m_fromLayer) {
    return false;
  }
  if (m_vtype == ROWVECTOR) {
    if (!(*(vadd->GetInputN(0)) <= BSSizeToSize(m_bs))
	&& !(*(vadd->GetInputN(1)) <= BSSizeToSize(m_bs))) {
      return true;
    } else {
      return false;
    }
  } else if (m_vtype == COLVECTOR) {
    if (!(*(vadd->GetInputM(0)) <= BSSizeToSize(m_bs))
	&& !(*(vadd->GetInputM(1)) <= BSSizeToSize(m_bs))) {
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

  Split *split1 = new Split(m_vtype == COLVECTOR ? PARTDOWN : PARTRIGHT, POSSTUNIN, true);
  split1->AddInput(vadd->Input(1), vadd->InputConnNum(1));

  Split *split0 = new Split(m_vtype == COLVECTOR ? PARTDOWN : PARTRIGHT, POSSTUNIN, false);
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

  Combine *com0 = split0->CreateMatchingCombine(0);
  Combine *com1 = split1->CreateMatchingCombine(1, 1, newVAdd, 0);

  Poss *loopPoss = new Poss(2, com0, com1);

  Loop *loop = new Loop(LLDLALOOP, loopPoss, USELLDLAMU);
  loop->SetDimName(m_vtype == COLVECTOR ? DIMM : DIMN);
  
  node->m_poss->AddLoop(loop);
  node->RedirectChildren(loop->OutTun(1), 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

bool VAddLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == VAdd::GetClass()) {
    const VAdd *vadd = (VAdd*) node;
    if (vadd->GetLayer() != m_fromLayer) {
      return false;
    }
    if (*(vadd->GetInputM(1)) <= m_bs &&
	*(vadd->GetInputN(1)) <= m_bs) {
      return true;
    } else {
      return false;
    }
  }
  else {
    throw;
  }
}

void VAddLowerLayer::Apply(Node *node) const
{
  VAdd *vadd = (VAdd*) node;
  vadd->SetLayer(m_toLayer);
}

string VAddLowerLayer::GetType() const
{
  return "VAdd lower layer " + LayerNumToStr(m_fromLayer)
    + " to " + LayerNumToStr(m_toLayer);
}

#endif // DOLLDLA
