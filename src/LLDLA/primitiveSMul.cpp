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

#endif //DOLLDLA
