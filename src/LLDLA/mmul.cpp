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
#include "layers.h"
#include "mmul.h"
#include "loopSupport.h"

#if DOLLDLA

MMul::MMul(Type type, Layer layer)
  : Gemm(layer, NORMAL, NORMAL, COEFONE, COEFONE, type)
{
  m_regWidth = arch->VecRegWidth(type);
}

Phase MMul::MaxPhase() const {
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

void MMul::PrintCode(IndStream &out)
{
  if (m_layer == ABSLAYER) {

    if (m_alpha.m_val == COEFVALONE && m_beta.m_val == COEFVALONE) {
      if (m_type == REAL_DOUBLE) {
	*out << "simple_mmul( ";
      } else if (m_type == REAL_SINGLE) {
	*out << "simple_mmul_float( ";
      }
      *out << InputDataType(2).m_numRowsVar << ", " <<
	InputDataType(2).m_numColsVar << ", " <<
	InputDataType(0).m_numColsVar << ", " <<
	GetInputName(0).str() << ", " <<
	InputDataType(0).m_rowStrideVar << ", " <<
	InputDataType(0).m_colStrideVar << ", " <<
	GetInputName(1).str() << ", " <<
	InputDataType(1).m_rowStrideVar << ", " <<
	InputDataType(1).m_colStrideVar << ", " <<
	GetInputName(2).str() << ", " <<
	InputDataType(2).m_rowStrideVar << ", " <<
	InputDataType(2).m_colStrideVar << ");\n";
      return;
    }
    else {
      LOG_FAIL("replacement for throw call");
    }
  }

  if (m_layer != LLDLAPRIMITIVELAYER) {
    cout << "ERROR: Attempt to generate code from non-primitive MMul\n";
    LOG_FAIL("replacement for throw call");
  }

  const DataTypeInfo &inInfo = InputDataType(2);
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

void MMul::PrintRowStride(IndStream &out)
{
  if (m_alpha.m_val == COEFVALONE && m_beta.m_val == COEFVALONE) {
    *out << "row_stride_mmul_2x2_2x2( " <<
      GetInputName(0).str() << ", " <<
      InputDataType(0).m_rowStrideVar << ", " << 
      GetInputName(1).str() << ", " <<
      InputDataType(1).m_rowStrideVar << ", " <<
      GetInputName(2).str() << ", " <<
      InputDataType(2).m_rowStrideVar << " );\n";
  }
  else
    LOG_FAIL("replacement for throw call");
}

void MMul::PrintColStride(IndStream &out)
{
  if (m_alpha.m_val == COEFVALONE && m_beta.m_val == COEFVALONE) {
    *out << "col_stride_mmul_2x2_2x2( " <<
      GetInputName(0).str() << ", " <<
      InputDataType(0).m_rowStrideVar << ", " << 
      GetInputName(1).str() << ", " <<
      InputDataType(1).m_rowStrideVar << ", " <<
      GetInputName(2).str() << ", " <<
      InputDataType(2).m_rowStrideVar << " );\n";
  }
  else
    LOG_FAIL("replacement for throw call");
}

// Currently only handles the case where alpha = beta = 1.0
void MMul::PrintGeneralStride(IndStream &out)
{
  if (m_alpha.m_val == COEFVALONE && m_beta.m_val == COEFVALONE) {
    *out << "gen_stride_mmul_2x2_2x2( " <<
      GetInputName(0).str() << ", " <<
      InputDataType(0).m_rowStrideVar << ", " <<
      InputDataType(0).m_colStrideVar << ", " <<
      GetInputName(1).str() << ", " <<
      InputDataType(1).m_rowStrideVar << ", " <<
      InputDataType(1).m_colStrideVar << ", " <<
      GetInputName(2).str() << ", " <<
      InputDataType(2).m_rowStrideVar << ", " <<
      InputDataType(2).m_colStrideVar << ");\n";
  }
  else
    LOG_FAIL("replacement for throw call");
}

void MMul::Prop()
{
  if (!IsValidCost(m_cost)) {
    Gemm::Prop();

    /*    if (m_layer != LLDLAPRIMITIVELAYER) {
      cout << "ERROR: MMul appears in layer " <<  LayerNumToStr(m_layer) << "\n" ;
      LOG_FAIL("replacement for throw call");
      }
    
    if (*GetInputM(0) != m_regWidth || *GetInputN(0) != m_regWidth) 
      cout << "ERROR1: MMul only operates on m_regWidth by m_regWidth inputs\n";

    if (*GetInputM(1) != m_regWidth || *GetInputN(1) != m_regWidth) {
      GetInputM(1)->Print();
      cout << endl;
      GetInputN(1)->Print();
      cout << "ERROR2: MMul only operates on m_regWidth by m_regWidth inputs\n";
    }

    if (*GetInputM(2) != m_regWidth || *GetInputN(2) != m_regWidth) 
    cout << "ERROR3: MMul only operates on m_regWidth by m_regWidth inputs\n";*/

    switch(m_layer) {
    case (ABSLAYER):
      m_cost = 2 * GetInputM(0)->SumProds111(*GetInputN(1), *GetInputM(1));
      break;
    default:
      m_cost = 0;
    }
  }
}

Node* MMul::BlankInst()
{
  return new MMul(REAL_SINGLE, ABSLAYER);
}

NodeType MMul::GetType() const
{
  return "MMul" + LayerNumToStr(GetLayer());
}

MMulToPrim::MMulToPrim(Layer fromLayer, Layer toLayer, Type type)
{
  m_fromLayer = fromLayer;
  m_toLayer = toLayer;
  m_type = type;
  m_regWidth = arch->VecRegWidth(m_type);
}

string MMulToPrim::GetType() const
{
  return "MMulToPrim";
}

bool MMulToPrim::CanApply(const Node *node) const
{
  if (CurrPhase == LLDLAPRIMPHASE) {
    if (node->GetNodeClass() == Gemm::GetClass()) {
      const Gemm *gemm = (Gemm*)node;
      if (gemm->m_transA != NORMAL || gemm->m_transB != NORMAL)
	return false;
      if (gemm->GetLayer() != m_fromLayer)
	return false;
      
      if ((*(gemm->GetInputM(2)) <= m_regWidth) &&
	  (*(gemm->GetInputN(2)) <= m_regWidth) &&
	  (*(gemm->GetInputN(0)) <= m_regWidth))
	{
	  return true;
	}
      else
	return false;
    }
  }
  return false;
}

void MMulToPrim::Apply(Node *node) const
{
  Gemm *gemm = (Gemm*)node;

  NodeConn *connA, *connB, *connC;
  connA = gemm->m_inputs[0];
  connB = gemm->m_inputs[1];
  connC = gemm->m_inputs[2];
  
  MMul *prim = new MMul(gemm->m_type,
			m_toLayer);

  prim->AddInputs(6, 
		  connA->m_n, connA->m_num,
		  connB->m_n, connB->m_num,
		  connC->m_n, connC->m_num);

  node->m_poss->AddNode(prim);

  node->RedirectChildren(prim,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

#endif //DOLLDLA
