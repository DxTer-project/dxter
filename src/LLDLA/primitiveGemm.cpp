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
#include "layers.h"
#include "primitiveGemm.h"
#include "loopSupport.h"

#if DOLLDLA

PrimitiveGemm::PrimitiveGemm(Coef alpha, Coef beta, Type type, Layer layer)
  : Gemm(layer, NORMAL, NORMAL, alpha, beta, type)
{
}

void PrimitiveGemm::PrintCode(IndStream &out)
{
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

void PrimitiveGemm::PrintRowStride(IndStream &out)
{
  *out << "ROW STRIDE not yet implemented\n";
}

void PrimitiveGemm::PrintColStride(IndStream &out)
{
  *out << "COL STRIDE not yet implemented\n";
}

// Currently only handles the case where alpha = beta = 1.0
void PrimitiveGemm::PrintGeneralStride(IndStream &out)
{
  if (m_alpha.m_val == COEFVALONE && m_beta.m_val == COEFVALONE) {
    *out << "gen_stride_mmul_2x2_2x2( " <<
      GetInputName(0).str() << ", " <<
      GetInputName(1).str() << ", " <<
      GetInputName(2).str() << ");\n";
  }
  else
    throw;
}

void PrimitiveGemm::Prop()
{
  if (!IsValidCost(m_cost)) {
    Gemm::Prop();

    if (m_layer != LLDLAPRIMITIVELAYER) {
      cout << "ERROR: PrimitiveGemm appears in layer " <<  LayerNumToStr(m_layer) << "\n" ;
      throw;
    }
    
    if (*GetInputM(0) != LLDLA_MU || *GetInputN(0) != LLDLA_MU) 
      cout << "ERROR1: PrimitiveGemm only operates on LLDLA_MU by LLDLA_MU inputs\n";

    if (*GetInputM(1) != LLDLA_MU || *GetInputN(1) != LLDLA_MU) {
      GetInputM(1)->Print();
      cout << endl;
      GetInputN(1)->Print();
      cout << "ERROR2: PrimitiveGemm only operates on LLDLA_MU by LLDLA_MU inputs\n";
    }

    if (*GetInputM(2) != LLDLA_MU || *GetInputN(2) != LLDLA_MU) 
      cout << "ERROR3: PrimitiveGemm only operates on LLDLA_MU by LLDLA_MU inputs\n";
    
    m_cost = ZERO;
  }
}

Node* PrimitiveGemm::BlankInst()
{
  return new PrimitiveGemm(COEFONE, COEFONE, REAL, ABSLAYER);
}

NodeType PrimitiveGemm::GetType() const
{
  return "PrimitiveGemm" + LayerNumToStr(GetLayer());
}

string LLDLAGemmToPrim::GetType() const
{
  return "LLDLAGemmToPrim";
}

bool LLDLAGemmToPrim::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Gemm::GetClass()) {
    const Gemm *gemm = (Gemm*)node;
    if (gemm->m_transA != NORMAL || gemm->m_transB != NORMAL)
      return false;
    if (gemm->GetLayer() != m_fromLayer)
      return false;

    if ((*(gemm->GetInputM(2)) <= LLDLA_MU) &&
	(*(gemm->GetInputN(2)) <= LLDLA_MU) &&
	(*(gemm->GetInputN(0)) <= LLDLA_MU))
      {
	return true;
      }
    else
      return false;
  }
  return false;
}

void LLDLAGemmToPrim::Apply(Node *node) const
{
  Gemm *gemm = (Gemm*)node;

  NodeConn *connA, *connB, *connC;
  connA = gemm->m_inputs[0];
  connB = gemm->m_inputs[1];
  connC = gemm->m_inputs[2];
  
  PrimitiveGemm *prim = new PrimitiveGemm(gemm->m_alpha,
					  gemm->m_beta,
					  gemm->m_type,
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
