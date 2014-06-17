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

#include "LLDLAGemmTransformations.h"
#include "loop.h"
#include "transpose.h"

#if DOLLDLA

LLDLAGemmLoopExp::LLDLAGemmLoopExp(Layer fromLayer, Layer toLayer, DimName dim, BSSize bsSize)
  : GemmLoopExp(fromLayer, toLayer, (dim==DIMM ? 0 : (dim==DIMK ? 1 : (dim==DIMN ? 2 : 5)))),
    m_bsSize(BSSizeToSize(bsSize))
{
}


string LLDLAGemmLoopExp::GetType() const
{
  return GemmLoopExp::GetType() + " for LLDLA";
}
  
bool LLDLAGemmLoopExp::CanApply(const Node *node) const
{
  if (!GemmLoopExp::CanApply(node))
    return false;
  const Gemm *gemm = (Gemm*)node;
  switch (m_dim) {
  case (0):
    {
      return !(*(gemm->GetInputM(2)) <= m_bsSize);
      //DIMM
      break;
    }
  case (1):
    {
      if (gemm->m_transA == NORMAL)
	return !(*(gemm->GetInputN(0)) <= m_bsSize);
      else if (gemm->m_transA != CONJ)
	return !(*(gemm->GetInputM(0)) <= m_bsSize);
      else
	throw;
      //DIMK
      break;
    }
  case (2):
    {
      return !(*(gemm->GetInputN(2)) <= m_bsSize);
      //DIMN
    }
  default:
    throw;
  }
}

bool GemmTransToNotTrans::CanApply(const Node *node) const
{
  const Gemm *gemm = (Gemm*)node;
  if (gemm->GetLayer() == m_layer) {
    return gemm->m_transA != NORMAL || gemm->m_transB != NORMAL;
  }
  return false;
}

void GemmTransToNotTrans::Apply(Node *node) const
{
  Gemm *gemm = (Gemm*)node;
  if (gemm->m_transA != NORMAL) {
    InsertTranspose(gemm->m_transA, gemm, 0, true);
    gemm->m_transA = NORMAL;
  }
  if (gemm->m_transB != NORMAL) {
    InsertTranspose(gemm->m_transB, gemm, 1, true);
    gemm->m_transB = NORMAL;
  }
}

bool LLDAGemmLowerLayer::CanApply(const Node *node) const
{
  if (node->GetNodeClass() == Gemm::GetClass()) {
    const Gemm *gemm = (Gemm*)node;
    if (gemm->GetLayer() != m_fromLayer)
      return false;
    if (*(gemm->GetInputM(0)) <= m_bs &&
	*(gemm->GetInputN(0)) <= m_bs &&
	*(gemm->GetInputN(1)) <= m_bs)
      return true;
    else
      return false;
  }
  throw;
}

void LLDAGemmLowerLayer::Apply(Node *node) const
{
  Gemm *gemm = (Gemm*)node;
  gemm->SetLayer(m_toLayer);
}

string LLDAGemmLowerLayer::GetType() const
{ 
  return "Gemm lower layer " + LayerNumToStr(m_fromLayer) 
  + " to " + LayerNumToStr(m_toLayer);
}


#endif
