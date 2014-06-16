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

#if DOLLDLA

LLDLAGemmLoopExp::LLDLAGemmLoopExp(Layer fromLayer, Layer toLayer, DimName dim, BSSize bsSize)
  : GemmLoopExp(fromLayer, toLayer, (dim==DIMM ? 0 : (dim==DIMN ? 2 : 5))),
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
  if (gemm->m_transA != NORMAL ||
      gemm->m_transB != NORMAL)
    return false;
  switch (m_dim) {
  case (0):
    {
      return !(*(gemm->GetInputM(0)) <= m_bsSize);
      //DIMM
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


#endif
