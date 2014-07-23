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

#include "regArith.h"

#if DOLLDLA

void FMAdd::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3, 1>::Prop();

    if (*GetInputM(0) != *GetInputM(1) || *GetInputM(0) != *GetInputM(2)) {
      throw;
    }

    if (*GetInputN(0) != *GetInputN(1) || *GetInputN(0) != *GetInputN(2)) {
      throw;
    }

    m_cost = 0;
  }
}

void FMAdd::PrintCode(IndStream &out)
{
  out.Indent();
  string aStr = GetInputNameStr(0);
  string bStr = GetInputNameStr(1);
  string cStr = GetInputNameStr(2);
  *out << "VEC_PD_FMA( " << aStr << ", " << bStr << ", " << cStr << " );\n";
  return;
}

void Add::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2, 1>::Prop();

    if (*GetInputM(0) != *GetInputM(1)) {
      throw;
    }

    if (*GetInputN(0) != *GetInputN(1)) {
      throw;
    }

    m_cost = 0;
  }
}

void Add::PrintCode(IndStream &out)
{
  out.Indent();
  string aStr = GetInputNameStr(0);
  string bStr = GetInputNameStr(1);
  *out << "VEC_PD_ADD( " << aStr << ", " << bStr << " );\n";
  return;
}

void Mul::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2, 1>::Prop();

    if (*GetInputM(0) != *GetInputM(1)) {
      throw;
    }

    if (*GetInputN(0) != *GetInputN(1)) {
      throw;
    }

    m_cost = 0;
  }
}

void Mul::PrintCode(IndStream &out)
{
  out.Indent();
  string aStr = GetInputNameStr(0);
  string bStr = GetInputNameStr(1);
  *out << "VEC_PD_MUL( " << aStr << ", " << bStr << " );\n";
  return;
}

#endif // DOLLDLA
