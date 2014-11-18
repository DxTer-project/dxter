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

#include "scalarArith.h"

#if DOLLDLA

void AddScalars::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2, 1>::Prop();
    m_cost = arch->ContigVecLoadCost() + arch->ContigVecStoreCost();
  }
}

void AddScalars::PrintCode(IndStream &out)
{
  out.Indent();
  string aStr = GetInputNameStr(0);
  string bStr = GetInputNameStr(1);
  *out << bStr << " = " << bStr << " + " << aStr << ";\n";
  return;
}

void MulScalars::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2, 1>::Prop();
    m_cost = arch->ContigVecLoadCost() + arch->ContigVecStoreCost();
  }
}

void MulScalars::PrintCode(IndStream &out)
{
  out.Indent();
  string aStr = GetInputNameStr(0);
  string bStr = GetInputNameStr(1);
  *out << bStr << " = " << bStr << " * " << aStr << ";\n";
  return;
}

void SetScalarToZero::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<1, 1>::Prop();
    m_cost = arch->ContigVecStoreCost();
  }
  return;
}

void SetScalarToZero::PrintCode(IndStream &out)
{
  out.Indent();
  string aStr = GetInputNameStr(0);
  *out << aStr << " = 0.0;\n";
  return;
}

#endif // DOLLDLA
