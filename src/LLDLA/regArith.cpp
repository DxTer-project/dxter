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

FMAdd::FMAdd(Type type)
{
  m_type = type;
}

void FMAdd::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<3, 1>::Prop();

// TODO: Delete shape checks and instead check that operands
// are registers

/*    if (*GetInputM(0) != *GetInputM(1) || *GetInputM(0) != *GetInputM(2)) {
      throw;
    }

    if (*GetInputN(0) != *GetInputN(1) || *GetInputN(0) != *GetInputN(2)) {
      throw;
      }*/

    m_cost = 0;
  }
}

void FMAdd::PrintCode(IndStream &out)
{
  out.Indent();
  string aStr = GetInputNameStr(0);
  string bStr = GetInputNameStr(1);
  string cStr = GetInputNameStr(2);
  *out << arch->FMACode(m_type, aStr, bStr, cStr, cStr);
  return;
}

Add::Add(Type type)
{
  m_type = type;
}

void Add::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2, 1>::Prop();

// TODO: Delete shape checks and instead check that
// inputs are registers

     m_cost = 0;
  }
}

void Add::PrintCode(IndStream &out)
{
  out.Indent();
  string aStr = GetInputNameStr(0);
  string bStr = GetInputNameStr(1);
  *out << arch->AddCode(m_type, aStr, bStr, bStr);
  return;
}

Mul::Mul(Type type)
{
  m_type = type;
}

void Mul::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2, 1>::Prop();

// TODO: Delete shape checks and instead check that
// inputs are registers

    m_cost = 0;
  }
}

void Mul::PrintCode(IndStream &out)
{
  out.Indent();
  string aStr = GetInputNameStr(0);
  string bStr = GetInputNameStr(1);
  *out << arch->MulCode(m_type, aStr, bStr, bStr);
  return;
}

ZeroReg::ZeroReg(Type type)
{
  m_type = type;
}

void ZeroReg::Prop()
{
// TODO: Add full prop method with register checks
  return;
}

void ZeroReg::PrintCode(IndStream &out)
{
  out.Indent();
  string aStr = GetInputNameStr(0);
  *out << arch->ZeroVar(m_type, aStr);
  return;
}

AccumReg::AccumReg(Type type)
{
  m_type = type;
}

void AccumReg::Prop()
{
// TODO: Add full prop method with register checks
  return;
}

void AccumReg::PrintCode(IndStream &out)
{
  out.Indent();
  string vecStr = GetInputNameStr(0);
  string accStr = GetInputNameStr(1);
  *out << arch->AccumCode(m_type, accStr, vecStr);
  return;
}

#endif // DOLLDLA