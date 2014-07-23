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

#include "regLoadStore.h"

void LoadToRegs::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLANode::Prop();
    if (m_inputs.size() != 1)
      throw;
    
    if (*(GetInputM(0)) != NUMREGSPERLOAD) {
      // this isn't 1 x NUMREGSPERLOAD
      if (*(GetInputM(0)) != 1 || *(GetInputN(0)) != NUMREGSPERLOAD)
	throw;
    }
    else if (*(GetInputN(0)) != 1) {
      // NUMREGSPERLOAD rows but not 1 column
      throw;
    }

    Input(0)->Prop();
    
    m_cost = 0;
  }
}

void LoadToRegs::PrintCode(IndStream &out)
{
  out.Indent();
  string toLoadName = GetInputNameStr(0);
  string loadStr = GetNameStr(0);
  *out << loadStr << " = VEC_PD_LOAD( " << toLoadName << " );\n";
  return;
}

const Sizes* LoadToRegs::GetM(ConnNum num) const
{
  if (num != 0)
    throw;
  return GetInputM(0);
}

const Sizes* LoadToRegs::GetN(ConnNum num) const
{
  if (num != 0)
    throw;
  return GetInputN(0);
}

Name LoadToRegs::GetName(ConnNum num) const
{
  if (num != 0)
    throw;
  Name name = GetInputName(0);
  name.m_name += "_regs";
  return name;
}

void LoadToRegs::AddVariables(VarSet &set) const
{
  string varDecl = "v2df_t " + GetInputNameStr(0)+ "_regs";
  Var var(DirectVarDeclType, varDecl);
  set.insert(var);
}

void StoreFromRegs::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();
    
    if (*(GetInputM(0)) != NUMREGSPERLOAD) {
      // this isn't 1 x NUMREGSPERLOAD
      if (*(GetInputM(0)) != 1 || *(GetInputN(0)) != NUMREGSPERLOAD)
	throw;
    }
    else if (*(GetInputN(0)) != 1) {
      // NUMREGSPERLOAD rows but not 1 column
      throw;
    }

    if (*GetInputM(0) != *GetInputM(1))
      throw;

    if (*GetInputN(0) != *GetInputN(1))
      throw;
    
    m_cost = 0;
  }
}

void StoreFromRegs::PrintCode(IndStream &out)
{
  out.Indent();
  string regVarName = GetInputNameStr(0);
  string memoryVarName = GetInputNameStr(1);
  *out << "VEC_PD_STORE( " << memoryVarName << ", " << regVarName << " );\n";
  return;
}
