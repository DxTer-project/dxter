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

#if DOLLDLA
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
  string toLoadName2, toLoadPair;
  string loadStr = GetNameStr(0);
  // Decide which load instruction is needed based on
  // dimension and stride of input vector
  Stride inputRowStride = InputDataType(0).m_rowStride;
  Stride inputColStride = InputDataType(0).m_colStride;

  if (IsInputColVector(0)) {
    if (IsUnitStride(inputRowStride)) {
      *out << "VEC_PTR_PD_LOAD( " << loadStr << ", " << toLoadName << " );\n";
      return;
    } else {
      toLoadName2 = toLoadName + " + " + InputDataType(0).m_rowStrideVar;
      toLoadPair = toLoadName + ", " + toLoadName2;
      *out << "VEC_PPTR_PD_LOAD( " << loadStr << ", " << toLoadPair << " );\n";
      return;
    }
  } else if (IsInputRowVector(0)) {
    if (IsUnitStride(inputColStride)) {
      *out << "VEC_PTR_PD_LOAD( " << loadStr << ", " << toLoadName << " );\n";
      return;
    } else {
      toLoadName2 = toLoadName + " + " + InputDataType(0).m_colStrideVar;
      toLoadPair = toLoadName + ", " + toLoadName2;
      *out << "VEC_PPTR_LOAD( " << loadStr << ", " << toLoadPair << " );\n";
      return;
    }
  } else {
    cout << "ERROR: Input to vector register load is neither row nor column vector\n";
    throw;
  }
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
  string varDecl = "v2df_t " + GetInputNameStr(0)+ "_regs;\n";
  Var var(DirectVarDeclType, varDecl);
  set.insert(var);
}

void StoreFromRegs::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();

    // TODO: Check that the correct input # is a register
    
    /*    if (*(GetInputM(0)) != NUMREGSPERLOAD) {
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
    throw;*/
    
    m_cost = 0;
  }
}

void StoreFromRegs::PrintCode(IndStream &out)
{
  out.Indent();
  string regVarName = GetInputNameStr(0);
  string storeLocation1 = GetInputNameStr(1);
  string storeLocation2;
  // Decide which store instruction is needed based on
  // dimension and stride of input vector
  Stride inputRowStride = InputDataType(1).m_rowStride;
  Stride inputColStride = InputDataType(1).m_colStride;

  if (IsInputColVector(1)) {
    if (IsUnitStride(inputRowStride)) {
      *out << "VEC_PTR_PD_STORE( " << regVarName << ", " << storeLocation1 << " );\n";
      return;
    } else {
      storeLocation2 = storeLocation1 + " + " + InputDataType(1).m_rowStrideVar;
      *out << "VEC_PTR_PD_SET( 0, " << regVarName << ", " << storeLocation1 << " );\n";
      out.Indent();
      *out << "VEC_PTR_PD_SET( 1, " << regVarName << ", " << storeLocation2 << " );\n";
      return;
    }
  } else if (IsInputRowVector(1)) {
    if (IsUnitStride(inputColStride)) {
      *out << "VEC_PTR_PD_STORE( " << regVarName << ", " << storeLocation1 << " );\n";
      return;
    } else {
      storeLocation2 = storeLocation1 + " + " + InputDataType(1).m_colStrideVar;
      *out << "VEC_PTR_PD_SET( 0, " << regVarName << ", " << storeLocation1 << " );\n";
      out.Indent();
      *out << "VEC_PTR_PD_SET( 1, " << regVarName << ", " << storeLocation2 << " );\n";
      return;
    }
  } else {
    cout << "ERROR: Input to vector register store is neither row nor column vector\n";
    throw;
  }
  return;
}

void DuplicateRegLoad::Prop()
{
  if (!IsValidCost(m_cost)) {
    if (m_inputs.size() != 1)
      throw;
    if ((*(GetInputM(0)) != 1) ||
	(*(GetInputN(0)) != 1))
      throw;
    Input(0)->Prop();
    m_cost = 0;
  }
}

void DuplicateRegLoad::PrintCode(IndStream &out)
{
  out.Indent();
  *out << "VEC_PTR_DUP_LOAD( " << GetNameStr(0) << ", " << GetInputNameStr(0) << " );\n";
}

void DuplicateRegLoad::ClearDataTypeCache()
{
  m_mSizes.ClearSizes();
  m_nSizes.ClearSizes();
}

void DuplicateRegLoad::BuildDataTypeCache()
{
  if (m_mSizes.m_entries.empty()) {
    m_info = InputDataType(0);
    m_info.m_numRowsVar = "vector register size";
    unsigned int num = GetInputM(0)->NumSizes();
    m_mSizes.AddRepeatedSizes(NUMREGSPERLOAD, num, 1);
    m_nSizes.AddRepeatedSizes(1, num, 1);
  }
}

const Sizes* DuplicateRegLoad::GetM(ConnNum num) const
{
  if (num != 0)
    throw;
  return &m_mSizes;
}

const Sizes* DuplicateRegLoad::GetN(ConnNum num) const
{
  if (num != 0)
    throw;
  return &m_nSizes;
}

Name DuplicateRegLoad::GetName(ConnNum num) const
{
  if (num != 0)
    throw;
  Name name = GetInputName(0);
  name.m_name += "_regDup";
  return name;
}

void DuplicateRegLoad::AddVariables(VarSet &set) const
{
  string varDecl = "v2df_t " + GetInputNameStr(0)+ "_regDup;";
  Var var(DirectVarDeclType, varDecl);
  set.insert(var);
}


void TempVecReg::Prop()
{
  if (!IsValidCost(m_cost)) {
    if (m_inputs.size() != 1)
      throw;
    if ((*(GetInputM(0)) != 1) ||
	(*(GetInputN(0)) != 1))
      throw;
    Input(0)->Prop();
    m_cost = 0;
  }
}

void TempVecReg::PrintCode(IndStream &out)
{
  out.Indent();
  *out << "VEC_SET_ZERO( " << GetNameStr(0) << " );\n";
}

void TempVecReg::ClearDataTypeCache()
{
  m_mSizes.ClearSizes();
  m_nSizes.ClearSizes();
}

void TempVecReg::BuildDataTypeCache()
{
  if (m_mSizes.m_entries.empty()) {
    m_info = InputDataType(0);
    m_info.m_numRowsVar = "vector register size";
    unsigned int num = GetInputM(0)->NumSizes();
    m_mSizes.AddRepeatedSizes(NUMREGSPERLOAD, num, 1);
    m_nSizes.AddRepeatedSizes(1, num, 1);
  }
}

const Sizes* TempVecReg::GetM(ConnNum num) const
{
  if (num != 0)
    throw;
  return &m_mSizes;
}

const Sizes* TempVecReg::GetN(ConnNum num) const
{
  if (num != 0)
    throw;
  return &m_nSizes;
}

Name TempVecReg::GetName(ConnNum num) const
{
  if (num != 0)
    throw;
  Name name = GetInputName(0);
  name.m_name += "_regTemp";
  return name;
}

void TempVecReg::AddVariables(VarSet &set) const
{
  string varDecl = "v2df_t " + GetInputNameStr(0)+ "_regTemp;";
  Var var(DirectVarDeclType, varDecl);
  set.insert(var);
}


#endif //DOLLDLA
