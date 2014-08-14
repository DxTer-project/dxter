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

    if (*(GetInputM(0)) != LLDLA_MU) {
      // this isn't 1 x LLDLA_MU
      cout << "*GetInputM(0) = " << endl;
      GetInputM(0)->Print();
      if (*(GetInputM(0)) != 1 || *(GetInputN(0)) != LLDLA_MU) {
	cout << "Error: Incorrect dimensions for register load\n";
	cout << "*GetInputM(0) != 1 ? " << std::to_string(*GetInputM(0) != 1) << endl;
	cout << "*GetInputN(0) != LLDLA_MU ? " << std::to_string(*GetInputN(0) != LLDLA_MU) << endl;
	throw;
      }
    } else if (*(GetInputN(0)) != 1) {
      GetInputN(0)->Print();
      // LLDLA_MU rows but not 1 column
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
  string toLoad;
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
#if USE_DOUBLE_PRECISION
      toLoad = toLoadName;
      *out << "tmp[0] = *" << toLoadName << endl;
      for (int i = 1; i < LLDLA_MU; i++) {
	out.Indent();
	toLoad += ", " + toLoadName + " + " + std::to_string((long long int) i) + " * " + InputDataType(0).m_rowStrideVar;
      }
      out.Indent();
      *out << "VEC_PPTR_PD_LOAD( " << loadStr << ", " << toLoad << " );\n";
#else
      toLoad = toLoadName;
      *out << "tmp[0] = *" << toLoadName << ";\n";
      for (int i = 1; i < LLDLA_MU; i++) {
	string valToLoad = toLoadName + " + " + std::to_string((long long int) i) + " * " + InputDataType(0).m_rowStrideVar;
	out.Indent();
	*out << "tmp[ " << std::to_string((long long int) i) << " ] = *(" << valToLoad << ");\n";
	toLoad += ", " + valToLoad;
      }
      out.Indent();
      *out << "VEC_PPTR_PD_LOAD( " << loadStr << ", " << toLoad << " );\n";
#endif // USE_DOUBLE PRECISION
      return;
    }
  } else if (IsInputRowVector(0)) {
    if (IsUnitStride(inputColStride)) {
      *out << "VEC_PTR_PD_LOAD( " << loadStr << ", " << toLoadName << " );\n";
      return;
    } else {
#if USE_DOUBLE_PRECISION
      toLoad = toLoadName;
      for (int i = 1; i < LLDLA_MU; i++) {
	toLoad += ", " + toLoadName + " + " + std::to_string((long long int) i) + " * " + InputDataType(0).m_colStrideVar;
      }

      *out << "VEC_PPTR_PD_LOAD( " << loadStr << ", " << toLoad << " );\n";
#else
      toLoad = toLoadName;
      *out << "tmp[0] = *" << toLoadName << ";\n";
      for (int i = 1; i < LLDLA_MU; i++) {
      string valToLoad = toLoadName + " + " + std::to_string((long long int) i) + " * " + InputDataType(0).m_colStrideVar;
	*out << "tmp[ " << std::to_string((long long int) i) << " ] = *(" << valToLoad << ");\n";
	toLoad += ", " + toLoadName + " + " + std::to_string((long long int) i) + " * " + InputDataType(0).m_colStrideVar;
      }

      *out << "VEC_PPTR_PD_LOAD( " << loadStr << ", " << toLoad << " );\n";
      
#endif // USE_DOUBLE_PRECISION
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
  string varDecl = "vec_reg " + GetInputNameStr(0)+ "_regs;\n";
  Var var(DirectVarDeclType, varDecl);
  set.insert(var);
}

void StoreFromRegs::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();

    // TODO: Check that the correct input # is a register
    
    /*    if (*(GetInputM(0)) != LLDLA_MU) {
      // this isn't 1 x LLDLA_MU
      if (*(GetInputM(0)) != 1 || *(GetInputN(0)) != LLDLA_MU)
	throw;
    }
    else if (*(GetInputN(0)) != 1) {
      // LLDLA_MU rows but not 1 column
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
  string regVarName = GetInputNameStr(0);
  string storeLocation = GetInputNameStr(1);
  // Decide which store instruction is needed based on
  // dimension and stride of input vector
  Stride inputRowStride = InputDataType(1).m_rowStride;
  Stride inputColStride = InputDataType(1).m_colStride;

  if (IsInputColVector(1)) {
    if (IsUnitStride(inputRowStride)) {
      out.Indent();
      *out << "VEC_PTR_PD_STORE( " << regVarName << ", " << storeLocation << " );\n";
      return;
    } else {
      StoreNonContigLocations(out, regVarName, storeLocation, InputDataType(1).m_rowStrideVar);
      return;
    }
  } else if (IsInputRowVector(1)) {
    if (IsUnitStride(inputColStride)) {
      out.Indent();
      *out << "VEC_PTR_PD_STORE( " << regVarName << ", " << storeLocation << " );\n";
      return;
    } else {
      StoreNonContigLocations(out, regVarName, storeLocation, InputDataType(1).m_colStrideVar);
      return;
    }
  } else {
    cout << "ERROR: Input to vector register store is neither row nor column vector\n";
    throw;
  }
  return;
}

void StoreFromRegs::StoreNonContigLocations(IndStream &out, string regVarName, string storePtr, string strideVar)
{
  out.Indent();
  *out << "VEC_PTR_PD_SET( 0, " + regVarName + ", " + storePtr + " );\n";
  for (int i = 1; i < LLDLA_MU; i++) {
    out.Indent();
    string storePtrStr = storePtr + " + " + std::to_string((long long int) i) + " * " + strideVar;
    *out << "VEC_PTR_PD_SET( " + std::to_string((long long int) i) + ", " + regVarName + ", " + storePtrStr + " );\n";
  }
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
    m_info.m_numColsVar = "1";
    unsigned int num = GetInputM(0)->NumSizes();
    m_mSizes.AddRepeatedSizes(LLDLA_MU, num, 1);
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
  string varDecl = "vec_reg " + GetInputNameStr(0)+ "_regDup;";
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
    m_info.m_numColsVar = "vector register size";
    unsigned int num = GetInputM(0)->NumSizes();
    m_mSizes.AddRepeatedSizes(LLDLA_MU, num, 1);
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
  string varDecl = "vec_reg " + GetInputNameStr(0)+ "_regTemp;";
  Var var(DirectVarDeclType, varDecl);
  set.insert(var);
}

#endif //DOLLDLA
