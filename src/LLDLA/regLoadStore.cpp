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

LoadToRegs::LoadToRegs(Type type)
{
  m_type = type;
  m_regWidth = arch->VecRegWidth(type);
}

void LoadToRegs::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLANode::Prop();
    if (m_inputs.size() != 1)
      throw;

    if (*(GetInputM(0)) != m_regWidth) {
      // this isn't 1 x m_regWidth
      if (*(GetInputM(0)) != 1 || *(GetInputN(0)) != m_regWidth) {
	cout << "Error: Incorrect dimensions for register load\n";
	GetInputN(0)->Print();
	cout << "m_regWidth = " << std::to_string((long long int) m_regWidth) << endl;
	cout << "*GetInputM(0) != 1 ? " << std::to_string((long long int) (*GetInputM(0) != 1)) << endl;
	cout << "*GetInputN(0) != m_regWidth ? " << std::to_string((long long int) (*GetInputN(0) != m_regWidth)) << endl;
	throw;
      }
    } else if (*(GetInputN(0)) != 1) {
      GetInputN(0)->Print();
      // m_regWidth rows but not 1 column
      throw;
    }

    Input(0)->Prop();
    if (IsInputColVector(0)) {
      if (IsUnitStride(InputDataType(0).m_rowStride)) {
	m_cost = CONTIG_VECTOR_LOAD_COST;
      } else {
	m_cost = m_regWidth * CONTIG_VECTOR_LOAD_COST;
      }
    } else {
      if (IsUnitStride(InputDataType(0).m_colStride)) {
	m_cost = CONTIG_VECTOR_LOAD_COST;
      } else {
	m_cost = m_regWidth * CONTIG_VECTOR_LOAD_COST;
      }
    }
  }
}

void LoadToRegs::PrintCode(IndStream &out)
{
  out.Indent();
  string toLoadName = GetInputNameStr(0);
  string loadStr = GetNameStr(0);
  // Decide which load instruction is needed based on
  // dimension and stride of input vector
  Stride inputRowStride = InputDataType(0).m_rowStride;
  Stride inputColStride = InputDataType(0).m_colStride;
  
  string strideVar = "ERROR: STRIDE NOT DEFINED\n";
  bool isStridedLoad;

  if (IsInputColVector(0)) {
    if (IsUnitStride(inputRowStride)) {
      isStridedLoad = false;
    } else {
      isStridedLoad = true;
      strideVar = InputDataType(0).m_rowStrideVar;
    }
  } else {
    if (IsUnitStride(inputColStride)) {
      isStridedLoad = false;
    } else {
      isStridedLoad = true;
      strideVar = InputDataType(0).m_colStrideVar;
    }    
  }

  if (isStridedLoad) {
    *out << arch->StridedLoad(m_type, toLoadName, loadStr, strideVar);
  } else {
    *out << arch->ContiguousLoad(m_type, toLoadName, loadStr);
  }

  /*  if (IsInputColVector(0)) {
    if (IsUnitStride(inputRowStride)) {
      *out << "VEC_PTR_PD_LOAD( " << loadStr << ", " << toLoadName << " );\n";
      return;
    } else {
#if USE_DOUBLE_PRECISION
      toLoad = toLoadName;
      *out << "tmp[0] = *" << toLoadName << ";\n";
      for (int i = 1; i < m_regWidth; i++) {
	out.Indent();
	toLoad += ", " + toLoadName + " + " + std::to_string((long long int) i) + " * " + InputDataType(0).m_rowStrideVar;
      }
      *out << "VEC_PPTR_PD_LOAD( " << loadStr << ", " << toLoad << " );\n";
#else
      toLoad = toLoadName;
      *out << "tmp[0] = *" << toLoadName << ";\n";
      for (int i = 1; i < m_regWidth; i++) {
	string valToLoad = toLoadName + " + " + std::to_string((long long int) i) + " * " + InputDataType(0).m_rowStrideVar;
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
      for (int i = 1; i < m_regWidth; i++) {
	toLoad += ", " + toLoadName + " + " + std::to_string((long long int) i) + " * " + InputDataType(0).m_colStrideVar;
      }

      *out << "VEC_PPTR_PD_LOAD( " << loadStr << ", " << toLoad << " );\n";
#else
      toLoad = toLoadName;
      *out << "tmp[0] = *" << toLoadName << ";\n";
      for (int i = 1; i < m_regWidth; i++) {
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
    }*/
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

StoreFromRegs::StoreFromRegs(Type type)
{
  m_type = type;
  m_regWidth = arch->VecRegWidth(m_type);
}

void StoreFromRegs::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();

    if (IsInputColVector(1)) {
      if (IsUnitStride(InputDataType(1).m_rowStride)) {
	m_cost = CONTIG_VECTOR_STORE_COST;
      } else {
	m_cost = m_regWidth * CONTIG_VECTOR_STORE_COST;
      }
    } else {
      if (IsUnitStride(InputDataType(1).m_colStride)) {
	m_cost = CONTIG_VECTOR_STORE_COST;
      } else {
	m_cost = m_regWidth * CONTIG_VECTOR_STORE_COST;
      }
    }
    
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

  string strideVar = "ERROR: STRIDE NOT DEFINED\n";
  bool isStridedLoad;

  if (IsInputColVector(1)) {
    if (IsUnitStride(inputRowStride)) {
      isStridedLoad = false;
    } else {
      isStridedLoad = true;
      strideVar = InputDataType(1).m_rowStrideVar;
    }
  } else {
    if (IsUnitStride(inputColStride)) {
      isStridedLoad = false;
    } else {
      isStridedLoad = true;
      strideVar = InputDataType(0).m_colStrideVar;
    }    
  }

  if (isStridedLoad) {
    *out << arch->StridedStore(m_type, storeLocation, regVarName, strideVar);
  } else {
    *out << arch->ContiguousStore(m_type, storeLocation, regVarName);
  }


  /*  if (IsInputColVector(1)) {
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
    }*/
  return;
}

void StoreFromRegs::StoreNonContigLocations(IndStream &out, string regVarName, string storePtr, string strideVar)
{
  out.Indent();
  *out << "VEC_PTR_PD_SET( 0, " + regVarName + ", " + storePtr + " );\n";
  for (int i = 1; i < m_regWidth; i++) {
    out.Indent();
    string storePtrStr = storePtr + " + " + std::to_string((long long int) i) + " * " + strideVar;
    *out << "VEC_PTR_PD_SET( " + std::to_string((long long int) i) + ", " + regVarName + ", " + storePtrStr + " );\n";
  }
}

DuplicateRegLoad::DuplicateRegLoad(Type type)
{
  m_type = type;
  m_regWidth = arch->VecRegWidth(m_type);
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
  *out << arch->DuplicateLoad(m_type, GetInputNameStr(0), GetNameStr(0));
  return;
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
    m_mSizes.AddRepeatedSizes(m_regWidth, num, 1);
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

TempVecReg::TempVecReg(Type type)
{
  m_type = type;
  m_regWidth = arch->VecRegWidth(type);
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
  *out << arch->ZeroVar(m_type, GetNameStr(0));
  return;
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
    m_mSizes.AddRepeatedSizes(m_regWidth, num, 1);
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
