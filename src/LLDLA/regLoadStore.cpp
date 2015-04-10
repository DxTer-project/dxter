/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2015, The University of Texas and Bryan Marker

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

#include "costModel.h"

void LoadToRegs::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLANode::Prop();
    if (m_inputs.size() != 1) {
      LOG_FAIL("replacement for throw call");
      throw;
    }

    if (*(GetInputM(0)) != GetVecRegWidth()) {
      // this isn't 1 x GetVecRegWidth()
      if (*(GetInputM(0)) != 1 || *(GetInputN(0)) != GetVecRegWidth()) {
	cout << "Error: Incorrect dimensions for register load\n";
	cout << "Input name: " << GetInputName(0).str() << endl;
	GetInputN(0)->Print();
	cout << "GetVecRegWidth() = " << std::to_string((long long int) GetVecRegWidth()) << endl;
	cout << "GetDataType() == REAL_DOUBLE ? " << std::to_string((long long int) (GetDataType() == REAL_DOUBLE)) << endl;
	cout << "*GetInputM(0) != 1 ? " << std::to_string((long long int) (*GetInputM(0) != 1)) << endl;
	cout << "*GetInputN(0) != GetVecRegWidth() ? " << std::to_string((long long int) (*GetInputN(0) != GetVecRegWidth())) << endl;
	LOG_FAIL("replacement for throw call");
      }
    } else if (*(GetInputN(0)) != 1) {
      GetInputN(0)->Print();
      // GetVecRegWidth() rows but not 1 column
      LOG_FAIL("replacement for throw call");
    }

    Input(0)->Prop();
    if (IsInputColVector(0)) {
      if (IsUnitStride(InputDataType(0).m_rowStride)) {
	m_cost = costModel->ContigVecLoadCost();
      } else {
	m_cost = GetVecRegWidth() * costModel->ContigVecLoadCost();
      }
    } else {
      if (IsUnitStride(InputDataType(0).m_colStride)) {
	m_cost = costModel->ContigVecLoadCost();
      } else {
	m_cost = GetVecRegWidth() * costModel->ContigVecLoadCost();
      }
    }
  }
}

void LoadToRegs::PrintCode(IndStream &out) {
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
    *out << arch->StridedLoad(GetDataType(), toLoadName, loadStr, strideVar);
  } else {
    *out << arch->ContiguousLoad(GetDataType(), toLoadName, loadStr);
  }

  return;
}

const Sizes* LoadToRegs::GetM(ConnNum num) const {
  if (num != 0) {
    LOG_FAIL("replacement for throw call");
  }
  return GetInputM(0);
}

const Sizes* LoadToRegs::GetN(ConnNum num) const {
  if (num != 0) {
    LOG_FAIL("replacement for throw call");
  }
  return GetInputN(0);
}

Name LoadToRegs::GetName(ConnNum num) const {
  if (num != 0) {
    LOG_FAIL("replacement for throw call");
  }
  Name name = GetInputName(0);
  name.m_name += "_regs";
  return name;
}

void LoadToRegs::AddVariables(VarSet &set) const {
  string varDecl = arch->TypeName(GetDataType()) + " " + GetInputNameStr(0)+ "_regs;\n";
  Var var(DirectVarDeclType, varDecl, GetDataType());
  set.insert(var);
}

void PackedLoadToRegs::Prop() {
  if (!IsValidCost(m_cost)) {
    if (!(IsInputRowVector(0) || IsInputColVector(0)) || !InputIsResidual(0)) {
      cout << "Input to PackedLoadToRegs is not a row or column vector" << endl;
      throw;
    }
    m_cost = GetVecRegWidth() * costModel->ContigVecLoadCost();
  }
}

int PackedLoadToRegs::ComputeResidual() {
  int residual;
  if (IsInputColVector(0)) {
    residual = GetInputNumRows(0);
  } else {
    residual = GetInputNumCols(0);
  }
  if (residual >= GetVecRegWidth()) {
    cout << "Error: PackedLoadToRegs::ComputeResidual gives bad residual" << endl;
  }
  return residual;
}

void PackedLoadToRegs::PrintCode(IndStream &out) {
  string toLoadName = GetInputNameStr(0);
  string loadStr = GetNameStr(0);

  Stride inputRowStride = InputDataType(0).m_rowStride;
  Stride inputColStride = InputDataType(0).m_colStride;
   
  string strideVar = "ERROR: STRIDE NOT DEFINED\n";

  if (IsInputColVector(0)) {
    if (IsUnitStride(inputRowStride)) {
      strideVar = "1";
    } else {
      strideVar = InputDataType(0).m_rowStrideVar;
    }
  } else {
    if (IsUnitStride(inputColStride)) {
      strideVar = "1";
    } else {
      strideVar = InputDataType(0).m_colStrideVar;
    }    
  }
  int residual = ComputeResidual();
  out.Indent();
  *out << arch->PackedLoad(GetDataType(), toLoadName, loadStr, strideVar, residual);
  return;
}

const Sizes* PackedLoadToRegs::GetM(ConnNum num) const {
  if (num != 0) {
    LOG_FAIL("replacement for throw call");
  }
  return GetInputM(0);
}

const Sizes* PackedLoadToRegs::GetN(ConnNum num) const {
  if (num != 0) {
    LOG_FAIL("replacement for throw call");
  }
  return GetInputN(0);
}

Name PackedLoadToRegs::GetName(ConnNum num) const {
  if (num != 0) {
    LOG_FAIL("replacement for throw call");
  }
  Name name = GetInputName(0);
  name.m_name += "_regs_packed";
  return name;
}

void PackedLoadToRegs::AddVariables(VarSet &set) const {
  string varDecl = arch->TypeName(GetDataType()) + " " + GetInputNameStr(0)+ "_regs_packed;\n";
  Var var(DirectVarDeclType, varDecl, GetDataType());
  set.insert(var);
}

void StoreFromRegs::Prop() {
  if (!IsValidCost(m_cost)) {
    DLAOp<2,1>::Prop();

    if (IsInputColVector(1)) {
      if (IsUnitStride(InputDataType(1).m_rowStride)) {
	m_cost = costModel->ContigVecStoreCost();
      } else {
	m_cost = GetVecRegWidth() * costModel->ContigVecStoreCost();
      }
    } else {
      if (IsUnitStride(InputDataType(1).m_colStride)) {
	m_cost = costModel->ContigVecStoreCost();
      } else {
	m_cost = GetVecRegWidth() * costModel->ContigVecStoreCost();
      }
    }
    
    m_cost = 0;
  }
}

void StoreFromRegs::PrintCode(IndStream &out) {
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

  out.Indent();
  if (isStridedLoad) {
    *out << arch->StridedStore(GetDataType(), storeLocation, regVarName, strideVar);
  } else {
    *out << arch->ContiguousStore(GetDataType(), storeLocation, regVarName);
  }

  return;
}

void StoreFromRegs::StoreNonContigLocations(IndStream &out, string regVarName, string storePtr, string strideVar)
{
  out.Indent();
  *out << "VEC_PTR_PD_SET( 0, " + regVarName + ", " + storePtr + " );\n";
  for (int i = 1; i < GetVecRegWidth(); i++) {
    out.Indent();
    string storePtrStr = storePtr + " + " + std::to_string((long long int) i) + " * " + strideVar;
    *out << "VEC_PTR_PD_SET( " + std::to_string((long long int) i) + ", " + regVarName + ", " + storePtrStr + " );\n";
  }
}

void UnpackStoreFromRegs::Prop() {
  if (!IsValidCost(m_cost)) {
    if (!(IsInputRowVector(1) || IsInputColVector(1)) || !InputIsResidual(1)) {
      m_cost = GetVecRegWidth() * costModel->ContigVecStoreCost();
    }
  }
}

int UnpackStoreFromRegs::ComputeResidual() {
  int residual;
  if (IsInputColVector(1)) {
    residual = GetInputNumRows(1);
  } else {
    residual = GetInputNumCols(1);
  }
  if (residual >= GetVecRegWidth()) {
    cout << "Error: PackedLoadToRegs::ComputeResidual gives bad residual" << endl;
    throw;
  }
  return residual;
}

void UnpackStoreFromRegs::PrintCode(IndStream &out)
{
  string regVarName = GetInputNameStr(0);
  string storeLocation = GetInputNameStr(1);

  Stride inputRowStride = InputDataType(1).m_rowStride;
  Stride inputColStride = InputDataType(1).m_colStride;

  string strideVar = "ERROR: STRIDE NOT DEFINED\n";
  //  bool isStridedLoad;

  if (IsInputColVector(1)) {
    if (IsUnitStride(inputRowStride)) {
      strideVar = "1";
    } else {
      strideVar = InputDataType(1).m_rowStrideVar;
    }
  } else {
    if (IsUnitStride(inputColStride)) {
      strideVar = "1";
    } else {
      strideVar = InputDataType(0).m_colStrideVar;
    }    
  }

  int residual = ComputeResidual();

  out.Indent();
  *out << arch->UnpackStore(GetDataType(), storeLocation, regVarName, strideVar, residual);
  return;
}

void DuplicateRegLoad::Prop() {
  if (!IsValidCost(m_cost)) {
    if (m_inputs.size() != 1) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
    if ((*(GetInputM(0)) != 1) ||
	(*(GetInputN(0)) != 1)) {
      cout << "Bad duplicate load arg, dimensions are: " << endl;
      GetInputM(0)->Print();
      GetInputN(0)->Print();
      LOG_FAIL("replacement for throw call");
      throw;
    }
    Input(0)->Prop();
    m_cost = 0;
  }
}

void DuplicateRegLoad::PrintCode(IndStream &out)
{
  out.Indent();
  *out << arch->DuplicateLoad(GetDataType(), GetInputNameStr(0), GetNameStr(0));
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
    m_mSizes.AddRepeatedSizes(GetVecRegWidth(), num);
    m_nSizes.AddRepeatedSizes(1, num);
  }
}

const Sizes* DuplicateRegLoad::GetM(ConnNum num) const
{
  if (num != 0) {
    LOG_FAIL("replacement for throw call");
  }
  return &m_mSizes;
}

const Sizes* DuplicateRegLoad::GetN(ConnNum num) const
{
  if (num != 0) {
    LOG_FAIL("replacement for throw call");
  }
  return &m_nSizes;
}

Name DuplicateRegLoad::GetName(ConnNum num) const
{
  if (num != 0) {
    LOG_FAIL("replacement for throw call");
  }
  Name name = GetInputName(0);
  name.m_name += "_regDup";
  return name;
}

void DuplicateRegLoad::AddVariables(VarSet &set) const {
  string varDecl = arch->TypeName(GetDataType()) + " " + GetInputNameStr(0)+ "_regDup;";
  Var var(DirectVarDeclType, varDecl, GetDataType());
  set.insert(var);
}

void TempVecReg::Prop() {
  if (!IsValidCost(m_cost)) {
    if (m_inputs.size() != 1) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
    if ((*(GetInputM(0)) != 1) ||
	(*(GetInputN(0)) != 1)) {
      cout << "TempVecReg bad Prop" << endl;
      cout << Input(0)->GetNodeClass() << endl;
      GetInputM(0)->Print();
      GetInputN(0)->Print();
      LOG_FAIL("replacement for throw call");
      throw;
    }
    Input(0)->Prop();
    m_cost = 0;
  }
}

void TempVecReg::PrintCode(IndStream &out) {
  out.Indent();
  *out << arch->ZeroVar(GetDataType(), GetNameStr(0));
  return;
}

void TempVecReg::ClearDataTypeCache() {
  m_mSizes.ClearSizes();
  m_nSizes.ClearSizes();
}

void TempVecReg::BuildDataTypeCache() {
  if (m_mSizes.m_entries.empty()) {
    m_info = InputDataType(0);
    m_info.m_numRowsVar = "vector register size";
    m_info.m_numColsVar = "vector register size";
    unsigned int num = GetInputM(0)->NumSizes();
    m_mSizes.AddRepeatedSizes(GetVecRegWidth(), num);
    m_nSizes.AddRepeatedSizes(1, num);
  }
}

const Sizes* TempVecReg::GetM(ConnNum num) const
{
  if (num != 0) {
    LOG_FAIL("replacement for throw call");
  }
  return &m_mSizes;
}

const Sizes* TempVecReg::GetN(ConnNum num) const
{
  if (num != 0) {
    LOG_FAIL("replacement for throw call");
  }
  return &m_nSizes;
}

Name TempVecReg::GetName(ConnNum num) const
{
  if (num != 0) {
    LOG_FAIL("replacement for throw call");
  }
  Name name = GetInputName(0);
  name.m_name += "_regTemp";
  return name;
}

void TempVecReg::AddVariables(VarSet &set) const
{
  string varDecl = arch->TypeName(GetDataType()) + " " +  GetInputNameStr(0)+ "_regTemp;";
  Var var(DirectVarDeclType, varDecl, GetDataType());
  set.insert(var);
}


#endif //DOLLDLA
