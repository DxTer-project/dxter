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

#include "avx.h"

#if DOLLDLA

#include "eliminateMaskedStoreLoad.h"
#include "maskedStore.h"
#include "packedLoadToMaskedLoad.h"
#include "unpackStoreToMaskedStore.h"
#include "regLoadStore.h"

AVX::AVX() {
  m_name = "AVX";
  auto maskedLoad = new pair<string, SingleTrans*>(PackedLoadToRegs::GetClass(), new PackedLoadToMaskedLoad());
  m_archTrans.push_back(*maskedLoad);
  auto maskedStore = new pair<string, SingleTrans*>(UnpackStoreFromRegs::GetClass(), new UnpackStoreToMaskedStore());
  m_archTrans.push_back(*maskedStore);
  auto maskedSLE = new pair<string, SingleTrans*>(MaskedStore::GetClass(), new EliminateMaskedStoreLoad());
  m_archTrans.push_back(*maskedSLE);
}

string AVX::MaskRegisterDeclaration(Type dataType, string varName, unsigned int residualSize) {
  string varDecl = "__m256i " + varName + " = ";
  int numEnts = arch->VecRegWidth(dataType);
  string opName;

  if (dataType == REAL_SINGLE) {
    opName = "_mm256_set_epi32(";
  } else if (dataType == REAL_DOUBLE) {
    opName = "_mm256_set_epi64x(";
  } else {
    throw;
  }
  varDecl += opName;
  for (int i = numEnts; i > 0; i--) {
    if (i <= residualSize) {
      varDecl += "-1";
    } else {
      varDecl += "0";
    }
    if (i > 1) {
      varDecl += ", ";
    }
  }
  varDecl += ");";
  return varDecl;
}

string AVX::GlobalDeclarations() {  
  return "";
}

string AVX::SetupFunc() {
  return "void AVX_setup() {}";
}

#endif // DOLLDLA
