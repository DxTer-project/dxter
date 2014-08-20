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

#pragma once

#include "layers.h"
#include "base.h"

#if DOLLDLA

#include <string>

using namespace std;

enum Stride { UNITSTRIDE,
	      NONUNITSTRIDE,
	      BADSTRIDE };

enum VecType { ROWVECTOR,
	       COLVECTOR };

inline bool IsUnitStride(const Stride &stride)
{
  switch (stride)
    {
    case (UNITSTRIDE):
      return true;
    case (NONUNITSTRIDE):
      return false;
    case (BADSTRIDE):
    default:
      throw;
    }
}

class DataTypeInfo
{
 public:
  Stride m_rowStride;
  Stride m_colStride;
  string m_numRowsVar;
  string m_numColsVar;
  string m_rowStrideVar;
  string m_colStrideVar;

  DataTypeInfo();
  DataTypeInfo(Stride rowStride, Stride colStride,
	       string numRowsVar, string numColsVar,
	       string rowStrideVar, string colStrideVar);

  DataTypeInfo& operator=(const DataTypeInfo &rhs);

};

class Architecture
{
 public:
  // Single precision
  virtual int SVecRegWidth();
  virtual string STypeName();
  virtual string SAddCode(string operand1, string operand2, string result);
  virtual string SMulCode(string operand1, string operand2, string result);
  virtual string SFMACode(string operand1, string operand2, string operand3, string result);
  virtual string SAccumCode(string memPtr, string startinLoc);
  virtual string SContiguousLoad(string memPtr, string receivingLoc);
  virtual string SStridedLoad(string memPtr, string receivingLoc, string stride);
  virtual string SDuplicateLoad(string memPtr, string receivingLoc);
  virtual string SContiguousStore(string memPtr, string startingLoc);
  virtual string SStridedStore(string memPtr, string startingLoc);

  // Double precision
  virtual int DVecRegWidth();
  virtual string DTypeName();
  virtual string DAddCode(string operand1, string operand2, string result);
  virtual string DMulCode(string operand1, string operand2, string result);
  virtual string DFMACode(string operand1, string operand2, string operand3, string result);
  virtual string DAccumCode(string memPtr, string startinLoc);
  virtual string DContiguousLoad(string memPtr, string receivingLoc);
  virtual string DStridedLoad(string memPtr, string receivingLoc, string stride);
  virtual string DDuplicateLoad(string memPtr, string receivingLoc);
  virtual string DContiguousStore(string memPtr, string startingLoc);
  virtual string DStridedStore(string memPtr, string startingLoc);

  // General
  virtual int VecRegWidth(Type type);
};

extern Architecture* arch;

#endif //DOLLDLA
