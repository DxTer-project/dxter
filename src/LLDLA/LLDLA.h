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

#ifndef LLDLA_H_
#define LLDLA_H_

#pragma once

#include "architecture.h"
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
      LOG_FAIL("replacement for throw call");
      throw;
    }
}

class DataTypeInfo
{
 public:
  Size m_rowStrideVal;
  Size m_colStrideVal;

  Stride m_rowStride;
  Stride m_colStride;

  string m_numRowsVar;
  string m_numColsVar;
  string m_rowStrideVar;
  string m_colStrideVar;

  Type m_type;

  DataTypeInfo();
  DataTypeInfo(Size rowStrideVal, Size colStrideVal,
	       string numRowsVar, string numColsVar,
	       string rowStrideVar, string colStrideVar,
	       Type type);

  DataTypeInfo& operator=(const DataTypeInfo& rhs);
  string ToString();
};

#endif // DOLLDLA

#endif // LLDLA_H_
