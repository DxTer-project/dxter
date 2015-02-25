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

#include "LLDLA.h"

#if DOLLDLA

DataTypeInfo::DataTypeInfo()
  : m_rowStride(BADSTRIDE),
    m_colStride(BADSTRIDE)
{
  m_rowStride = BADSTRIDE;
  m_colStride = BADSTRIDE;
}

DataTypeInfo::DataTypeInfo(Size numRows, Size numCols,
			   Size rowStrideVal, Size colStrideVal,
			   string numRowsVar, string numColsVar,
			   string rowStrideVar, string colStrideVar,
			   Type type)
  : //m_numRows(numRows),
    //m_numCols(numCols),
    //    m_rowStrideVal(rowStrideVal),
    //    m_colStrideVal(colStrideVal),
    m_numRowsVar(numRowsVar),
    m_numColsVar(numColsVar),
    m_rowStrideVar(rowStrideVar),
    m_colStrideVar(colStrideVar),
    m_type(type)
{
  /*  if (m_rowStrideVal == 1) {
    m_rowStride = UNITSTRIDE;
  } else {
    m_rowStride = NONUNITSTRIDE;
  }

  if (m_colStrideVal == 1) {
    m_colStride = UNITSTRIDE;
  } else {
    m_colStride = NONUNITSTRIDE;
    }*/

}

DataTypeInfo& DataTypeInfo::operator=(const DataTypeInfo &rhs)
{
  /*  m_numRows = rhs.m_numRows;
  m_numCols = rhs.m_numCols;
  m_rowStrideVal = rhs.m_rowStrideVal;
  m_colStrideVal = rhs.m_colStrideVal;*/
  m_rowStride = rhs.m_rowStride;
  m_colStride = rhs.m_colStride;
  m_numRowsVar = rhs.m_numRowsVar;
  m_numColsVar = rhs.m_numColsVar;
  m_rowStrideVar = rhs.m_rowStrideVar;
  m_colStrideVar = rhs.m_colStrideVar;
  m_type = rhs.m_type;
  return *this;
}

bool DataTypeInfo::IsGenStride() const {
  /*  if (m_rowStrideVal != 1 && m_colStrideVal != 1) {
    cout << "NOT GENERAL STRIDE" << endl;
    cout << "m_rowStrideVal == 1 ? " << (m_rowStrideVal == 1) << endl;
    cout << "m_colStrideVal == 1 ? " << (m_colStrideVal == 1) << endl;
    return true;
    }*/
  //  return false;
  throw;
}

bool DataTypeInfo::IsContiguous() const {
  /*  if (m_rowStrideVal == 1 && m_colStrideVal == m_numRows) {
    return true;
  }

  if (m_colStrideVal == 1 && m_rowStrideVal == m_numCols) {
    return true;
  }

  if (IsRowVector() && m_rowStrideVal > m_numRows) {
    return true;
  }

  if (IsColVector() && m_colStrideVal > m_numCols) {
    return true;
    }*/

  //  return false;
  throw;
}

bool DataTypeInfo::IsRowVector() const {
  //  return m_numRows == 1;
  throw;
}

bool DataTypeInfo::IsColVector() const {
  //  return m_numCols == 1;
  throw;
}

bool DataTypeInfo::IsScalar() const {
  return IsRowVector() && IsColVector();
}

bool DataTypeInfo::IsSameSizeAs(const DataTypeInfo& other) const {
  //  return m_numRows == other.m_numRows && m_numCols == other.m_numCols;
  throw;
}

string DataTypeInfo::ToString() {
  string dataStr = "Num rows var name: " + m_numRowsVar + "\n";
  dataStr += "Num cols var name: " + m_numColsVar + "\n";
  dataStr += "Row stride var name: " + m_rowStrideVar + "\n";
  dataStr += "Col stride var name: " + m_colStrideVar + "\n";
  /*  dataStr += "Num rows = " + std::to_string(m_numRows);
  dataStr += "\nNum cols = " + std::to_string(m_numCols);
  dataStr += "\nRow stride = " + std::to_string(m_rowStrideVal);
  dataStr += "\nColStride = " + std::to_string(m_colStrideVal) + "\n";*/
  return dataStr;
}

#endif
