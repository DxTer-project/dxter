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

//declared in base.h
#include "base.h"

#if DOTENSORS


Permutation::Permutation(string start, string end)
{
  if (start == end)
    throw;
  if (start.length() != end.length())
    throw;
  string::iterator iter = start.begin();
  for(; iter != start.end(); ++iter) {
    m_permutation.push_back(end.find(*iter));
  }
}

bool Permutation::operator==(const Permutation &rhs) const
{
  return m_permutation == rhs.m_permutation;
}

void Permutation::SetToDefault(Dim numDims)
{
  for(Dim dim = 0; dim < numDims; ++dim)
    m_permutation.push_back(dim);
}

Permutation& Permutation::operator=(const Permutation &rhs)
{
  m_permutation = rhs.m_permutation;
  return *this;
}

string Permutation::Str() const
{
  string str;
  DimVecConstIter iter = m_permutation.begin();
  for(; iter != m_permutation.end(); ++iter) {
    str += std::to_string(*iter);
  }
  return str;
}

Permutation Permutation::ComposeWith(const Permutation &perm) const
{
  if (perm.Size() != Size())
    throw;
  Permutation newPerm;
  DimVecConstIter iter = m_permutation.begin();
  for(; iter != m_permutation.end(); ++iter)
    newPerm.m_permutation.push_back(perm.Map(*iter));
  return newPerm;
}
#endif
