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



#pragma once

#include "base.h"
#include "linElem.h"

typedef std::vector<LinElem*> LinElemVec;
typedef LinElemVec::iterator LinElemVecIter;

class Linearization
{
 public:
  LinElemVec m_order;
  LinElemVec m_clears;

  Cost m_cost;

 Linearization() : m_cost(-1) {}
  ~Linearization();
  void Clear();

  void InsertVecClearing(const StrSet &stillLive, const StrSet &alwaysLive);
  Cost GetCostNoRecursion(const StrSet &stillLive, const StrSet &alwaysLive);

  void Push(LinElem *elem);
  void PopNonClear();
  void operator=(const Linearization &rhs);
  void Print(IndStream &out);
  
  bool LiveAfter(unsigned int loc, const string &name, const StrSet &alwaysLive) const;

  bool EnforceMemConstraint(Cost maxMem, const StrSet &stillLive, const StrSet &alwaysLive);
};
