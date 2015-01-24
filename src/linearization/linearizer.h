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
#include "linearization.h"

typedef std::map<const void*,LinElem*> PtrToLinElemMap;
typedef PtrToLinElemMap::iterator PtrToLinElemMapIter;
typedef std::set<LinElem*> LinElemSet;
typedef LinElemSet::iterator LinElemSetIter;

class Linearizer
{
 public:
  LinElemVec m_elems;
  Linearization m_lin;

  Linearizer(const Poss *poss);
  ~Linearizer();

  void ClearCurrLinearization();
  bool HasCurrLinearization() const;

  void FindAnyLinearization();
  void FindOptimalLinearization(const StrSet &stillLive);

  void InsertVecClearing(const StrSet &stillLive);

  LinElem* FindOrAdd(Node *node, PtrToLinElemMap &map);
  LinElem* FindOrAdd(BasePSet *set, PtrToLinElemMap &map);

  void PrintConnections() const;
};
