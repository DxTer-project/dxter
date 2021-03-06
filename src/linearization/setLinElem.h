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
#include "graphIter.h"

class BasePSet;


class SetLinElem : public LinElem
{
 public:
  BasePSet *m_set;
  StrSet m_live;

 SetLinElem(BasePSet *set) : LinElem(), m_set(set) {}
  
  virtual void Print(IndStream &out);
  void Print(IndStream &out, GraphIter *graphIter, Poss *poss);
  virtual StrVec PossiblyDyingVars() const;
  virtual StrSet NewVars() const;
  virtual VarCostMap NewVarsAndCosts() const;
  virtual bool UsesInputVar(const string &var) const;
  virtual void CacheLiveVars(const StrSet &stillLive);
  virtual void ClearCache();
  virtual bool IsSet() const {return true;}
};
