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

#include "linElem.h"

class ClearLinElem : public LinElem
{
 public:
  string m_name;

 ClearLinElem(string &name) : m_name(name) {}
 ClearLinElem(const ClearLinElem &clear) : m_name(clear.m_name) {}
  
  virtual bool IsClear() const {return true;}

  virtual void Print();
  virtual StrVec PossiblyDyingVars() const {throw;}
  virtual StrSet NewVars() const {throw;}
  virtual VarCostMap NewVarsAndCosts() const {throw;}
  virtual bool UsesInputVar(const string &var) const {throw;}
  virtual void CacheLiveVars(const StrSet &stillLive) {throw;}
  virtual void ClearCache() {}
};
