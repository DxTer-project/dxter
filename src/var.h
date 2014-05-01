/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2013, The University of Texas and Bryan Marker

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


#include "base.h"

enum VarType {
#if DOTENSORS
  TensorVarType,
  ModeArrayVarType,
#endif
  InvalidType
};

class Var
{
 public:
  VarType m_type;
  union {
    Name *m_name;
    DimVec *m_vec;
  };
  string m_compStr;
 Var() : m_type(InvalidType) {}
  Var(const Name &name);
  Var(const DimVec &vec);
  Var(const Var &var);
  ~Var();
    Var& operator=(const Var &rhs);
  virtual string CompStr() const {return m_compStr;}
  virtual void PrintDecl(IndStream &out) const;
};

struct VarCompare {
  bool operator() (const Var& lhs, const Var& rhs) const{
    return lhs.CompStr() < rhs.CompStr();
  }
};

typedef set<Var,VarCompare> VarSet; // inefficient compare, so only do with small sizes
typedef VarSet::iterator VarSetIter;
typedef VarSet::const_iterator VarSetConstIter;

