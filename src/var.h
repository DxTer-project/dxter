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

#if DOTENSORS

string ModeArrayVarName(const DimVec &vec);
string IndexPairVarName(Dim dim1, Dim dim2);
string ModeArrayPairVarName(const DimVec &arr1, const DimVec &arr2);
string TensorDistVarName(const DistType &type);
string IndexArrayVarName(const string &indices);
#endif



enum VarType {
#if DOTENSORS
  TensorVarType,
  ModeArrayVarType,
  IndexPairType,
  ModeArrayPairVarType,
  TensorDistVarType,
  IndexArrayType,
#endif
  InvalidType
};

class Var
{
 public:
  VarType m_type;
  union {
    Name *m_name;
#if DOTENSORS
    DimVec *m_vec;
    std::pair<DimVec, DimVec> *m_arrPair;
    std::pair<Dim,Dim> *m_pair;
    string *m_indices;
#endif
#if DODM
    DistType *m_distType;
#endif
  };
  string m_compStr;
 Var() : m_type(InvalidType) {}
#if DOTENSORS
  Var(const Name &name);
  Var(const DimVec &vec);
  Var(const DimVec &vec1, const DimVec &vec2);
  Var(Dim dim1, Dim dim2);
  Var(const string &indices);
#endif
  Var(const Var &var);
#if DODM
  Var(const DistType &type);
#endif
  ~Var();
    Var& operator=(const Var &rhs);
  string CompStr() const {return m_compStr;}
  void PrintDecl(IndStream &out) const;
  string GetVarName() const;
};

struct VarCompare {
  bool operator() (const Var& lhs, const Var& rhs) const{
    return lhs.CompStr() < rhs.CompStr();
  }
};

typedef set<Var,VarCompare> VarSet; // inefficient compare, so only do with small sizes
typedef VarSet::iterator VarSetIter;
typedef VarSet::const_iterator VarSetConstIter;

