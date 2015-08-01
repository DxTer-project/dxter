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

#include "layers.h"
#include "sortable.h"
#include "rqoBasis.h"
#include "transform.h"

#if DORQO



class Union : public Sortable
{
 protected:
  DataTypeInfo m_dataTypeInfo;
  

 public:
  vector<string> m_fields;

  Union();
  Union(string sortBy, vector<string> inFields);
  virtual NodeType GetType() const;
  static Node* BlankInst() { return  new Union; }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual const DataTypeInfo& DataType(ConnNum num) const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual Cost GetCost() {return 0;}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "union";}
  virtual void ClearDataTypeCache();
  virtual void BuildDataTypeCache();
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
};


#endif
