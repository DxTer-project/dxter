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
#include "rqoHelperNodes.h"
#include "rqoBasis.h"

#if DORQO


class NIndexedNode : public InputNode
{

 protected:
  DataTypeInfo m_dataTypeInfo;

 public:
  set<int> m_indeces;

  NIndexedNode();
  NIndexedNode(string name, string sortBy, set<string> fields, string fileName, string query, set<int> indeces);
  virtual NodeType GetType() const {return m_type;}
  static Node* BlankInst() { return  new NIndexedNode; }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual const DataTypeInfo& DataType(ConnNum num) const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual Cost GetCost() {return 0;}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "nindexednode";}
  virtual Name GetName(ConnNum num) const;
  virtual void ClearDataTypeCache() {}
  virtual void BuildDataTypeCache() {}
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
};

#endif