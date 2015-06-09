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
#include "node.h"

#if DORQO

class DataTypeInfo
{
 public:
  DataTypeInfo& operator=(const DataTypeInfo &rhs) {return *this;}
};


class InputNode : public Node
{
  NodeType m_type;
 protected:
  DataTypeInfo m_dataTypeInfo;
  Name m_varName;

 public:
  InputNode();
  InputNode(string str);
  virtual NodeType GetType() const {return m_type;}
  static Node* BlankInst() { return  new InputNode; }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual const DataTypeInfo& DataType(ConnNum num) const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual Cost GetCost() {return 0;}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "inputNode";}
  virtual Name GetName(ConnNum num) const;
  virtual void ClearDataTypeCache() {}
  virtual void BuildDataTypeCache() {}
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
};

class OutputNode : public Node
{
 public:
  virtual NodeType GetType() const {return "out";}
  static Node* BlankInst() { return  new OutputNode; }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual const DataTypeInfo& DataType(ConnNum num) const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "outputNode";}
  virtual Name GetName(ConnNum num) const {throw;}
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
  virtual Cost GetCost() {return 0;}
  virtual void ClearDataTypeCache() {}
  virtual void BuildDataTypeCache() {}
};

class TempVarNode : public Node
{
 public:
  DataTypeInfo m_info;
  string m_name;

  TempVarNode();
  TempVarNode(string name);

  virtual NodeType GetType() const;
  static Node* BlankInst() { return  new TempVarNode; }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);

  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "tempVar";}
  virtual const DataTypeInfo& DataType(ConnNum num) const;
  virtual Name GetName(ConnNum num) const;
  virtual void ClearDataTypeCache();
  virtual void BuildDataTypeCache();
  virtual bool IsReadOnly() const {return true;}
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
  virtual bool IsDataDependencyOfInput() const {return false;}
  virtual Cost GetCost() {return 0;}
};

#endif //DORQO
