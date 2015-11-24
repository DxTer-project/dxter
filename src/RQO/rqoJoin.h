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
#include <climits>

#if DORQO



class Join : public Sortable
{
 protected:
  DataTypeInfo m_dataTypeInfo;
  

 public:
  vector<string> m_in0Fields;
  vector<string> m_in1Fields;

  Join();
  Join(string sortBy, vector<string> in0Fields, vector<string> in1Fields);
  virtual NodeType GetType() const;
  static Node* BlankInst() { return  new Join; }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual const DataTypeInfo& DataType(ConnNum num) const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual Cost GetCost() {return INT_MAX;}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "join";}
  virtual void ClearDataTypeCache();
  virtual void BuildDataTypeCache();
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
  virtual bool IsJoin() const {return true;}
  virtual Join* CreateCopyOfJoin() const;
  virtual int Outputs() {return (Input(0)->Outputs() > Input(1)->Outputs()) ? Input(0)-> Outputs() : Input(1)->Outputs();}
};

class SwapNodes : public SingleTrans
{
 public:
  unsigned int m_inNum;
  ClassType m_type;
  SwapNodes(unsigned int inNum, ClassType type);
  virtual string GetType() const {return "Switch two node's positions " + std::to_string(m_inNum) + " on " + m_type;}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

class JoinToHash : public SingleTrans
{
public:
  virtual string GetType() const {return "Turn Join to HashJoin";}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

class JoinToNested : public SingleTrans
{
public:
  virtual string GetType() const {return "Turn Join to NestedJoin";}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

class JoinToMerge : public SingleTrans
{
public:
  virtual string GetType() const {return "Turn Join to MergeJoin";}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

#endif
