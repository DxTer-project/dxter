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
#if DOBOOL

#include "DLAOp.h"

class Not : public Node
{
  Name m_name;
  static DataTypeInfo m_info;

  public:

  Not();
  string GetUniqueName();
  static Node* BlankInst() { return new Not; }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "not";}
  virtual Cost GetCost() {return 1;}
  virtual Name GetName(ConnNum num) const;
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
  virtual const DataTypeInfo& DataType(ConnNum num) const { return m_info;}
  virtual NodeType GetType() const {return "not";}
  virtual Phase MaxPhase() const {return NUMPHASES;}
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  //  virtual void AddVariables(VarSet &set) const;
};

class NotTrue : public SingleTrans
{
 public:
  virtual string GetType() const {return "Not with True";}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

class NotFalse : public SingleTrans
{
 public:
  virtual string GetType() const {return "Not with False";}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

class NotNot : public SingleTrans
{
 public:
  virtual string GetType() const {return "Not Not";}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

#endif //DOBOOL
