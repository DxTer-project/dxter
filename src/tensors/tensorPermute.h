/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2014, The University of Texas and Bryan Marker

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

#if DOTENSORS

#include "DLANode.h"
#include "transform.h"
#include "DLAOp.h"

class Permute : public DLANode
{
 public:
  DataTypeInfo m_info;
  Permutation m_permutation;
  Permute(string start, string end, Layer layer);
  Permute(const Permutation &permutation, Layer layer);
  virtual const DataTypeInfo& DataType(ConnNum num) const {return m_info;}
  static Node* BlankInst() { return new Permute("ab","ba",ABSLAYER); }
  bool KeepsInputVarLive(Node *input, ConnNum numIn, ConnNum &numOut) const {return false;}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual NodeType GetType() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "permute";}
  virtual void BuildDataTypeCache();
  virtual const Dim NumDims(ConnNum num) const;
  virtual const Sizes* Len(ConnNum num, Dim dim) const;
  virtual const Sizes* LocalLen(ConnNum num, Dim dim) const;
  virtual Name GetName(ConnNum num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
  virtual Phase MaxPhase() const;
  virtual void AddVariables(VarSet &set) const;
};

class LowerPermute : public SingleTrans
{
 public:
  virtual string GetType() const { return "LowerPermute"; }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

class PermuteLoopHoist : public SingleTrans
{
 public:
  virtual string GetType() const { return "PermuteLoopHoist"; }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

class CombinePermutations : public SingleTrans
{
 public:
  virtual string GetType() const { return "CombinePermutations"; }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};



#endif //DOTENSORS
