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
#if DOBLIS||DOLLDLA

#include "DLANode.h"
#include "DLAOp.h"
#include "pack.h"
#include "comm.h"
#include "loopSupport.h"
#include "LLDLA.h"

class Transpose : public DLAOp<1,1>
{
 public:
  Trans m_trans;
#if DOLLDLA
  DataTypeInfo m_info;
#endif
#if DOBLIS
  bool m_objTrans;
  Transpose(Trans trans, bool objectTrans);
#elif DOLLDLA
  Transpose(Trans trans);
#endif
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void PrintCode(IndStream &out);
#if DOBLIS
  static Node* BlankInst() { return new Transpose(NORMAL, false); }
#elif DOLLDLA
  static Node* BlankInst() { return new Transpose(NORMAL); }
#endif
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "Transpose";}
  virtual NodeType GetType() const {return "Transpose";}
  void Flatten(ofstream &out) const;
  void Unflatten(ifstream &in, SaveInfo &info);
  virtual Name GetName(ConnNum num) const;
  virtual void Prop();
  virtual unsigned int NumOutputs() const {return 1;}
  virtual const SizeList* GetM(ConnNum num) const;
  virtual const SizeList* GetN(ConnNum num) const;
#if DODM
  virtual const SizeList* LocalM(ConnNum num) const;
  virtual const SizeList* LocalN(ConnNum num) const;
#endif
  virtual bool Overwrites(const Node *input, ConnNum num) const;
#if DOLLDLA
  virtual const DataTypeInfo& DataType(ConnNum num) const;
  virtual void ClearDataTypeCache();
  virtual void BuildDataTypeCache();
  virtual void AddVariables(VarSet &set) const;
#endif
};

class CombineTranspose : public SingleTrans
{
 public:
  virtual string GetType() const {return "Combine transpose";}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

#if DOBLIS
Transpose* AddTranspose(Trans trans, bool objTrans, Node *input, ConnNum num, bool addToPoss);
//Node* AddTrans(Trans trans, bool objTrans, Node *input, ConnNum num, bool addToPoss);
Transpose* InsertTranspose(Trans trans, bool objTrans, Node *node, ConnNum inNum, bool addToPoss);
#elif DOLLDLA
Transpose* AddTranspose(Trans trans, Node *input, ConnNum num, bool addToPoss);
//Node* AddTrans(Trans trans, Node *input, ConnNum num, bool addToPoss);
Transpose* InsertTranspose(Trans trans, Node *node, ConnNum inNum, bool addToPoss);
#endif

#endif
