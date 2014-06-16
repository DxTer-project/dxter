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
#if DOBLIS||DOLLDLA

#include "DLANode.h"
#include "DLAOp.h"
#include "pack.h"
#include "comm.h"
#include "loopSupport.h"

class Transpose : public DLAOp<1,1>
{
 public:
  Trans m_trans;
  bool m_objTrans;
  Transpose(Trans trans, bool objectTrans);
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void PrintCode(IndStream &out);
  static Node* BlankInst() {return new Transpose(NORMAL, false); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "Transpose";}
  virtual NodeType GetType() const {return "Transpose";}
  void Flatten(ofstream &out) const;
  void Unflatten(ifstream &in, SaveInfo &info);
  virtual Name GetName(unsigned int num) const;
  virtual void Prop();
  virtual unsigned int NumOutputs() const {return 1;}
  virtual const Sizes* GetM(unsigned int num) const;
  virtual const Sizes* GetN(unsigned int num) const;
  virtual const Sizes* LocalM(unsigned int num) const;
  virtual const Sizes* LocalN(unsigned int num) const;
  virtual bool Overwrites(const Node *input, unsigned int num) const;
};

class CombineTranspose : public SingleTrans
{
 public:
  virtual string GetType() const {return "Combine transpose";}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

Transpose* AddTranspose(Trans trans, bool objTrans, Node *input, unsigned int num, bool addToPoss);

Node* AddTrans(Trans trans, bool objTrans, Node *input, unsigned int num, bool addToPoss);

Transpose* InsertTranspose(Trans trans, bool objTrans, Node *node, unsigned int inNum, bool addToPoss);


#endif
