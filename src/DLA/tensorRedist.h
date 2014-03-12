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


#pragma once

#include "layers.h"

#if DOTENSORS

#include "DLANode.h"
#include "transform.h"

class RedistNode : public DLANode
{
 public:
  DistType m_destType;
  SizesArray m_lsizes;
  string m_name;
  RedistNode();
  RedistNode(const DistType &destType);
  virtual ~RedistNode();
  virtual const DistType& GetDistType(unsigned int num) const { return m_destType; }
  static Node* BlankInst() { return  new RedistNode; }
  bool KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const {return false;}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual bool IsRedistNode() const {return true;}
  virtual NodeType GetType() const;
  virtual void SanityCheck();
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  RedistNode* CreateTrans(Trans trans);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "redist";}
  virtual void ClearSizeCache();
  virtual void BuildSizeCache();
  virtual const Dim NumDims(unsigned int num) const;
  virtual const Sizes* Len(unsigned int num, Dim dim) const;
  virtual const Sizes* LocalLen(unsigned int num, Dim dim) const;
  virtual Name GetName(unsigned int num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool Overwrites(const Node *input, unsigned int num) const {return false;}
};

class RemoveWastedRedist : public SingleTrans
{
 public:
  virtual string GetType() const { return "RemoveWastedRedist"; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

class RemoveNOPRedistribs : public SingleTrans
{
 public:
  virtual string GetType() const { return "RemoveNOPRedist"; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

class CombineRedistribs : public SingleTrans
{
 public:
  virtual string GetType() const { return "CombineRedist"; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

#endif
