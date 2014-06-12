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
#include "tensorSumScatter.h"

class RedistNode : public DLANode
{
 public:
  DistType m_destType;
  SizesArray m_lsizes;
  bool m_isArray;
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
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "redist";}
  virtual void ClearDataTypeCache();
  virtual void BuildDataTypeCache();
  virtual const Dim NumDims(unsigned int num) const;
  virtual const Sizes* Len(unsigned int num, Dim dim) const;
  virtual const Sizes* LocalLen(unsigned int num, Dim dim) const;
  virtual Name GetName(unsigned int num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool Overwrites(const Node *input, unsigned int num) const {return false;}
  virtual Phase MaxPhase() const;
  virtual void AddVariables(VarSet &set) const;
};

class AllReduceNode : public DLAOp<1,1>
{
 public:
  DimVec m_sumDims;
  string m_sumIndices;
 AllReduceNode() : DLAOp<1,1>() {}
  AllReduceNode(const DimVec &sumDims, const string &sumIndices);
  static Node* BlankInst() { return  new AllReduceNode; }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual bool IsRedistNode() const {return true;}
  virtual NodeType GetType() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  //  RedistNode* CreateTrans(Trans trans);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "AllGather";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual const DistType& GetDistType(unsigned int num) const { return InputDistType(0); }
};



class RemoveWastedRedist : public SingleTrans
{
 public:
  virtual string GetType() const { return "RemoveWastedRedist"; }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

class RemoveNOPRedistribs : public SingleTrans
{
 public:
  virtual string GetType() const { return "RemoveNOPRedist"; }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

class CombineRedistribs : public SingleTrans
{
 public:
  virtual string GetType() const { return "CombineRedist"; }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

class SplitRedistribs : public SingleTrans
{
 public:
  Dim m_dim;
 SplitRedistribs(Dim dim) : m_dim(dim) {}
  virtual string GetType() const { return (string)"SplitRedist" + (char)(m_dim+48); }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};


class SingleIndexAllToAll : public SingleTrans
{
 public:
  Dim m_dim;
 SingleIndexAllToAll(Dim dim) : m_dim(dim) {}
  virtual string GetType() const { return (string)"SingleIndexAllToAll" + (char)(m_dim+48); }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};



class SplitAllGathers : public SingleTrans
{
 public:
  Dim m_dim;
 SplitAllGathers(Dim dim) : m_dim(dim) {}
  virtual string GetType() const { return (string)"SplitAllGathers" + (char)(m_dim+48); }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

#endif
