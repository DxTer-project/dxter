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

void GetSuffix(const DimVec &dims1, const DimVec &dims2, 
	       DimVec &suff);


class RedistNode : public DLANode
{
 public:
  DataTypeInfo m_info;
  SizesArray m_lsizes;
  bool m_isArray;
  string m_name;
  RedistNode();
  RedistNode(const DistType &destType);
  virtual ~RedistNode();
  virtual const DataTypeInfo& DataType(ConnNum num) const {return m_info; }
  static Node* BlankInst() { return  new RedistNode; }
  bool KeepsInputVarLive(Node *input, ConnNum numIn, ConnNum &numOut) const {return false;}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual NodeType GetType() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "redist";}
  virtual void ClearDataTypeCache();
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
  virtual NodeType GetType() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  //  RedistNode* CreateTrans(Trans trans);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "AllGather";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual const DataTypeInfo& DataType(ConnNum num) const {return InputDataType(0); }
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



class CombineDisappearingModes : public SingleTrans
{
 public:
  Dim m_srcDim, m_destDim;
 CombineDisappearingModes(Dim srcDim, Dim destDim) : m_srcDim(srcDim), m_destDim(destDim) {}
  virtual string GetType() const { return (string)"CombineDisappearingModes " + (char)(m_srcDim+48) + (char)(m_destDim+48); }
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


//Separating legal AllToAll like the following
//[(x|_|y),(z|_|w)] <- [(x|_|w),(z|_|y)]
class DoubleIndexAllToAll : public SingleTrans
{
 public:
  Dim m_dim;
 DoubleIndexAllToAll(Dim dim) : m_dim(dim) {}
  virtual string GetType() const { return (string)"DoubleIndexAllToAll" + (char)(m_dim+48); }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};

/*
[([x|_|p(y)|_|z),(w)]
[([x|_|z),(w|_|y)]

into

[([x|_|p(y)|_|z),(w)]
[([x|_|z|_|y),(w)]
[([x|_|z),(w|_|y)]
*/
class PermuteDistribution : public SingleTrans
{
 public:
  Dim m_srcDim;
  Dim m_destDim;
 PermuteDistribution(Dim srcDim, Dim destDim) : m_srcDim(srcDim), m_destDim(destDim) {}
  virtual string GetType() const { return (string)"PermuteDistribution " + (char)(m_srcDim+48) + (char)(m_destDim+48); }
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

class SplitAllAllGathers : public SingleTrans
{
 public:
  virtual string GetType() const { return (string)"SplitAllAllGathers";}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

#endif
