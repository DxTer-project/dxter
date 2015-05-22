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

#if DOTENSORS

#define ALLMULTIMODEALLGATHER 1

#include "DLANode.h"
#include "transform.h"
#include "DLAOp.h"
#include "tensorSumScatter.h"
#include "base.h"

void GetSuffix(const DimVec &dims1, const DimVec &dims2, 
	       DimVec &suff);

string GetAlignmentSource(Node *node, ConnNum inNum);

class RedistNode : public DLANode
{
 public:
  DataTypeInfo m_info;
  SizesVec m_lsizes;
  bool m_isArray;
  string m_name;
  string m_align;
  DimVec m_alignModes;
  DimVec m_alignModesSrc;
  
  RedistNode();
  RedistNode(const DistType &destType);
  RedistNode(const DistType &destType, const string &align, const DimVec &alignModes, const DimVec &alignModesSrc);
  RedistNode(const DistType &destType, const Permutation &perm, const string &align, 
	     const DimVec &alignModes, const DimVec &alignModesSrc);
  virtual const DataTypeInfo& DataType(ConnNum num) const;
  static Node* BlankInst() { return  new RedistNode; }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual NodeType GetType() const;
  virtual void Prop();
  bool IsPrimitive() const;
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "redist";}
  virtual void ClearDataTypeCache();
  virtual void BuildDataTypeCache();
  virtual const Dim NumDims(ConnNum num) const;
  virtual const SizeList* Len(ConnNum num, Dim dim) const;
  virtual const SizeList* LocalLen(ConnNum num, Dim dim) const;
  virtual Name GetName(ConnNum num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
  virtual Phase MaxPhase() const;
  virtual void AddVariables(VarSet &set) const;
  bool HasSamePermAndAlign(const RedistNode *redist) const;
  bool HasSameAlign(const RedistNode *redist) const;
  bool HasSamePerm(const RedistNode *redist) const;
  bool HasNoAlign() const;
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
  Dim m_dim;
 CombineDisappearingModes(Dim dim) : m_dim(dim) {}
  virtual string GetType() const { return (string)"CombineDisappearingModes " + (char)(m_dim+48); }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};


//Implements [D0] <- [D1] as
// [D10] <- [D1]    (Local copy)
// [D01] <- [D10]   (Permutation)
// [D0] <- [D01]    (AllGather)
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
//[p(x|_|y),p(z|_|w)] <- [(p(x)|_|w),(p(z)|_|y)]
// into
//[(p(x)|_|y),(p(z)|_|w)] <- [(p(x)|_|w),(p(z)|_|y)]
//[p(x|_|y),p(z|_|w)] <- [(p(x)|_|y),(p(z)|_|w)]
// (Also, grid modes can disappear in the final redistribution)
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

class DoubleIndexAllToAll2 : public SingleTrans
{
 public:
  Dim m_dim;
 DoubleIndexAllToAll2(Dim dim) : m_dim(dim) {}
  virtual string GetType() const { return (string)"DoubleIndexAllToAll2" + (char)(m_dim+48); }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};


//E.g., [21,03] -> [01,23]
// as [21,03] -> [12,03]
//    [12,03] -> [12,30]
//    [12,30] -> [21,03]
//Separating legal AllToAll like the following
//[(x|_|y),(z|_|w)] <- [(z),(x)]
// to
//[(x),(z)] <- [(z),(x)]
//[(x|_|y),(z|_|w)] <- [(x),(z)]
class DoubleIndexAllToAllPrefix : public SingleTrans
{
 public:
  Dim m_dim;
 DoubleIndexAllToAllPrefix(Dim dim) : m_dim(dim) {}
  virtual string GetType() const { return (string)"DoubleIndexAllToAllPrefix" + (char)(m_dim+48); }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};

class MultiIndexAllToAll : public SingleTrans
{
 public:
  Dim m_dim;
 MultiIndexAllToAll(Dim dim) : m_dim(dim) {}
  virtual string GetType() const { return (string)"MultiIndexAllToAll" + (char)(m_dim+48); }
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

class CombinePermuteRedists : public SingleTrans
{
 public:
  virtual string GetType() const { return (string)"CombinePermuteRedists";}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

#if ALLMULTIMODEALLGATHER
class CombineAllGathers : public SingleTrans
{
 public:
  virtual string GetType() const { return (string)"CombineAllGathers";}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};
#endif //ALLMULTIMODEALLGATHER

class CombineMovingModes : public SingleTrans
{
 public:
  Dim m_dim;
 CombineMovingModes(Dim dim) : m_dim(dim) {}
  virtual string GetType() const { return (string)"CombineMovingModes" + (char)(m_dim+48); }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};


#endif
