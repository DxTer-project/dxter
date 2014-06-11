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

#include "transform.h"

#if DOELEM
enum OneDDistType {
  STAR,
  MC,
  MR,
  VC,
  VR,
  LASTONEDDISTTYPE };

bool CanUseUp(const Node *root);
bool CanUseDown(Node *root);
bool CanTrans(DistType src, DistType dest, bool tryOpposite = false);

Trans UpdateTrans(Trans trans, DistType dist);
Trans DistTransType(DistType dist);
Trans SwitchTrans(Trans trans, Type type);
DistType GetNonTrans(DistType dist);
DistType TransType(DistType dist, Trans trans);
bool IsTransType(DistType dist);

OneDDistType GetColDist( DistType dist);
OneDDistType GetRowDist( DistType dist);

//void GetLocalSizes(DistType dist, const Sizes *m, const Sizes *n, Sizes &localM, Sizes &localN);
//void GetLocalSize(DistType dist, Size m, Size n, Size &localM, Size &localN);

class RemoveWastedRedist : public SingleTrans
{
  string m_type;
  DistType m_destType;
 public:
  RemoveWastedRedist(DistType destType);
  virtual string GetType() const { return m_type; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

template<DistType SrcType, DistType DestType>
class ExpandRedistribution : public SingleTrans
{
  string m_type;
 public:
  ExpandRedistribution();
  virtual string GetType() const { return m_type; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool WorthApplying(const Node *node) const;
  virtual bool ExpansionHasPossibleTrans() const;
};

class CombineRedistribs : public SingleTrans
{
  string m_type;
  DistType m_srcType;
  DistType m_destType;
 public:
  CombineRedistribs(DistType srcType, DistType destType);
  virtual string GetType() const { return m_type; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

class RemoveNOPRedistribs : public SingleTrans
{
  string m_type;
  DistType m_distribType;
 public:
  RemoveNOPRedistribs(DistType type);
  virtual string GetType() const { return m_type; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

class FindMidDistributions : public SingleTrans
{
  string m_type;
  DistType m_srcType;
  DistType m_midType;
  DistType m_destType;
 public:
  FindMidDistributions(DistType srcType, DistType midType, DistType destType);
  virtual string GetType() const { return m_type; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

class ReplaceWithTrans : public SingleTrans
{
  string m_type;
  DistType m_srcType;
  DistType m_origDestType;
  DistType m_newDestType;
 public:
  ReplaceWithTrans(DistType srcType, DistType origDestType, DistType newDestType);
  virtual string GetType() const { return m_type; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

class UseTransposedRedist : public SingleTrans
{
  DistType m_destType;
  DistType m_srcType1;
  DistType m_srcType2;
 public:
  UseTransposedRedist(DistType destType, DistType srcType1, DistType srcType2); 
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

class TransTransformation : public SingleTrans
{
 public:
  unsigned int m_argNum;
  Trans m_trans;
  TransTransformation(unsigned int argNum, Trans trans) : m_argNum(argNum), m_trans(trans) {}
  void Apply(Poss *poss, Node *node) const;
  virtual void PreApply(Node *node) const {}
  virtual void PostApply(Node *node) const {}
  virtual string GetType() const;
  virtual string GetTransType() const = 0;
};


class RedistTrans : public TransTransformation
{
 public:
  RedistTrans(Trans trans) : TransTransformation(0,trans) {}
  virtual string GetTransType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual bool WorthApplying(const Node *node) const;
};


class RedistNode : public DLANode
{
 public:
  DistType m_destType;
  Sizes *m_mSizes, *m_nSizes;
 RedistNode() : m_destType(UNKNOWN), m_mSizes(NULL), m_nSizes(NULL) {}
  RedistNode(DistType destType);
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
  virtual bool CanTrans() const;
  RedistNode* CreateTrans(Trans trans);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "redist";}
  virtual void ClearSizeCache();
  virtual void BuildSizeCache();
  virtual const Sizes* GetM(unsigned int num) const;
  virtual const Sizes* GetN(unsigned int num) const;
  virtual const Sizes* LocalM(unsigned int num) const;
  virtual const Sizes* LocalN(unsigned int num) const;
  virtual Name GetName(unsigned int num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  static Cost GetCost(DistType srcType, DistType destType, const Sizes *ms, const Sizes *ns);
  virtual bool Overwrites(const Node *input, unsigned int num) const {return false;}
};

class SumScatterNode : public DLANode
{
 public:
  Coef m_coeff;
 SumScatterNode() : m_coeff(COEFZERO) {}
 SumScatterNode(Coef coeff) : m_coeff(coeff) {}
  virtual const DistType& GetDistType(unsigned int num) const { return InputDistType(1); }
  static Node* BlankInst() { return  new SumScatterNode; }
  bool KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const;
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual NodeType GetType() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "SumScatter";}
  virtual const Sizes* GetM(unsigned int num) const;
  virtual const Sizes* GetN(unsigned int num) const;
  virtual const Sizes* LocalM(unsigned int num) const;
  virtual const Sizes* LocalN(unsigned int num) const;
  virtual Name GetName(unsigned int num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  static Cost GetCost(const Sizes *localM, const Sizes *localN, DistType destType, DistType srcType);
  virtual bool Overwrites(const Node *input, unsigned int num) const;
};

class SumScatterFrom : public DLANode
{
 public:
  SumScatterFrom() {}
  virtual const DistType& GetDistType(unsigned int num) const { return InputDistType(1); }
  static Node* BlankInst() { return  new SumScatterFrom; }
  bool KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const;
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual NodeType GetType() const {return "sumScatterFrom";}
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "SumScatterFrom";}
  virtual const Sizes* GetM(unsigned int num) const;
  virtual const Sizes* GetN(unsigned int num) const;
  virtual const Sizes* LocalM(unsigned int num) const;
  virtual const Sizes* LocalN(unsigned int num) const;
  virtual Name GetName(unsigned int num) const;
  static Cost GetCost(DistType destType, DistType srcType, const Sizes *localMs, const Sizes *localNs);
  virtual bool Overwrites(const Node *input, unsigned int num) const;
};


class SumOverCommNode : public DLANode
{
 public:
  SumOverCommNode() {} 
  virtual const DistType& GetDistType(unsigned int num) const { return InputDistType(0); }
  static Node* BlankInst() { return  new SumOverCommNode; }
  bool KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const;
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual NodeType GetType() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "SumOverComm";}
  virtual const Sizes* GetM(unsigned int num) const;
  virtual const Sizes* GetN(unsigned int num) const;
  virtual const Sizes* LocalM(unsigned int num) const;
  virtual const Sizes* LocalN(unsigned int num) const;
  virtual Name GetName(unsigned int num) const;
  static Cost GetCost(DistType destType, const Sizes *localMs, const Sizes *localNs);
  virtual bool Overwrites(const Node *input, unsigned int num) const;
};

class UniqueTransTrans : public SingleTrans
{
 public:
  virtual string GetType() const {return "Unique Trans Trans";}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};
  

template<DistType destType>
Cost GetRedistCost(DistType srcType, Size m, Size n);


template<DistType destType>
inline 
Cost GetRedistCost(DistType srcType, const Sizes *m, const Sizes *n) 
{
  Cost cost = 0;
  unsigned int length = m->NumSizes();
  unsigned int length2 = n->NumSizes();
  if (length != length2)
    throw;
  for(unsigned int i = 0; i < length; ++i)
    cost += GetRedistCost<destType>(srcType, (*m)[i], (*n)[i]);
  return cost;
}


#endif //DOELEM
