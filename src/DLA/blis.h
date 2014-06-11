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
#if DOBLIS

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
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

Transpose* AddTranspose(Trans trans, bool objTrans, Node *input, unsigned int num, bool addToPoss);

Node* AddTrans(Trans trans, bool objTrans, Node *input, unsigned int num, bool addToPoss);

Transpose* InsertTranspose(Trans trans, bool objTrans, Node *node, unsigned int inNum, bool addToPoss);

//Input 0 is the off-diagonal partition of the triangular matrix
//Input 1 is the on-diagonal partition
//Input 2 is the on-diagonal partition of the general matrix
class GetUpToDiag : public DLANode
{
 public:
  PartDir m_dir;
  Tri m_tri;
  Sizes *m_sizes, *m_lsizes;
  GetUpToDiag(Tri tri, PartDir dir);
  virtual const Sizes* GetM(unsigned int num) const;
  virtual const Sizes* GetN(unsigned int num) const;
  virtual const Sizes* LocalM(unsigned int num) const;
  virtual const Sizes* LocalN(unsigned int num) const;
  virtual Name GetName(unsigned int num) const;
  virtual void Prop();
  virtual unsigned int NumOutputs() const;
  virtual void ClearSizeCache();
  virtual void BuildSizeCache();

  virtual void UpdateInnerPackingMultiple(PackSize size) {}

  virtual void PrintCode(IndStream &out);
  static Node* BlankInst() {return new GetUpToDiag(NOTTRI,PARTDOWN); }
  bool KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const {return false;}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "GetUpToDiag";}
  virtual NodeType GetType() const {return "GetUpToDiag";}
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  void Flatten(ofstream &out) const;
  void Unflatten(ifstream &in, SaveInfo &info);
  virtual bool Overwrites(const Node *input, unsigned int num) const { return false;}
};

class CombineDiag : public DLAOp<2,1>
{
 public:
  static Node* BlankInst() {return new CombineDiag; }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "CombineDiag";}
  virtual NodeType GetType() const {return "CombineDiag";}
  virtual void Prop();
  virtual void PrintCode(IndStream &out) {}
  virtual bool Overwrites(const Node *input, unsigned int num) const {return false;}
};


class SetObjProps : public DLAOp<1,1>
{
 public:
  Tri m_tri;
  TriStruct m_struct;
  Diag m_diag;
  SetObjProps(Tri tri, Diag diag, TriStruct triStruct);
  static Node* BlankInst() {return new SetObjProps(NOTTRI,NOTTRIDIAG,GEN); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "SetObjProps";}
  virtual NodeType GetType() const {return "SetObjProps";}
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  void Flatten(ofstream &out) const;
  void Unflatten(ifstream &in, SaveInfo &info);
};

class Copy : public DLAOp<2,1>
{
 public:
  static Node* BlankInst() {return new Copy; }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "Copy";}
  virtual NodeType GetType() const {return "Copy";}
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
};

class ParallelizeTrans : public SingleTrans
{
 public:
  Comm m_comm;
 ParallelizeTrans(Comm comm) : m_comm(comm) {}
};

class ParallelizeMDim : public ParallelizeTrans
{
 public:
  ParallelizeMDim(Comm comm)
    : ParallelizeTrans(comm) {}
  virtual string GetType() const {return "ParallelizeMDim"+CommToStr(m_comm);}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

class ParallelizeInnerNDim : public ParallelizeTrans
{
 public:
  ParallelizeInnerNDim(Comm comm)
    : ParallelizeTrans(comm) {}
  virtual string GetType() const {return "ParallelizeInnerNDim"+CommToStr(m_comm);}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

class ParallelizeOuterNDim : public ParallelizeTrans
{
 public:
  ParallelizeOuterNDim(Comm comm)
    : ParallelizeTrans(comm) {}
  virtual string GetType() const {return "ParallelizeOuterNDim"+CommToStr(m_comm);}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

class ParallelizeK : public ParallelizeTrans
{
 public:
  ParallelizeK(Comm comm)
    : ParallelizeTrans(comm) {}
  virtual string GetType() const {return "ParallelizeK"+CommToStr(m_comm);}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};
/*
class IncreaseParallelizedLoop : public SingleTrans
{
  virtual string GetType() const {return "IncreaseParallelizedLoop";}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};
*/
bool LegalParallelizationNestingUp(const Node *node, Comm comm);
bool LegalParallelizationNestingDown(const PSet *pset, Comm comm);

bool FoundBarrier(const Node *node, unsigned int input, Comm comm);

Cost AdditionalCostForBringingIntoL2(Node *node, unsigned int num, Size numAElems, Comm comm);
Cost AdditionalCostOfBarrier(Comm comm, unsigned int numHits);
#endif
