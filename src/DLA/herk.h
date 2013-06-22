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

#include "transform.h"
#include "DLAOp.h"
#include "distributions.h"
#include "lowerLayer.h"

Loop* HerkLoopVar1(Node *Ain, unsigned int Anum, 
		   Node *Cin, unsigned int Cnum,
		   Tri tri,
		   Trans trans,
		   Coef alpha, Coef beta, Type type,
		   Layer layer);

Loop* HerkLoopVar2(Node *Ain, unsigned int Anum, 
		   Node *Cin, unsigned int Cnum,
		   Tri tri,
		   Trans trans,
		   Coef alpha, Coef beta, Type type,
		   Layer layer);

Loop* HerkLoopVar5(Node *Ain, unsigned int Anum, 
		   Node *Cin, unsigned int Cnum,
		   Tri tri,
		   Trans trans,
		   Coef alpha, Coef beta, Type type,
		   Layer layer);

Loop* BLISHerkLoop(Node *Ain, unsigned int Anum, 
		   Node *Bin, unsigned int Bnum,
		   Node *Cin, unsigned int Cnum,
		   Tri tri,
		   Coef alpha, Type type,
		   Layer layer);

class HerkProps
{
 public:
  Trans m_transA, m_transB;
  Tri m_tri;
  Coef m_alpha;
  Coef m_beta;
  Type m_type;
  HerkProps(Tri tri, Trans transA, Trans transB, Coef alpha, Coef beta, Type type);
  void Duplicate(const HerkProps *orig);
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
};

class Herk : public DLAOp<2,1>, public HerkProps
{
 public:
  Herk(Layer layer, Tri tri, Trans trans, 
       Coef alpha, Coef beta, Type type);
  virtual NodeType GetType() const;
  static Node* BlankInst() { return new Herk(ABSLAYER, LOWER, NORMAL, COEFONE, COEFONE, REAL); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual DistType GetDistType(unsigned int num) const;
  virtual void SanityCheck();
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "Herk";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual Phase MaxPhase() const;
  virtual bool ShouldCullDP() const;
};

class HerkLoopExp : public SingleTrans
{
 public:  
  unsigned int m_var;
  Layer m_fromLayer, m_toLayer;
 HerkLoopExp(Layer fromLayer, Layer toLayer, unsigned int variant) 
   : m_var(variant), m_fromLayer(fromLayer), m_toLayer(toLayer) {}
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};

//This is a triangular rank-k update used in Elemental
class TriRK : public DLAOp<3,1>, public HerkProps
{
 public:
  TriRK(Layer layer, Tri tri, Trans transA, Trans transB, Coef alpha, Coef beta, Type type);
  virtual NodeType GetType() const;
  static Node* BlankInst() { return new TriRK(SMLAYER, LOWER, NORMAL, NORMAL, COEFONE, COEFONE, REAL); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual DistType GetDistType(unsigned int num) const;
  virtual void SanityCheck();
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "TriRK";}
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual bool CanTransposeInputs() const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
};

class TriRKTrans : public TransTransformation
{
 public:
  TriRKTrans(unsigned int argNum, Trans trans) : TransTransformation(argNum, trans) {}
  virtual string GetTransType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
};

class DistHerkToLocalTriRK : public SingleTrans
{
 public:
  virtual string GetType() const {return "Distributed Herk to Local TriRK";}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};

class BLISHerkLoopExp : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  BLISHerkLoopExp(Layer from, Layer to)
    : m_fromLayer(from), m_toLayer(to) {}
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
  virtual Cost RHSCostEstimate(const Node *node) const {throw;}
};

class HerkBP : public DLAOp<3,1>
{
 public:
  Tri m_tri;
  Coef m_alpha;
  HerkBP(Layer layer, Tri tri, Coef alpha);
  static Node* BlankInst() { return new HerkBP(S3LAYER,UPPER,COEFONE); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "HerkBP";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual DistType GetDistType(unsigned int num) const { return InputDistType(2);}
  virtual Phase MaxPhase() const {return NUMPHASES;}
  virtual NodeType GetType() const {return "HerkBP";}
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
};

class HerkLowerLayer : public LowerLayer
{
 public:
 HerkLowerLayer(Layer fromLayer, Layer toLayer, Dim dim, Size bs)
   : LowerLayer(fromLayer, toLayer, dim, bs) {}
  virtual string GetType() const { return "Herk lower layer"; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};
