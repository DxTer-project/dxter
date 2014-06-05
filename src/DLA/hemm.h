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
#if DOELEM||DOBLIS
#include "transform.h"
#include "DLAOp.h"
#include "elemRedist.h"
#include "blis.h"
#include "lowerLayer.h"

Loop* HemmLoopVar4(Node *Ain, unsigned int Anum,
		   Node *Bin, unsigned int Bnum,
		   Node *Cin, unsigned int Cnum,
		   Side side, Tri tri,
		   Coef alpha, Coef beta,
		   Type type,
		   Layer layer);

Loop* HemmLoopVar8(Node *Ain, unsigned int Anum,
		   Node *Bin, unsigned int Bnum,
		   Node *Cin, unsigned int Cnum,
		   Side side, Tri tri,
		   Coef alpha, Coef beta,
		   Type type,
		   Layer layer);

#if DOELEM
//This has an optimization found in Elemental hardcoded
// into the algorithm.  Think about it as a simplifier
// I manually applied.
Loop* HemmLoopVar8Altered(Node *Ain, unsigned int Anum,
			  Node *Bin, unsigned int Bnum,
			  Node *Cin, unsigned int Cnum,
			  Side side, Tri tri,
			  Coef alpha, Coef beta,
			  Type type,
			  Layer layer);
#endif



class Hemm : public DLAOp<3,1>
{
 public:
  Side m_side;
  Tri m_tri;
  Coef m_alpha, m_beta;
  Type m_type;
  Hemm(Layer layer, Side side, Tri tri, Coef alpha, Coef beta, Type type);
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "Hemm";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual NodeType GetType() const;
  static Node* BlankInst() { return new Hemm(ABSLAYER, LEFT, LOWER, COEFONE, COEFONE, REAL); }
  virtual Node* GetNewInst() { return BlankInst(); }
#if DOELEM
  virtual const DistType& GetDistType(unsigned int num) const;
#endif
  virtual void SanityCheck();
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual Phase MaxPhase() const;
  static Cost GetCost(Layer layer, Side side, const Sizes *localMs, const Sizes *localNs);
  virtual bool ShouldCullDP() const;
};

class HemmLoopExp : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  int m_var;
 HemmLoopExp(Layer fromLayer, Layer toLayer, int var) 
   : m_fromLayer(fromLayer), m_toLayer(toLayer), m_var(var) {}
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};

#if DOELEM
class DistHemmToLocalHemm : public SingleTrans
{
  DistType m_leftType, m_rightType;
 public:
 DistHemmToLocalHemm(DistType leftType, DistType rightType) :m_leftType(leftType),m_rightType(rightType) {}
  virtual string GetType() const {return "Distributed Hemm to Local Hemm " + DistTypeToStr(m_leftType) + "," + DistTypeToStr(m_rightType);}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
  virtual Cost RHSCostEstimate(const Node *node) const;
};
#endif

#if DOELEM
//Special computation kernel in Elemental
class LocalSymmAcc : public DLAOp<5,2>
{
 public:
  Side m_side;
  Tri m_tri;
  Type m_type;
  Coef m_alpha;
  LocalSymmAcc(Side side, Tri tri, Type type, Coef alpha);
  virtual NodeType GetType() const { return "LocalSymmAcc"; }
  static Node* BlankInst() { return  new LocalSymmAcc(LEFT, UPPER, COMPLEX, COEFONE); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual const DistType& GetDistType(unsigned int num) const { return InputDistType(3+num); }
  virtual void SanityCheck();
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "localSymmAcc";}
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  static Cost GetCost(Side side, const Sizes *localMs, const Sizes *localNs);
};


class DistHemmToLocalHemmStatA : public SingleTrans
{
 public:
  virtual string GetType() const {return "Distributed Hemm to Local Hemm stat A";}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
  virtual Cost RHSCostEstimate(const Node *node) const;
};
#endif

#if DOBLIS
class BLISHemmLoopExp : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
 BLISHemmLoopExp(Layer fromLayer, Layer toLayer) 
   : m_fromLayer(fromLayer), m_toLayer(toLayer) {}
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  virtual bool IsRef() const {return true;}
};

class BLISHemmToGemm: public SingleTrans
{
 public:
  Layer m_layer;
  BLISHemmToGemm(Layer layer)
    : m_layer(layer) {}
  virtual string GetType() const { return "Hemm to Gemm"; }
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};
#endif

class HemmLowerLayer : public LowerLayer
{
 public:
 HemmLowerLayer(Layer fromLayer, Layer toLayer, DimName dim, Size bs)
   : LowerLayer(fromLayer, toLayer, dim, bs) {}
  virtual string GetType() const;
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};
#endif
